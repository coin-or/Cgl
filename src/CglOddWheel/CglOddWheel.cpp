/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class for separating violated odd-cycles. It contains
 * a lifting module that tries to transform the odd-cycles
 * into odd-wheels.
 *
 * @file CglOddWheel.cpp
 * @brief Odd-wheel cut separator
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "CglOddWheel.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "CoinTime.hpp"
#include "CoinCutPool.hpp"
#include "CoinOddWheelSeparator.hpp"

#define ODDHWC_EPS 1e-6

size_t CglOddWheel::sepCuts = 0;
double CglOddWheel::sepTime = 0.0;

static void *xmalloc( const size_t size );

CglOddWheel::CglOddWheel(size_t extMethod) : cap_(0), extMethod_(extMethod), stats_(Stats()) {
    idxs_ = NULL;
    idxMap_ = NULL;
    coefs_ = NULL;
    x_ = NULL;
    rc_ = NULL;
}

CglOddWheel::CglOddWheel(const CglOddWheel& rhs) {
    this->cap_ = rhs.cap_;
    this->extMethod_ = rhs.extMethod_;
    // Not copied: a clone has made no call of its own yet.
    this->stats_ = Stats();

    if (this->cap_ > 0) {
        this->idxs_ = (int*)xmalloc(sizeof(int) * this->cap_);
        this->idxMap_ = (int*)xmalloc(sizeof(int) * this->cap_);
        this->coefs_ = (double*)xmalloc(sizeof(double) * this->cap_);
        this->x_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
        this->rc_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    } else {
        this->idxs_ = NULL;
        this->idxMap_ = NULL;
        this->coefs_ = NULL;
        this->x_ = NULL;
        this->rc_ = NULL;
    }
}

CglOddWheel::~CglOddWheel() {
    if (this->idxs_) {
        free(this->idxs_);
    }
    if (this->idxMap_) {
        free(this->idxMap_);
    }
    if (this->coefs_) {
        free(this->coefs_);
    }
    if (this->x_) {
        free(this->x_);
    }
    if (this->rc_) {
        free(this->rc_);
    }
}

void CglOddWheel::refreshSolver(OsiSolverInterface *solver) {
    solver->checkCGraph();
  // Get integer information
    solver->getColType(true);
}

CglCutGenerator * CglOddWheel::clone() const {
    return new CglOddWheel(*this);
}

void CglOddWheel::generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info ) {
    if (si.getNumCols() == 0 || si.getNumRows() == 0) {
        return;
    }

    stats_ = Stats();
    double startSep = CoinCpuTime();
    const size_t numCols = si.getNumCols();
    const CoinConflictGraph *cgraph = si.getCGraph();

	if(numCols != cgraph->size() / 2) {
        fprintf(stderr, "Invalid conflict graph! Number of columns %ld ... in graph %ld\n", numCols, cgraph->size() / 2);
        exit(EXIT_FAILURE);
    }

    checkMemory(numCols);

    const double startSetup = CoinGetTimeOfDay();
    const double *colSol = si.getColSolution();
    const double *rCost = si.getReducedCost();
    for(size_t i = 0; i < numCols; i++) {
        x_[i] = colSol[i];
        rc_[i] = rCost[i];
        x_[i + numCols] = 1.0 - x_[i];
        rc_[i + numCols] = -rc_[i];
    }

    // tSeparator covers the separator's construction (which runs
    // fillActiveColumns) and its search; the cut pool's allocation sits in
    // between and is counted with it rather than reordering the code.
    const double startSeparator = CoinGetTimeOfDay();
    stats_.tSetup = startSeparator - startSetup;

    CoinOddWheelSeparator oddH(cgraph, x_, rc_, extMethod_);
    if (maxSeconds_ > 0.0)
        oddH.setMaxSeconds(maxSeconds_);
    CoinCutPool cutPool(x_, numCols, "OddWheel");

    oddH.searchOddWheels();

    const double startCutPool = CoinGetTimeOfDay();
    stats_.tSeparator = startCutPool - startSeparator;
    stats_.sep = oddH.stats();

    // Same rationale as CglBKClique::insertCuts(): only pay for the
    // per-column best-score filtering when there are enough candidates
    // for it to actually matter, and exempt small models outright.
  const char *minColsEnv = getenv("CBC_CLIQUE_POOL_MIN_COLS");
  const size_t minCols = minColsEnv ? (size_t)strtol(minColsEnv, nullptr, 10) : 500;
  const bool smallModel = numCols < minCols;

  const char *alwaysFilterEnv = getenv("CBC_CLIQUE_POOL_ALWAYS_FILTER");
  if (alwaysFilterEnv && atoi(alwaysFilterEnv) != 0) {
    cutPool.setFilteringEnabled(true);
  } else {
    const char *minCandEnv = getenv("CBC_CLIQUE_POOL_MIN_CANDIDATES");
    const size_t minCandidates = minCandEnv ? (size_t)strtol(minCandEnv, nullptr, 10) : 20;
    cutPool.setFilteringEnabled(!smallModel && oddH.numOddWheels() >= minCandidates);
  }

  // Orthogonality/parallelism-based cut selection (see CglBKClique for
  // rationale and A/B benchmark result). Disabled (1.0) by default; opt
  // in via CBC_CLIQUE_POOL_MAX_PARALLELISM for further experimentation.
  const bool forceFilter = alwaysFilterEnv && atoi(alwaysFilterEnv) != 0;
  const char *maxParEnv = getenv("CBC_CLIQUE_POOL_MAX_PARALLELISM");
  const double maxPar = maxParEnv ? atof(maxParEnv) : 1.0;
  cutPool.setMaxParallelism((smallModel && !forceFilter) ? 1.0 : maxPar);

    /* adding odd holes */
    for(size_t j = 0; j < oddH.numOddWheels(); j++) {
        const size_t *oddEl = oddH.oddHole(j);
        const size_t oddSize = oddH.oddHoleSize(j);
        double rhs = oddH.oddWheelRHS(j);

        if(oddSize < 5) {
            fprintf(stderr, "Invalid size of cut: %lu\n", oddSize);
            exit(EXIT_FAILURE);
        }

        int realSize = 0;
        size_t duplicated = 0;
        std::fill(idxMap_, idxMap_ + numCols, -1);

        /* Translating conflict graph nodes into columns.
         *
         * A node j means x_j = 1 and a node j + numCols means x_j = 0, and
         * addVariableComplementConflicts() makes those two nodes adjacent, so a
         * cycle may legitimately pass through a variable *and* its own
         * complement. Both then map to column j. Summing the coefficients is
         * exact algebra on x_j + (1 - x_j): the pair cancels to zero and the
         * RHS keeps its -1, leaving an inequality that is still valid and still
         * violated by exactly the amount the separator measured on the doubled
         * value vector. Rejecting the cut instead discarded every odd hole
         * found on all 237 replay fixtures (2528 of 2528).
         *
         * With no repeated column there is no accumulation and no compaction,
         * so cuts that were emitted before come out unchanged. */
        for(size_t k = 0; k < oddSize; k++) {
            const bool complement = (oddEl[k] >= numCols);
            const int col = complement ? ((int)(oddEl[k] - numCols)) : ((int)oddEl[k]);
            const double coef = complement ? -1.0 : 1.0;

            if (complement) {
                rhs -= 1.0;
            }

            if(idxMap_[col] == -1) {
                idxMap_[col] = realSize;
                idxs_[realSize] = col;
                coefs_[realSize] = coef;
                realSize++;
            } else {
                coefs_[idxMap_[col]] += coef;
                duplicated++;
            }
        }

        const size_t centerSize = oddH.wheelCenterSize(j);
        const size_t *centerIdx = oddH.wheelCenter(j);
        if (centerSize && fabs(rhs) >= ODDHWC_EPS) {
            const double oldRhs = rhs;
            /* The wheel centres form a clique, so at most one of them is 1 and
             * lifting them all with coefficient oldRhs stays valid. A centre
             * whose column already appears in the cycle accumulates for the
             * same reason as above; dropping the whole wheel over it would
             * throw away the plain odd-hole cut as well. */
            for (size_t k = 0; k < centerSize; k++) {
                const bool complement = (centerIdx[k] >= numCols);
                const int col = complement ? ((int)(centerIdx[k] - numCols)) : ((int)centerIdx[k]);
                const double coef = complement ? (-1.0 * oldRhs) : oldRhs;

                if (complement) {
                    rhs -= oldRhs;
                }

                if (idxMap_[col] == -1) {
                    idxMap_[col] = realSize;
                    idxs_[realSize] = col;
                    coefs_[realSize] = coef;
                    realSize++;
                } else {
                    coefs_[idxMap_[col]] += coef;
                    duplicated++;
                }
            }
        }

        if (duplicated) {
            /* Every coefficient here is integral -- +-1 from the cycle and
             * +-floor(|C|/2) from the centres -- so testing against ODDHWC_EPS
             * separates an exact zero from the smallest survivor. */
            stats_.cutsDuplicatedIdx++;
            int keep = 0;
            for (int k = 0; k < realSize; k++) {
                if (fabs(coefs_[k]) >= ODDHWC_EPS) {
                    idxs_[keep] = idxs_[k];
                    coefs_[keep] = coefs_[k];
                    keep++;
                }
            }
            stats_.cutsZeroCoefs += (size_t)(realSize - keep);
            realSize = keep;
        }

        if (realSize == 0) {
            /* Nothing left on the left-hand side. Since the row was violated,
             * `0 <= rhs` with rhs < 0 would be an infeasibility claim, and that
             * is not a cut separator's call to make. */
            stats_.cutsEmpty++;
            continue;
        }

        stats_.cutsBeforePool++;
        cutPool.add(idxs_, coefs_, realSize, rhs);
    }

    cutPool.removeNullCuts();
    cutPool.filterByParallelism();
    stats_.cutsAfterPool = cutPool.numCuts();

    const size_t numberRowCutsBefore = cs.sizeRowCuts();
    for(size_t i = 0; i < cutPool.numCuts(); i++) {
        osrc_.setRow(cutPool.cutSize(i) , cutPool.cutIdxs(i), cutPool.cutCoefs(i));
        osrc_.setUb(cutPool.cutRHS(i));
        cs.insertIfNotDuplicate(osrc_);
    }

    int numberRowCutsAfter = cs.sizeRowCuts();
    CglOddWheel::sepCuts += numberRowCutsAfter - numberRowCutsBefore;
    stats_.rowCutsAdded = numberRowCutsAfter - numberRowCutsBefore;
    stats_.tCutPool = CoinGetTimeOfDay() - startCutPool;

    if(!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
        numberRowCutsAfter = cs.sizeRowCuts();
        for(int i = numberRowCutsBefore; i < numberRowCutsAfter; i++) {
            cs.rowCutPtr(i)->setGloballyValid();
        }
    }

	CglOddWheel::sepTime += (CoinCpuTime() - startSep);
}

void CglOddWheel::checkMemory(const size_t newNumCols) {
    if (cap_ < newNumCols) {
        if (cap_ > 0) {
#ifdef DEBUGCG
            assert(idxs_);
            assert(idxMap_);
            assert(coefs_);
            assert(x_);
            assert(rc_);
#endif
            free(idxs_);
            free(idxMap_);
            free(coefs_);
            free(x_);
            free(rc_);
        }

        idxs_ = (int*)xmalloc(sizeof(int) * newNumCols);
        idxMap_ = (int*)xmalloc(sizeof(int) * newNumCols);
        coefs_ = (double*)xmalloc(sizeof(double) * newNumCols);
        x_ = (double*)xmalloc(sizeof(double) * newNumCols * 2);
        rc_ = (double*)xmalloc(sizeof(double) * newNumCols * 2);
        cap_ = newNumCols;
    }
}

static void *xmalloc( const size_t size ) {
    void *result = malloc( size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
        abort();
    }

    return result;
}

void CglOddWheel::setExtendingMethod(size_t extMethod) {
    if(extMethod > 2) {
        fprintf(stderr, "Invalid value for parameter extMethod (%lu).\n", extMethod);
        abort();
    }

    this->extMethod_ = extMethod;
}
