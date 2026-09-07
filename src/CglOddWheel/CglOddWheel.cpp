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
#include <algorithm>
#include <utility>
#include <vector>

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

/** Dump one odd wheel's structure as a JSON line, for offline inspection and
 *  drawing. Off unless CGL_ODDWHEEL_DUMP names a file, and deliberately driven
 *  by the environment rather than by a setter: a new Stats field or accessor
 *  changes CglOddWheel's ABI, and Cgl and Cbc are built and installed
 *  separately here, so a mismatched pair segfaults.
 *
 *  CGL_ODDWHEEL_DUMP_CYCLE=n   only wheels whose cycle has exactly n nodes
 *  CGL_ODDWHEEL_DUMP_CENTRE=1  only wheels that actually received a centre
 *
 *  Everything written is read straight off the graph and the value vector, so
 *  the drawing is a picture of what the separator saw, not a reconstruction:
 *  `conflicts` is the induced subgraph on cycle+centre as cgraph reports it,
 *  and `lhs`/`rhs` are the emitted row evaluated at x. */
static void dumpOddWheel(const CoinConflictGraph *cgraph, size_t numCols,
  const double *x, const size_t *cycle, size_t cycleSize,
  const size_t *centre, size_t centreSize, double alpha, double cycleRhs,
  const int *idxs, const double *coefs, int cutSize, double rhs)
{
    const char *path = getenv("CGL_ODDWHEEL_DUMP");
    if (!path)
        return;
    const char *wantCycle = getenv("CGL_ODDWHEEL_DUMP_CYCLE");
    if (wantCycle && (size_t)atoi(wantCycle) != cycleSize)
        return;
    const char *wantCentre = getenv("CGL_ODDWHEEL_DUMP_CENTRE");
    if (wantCentre && atoi(wantCentre) != 0 && centreSize == 0)
        return;

    FILE *f = fopen(path, "a");
    if (!f)
        return;

    /* cycle then centre, in one array, so the conflict block below can index
     * both halves uniformly. */
    std::vector< size_t > nodes;
    nodes.reserve(cycleSize + centreSize);
    for (size_t k = 0; k < cycleSize; k++)
        nodes.push_back(cycle[k]);
    for (size_t k = 0; k < centreSize; k++)
        nodes.push_back(centre[k]);

    fprintf(f, "{\"cycleSize\":%lu,\"centreSize\":%lu,\"alpha\":%g,"
               "\"cycleRhs\":%g,\"rhs\":%g,\"numCols\":%lu",
      (unsigned long)cycleSize, (unsigned long)centreSize, alpha, cycleRhs,
      rhs, (unsigned long)numCols);

    fprintf(f, ",\"nodes\":[");
    for (size_t k = 0; k < nodes.size(); k++) {
        const bool compl_ = (nodes[k] >= numCols);
        const size_t col = compl_ ? (nodes[k] - numCols) : nodes[k];
        /* z is the doubled-graph value the separator ranked on: x_col for a
         * plain node, 1 - x_col for a complemented one. */
        const double z = compl_ ? (1.0 - x[col]) : x[col];
        fprintf(f, "%s{\"node\":%lu,\"col\":%lu,\"compl\":%d,\"x\":%.17g,"
                   "\"z\":%.17g,\"role\":\"%s\"}",
          k ? "," : "", (unsigned long)nodes[k], (unsigned long)col,
          compl_ ? 1 : 0, x[col], z, k < cycleSize ? "cycle" : "centre");
    }
    fprintf(f, "]");

    fprintf(f, ",\"conflicts\":[");
    bool first = true;
    for (size_t a = 0; a < nodes.size(); a++)
        for (size_t b = a + 1; b < nodes.size(); b++)
            if (cgraph->conflicting(nodes[a], nodes[b])) {
                fprintf(f, "%s[%lu,%lu]", first ? "" : ",",
                  (unsigned long)a, (unsigned long)b);
                first = false;
            }
    fprintf(f, "]");

    double lhs = 0.0;
    fprintf(f, ",\"cut\":[");
    for (int k = 0; k < cutSize; k++) {
        fprintf(f, "%s{\"col\":%d,\"coef\":%g,\"x\":%.17g}",
          k ? "," : "", idxs[k], coefs[k], x[idxs[k]]);
        lhs += coefs[k] * x[idxs[k]];
    }
    fprintf(f, "],\"lhs\":%.17g,\"viol\":%.17g}\n", lhs, lhs - rhs);

    fclose(f);
}

CglOddWheel::CglOddWheel(size_t extMethod) : cap_(0), extMethod_(extMethod), verifyPrepare_(false), useGate_(true), checkValidity_(false), stats_(Stats()) {
    idxs_ = NULL;
    idxMap_ = NULL;
    coefs_ = NULL;
    x_ = NULL;
    rc_ = NULL;
}

CglOddWheel::CglOddWheel(const CglOddWheel& rhs) {
    this->cap_ = rhs.cap_;
    this->extMethod_ = rhs.extMethod_;
    this->verifyPrepare_ = rhs.verifyPrepare_;
    this->useGate_ = rhs.useGate_;
    this->checkValidity_ = rhs.checkValidity_;
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
    if (verifyPrepare_)
        oddH.setVerifyPrepare(true);
    if (!useGate_)
        oddH.setUseFutilityGate(false);
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
        /* k = floor(|C|/2): the right-hand side of the odd-cycle inequality in
         * the conflict graph's own z-space, where z_j = x_j and
         * z_{j+numCols} = 1 - x_j. Kept separate from rhs, which the loops below
         * translate into x-space as they go. */
        const double cycleRhs = oddH.oddWheelRHS(j);
        double rhs = cycleRhs;

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
         * value vector.
         *
         * Rejecting the cut instead is expensive, and this is the current
         * measurement of how much: over the 237 replay fixtures at --rounds=1,
         * 331 of 1228 odd wheels (27%) have a column appearing both plain and
         * complemented, and none of the 1228 cancels away entirely (cutsEmpty
         * 0). setCheckValidity() re-derives all 1228 from the graph and agrees
         * on the count (certComplPair 331 == cutsDuplicatedIdx 331) and on every
         * coefficient. Complemented nodes are not a corner case at all here:
         * 1170 of the 1228 cycles (95%) traverse at least one. At --rounds=4 the
         * share only grows -- 512 of 1633 (31%) repeated, 1567 (96%) with a
         * complemented node -- so this is not an artifact of one cut pass.
         *
         * An earlier version of this comment claimed 2528 of 2528, which was a
         * *pre-fix* figure -- CoinOddWheelSeparator::addOddHole() used to store
         * the whole scratch buffer, so every hole tripped this guard for a
         * reason unrelated to complements. Post-fix, 2528 splits as 1227 holes
         * plus 1301 duplicates.
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
        const double alpha = centerSize ? cycleRhs : 0.0;
        if (centerSize) {
            /* The odd-wheel inequality in z-space is
             *
             *     sum_{v in C} z_v  +  alpha * sum_{w in W} z_w  <=  k
             *
             * and it is valid for every alpha <= k. W is a clique, so at most
             * one z_w is 1; each w conflicts with all of C, so z_w = 1 forces
             * every z_v to 0 and the left-hand side is alpha <= k. Otherwise
             * every z_w is 0 and the left-hand side is the odd-cycle sum, at
             * most k. alpha = k is therefore the strongest valid choice, and
             * alpha > k is cut off by exactly the z_w = 1 point.
             *
             * The lifting coefficient has to be that z-space k, which is why
             * cycleRhs is captured before the cycle loop runs. rhs has by now
             * absorbed one -1 per complement in the cycle, so lifting with it
             * used alpha = k - |C-|: still valid, but weaker on every cycle that
             * traverses a complement -- 95% of them here -- and for |C-| > k
             * actually negative, which makes the wheel weaker than the plain odd
             * cycle it came from. The `fabs(rhs) >= ODDHWC_EPS` guard that used
             * to sit on this branch existed to keep the degenerate alpha = 0 case
             * out, and dropped the whole centre when it fired; alpha = k >= 2 for
             * any cycle of length >= 5, so it can no longer fire and is gone.
             *
             * Raising alpha cannot lose a cut: in z-space the right-hand side
             * does not move and the left-hand side grows by
             * sum_{w in W} z*_w >= 0, so the violation the cut pool measures is
             * non-decreasing.
             *
             * A centre whose column already appears in the cycle accumulates for
             * the same reason as above -- and with alpha = k >= 2 such a pair can
             * no longer cancel to zero, since the cycle contributes +-1 and the
             * centre +-k. */
            for (size_t k = 0; k < centerSize; k++) {
                const bool complement = (centerIdx[k] >= numCols);
                const int col = complement ? ((int)(centerIdx[k] - numCols)) : ((int)centerIdx[k]);
                const double coef = complement ? (-1.0 * alpha) : alpha;

                if (complement) {
                    rhs -= alpha;
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

        if (checkValidity_)
            certifyOddWheel(cgraph, numCols, oddEl, oddSize, centerIdx,
              centerSize, alpha, idxs_, coefs_, realSize, rhs);

        dumpOddWheel(cgraph, numCols, x_, oddEl, oddSize, centerIdx, centerSize,
          alpha, cycleRhs, idxs_, coefs_, realSize, rhs);

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

void CglOddWheel::certifyOddWheel(const CoinConflictGraph *cgraph, size_t numCols,
  const size_t *cycle, size_t cycleSize,
  const size_t *center, size_t centerSize, double alpha,
  const int *idxs, const double *coefs, int nz, double rhs)
{
    stats_.certChecked++;

    /* k is what makes the inequality valid, and it is a property of the cycle
     * alone -- not of `rhs`, which the complement translation has already moved
     * by -1 per complemented node. */
    const double k = floor(cycleSize / 2.0);

    /* --- the cycle: odd, distinct nodes, consecutive ones in conflict ------ */
    bool badCycle = (cycleSize < 5) || ((cycleSize % 2) == 0);
    size_t nCompCycle = 0;
    for (size_t i = 0; i < cycleSize && !badCycle; i++) {
        if (cycle[i] >= numCols)
            nCompCycle++;
        for (size_t j = i + 1; j < cycleSize; j++) {
            if (cycle[i] == cycle[j]) {
                badCycle = true;
                break;
            }
        }
        if (!badCycle && !cgraph->conflicting(cycle[i], cycle[(i + 1) % cycleSize]))
            badCycle = true;
    }
    if (badCycle) {
        stats_.certBadCycle++;
        return;
    }
    if (nCompCycle)
        stats_.certComplCycle++;
    if (nCompCycle && centerSize)
        stats_.certCenterOnComplCycle++;
    if ((double)nCompCycle >= k)
        stats_.certComplAtLeastK++;

    /* --- the centres: adjacent to all of C, and a clique among themselves -- */
    size_t nCompCenter = 0;
    for (size_t i = 0; i < centerSize; i++) {
        if (center[i] >= numCols)
            nCompCenter++;
        for (size_t v = 0; v < cycleSize; v++) {
            if (!cgraph->conflicting(center[i], cycle[v])) {
                stats_.certBadCenterAdj++;
                return;
            }
        }
        for (size_t j = i + 1; j < centerSize; j++) {
            if (!cgraph->conflicting(center[i], center[j])) {
                stats_.certBadCenterClq++;
                return;
            }
        }
    }
    if (nCompCenter)
        stats_.certComplCenter++;
    if (centerSize && alpha > k + ODDHWC_EPS) {
        stats_.certBadAlpha++;
        return;
    }

    /* --- the translation: re-derive the column-space cut from scratch ------
     *
     * Node j means x_j = 1 and node j + numCols means x_j = 0, i.e. z_{j+numCols}
     * = 1 - x_j, so a complemented node contributes -coef to column j and -coef
     * to the right-hand side. Accumulating per column is what makes a wheel that
     * uses both a variable and its own complement come out right: the two
     * coefficients cancel and the rhs keeps both shifts. Rebuilding it here
     * rather than re-reading idxMap_ is the point -- it is an independent
     * computation of the same thing. */
    const size_t nNodes = cycleSize + centerSize;
    std::vector<std::pair<size_t, double> > exp;
    exp.reserve(nNodes);
    double expRhs = k;
    for (size_t i = 0; i < nNodes; i++) {
        const size_t nd = (i < cycleSize) ? cycle[i] : center[i - cycleSize];
        const double w = (i < cycleSize) ? 1.0 : alpha;
        const bool comp = (nd >= numCols);
        exp.push_back(std::make_pair(comp ? (nd - numCols) : nd, comp ? -w : w));
        if (comp)
            expRhs -= w;
    }
    std::sort(exp.begin(), exp.end());

    /* A column present both plain and complemented is the edge case that used to
     * be rejected outright; count it so its frequency is a measured number
     * rather than an argument. After the sort it is simply a repeated column. */
    for (size_t i = 1; i < exp.size(); i++) {
        if (exp[i].first == exp[i - 1].first) {
            stats_.certComplPair++;
            break;
        }
    }

    /* Merge repeats, drop exact zeros -- the same two steps the emitted cut has
     * already taken -- then compare column by column. */
    std::vector<std::pair<size_t, double> > merged;
    for (size_t i = 0; i < exp.size(); i++) {
        if (!merged.empty() && merged.back().first == exp[i].first)
            merged.back().second += exp[i].second;
        else
            merged.push_back(exp[i]);
    }
    bool bad = fabs(expRhs - rhs) > ODDHWC_EPS;
    std::vector<std::pair<size_t, double> > got;
    got.reserve((size_t)nz);
    for (int i = 0; i < nz; i++)
        got.push_back(std::make_pair((size_t)idxs[i], coefs[i]));
    std::sort(got.begin(), got.end());

    size_t g = 0;
    for (size_t i = 0; i < merged.size() && !bad; i++) {
        if (fabs(merged[i].second) < ODDHWC_EPS)
            continue;
        if (g >= got.size() || got[g].first != merged[i].first
          || fabs(got[g].second - merged[i].second) > ODDHWC_EPS)
            bad = true;
        g++;
    }
    if (!bad && g != got.size())
        bad = true;
    if (bad)
        stats_.certBadTranslate++;
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
