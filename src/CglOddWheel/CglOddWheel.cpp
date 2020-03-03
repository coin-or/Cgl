#include <cstdio>
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

CglOddWheel::CglOddWheel(size_t extMethod) : cap_(0), extMethod_(extMethod) {
    idxs_ = NULL;
    idxMap_ = NULL;
    coefs_ = NULL;
    x_ = NULL;
    rc_ = NULL;
}

CglOddWheel::CglOddWheel(const CglOddWheel& rhs) {
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

    this->cgraph_ = rhs.cgraph_;
    this->cap_ = rhs.cap_;
    this->idxs_ = (int*)xmalloc(sizeof(int) * this->cap_);
    this->idxMap_ = (int*)xmalloc(sizeof(int) * this->cap_);
    this->coefs_ = (double*)xmalloc(sizeof(double) * this->cap_);
    this->x_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    this->rc_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
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

CglCutGenerator * CglOddWheel::clone() const {
    CglOddWheel *oddHWC = new CglOddWheel();

    oddHWC->cgraph_ = this->cgraph_;
    oddHWC->cap_ = this->cap_;
    oddHWC->idxs_ = (int*)xmalloc(sizeof(int) * this->cap_);
    std::copy(this->idxs_, this->idxs_ + this->cap_, oddHWC->idxs_);
    oddHWC->idxMap_ = (int*)xmalloc(sizeof(int) * this->cap_);
    std::copy(this->idxMap_, this->idxMap_ + this->cap_, oddHWC->idxMap_);
    oddHWC->coefs_ = (double*)xmalloc(sizeof(double) * this->cap_);
    std::copy(this->coefs_, this->coefs_ + this->cap_, oddHWC->coefs_);

    oddHWC->x_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    std::copy(this->x_, this->x_ + (this->cap_ * 2), oddHWC->x_);
    oddHWC->rc_ = (double*)xmalloc(sizeof(double) * this->cap_);
    std::copy(this->rc_, this->rc_ + (this->cap_ * 2), oddHWC->rc_);

    return static_cast<CglCutGenerator*>(oddHWC);
}

void CglOddWheel::generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info ) {
    double startSep = CoinCpuTime();
    const size_t numCols = si.getNumCols();

	if(numCols != cgraph_->size() / 2) {
        fprintf(stderr, "Invalid conflict graph! Number of columns %ld ... in graph %ld\n", numCols, cgraph_->size() / 2);
        exit(EXIT_FAILURE);
    }

    checkMemory(numCols);

    const double *colSol = si.getColSolution();
    const double *rCost = si.getReducedCost();
    for(size_t i = 0; i < numCols; i++) {
        x_[i] = colSol[i];
        rc_[i] = rCost[i];
        x_[i + numCols] = 1.0 - x_[i];
        rc_[i + numCols] = -rc_[i];
    }

    CoinOddWheelSeparator oddH(cgraph_, x_, rc_, extMethod_);
    CoinCutPool cutPool(x_, numCols);

    oddH.searchOddWheels();

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
        size_t dup = 0;

        std::fill(idxMap_, idxMap_ + numCols, -1);

        for(size_t k = 0; k < oddSize; k++) {
            if(oddEl[k] < numCols) {
#ifdef DEBUGCG
                assert(oddEl[k] < std::numeric_limits<int>::max());
#endif
                if(idxMap_[oddEl[k]] == -1) {
                    idxMap_[oddEl[k]] = realSize;
                    idxs_[realSize] = oddEl[k];
                    coefs_[realSize] = 1.0;
                    realSize++;
                }
                else {
                    coefs_[idxMap_[oddEl[k]]] += 1.0;
                    dup++;
                }
            }
            else {
#ifdef DEBUGCG
                assert(oddEl[k] - numCols < std::numeric_limits<int>::max());
#endif
                if(idxMap_[oddEl[k]-numCols] == -1) {
                    idxMap_[oddEl[k]-numCols] = realSize;
                    idxs_[realSize] = ((int)(oddEl[k] - numCols));
                    coefs_[realSize] = -1.0;
                    realSize++;
                }
                else {
                    coefs_[idxMap_[oddEl[k]-numCols]] -= 1.0;
                    dup++;
                }
                rhs = rhs - 1.0;
            }
        }

        const size_t centerSize = oddH.wheelCenterSize(j);
        const size_t *centerIdx = oddH.wheelCenter(j);
        if (centerSize && fabs(rhs) >= ODDHWC_EPS) {
            const double oldRhs = rhs;
            for (size_t k = 0; k < centerSize; k++) {
                if (centerIdx[k] < numCols) {
#ifdef DEBUGCG
                    assert(centerIdx[k] < std::numeric_limits<int>::max());
#endif
                    if (idxMap_[centerIdx[k]] == -1) {
                        idxMap_[centerIdx[k]] = realSize;
                        idxs_[realSize] = centerIdx[k];
                        coefs_[realSize] = oldRhs;
                        realSize++;
                    } else {
                        coefs_[idxMap_[centerIdx[k]]] += oldRhs;
                        dup++;
                    }
                } else {
#ifdef DEBUGCG
                    assert(centerIdx[k] - numCols < std::numeric_limits<int>::max());
#endif
                    if (idxMap_[centerIdx[k] - numCols] == -1) {
                        idxMap_[centerIdx[k] - numCols] = realSize;
                        idxs_[realSize] = ((int) (centerIdx[k] - numCols));
                        coefs_[realSize] = -1.0 * oldRhs;
                        realSize++;
                    } else {
                        coefs_[idxMap_[centerIdx[k] - numCols]] -= oldRhs;
                        dup++;
                    }
                    rhs = rhs - oldRhs;
                }
            }
        }

        if(dup) {
            int last = 0;
            for(int k = 0; k < realSize; k++) {
                if(fabs(coefs_[k]) >= ODDHWC_EPS) {
                    idxs_[last] = idxs_[k];
                    coefs_[last] = coefs_[k];
                    last++;
                }
            }
            realSize = last;
        }

        cutPool.add(idxs_, coefs_, realSize, rhs);
    }

    cutPool.removeNullCuts();

    const size_t numberRowCutsBefore = cs.sizeRowCuts();
    for(size_t i = 0; i < cutPool.numCuts(); i++) {
        osrc_.setRow(cutPool.cutSize(i) , cutPool.cutIdxs(i), cutPool.cutCoefs(i));
        osrc_.setUb(cutPool.cutRHS(i));
        cs.insertIfNotDuplicate(osrc_);
    }

    int numberRowCutsAfter = cs.sizeRowCuts();
    CglOddWheel::sepCuts += numberRowCutsAfter - numberRowCutsBefore;

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