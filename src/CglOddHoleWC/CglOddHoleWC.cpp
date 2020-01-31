#include <cstdio>
#include <cassert>

#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "CoinTime.hpp"
#include "CglOddHoleWC.hpp"
#include "cut.h"
#include "oddhs.h"

size_t CglOddHoleWC::sepCuts = 0;
double CglOddHoleWC::sepTime = 0.0;

CglOddHoleWC::CglOddHoleWC() {
}

CglOddHoleWC::CglOddHoleWC(const CglOddHoleWC& rhs) {
}

CglCutGenerator * CglOddHoleWC::clone() const {
    CglOddHoleWC *aClq = new CglOddHoleWC();

    return static_cast<CglCutGenerator*>(aClq);
}

void CglOddHoleWC::generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info ) {
    double startSep = CoinCpuTime();
    const size_t numCols = si.getNumCols();
    const CGraph *cg = si.getCGraph();
    double *x = new double[numCols*2];
    double *rc = new double[numCols*2];
    int *idx = new int[numCols*2];
    int *idxMap = new int[numCols*2];
    double *coef = new double[numCols*2];
    OddHoleSep *oddhs = oddhs_create(cg);
    const CoinAbsFltEq equal(1.0e-12);
    OsiRowCut osrc;

	if(numCols != cgraph_size(cg)/2) {
        fprintf(stderr, "Invalid conflict graph! Number of columns %ld ... in graph %ld\n",
                numCols, cgraph_size(cg)/2);
        exit(EXIT_FAILURE);
    }

    const double *colSol = si.getColSolution();
    const double *rCost = si.getReducedCost();
    for(size_t i = 0; i < numCols; i++) {
        x[i] = colSol[i];
        rc[i] = rCost[i];
        x[i + numCols] = 1.0 - x[i];
        rc[i + numCols] = -rc[i];
    }

    CutPool *cutPool = cut_pool_create(numCols);
    oddhs_search_odd_holes(oddhs, x, rc);

    /* adding odd holes */
    for(size_t j = 0; j < oddhs_get_odd_hole_count(oddhs); j++) {
        const std::vector<size_t> &oddEl = oddhs_get_odd_hole(oddhs, j);
        const size_t oddSize = oddEl.size();
        double rhs = oddhs_get_odd_hole_rhs(oddhs, j);

        if(oddSize < 5) {
            fprintf(stderr, "Invalid size of cut: %lu\n", oddSize);
            exit(EXIT_FAILURE);
        }

        const size_t centerSize = oddhs_get_nwc_doh(oddhs, j);
        const std::vector<size_t> &centerIdx = oddhs_get_wc_doh(oddhs, j);

        const size_t cutSize = oddSize + centerSize;
        int realSize = 0;
        size_t dup = 0;

        std::fill(idxMap, idxMap + numCols, -1);

        for(size_t k = 0; k < oddSize; k++) {
            if(oddEl[k] < numCols) {
#ifdef DEBUG
                assert(oddEl[k] < std::numeric_limits<int>::max());
#endif
                if(idxMap[oddEl[k]] == -1) {
                    idxMap[oddEl[k]] = realSize;
                    idx[realSize] = oddEl[k];
                    coef[realSize] = 1.0;
                    realSize++;
                }
                else {
                    coef[idxMap[oddEl[k]]] += 1.0;
                    dup++;
                }
            }
            else {
#ifdef DEBUG
                assert(oddEl[k] - numCols < std::numeric_limits<int>::max());
#endif
                if(idxMap[oddEl[k]-numCols] == -1) {
                    idxMap[oddEl[k]-numCols] = realSize;
                    idx[realSize] = oddEl[k] - numCols;
                    coef[realSize] = -1.0;
                    realSize++;
                }
                else {
                    coef[idxMap[oddEl[k]-numCols]] -= 1.0;
                    dup++;
                }
                rhs = rhs - 1.0;
            }
        }

        if(centerSize && fabs(rhs) >= 1e-6) {
            const double oldRhs = rhs;
            for(size_t k = 0; k < centerSize; k++) {
                if(centerIdx[k] < numCols) {
#ifdef DEBUG
                    assert(centerIdx[k] < std::numeric_limits<int>::max());
#endif
                    if(idxMap[centerIdx[k]] == -1) {
                        idxMap[centerIdx[k]] = realSize;
                        idx[realSize] = centerIdx[k];
                        coef[realSize] = oldRhs;
                        realSize++;
                    }
                    else {
                        coef[idxMap[centerIdx[k]]] += oldRhs;
                        dup++;
                    }
                }
                else {
#ifdef DEBUG
                    assert(centerIdx[k] - numCols < std::numeric_limits<int>::max());
#endif
                    if(idxMap[centerIdx[k]-numCols] == -1) {
                        idxMap[centerIdx[k]-numCols] = realSize;
                        idx[realSize] = centerIdx[k] - numCols;
                        coef[realSize] = -1.0 * oldRhs;
                        realSize++;
                    }
                    else {
                        coef[idxMap[centerIdx[k]-numCols]] -= oldRhs;
                        dup++;
                    }
                    rhs = rhs - oldRhs;
                }
            }
        }

        if(dup) {
            int last = 0;
            for(int k = 0; k < realSize; k++) {
                if(fabs(coef[k]) >= 1e-6) {
                    idx[last] = idx[k];
                    coef[last] = coef[k];
                    last++;
                }
            }
            realSize = last;
        }

        cut_pool_insert(cutPool, idx, coef, realSize, rhs, x);
    }

    cut_pool_update(cutPool);

    const size_t numberRowCutsBefore = cs.sizeRowCuts();
    for(size_t i = 0; i < cut_pool_size(cutPool); i++) {
        const Cut *cut = cut_pool_get_cut(cutPool, i);
        osrc.setRow(cut_size(cut) , cut_get_idxs(cut), cut_get_coefs(cut));
        osrc.setUb(cut_get_rhs(cut));
        cs.insertIfNotDuplicate(osrc, equal);
    }

    size_t numberRowCutsAfter = cs.sizeRowCuts();
    CglOddHoleWC::sepCuts += numberRowCutsAfter - numberRowCutsBefore;

    if(!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
        numberRowCutsAfter = cs.sizeRowCuts();
        for(int i = numberRowCutsBefore; i < numberRowCutsAfter; i++) {
            cs.rowCutPtr(i)->setGloballyValid();
        }
    }

    cut_pool_free(&cutPool);
    oddhs_free( &oddhs );
    delete[] x;
    delete[] rc;
    delete[] idx;
    delete[] idxMap;
    delete[] coef;

	CglOddHoleWC::sepTime += (CoinCpuTime() - startSep);
}

CglOddHoleWC::~CglOddHoleWC() {

}
