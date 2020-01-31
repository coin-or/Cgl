#include <cstdio>
#include <cassert>
#include <limits>

#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "CoinTime.hpp"
#include "CglBKClique.hpp"
#include "bron_kerbosch.h"
#include "clique_separation.h"
#include "cut.h"

size_t CglBKClique::sepCuts = 0;
double CglBKClique::sepTime = 0.0;

CglBKClique::CglBKClique() : maxItBK(1000), maxItBKExt(100), extMethod(3) {
}

CglBKClique::CglBKClique(const CglBKClique& rhs) {
    this->maxItBK = rhs.maxItBK;
    this->maxItBKExt = rhs.maxItBKExt;
    this->extMethod = rhs.extMethod;
}

CglCutGenerator * CglBKClique::clone() const {
    CglBKClique *aClq = new CglBKClique();

    aClq->maxItBK = this->maxItBK;
    aClq->maxItBKExt = this->maxItBKExt;
    aClq->extMethod = this->extMethod;

    return static_cast<CglCutGenerator*>(aClq);
}

void CglBKClique::generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info) {
    double startSep = CoinCpuTime();

    const CoinAbsFltEq equal(1.0e-12);
    OsiRowCut osrc;

    const size_t numCols = si.getNumCols();

    const CGraph *cg = si.getCGraph();
    if(numCols != cgraph_size(cg)/2) {
        fprintf(stderr, "Invalid conflict graph! Number of columns %ld ... in graph %ld\n",
                numCols, cgraph_size(cg)/2);
        exit(EXIT_FAILURE);
    }

    double *x = new double[numCols*2];
    double *rc = new double[numCols*2];
    int *idxs = new int[numCols];
    double *coefs = new double[numCols];
    int *idxMap = new int[numCols];
    CliqueSeparation *sep = clq_sep_create(cg);

    const double *colSol = si.getColSolution();
    const double *rCost = si.getReducedCost();
    for(size_t i = 0; i < numCols; i++) {
        x[i] = colSol[i];
        rc[i] = rCost[i];
        x[i + numCols] = 1.0 - x[i];
        rc[i + numCols] = -rc[i];
    }

    clq_sep_set_rc(sep, rc);
    CutPool *cutPool = cut_pool_create(numCols);
    clq_sep_set_extend_method(sep, this->extMethod);
    clq_sep_set_max_it_bk(sep, this->maxItBK);
    clq_sep_set_max_it_bk_ext(sep, this->maxItBKExt);

    clq_sep_separate(sep, x);

    /* inserting cliques */
    const CliqueSet *clqSet = clq_sep_get_cliques(sep);
    for(size_t i = 0; i < clq_set_number_of_cliques(clqSet); i++) {
        const size_t clqSize = clq_set_clique_size(clqSet, i);
        const size_t *el = clq_set_clique_elements(clqSet, i);
        double rhs = 1.0;
        int cutSize = 0;
        size_t dup = 0;

        std::fill(idxMap, idxMap + numCols, -1);

        for (size_t j = 0; j < clqSize; j++) {
            if (el[j] < numCols) {
#ifdef DEBUG
                assert(el[j] < std::numeric_limits<int>::max());
#endif
                if(idxMap[el[j]] == -1) {
                    idxMap[el[j]] = cutSize;
                    idxs[cutSize] = el[j];
                    coefs[cutSize] = 1.0;
                    cutSize++;
                } else {
                    coefs[idxMap[el[j]]] += 1.0;
                    dup++;
                }
            } else {
#ifdef DEBUG
                assert(el[j] - numCols < std::numeric_limits<int>::max());
#endif
                rhs -= 1.0;
                if(idxMap[el[j]-numCols] == -1) {
                    idxMap[el[j]-numCols] = cutSize;
                    idxs[cutSize] = el[j] - numCols;
                    coefs[cutSize] = -1.0;
                    cutSize++;
                } else {
                    coefs[idxMap[el[j]-numCols]] -= 1.0;
                    dup++;
                }
            }
        }
#ifdef DEBUG
    assert(dup == 0 || dup == 1);
#endif
        if(dup) {
            int last = 0;

            for(int k = 0; k < cutSize; k++) {
                if(fabs(coefs[k]) >= 1e-6) {
                    idxs[last] = idxs[k];
                    coefs[last] = coefs[k];
                    last++;
                }
            }
            cutSize = last;
        }

        cut_pool_insert(cutPool, idxs, coefs, cutSize, rhs, x);
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
    CglBKClique::sepCuts += numberRowCutsAfter - numberRowCutsBefore;

    if(!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
        numberRowCutsAfter = cs.sizeRowCuts();
        for(size_t i = numberRowCutsBefore; i < numberRowCutsAfter; i++) {
            cs.rowCutPtr(i)->setGloballyValid();
        }
    }

    clq_sep_free(&sep);
    cut_pool_free(&cutPool);
    delete[] x;
    delete[] rc;
    delete[] idxs;
    delete[] coefs;
    delete[] idxMap;

    CglBKClique::sepTime += (CoinCpuTime() - startSep);
}

CglBKClique::~CglBKClique() {

}

void CglBKClique::setMaxItBK(size_t _maxItBK) {
    if(_maxItBK <= 0) {
        fprintf(stderr, "Invalid value for parameter maxItBK (%ld).\n", _maxItBK);
        exit(1);   
    }

    this->maxItBK = _maxItBK;
}

void CglBKClique::setMaxItBKExt(size_t _maxItBKExt) {
    if(_maxItBKExt <= 0) {
        fprintf(stderr, "Invalid value for parameter maxItBKExt (%ld).\n", _maxItBKExt);
        exit(1);   
    }
    
    this->maxItBKExt = _maxItBKExt;
}

void CglBKClique::setExtendingMethod(size_t _extMethod) {
    if(_extMethod < 0 || _extMethod > 4) {
        fprintf(stderr, "Invalid value for parameter extMethod (%ld).\n", _extMethod);
        exit(1);   
    }
    
    this->extMethod = _extMethod;
}
