#include <cstdio>
#include <cassert>
#include <OsiCuts.hpp>
#include <OsiRowCut.hpp>
#include <CoinTime.hpp>

#include "CglBKClique.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinBronKerbosch.hpp"
#include "CoinCliqueList.hpp"
#include "CoinCliqueExtender.hpp"
#include "CoinCutPool.hpp"

#define BKCLQ_MULTIPLIER 1000.0
#define BKCLQ_EPS 1e-6

size_t CglBKClique::sepCuts_ = 0;
double CglBKClique::sepTime_ = 0.0;

static void *xmalloc( const size_t size );

CglBKClique::CglBKClique() : cap_(0), maxCallsBK_(1000),
extMethod_(4), minFrac_(0.001), minViol_(0.02), pivotingStrategy_(3)
{
    minWeight_ = floor(BKCLQ_MULTIPLIER + (minViol_ * BKCLQ_MULTIPLIER));
    vertexWeight_ = NULL;
    rc_ = NULL;
    idxs_ = NULL;
    idxMap_ = NULL;
    coefs_ = NULL;
    inducedVert_ = NULL;
    currClq_ = NULL;
    maxCallsBK_ = 0;
    completeBK_ = false;
}

CglBKClique::CglBKClique(const CglBKClique& rhs) {
    if (this->vertexWeight_) {
        free(this->vertexWeight_);
    }
    if (this->rc_) {
        free(this->rc_);
    }
    if (this->idxs_) {
        free(this->idxs_);
    }
    if (this->idxMap_) {
        free(this->idxMap_);
    }
    if (this->coefs_) {
        free(this->coefs_);
    }
    if (this->inducedVert_) {
        free(this->inducedVert_);
    }
    if (this->currClq_) {
        free(this->currClq_);
    }
    this->cap_ = rhs.cap_;
    this->vertexWeight_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    this->rc_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    this->idxs_ = (int*)xmalloc(sizeof(int) * this->cap_);
    this->idxMap_ = (int*)xmalloc(sizeof(int) * this->cap_);
    this->coefs_ = (double*)xmalloc(sizeof(double) * this->cap_);
    this->inducedVert_ = (size_t*)xmalloc(sizeof(size_t) * this->cap_ * 2);
    this->currClq_ = (size_t*)xmalloc(sizeof(size_t) * this->cap_ * 2);

    this->pivotingStrategy_ = rhs.pivotingStrategy_;
    this->minFrac_ = rhs.minFrac_;
    this->minViol_ = rhs.minViol_;
    this->minWeight_ = rhs.minWeight_;
    this->maxCallsBK_ = rhs.maxCallsBK_;
    this->extMethod_ = rhs.extMethod_;
}

CglBKClique::~CglBKClique() {
    if (this->vertexWeight_) {
        free(this->vertexWeight_);
    }
    if (this->rc_) {
        free(this->rc_);
    }
    if (this->idxs_) {
        free(this->idxs_);
    }
    if (this->idxMap_) {
        free(this->idxMap_);
    }
    if (this->coefs_) {
        free(this->coefs_);
    }
    if (this->inducedVert_) {
        free(this->inducedVert_);
    }
    if (this->currClq_) {
        free(this->currClq_);
    }
}

CglCutGenerator * CglBKClique::clone() const {
    CglBKClique *aClq = new CglBKClique();

    aClq->cap_ = this->cap_;
    aClq->vertexWeight_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    std::copy(this->vertexWeight_, this->vertexWeight_ + (this->cap_ * 2), aClq->vertexWeight_);
    aClq->rc_ = (double*)xmalloc(sizeof(double) * this->cap_ * 2);
    std::copy(this->rc_, this->rc_ + (this->cap_ * 2), aClq->rc_);
    aClq->idxs_ = (int*)xmalloc(sizeof(int) * this->cap_);
    std::copy(this->idxs_, this->idxs_ + this->cap_, aClq->idxs_);
    aClq->idxMap_ = (int*)xmalloc(sizeof(int) * this->cap_);
    std::copy(this->idxMap_, this->idxMap_ + this->cap_, aClq->idxMap_);
    aClq->coefs_ = (double*)xmalloc(sizeof(double) * this->cap_);
    std::copy(this->coefs_, this->coefs_ + this->cap_, aClq->coefs_);
    aClq->inducedVert_ = (size_t*)xmalloc(sizeof(size_t) * this->cap_ * 2);
    std::copy(this->inducedVert_, this->inducedVert_ + (this->cap_ * 2), aClq->inducedVert_);
    aClq->currClq_ = (size_t*)xmalloc(sizeof(size_t) * this->cap_ * 2);
    std::copy(this->currClq_, this->currClq_ + (this->cap_ * 2), aClq->currClq_);

    aClq->pivotingStrategy_ = this->pivotingStrategy_;
    aClq->minFrac_ = this->minFrac_;
    aClq->minViol_ = this->minViol_;
    aClq->minWeight_ = this->minWeight_;
    aClq->maxCallsBK_ = this->maxCallsBK_;
    aClq->extMethod_ = this->extMethod_;
    aClq->maxCallsBK_ = this->maxCallsBK_;
    aClq->completeBK_ = this->completeBK_;

    return static_cast<CglCutGenerator*>(aClq);
}

void CglBKClique::generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info) {
    double startSep = CoinCpuTime();
    const CoinConflictGraph *cgraph = si.getCGraph();

#ifdef DEBUGCG
    if(si.getNumCols() != cgraph->size() / 2) {
        fprintf(stderr, "Invalid conflict graph! Number of columns %d ... in graph %lu\n",
                si.getNumCols(), cgraph->size() / 2);
        abort();
    }
#endif

    checkMemory(si.getNumCols());

    CoinCliqueList *initialCliques = separateCliques(si);

    if (initialCliques->nCliques() > 0) {
        if (!extMethod_) {
            insertCuts(si, info, initialCliques, cs);
        } else {
            CoinCliqueList *extCliques = extendCliques(si, initialCliques);
            insertCuts(si, info, extCliques, cs);
            delete extCliques;
        }
    }

    delete initialCliques;
    CglBKClique::sepTime_ += (CoinCpuTime() - startSep);
}

void CglBKClique::checkMemory(const size_t newNumCols) {
    if (cap_ < newNumCols) {
        if (cap_ > 0) {
#ifdef DEBUGCG
            assert(vertexWeight_);
            assert(rc_);
            assert(idxs_);
            assert(idxMap_);
            assert(coefs_);
            assert(inducedVert_);
            assert(currClq_);
#endif
            free(vertexWeight_);
            free(rc_);
            free(idxs_);
            free(idxMap_);
            free(coefs_);
            free(inducedVert_);
            free(currClq_);
        }

        vertexWeight_ = (double*)xmalloc(sizeof(double) * newNumCols * 2);
        rc_ = (double*)xmalloc(sizeof(double) * newNumCols * 2);
        idxs_ = (int*)xmalloc(sizeof(int) * newNumCols);
        idxMap_ = (int*)xmalloc(sizeof(int) * newNumCols);
        coefs_ = (double*)xmalloc(sizeof(double) * newNumCols);
        inducedVert_ = (size_t*)xmalloc(sizeof(size_t) * newNumCols * 2);
        currClq_ = (size_t*)xmalloc(sizeof(size_t) * newNumCols * 2);
        cap_ = newNumCols;
    }
}

CoinCliqueList* CglBKClique::separateCliques(const OsiSolverInterface &si) {
    const size_t numCols = si.getNumCols();
    const double *x = si.getColSolution();
    CoinCliqueList *initialCliques = new CoinCliqueList(4096, 32768);
    const CoinConflictGraph *cgraph = si.getCGraph();
    size_t n = 0;

    //generating the subgraph induced by the fractional variables and variables at one
    for (size_t i = 0; i < numCols; i++) {//variables
        const size_t degree = cgraph->degree(i);
        if (degree < 2) {
            //disconsidering variables that have no conflicts or just conflicts involving their complements
            continue;
        } else if (x[i] + BKCLQ_EPS <= minFrac_) {
            //variables at zero in x are disconsidered
            continue;
        } else {
            inducedVert_[n] = i;
            vertexWeight_[n] = x[i] * BKCLQ_MULTIPLIER;
            n++;
        }
    }
    for (size_t i = numCols; i < numCols * 2; i++) {//complements
        const size_t degree = cgraph->degree(i);
        const double xc = 1.0 - x[i - numCols];
        if (degree < 2) {
            //disconsidering variables that have no conflicts or just conflicts involving their complements
            continue;
        } else if (xc + BKCLQ_EPS <= minFrac_) {
            //variables at zero in x are disconsidered
            continue;
        } else {
            inducedVert_[n] = i;
            vertexWeight_[n] = xc * BKCLQ_MULTIPLIER;
            n++;
        }
    }

    CoinConflictGraph *ppcg = new CoinStaticConflictGraph(cgraph, n, inducedVert_);
#ifdef DEBUGCG
    assert(ppcg->size() == n);
#endif

    completeBK_ = true;

    if (ppcg->size() >= 2) {
        ppcg->computeModifiedDegree();
        CoinBronKerbosch *bk = new CoinBronKerbosch(ppcg, vertexWeight_, pivotingStrategy_);
        bk->setMaxCalls(maxCallsBK_);
        bk->setMinWeight(minWeight_);
        bk->findCliques();
        completeBK_ = bk->completedSearch();
        maxCallsBK_ = bk->numCalls();

        if (bk->nCliques() > 0) {
            for(size_t i = 0; i < bk->nCliques(); i++) {
                const size_t cliqueSize = bk->getCliqueSize(i);
                const size_t *cliqueIdxs = bk->getClique(i);
                for(size_t j = 0; j < cliqueSize; j++) {
                    currClq_[j] = inducedVert_[cliqueIdxs[j]];
                }
                initialCliques->addClique(cliqueSize, currClq_);
            }
        }

        delete bk;
    }

    delete ppcg;

    return initialCliques;
}

CoinCliqueList* CglBKClique::extendCliques(const OsiSolverInterface &si, const CoinCliqueList *initialCliques) {
    CoinCliqueList *extCliques = new CoinCliqueList(4096, 32768);
    const double *rCost = si.getReducedCost();
    const size_t numCols = si.getNumCols();
    const CoinConflictGraph *cgraph = si.getCGraph();

    //setting reduced costs
    for (size_t i = 0; i < numCols; i++) {
        rc_[i] = rCost[i];
        rc_[i + numCols] = -rc_[i];
    }

    CoinCliqueExtender clqe(cgraph, extMethod_, rc_, 100.0);

    for (size_t i = 0; i < initialCliques->nCliques(); i++) {
        const size_t *clqEl = initialCliques->cliqueElements(i);
        const size_t nClqEl = initialCliques->cliqueSize(i);
        const size_t nNewClqs = clqe.extendClique(clqEl, nClqEl);

        /* adds clique if it is not extended */
        if (!nNewClqs) {
            extCliques->addClique(nClqEl, clqEl);
        }
    }

    /* adding all extended cliques */
    for (size_t i = 0; i < clqe.nCliques(); i++) {
        extCliques->addClique(clqe.getCliqueSize(i), clqe.getClique(i));
    }

    return extCliques;
}

void CglBKClique::insertCuts(const OsiSolverInterface &si, const CglTreeInfo &info, const CoinCliqueList *cliques, OsiCuts &cs) {
    const double *x = si.getColSolution();
    const size_t numCols = si.getNumCols();
    CoinCutPool cutpool(x, numCols);

    for(size_t i = 0; i < cliques->nCliques(); i++) {
        const size_t clqSize = cliques->cliqueSize(i);
        const size_t *el = cliques->cliqueElements(i);
        double rhs = 1.0;
        int cutSize = 0;
        size_t dup = 0;

        std::fill(idxMap_, idxMap_ + numCols, -1);

        for (size_t j = 0; j < clqSize; j++) {
            if (el[j] < numCols) {
                if(idxMap_[el[j]] == -1) {
                    idxMap_[el[j]] = cutSize;
                    idxs_[cutSize] = ((int)el[j]);
                    coefs_[cutSize] = 1.0;
                    cutSize++;
                } else {
                    coefs_[idxMap_[el[j]]] += 1.0;
                    dup++;
                }
            } else {
                rhs -= 1.0;
                if(idxMap_[el[j]-numCols] == -1) {
                    idxMap_[el[j]-numCols] = cutSize;
                    idxs_[cutSize] = ((int)(el[j] - numCols));
                    coefs_[cutSize] = -1.0;
                    cutSize++;
                } else {
                    coefs_[idxMap_[el[j]-numCols]] -= 1.0;
                    dup++;
                }
            }
        }

#ifdef DEBUGCG
        assert(dup == 0 || dup == 1);
#endif

        if(dup) {
            int last = 0;

            for(int k = 0; k < cutSize; k++) {
                if(fabs(coefs_[k]) >= BKCLQ_EPS) {
                    idxs_[last] = idxs_[k];
                    coefs_[last] = coefs_[k];
                    last++;
                }
            }
            cutSize = last;
        }

        cutpool.add(idxs_, coefs_, cutSize, rhs);
    }

    cutpool.removeNullCuts();

    const size_t numberRowCutsBefore = cs.sizeRowCuts();
    for(size_t i = 0; i < cutpool.numCuts(); i++) {
        osrc_.setRow(cutpool.cutSize(i) , cutpool.cutIdxs(i), cutpool.cutCoefs(i));
        osrc_.setUb(cutpool.cutRHS(i));
        cs.insertIfNotDuplicate(osrc_);
    }

    size_t numberRowCutsAfter = cs.sizeRowCuts();
    CglBKClique::sepCuts_ += numberRowCutsAfter - numberRowCutsBefore;

    if(!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
        numberRowCutsAfter = cs.sizeRowCuts();
        for(size_t i = numberRowCutsBefore; i < numberRowCutsAfter; i++) {
            cs.rowCutPtr(i)->setGloballyValid();
        }
    }
}

void CglBKClique::setMaxCallsBK(size_t maxCallsBK) {
    this->maxCallsBK_ = maxCallsBK;
}

void CglBKClique::setExtendingMethod(size_t extMethod) {
    if(extMethod > 6) {
        fprintf(stderr, "Invalid value for parameter extMethod (%ld).\n", (size_t)extMethod);
        abort();
    }

    this->extMethod_ = extMethod;
}

void CglBKClique::setMinFrac(const double minFrac) {
    minFrac_ = minFrac;
}

void CglBKClique::setMinViol(const double minViol){
    minViol_ = minViol;
    minWeight_ = BKCLQ_MULTIPLIER + (minViol_ * BKCLQ_MULTIPLIER);
}

void CglBKClique::setPivotingStrategy(const size_t pivotingStrategy) {
    pivotingStrategy_ = pivotingStrategy;
}

static void *xmalloc( const size_t size ) {
    void *result = malloc( size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
        abort();
    }

    return result;
}
