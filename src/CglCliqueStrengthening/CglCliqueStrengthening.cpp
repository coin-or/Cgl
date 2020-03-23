#include <cfloat>
#include <cassert>
#include <algorithm>
#include "CglCliqueStrengthening.hpp"
#include "CoinStaticConflictGraph.hpp"
#include "CoinCliqueExtender.hpp"
#include "CoinCliqueSet.hpp"
#include "OsiSolverInterface.hpp"
#include "CglMessage.hpp"

#define CLQ_STR_EPS 1e-6
#define MAX_SIZE_CLIQUE_TO_BE_EXTENDED 256

static void *xmalloc( const size_t size );
static void *xcalloc( const size_t elements, const size_t size );

enum CliqueStatus {
    NotDominated = 1,
    Dominated = 2
};

class CliqueRows {
public:
    explicit CliqueRows(size_t linesToReserve, size_t nzsToReserve) {
        nRows_ = 0;
        starts_ = (size_t*)xmalloc(sizeof(size_t) * (linesToReserve + 1));
        starts_[0] = 0;
        rowIdx_ = (size_t*)xmalloc(sizeof(size_t) * linesToReserve);
        rowStatus_ = (CliqueStatus*)xmalloc(sizeof(CliqueStatus) * linesToReserve);
        elements_ = (size_t*)xmalloc(sizeof(size_t) * nzsToReserve);
    }

    ~CliqueRows() {
        free(starts_);
        free(elements_);
        free(rowIdx_);
        free(rowStatus_);
    }

    void addRow(size_t nz, const size_t els[], size_t rowIdx, CliqueStatus status) {
        std::copy(els, els + nz, elements_ + starts_[nRows_]);
        rowIdx_[nRows_] = rowIdx;
        rowStatus_[nRows_] = status;
        nRows_++;
        starts_[nRows_] = starts_[nRows_ - 1] + nz;
    }

    const size_t* row(size_t idxRow) const {
#ifdef DEBUGCG
        assert(idxRow < nRows_);
#endif
        return elements_ + starts_[idxRow];
    }

    size_t origIdxRow(size_t idxRow) const {
#ifdef DEBUGCG
        assert(idxRow < nRows_);
#endif
        return rowIdx_[idxRow];
    }

    size_t nz(size_t idxRow) const {
#ifdef DEBUGCG
        assert(idxRow < nRows_);
#endif
        return starts_[idxRow + 1] - starts_[idxRow];
    }

    CliqueStatus status(size_t idxRow) const {
#ifdef DEBUGCG
        assert(idxRow < nRows_);
#endif
        return rowStatus_[idxRow];
    }

    void setStatus(size_t idxRow, CliqueStatus status) const {
#ifdef DEBUGCG
        assert(idxRow < nRows_);
#endif
        rowStatus_[idxRow] = status;
    }

    size_t rows() const {
        return nRows_;
    }

private:
    size_t *starts_, *elements_, *rowIdx_;
    size_t nRows_;
    CliqueStatus *rowStatus_;
};

CglCliqueStrengthening::CglCliqueStrengthening() : nExtended_(0), nDominated_(0), handler_(NULL),
defaultHandler_(true) {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(2);
    messages_ = CglMessage();
}

CglCliqueStrengthening::CglCliqueStrengthening(const CglCliqueStrengthening &rhs) :
nExtended_(rhs.nExtended_), nDominated_(rhs.nDominated_),
defaultHandler_(rhs.defaultHandler_) {

    if (defaultHandler_) {
        handler_ = new CoinMessageHandler();
        handler_->setLogLevel(rhs.handler_->logLevel());
    } else {
        handler_ = rhs.handler_;
    }
    messages_ = rhs.messages_;
}

CglCliqueStrengthening &CglCliqueStrengthening::operator=(const CglCliqueStrengthening &rhs) {
    if (this != &rhs) {
        gutsOfDestructor();
        nExtended_ = rhs.nExtended_;
        nDominated_ = rhs.nDominated_;
        defaultHandler_ = rhs.defaultHandler_;
        if (defaultHandler_) {
            handler_ = new CoinMessageHandler();
            handler_->setLogLevel(rhs.handler_->logLevel());
        } else {
            handler_ = rhs.handler_;
        }
        messages_ = rhs.messages_;
    }

    return *this;
}

CglCliqueStrengthening::~CglCliqueStrengthening() {
    gutsOfDestructor();
}

void CglCliqueStrengthening::gutsOfDestructor() {
    if (defaultHandler_) {
        delete handler_;
        handler_ = NULL;
    }
}

// Pass in Message handler (not deleted at end)
void CglCliqueStrengthening::passInMessageHandler(CoinMessageHandler *handler) {
    if (defaultHandler_)
        delete handler_;
    defaultHandler_ = false;
    handler_ = handler;
}

// Set language
void CglCliqueStrengthening::newLanguage(CoinMessages::Language language) {
    messages_ = CglMessage(language);
}

void detectCliques(const OsiSolverInterface &mip, CliqueRows *cliques) {
    const int numRows = mip.getNumRows();
    const int numCols = mip.getNumCols();
    const CoinPackedMatrix *cpmRow = mip.getMatrixByRow();
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const double *Arhs = mip.getRightHandSide();
    const char *Asense = mip.getRowSense();
    const double *colLB = mip.getColLower();
    const double *colUB = mip.getColUpper();
    const char *colType = mip.getColType();
    const int *lengths = cpmRow->getVectorLengths();
    size_t *tmpRow = (size_t*)xmalloc(sizeof(size_t) * numCols);

    for (size_t i = 0; i < numRows; i++) {
        const size_t nz = lengths[i];
        const char sense = Asense[i];

        if (nz <= 1 || nz > MAX_SIZE_CLIQUE_TO_BE_EXTENDED) {
            continue;
        }

        if (sense != 'L' && sense != 'G') {
            continue;
        }

        const int *idxs = cpmRow->getIndices() + starts[i];
        const double *coefs = cpmRow->getElements() + starts[i];
        const double mult = (sense == 'G') ? -1.0 : 1.0;
        double rhs = mult * Arhs[i];

        bool testRow = true;
        double minCoef1 = std::numeric_limits< double >::max();
        double minCoef2 = std::numeric_limits< double >::max();
        for (size_t j = 0; j < nz; j++) {
            tmpRow[j] = idxs[j];

            double coefCol = coefs[j] * mult;
            const bool isBinary = (colType[tmpRow[j]] != 0) && (colLB[tmpRow[j]] == 1.0 || colLB[tmpRow[j]] == 0.0)
                                  && (colUB[tmpRow[j]] == 0.0 || colUB[tmpRow[j]] == 1.0);

            if (!isBinary) {
                testRow = false;
                break;
            }

            if (coefCol <=- CLQ_STR_EPS) {
                tmpRow[j] += numCols;
                coefCol = -coefCol;
                rhs = rhs + coefCol;
            }

            if (coefCol + CLQ_STR_EPS <= minCoef1) {
                minCoef2 = minCoef1;
                minCoef1 = coefCol;
            } else if (coefCol + CLQ_STR_EPS <= minCoef2) {
                minCoef2 = coefCol;
            }
        }

        if (!testRow) {
            continue;
        }

        if (minCoef1 + minCoef2 >= rhs + CLQ_STR_EPS) {
            cliques->addRow(nz, tmpRow, i, NotDominated);
        }
    }

    free(tmpRow);
}

void checkDominance(const size_t *extClqEl, size_t extClqSize, CliqueRows *cliques, const size_t **colClqs,
        const size_t *nColClqs, bool *ivRow, bool *ivCol, size_t numCols) {
#ifdef DEBUGCG
    for (size_t i = 0; i < cliques->rows(); i++) {
        assert(!ivRow[i]);
    }
    for (size_t i = 0; i < numCols * 2; i++) {
        assert(!ivCol[i]);
    }
#endif

    for (size_t i = 0; i < extClqSize; i++) {
        ivCol[extClqEl[i]] = true;
    }

    for (size_t i = 0; i < extClqSize; i++) {
        size_t col = extClqEl[i];

        // checking cliques where column col appear
        for (size_t j = 0; j < nColClqs[col]; j++) {
            const size_t clqRowIdx = colClqs[col][j];

            // skipping already dominated or already tested rows
            if (cliques->status(clqRowIdx) == Dominated || ivRow[clqRowIdx]) {
                continue;
            }

            ivRow[clqRowIdx] = true;

            const size_t *clqEl = cliques->row(clqRowIdx);
            const size_t clqNZ = cliques->nz(clqRowIdx);
            bool dominates = true;

            for (size_t k = 0; k < clqNZ; k++) {
                if (!ivCol[clqEl[k]]) {
                    dominates = false;
                    break;
                }
            }

            if (dominates) {
                cliques->setStatus(clqRowIdx, Dominated);
            }
        }
    }

    //clearing ivCol
    for (size_t i = 0; i < extClqSize; i++) {
        ivCol[extClqEl[i]] = false;
    }

    //clearing ivRow
    for (size_t i = 0; i < cliques->rows(); i++) {
        ivRow[i] = false;
    }
}

void CglCliqueStrengthening::strengthenCliques(OsiSolverInterface &model, size_t extMethod) {
    const int numCols = model.getNumCols();
    const int numRows = model.getNumRows();

    if (numCols == 0 || numRows == 0) {
    	return;
    }

    model.checkCGraph();

    const CoinConflictGraph *cgraph = model.getCGraph();
    
    CliqueRows cliques(numRows, model.getNumElements());

    nExtended_ = nDominated_ = 0;
    detectCliques(model, &cliques);

    if (cliques.rows() == 0) { //no cliques
        return;
    }

    //filling cliques per column
    size_t *nColClqs = (size_t *) xcalloc(numCols * 2, sizeof(size_t));
    size_t nElements = 0;
    for (size_t i = 0; i < cliques.rows(); i++) {
        const size_t *clqEl = cliques.row(i);
        const size_t clqSize = cliques.nz(i);
#ifdef DEBUGCG
        assert(clqSize >= 2);
#endif
        for (size_t j = 0; j < clqSize; j++) {
            nColClqs[clqEl[j]]++;
            nElements++;
        }
    }
    size_t **colClqs = (size_t **) xmalloc(sizeof(size_t *) * numCols * 2);
    colClqs[0] = (size_t *) xmalloc(sizeof(size_t) * nElements);
    for (size_t i = 1; i < numCols * 2; i++) {
        colClqs[i] = colClqs[i - 1] + nColClqs[i - 1];
        nColClqs[i - 1] = 0;
    }
    nColClqs[(2 * numCols) - 1] = 0;
    for (size_t i = 0; i < cliques.rows(); i++) {
        const size_t *clqEl = cliques.row(i);
        for (size_t j = 0; j < cliques.nz(i); j++) {
            const size_t col = clqEl[j];
            colClqs[col][nColClqs[col]++] = i;
        }
    }

    //filling reduced cost
    const double *redCost = model.getReducedCost();
    double *rc = (double*)xmalloc(sizeof(double) * numCols * 2);
    for (size_t i = 0; i < numCols; i++) {
        rc[i] = redCost[i];
        rc[i + numCols] = -rc[i];
    }

    const double *Arhs = model.getRightHandSide();
    const char *Asense = model.getRowSense();
    CoinCliqueSet newCliques(4096, 32768);
    OsiSolverInterface::OsiNameVec clqNames;
    bool *ivRow = (bool*)xcalloc(numRows, sizeof(bool));
    bool *ivCol = (bool*)xcalloc(numCols * 2, sizeof(bool));
    size_t last = 0;
    char name[256];

    CoinCliqueExtender clqe(cgraph, extMethod, rc);
    clqe.setMaxCandidates(512);

    for (size_t i = 0; i < cliques.rows(); i++) {
        const size_t clqOrigRowIdx = cliques.origIdxRow(i);
        const size_t *clqIdx = cliques.row(i);
        const size_t clqSize = cliques.nz(i);

#ifdef DEBUGCG
        CoinCliqueList::validateClique(cgraph, clqIdx, clqSize);
        assert(clqe.nCliques() >= last);
#endif

        if (cliques.status(i) == Dominated) {
            continue;
        }

        clqe.extendClique(clqIdx, clqSize);
        const size_t numClqs = clqe.nCliques() - last;

        if (numClqs > 0) {
            cliques.setStatus(i, Dominated);
            nExtended_++;

#ifdef DEBUGCG
            assert(numClqs == 1);
#endif
            const bool inserted = newCliques.insertIfNotDuplicate(clqe.getCliqueSize(last), clqe.getClique(last));
            if (inserted) {
                checkDominance(clqe.getClique(last), clqe.getCliqueSize(last), &cliques,
                        (const size_t**)colClqs, (const size_t*)nColClqs, ivRow, ivCol, numCols);
                sprintf(name, "%s_ext", model.getRowName(clqOrigRowIdx).c_str());
                clqNames.push_back(name);
            }
            last = clqe.nCliques();
        }
    }

#ifdef DEBUGCG
    assert(clqNames.size() == nCliques);
#endif

    free(rc);
    free(colClqs[0]);
    free(colClqs);
    free(nColClqs);
    free(ivRow);
    free(ivCol);

    if (newCliques.nCliques() == 0) {
        return;
    }

    /* removing dominated cliques */
    int *toRemove = (int*)xmalloc(sizeof(int) * numRows);
    int nToRemove = 0;
    for (size_t i = 0; i < cliques.rows(); i++) {
        if (cliques.status(i) == Dominated) {
            const size_t origRowIDx = cliques.origIdxRow(i);
            nDominated_++;
            toRemove[nToRemove++] = origRowIDx;
        }
    }
    if (nToRemove > 0) {
        model.deleteRows(nToRemove, toRemove);
    }
    free(toRemove);

    /* adding stronger cliques */
    const size_t nCliques = newCliques.nCliques();
    int *nrIdx = (int*)xmalloc(sizeof(int) * newCliques.totalElements());
    int *idxMap = (int*)xmalloc(sizeof(int) * numCols);//controls duplicated indexes (var and complement)
    double *nrCoef = (double*)xmalloc(sizeof(double) * newCliques.totalElements());
    int *nrStart = (int*)xmalloc(sizeof(int) * (nCliques + 1)); nrStart[0] = 0;
    double *nrLB = (double*)xmalloc(sizeof(double) * nCliques);
    double *nrUB = (double*)xmalloc(sizeof(double) * nCliques);
    size_t numVars = 0;
    for (size_t ic = 0; ic < nCliques; ic++) {
        const size_t extClqSize = newCliques.cliqueSize(ic);
        const size_t *extClqEl = newCliques.cliqueElements(ic);
        double rhs = 1.0;
        size_t duplicated = 0;

        std::fill(idxMap, idxMap + numCols, -1);

        for (size_t i = 0; i < extClqSize; i++) {
            if (extClqEl[i] < numCols) {
                if(idxMap[extClqEl[i]] == -1) {
                    idxMap[extClqEl[i]] = numVars;
                    nrIdx[numVars] = (int)extClqEl[i];
                	nrCoef[numVars] = 1.0;
                    numVars++;
                } else {
                    nrCoef[idxMap[extClqEl[i]]] += 1.0;
                    assert(nrCoef[idxMap[extClqEl[i]]] == 0.0);
                    duplicated++;
                }
            } else {
            	rhs -= 1.0;
            	if(idxMap[extClqEl[i]-numCols] == -1) {
            		idxMap[extClqEl[i]-numCols] = numVars;
            		nrIdx[numVars] = ((int)extClqEl[i] - numCols);
                	nrCoef[numVars] = -1.0;
                	numVars++;
                } else {
                	nrCoef[idxMap[extClqEl[i]-numCols]] -= 1.0;
                	assert(nrCoef[idxMap[extClqEl[i]-numCols]] == 0.0);
                    duplicated++;
                }
            }
        }

        assert(duplicated == 0 || duplicated == 1);
        if(duplicated == 1) {
            int last = nrStart[ic];
            rhs = 0.0;
            for(int k = nrStart[ic]; k < numVars; k++) {
    			assert(nrCoef[k] == -1.0 || nrCoef[k] == 0.0 || nrCoef[k] == 1.0);
                if(nrCoef[k] == -1.0 || nrCoef[k] == 1.0) {
                    nrIdx[last] = nrIdx[k];
                    nrCoef[last] = nrCoef[k];
                    last++;
                    if (nrCoef[k] == -1.0) {
                		rhs -= 1.0;
                	}
                }
            }
            numVars = last;

            sprintf(name, "%s_dup", clqNames[ic].c_str());
            clqNames[ic] = name;
        }

        nrLB[ic] = -DBL_MAX;
        nrUB[ic] = rhs;
        nrStart[ic + 1] = numVars;
    }

    const int lastOrigIdx = model.getNumRows();
    model.addRows(nCliques, nrStart, nrIdx, nrCoef, nrLB, nrUB);
    for (int i = 0; i < (int)nCliques; i++) {
        model.setRowName(lastOrigIdx + i, clqNames[i]);
    }

    free(nrIdx);
    free(nrCoef);
    free(nrLB);
    free(nrUB);
    free(nrStart);
    free(idxMap);

    handler_->message(CGL_PROCESS_CLQSTR, messages_) << nExtended_ << nDominated_ << CoinMessageEol;
}

static void *xmalloc( const size_t size ) {
    void *result = malloc( size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to allocate %zu bytes.", size);
        abort();
    }

    return result;
}

static void *xcalloc( const size_t elements, const size_t size ) {
    void *result = calloc( elements, size );
    if (!result) {
        fprintf(stderr, "No more memory available. Trying to callocate %zu bytes.", size * elements);
        abort();
    }

    return result;
}