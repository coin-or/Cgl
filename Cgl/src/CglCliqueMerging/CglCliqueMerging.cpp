#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <algorithm>
#include "CglMessage.hpp"
#include "OsiSolverInterface.hpp"
#include "CglCliqueMerging.hpp"
#include "cgraph.h"
#include "clique.h"
#include "clique_extender.h"
#include "vint_set.h"

enum CliqueType {
    NotAClique = 0,
    NotDominated = 1,
    Dominated = 2,
    MoreThanAClique = 3 // do not try to dominate this constraint
};

// simple sparse matrix
class SpsMtx {
public:
    explicit SpsMtx(size_t lines, size_t nzsToReserve = std::numeric_limits<int>::max()) :
    starts(lines + 1, std::numeric_limits<int>::max()), nzs(lines, 0) {
        if (nzsToReserve != (size_t)std::numeric_limits<int>::max()) {
            elements.reserve(nzsToReserve);
        }
    }

    void addRow(size_t line, size_t nz, const int els[], bool _sort = false) {
#ifdef DEBUG
        assert(line >= 0 && line < nzs.size());
#endif
        starts[line] = elements.size();
        nzs[line] = nz;
        elements.insert(elements.end(), els, els + nz);

        if (_sort) {
            std::sort(elements.begin() + starts[line], elements.begin() + starts[line] + nz);
        }
    }

    const size_t *row(size_t line) const {
#ifdef DEBUG
        assert(line >= 0 && line < nzs.size());
#endif

        if (nzs[line] == 0) {
            return NULL;
        }

        return &(elements[starts[line]]);
    }

    size_t nz(size_t line) const {
#ifdef DEBUG
        assert(line >= 0 && line < nzs.size());
#endif
        return nzs[line];
    }

private:
    std::vector<size_t> starts;
    std::vector<size_t> nzs;
    std::vector<size_t> elements;
};

#define MAX_SIZE_CLIQUE_TO_BE_EXTENDED 256

// global configurations and execution stats
size_t clqMergeVerbose = 0;

/*
double clqMergeSecsCheckClique = 0.0;
double clqMergeSecsExtendAndDominate = 0.0;
double clqMergeSecsAddAndRemove = 0.0;
double clqMergeSecsExtend = 0.0;
size_t clqMergeNExtended = 0;
size_t clqMergeNDominatedFull = 0;
size_t clqMergeNDominatedEqPart = 0;
*/

static bool dominates(size_t nClq1, const size_t clq1[], size_t nClq2, const size_t clq2[], bool *iv) {
    if (nClq1 < nClq2) {
        return false;
    }

    // filling incidence vector
    for (size_t i = 0; i < nClq1; i++) {
#ifdef DEBUG
        assert(!iv[clq1[i]]);
#endif

        iv[clq1[i]] = true;
    }

    bool res = true;
    // checking if clq2 is contained in clq1
    for (size_t i = 0; i < nClq2; i++) {
        if (!iv[clq2[i]]) {
            res = false;
            break;
        }
    }

    // clearing iv
    for (size_t i = 0; i < nClq1; i++) {
        iv[clq1[i]] = false;
    }

    return res;
}

static void add_clique(
        OsiSolverInterface &mip,
        CliqueSet *newCliques, // new cliques
        size_t size, const size_t el[], // clique being added
        std::vector<enum CliqueType> &cliqueState, // which are the clique rows and their state
        const SpsMtx &origClqs, // original cliques
        bool *iv, // temporary incidence vector
        bool *ivrt, // temporary incidence vector
        size_t *rtc, // list of rows checked
        const size_t **colClqs, // column cliques
        const size_t nColClqs[] // number of cliques of a column
) {
#ifdef DEBUG
    for (size_t i = 0; i < mip.getNumCols(); i++) {
        assert(!iv[i]);
    }
    for (size_t i = 0; i < mip.getNumRows(); i++) {
        assert(!ivrt[i]);
    }
#endif

    size_t add = 0;
    size_t nc1 = clq_set_number_of_cliques(newCliques);
    clq_set_add(newCliques, el, size, size);
    size_t nc2 = clq_set_number_of_cliques(newCliques);
    add = nc2 - nc1;

#ifdef DEBUG
    assert(nc2 >= nc1);
#endif
    if (add == 0) { // ignoring repeated cliques
        return;
    }

    size_t minColClq = std::numeric_limits<int>::max(), maxColClq = 0;

    // checking which rows should be checked considering
    // columns that appear on this clique
    size_t nrtc = 0;
    for (size_t i = 0; i < size; i++) {
        size_t col = el[i];

        minColClq = std::min(minColClq, col);
        maxColClq = std::max(maxColClq, col);

        if (col >= mip.getNumCols()) // complimentary variable
            col -= mip.getNumCols();
#ifdef DEBUG
        assert(col >= 0 && col < mip.getNumCols());
#endif
        // checking non-dominated cliques that column appear
        for (size_t j = 0; j < nColClqs[col]; j++) {
            size_t ir = colClqs[col][j];
#ifdef DEBUG
            assert(cliqueState[ir] != NotAClique);
#endif
            // skipping already dominated, already included or larger rows (cannot be dominated by this one )
            if (cliqueState[ir] == Dominated || ivrt[ir] || origClqs.nz(ir) > size) {
                continue;
            }

            ivrt[ir] = true;
            rtc[nrtc++] = ir;
        }
    }

    for (size_t i = 0; i < nrtc; i++) {
        size_t rowClique = rtc[i];
        ivrt[rowClique] = false;

        // quick check before going to the slow check
        size_t nzRow = origClqs.nz(rowClique);
        const size_t *clqRow = origClqs.row(rowClique);
        if (size == nzRow) {
            if ((clqRow[0] != minColClq) || (clqRow[nzRow - 1] != maxColClq)) {
                continue;
            }
        } else {
#ifdef DEBUG
            assert(nzRow < size);
#endif
            if ((clqRow[0] < minColClq) || (clqRow[nzRow - 1] > maxColClq)) {
                continue;
            }
        }

        if (dominates(size, el, nzRow, clqRow, iv)) {
            cliqueState[rowClique] = Dominated;
            if (clqMergeVerbose >= 2) {
                printf("\t\tdominates %s\n", mip.getRowName(rowClique).c_str());
            }
        } // dominates
    } // all rows to check
}

/* detect cliques, including those involving
 * complimentary variables and stores those cliques
 * in origClqs */
static size_t detect_cliques(
        OsiSolverInterface &mip,
        std::vector<enum CliqueType> &cliqueState,
        std::vector<size_t> &cliques,
        SpsMtx &origClqs) {
    const CoinPackedMatrix *cpmRow = mip.getMatrixByRow();
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const double *Arhs = mip.getRightHandSide();
    const char *Asense = mip.getRowSense();
    const double *colLB = mip.getColLower();
    const double *colUB = mip.getColUpper();
    cliques.clear();

    // used to fix column indexes
    int *cidx = new int[mip.getNumCols()];

    if (clqMergeVerbose >= 1) {
        printf("checking candidates for clique extension.\n");
    }

    size_t rows = mip.getNumRows();

    size_t minClqSize = std::numeric_limits<int>::max(), maxClqSize = 0;
    double avClgSize = 0.0;

    for (size_t i = 0; i < rows; i++) {
        size_t nz = cpmRow->getVectorLengths()[i];
        const char sense = Asense[i];

        if (nz <= 1 || nz > MAX_SIZE_CLIQUE_TO_BE_EXTENDED || sense == 'R') {
            continue;
        }

        const int *idx = cpmRow->getIndices() + starts[i];
        const double *coef = cpmRow->getElements() + starts[i];

        double rhs = Arhs[i];
        // all constraints as <=
        double mult = (sense == 'G') ? -1.0 : 1.0;

        /* all variables should be positive, integers */
        char varsOk = true;
        double minCoef = DBL_MAX;

        // to check if original and complimentary variables are involved
        size_t nOnes = 0, nMinusOne = 0;

        for (size_t j = 0; j < nz; j++) {
            if ((!mip.isInteger(idx[j])) || colLB[idx[j]] <= -1e-5 ||
                (coef[j] <= -1e8 && colUB[idx[j]] >= 1.0 + 1e-8)) {
                varsOk = false;
                break;
            }

            const double realCoef = mult * coef[j];

            // binaries, checking for original and complimentary variables
            if (colLB[idx[j]] >= -1e-8 && colUB[idx[j]] <= 1.0 + 1e-8) {
                if (fabs(realCoef - 1.0) <= 1e-8) {
                    ++nOnes;
                } else if (fabs(realCoef + 1.0) <= 1e-8) {
                    ++nMinusOne;
                }
            }

            minCoef = std::min(realCoef, minCoef);
        }

        if (!varsOk) {
            continue;
        }

        // first detect clique only with normal variables
        cliqueState[i] = (2 * minCoef >= mult * rhs + 1e-5 && rhs >= 1e-5) ? NotDominated : NotAClique;

        if (cliqueState[i] == NotDominated) {
            cliques.push_back(i);

            origClqs.addRow(i, nz, idx, true);

            minClqSize = std::min(minClqSize, nz);
            maxClqSize = std::max(maxClqSize, nz);
            avClgSize += nz;
        } // found a clique candidate
        else {
            // checking for clique involving normal variables and complementary variables
            if ((nOnes + nMinusOne == nz) && (fabs(rhs - (1.0 - nMinusOne)) <= 1e-8)) {
                std::copy(idx, idx + nz, cidx);
                cliqueState[i] = NotDominated;

                cliques.push_back(i);

                for (size_t j = 0; j < nz; j++) {
                    if (coef[j] < -1e-8) {
                        cidx[j] += mip.getNumCols();
                    }
                }

                origClqs.addRow(i, nz, cidx, true);

                minClqSize = std::min(minClqSize, nz);
                maxClqSize = std::max(maxClqSize, nz);
                avClgSize += nz;
            }
        }

        if (clqMergeVerbose >= 2 && cliqueState[i] == NotDominated) {
            printf("\trow %s: ", mip.getRowName(i).c_str());
            for (size_t j = 0; j < nz; j++) {
                printf(" %+g %s ", coef[j], mip.getColName(idx[j]).c_str());
            }
            char strSense[3] = "";
            switch (sense) {
                case 'E': {
                    strcpy(strSense, "=");
                    break;
                }
                case 'L': {
                    strcpy(strSense, "<=");
                    break;
                }
                case 'G': {
                    strcpy(strSense, ">=");
                    break;
                }
            } // sense

            printf("%s %g\n", strSense, rhs);
        } // if verbose

    } // all rows

    //clqMergeSecsCheckClique = ((double)clock()-startcq) / ((double)CLOCKS_PER_SEC);

    //if (clqMergeVerbose>=1)
    //  printf("model checked in %.4f seconds. %d candidate cliques for extension/merging. clique sizes range:[%d...%d], av %.2f.\n", clqMergeSecsCheckClique, nCliques, minClqSize, maxClqSize, avClgSize/((double)nCliques) );

    delete[] cidx;

    return cliques.size();
}

/* to sort cliques per size */
class CliqueSize {
public:
    CliqueSize()
            : clqIdx(0), size(0) {
    }

    CliqueSize(size_t clqIdx_, size_t size_)
            : clqIdx(clqIdx_), size(size_) {
    }

    size_t clqIdx;
    size_t size;

    bool operator<(const CliqueSize &other) const {
        return other.size < this->size;
    }

    CliqueSize &operator=(const CliqueSize &other) {
        this->clqIdx = other.clqIdx;
        this->size = other.size;
        return *this;
    }
};

CglCliqueMerging::CglCliqueMerging()
        : maxExtensions_(2), maxItBK_(4096), nExtended_(0), nDominated_(0), handler_(NULL), defaultHandler_(true) {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(2);
    messages_ = CglMessage();
}

CglCliqueMerging::CglCliqueMerging(const CglCliqueMerging &rhs)
        : maxExtensions_(rhs.maxExtensions_), maxItBK_(rhs.maxItBK_), nExtended_(rhs.nExtended_),
          nDominated_(rhs.nDominated_), defaultHandler_(rhs.defaultHandler_) {
    if (defaultHandler_) {
        handler_ = new CoinMessageHandler();
        handler_->setLogLevel(rhs.handler_->logLevel());
    } else {
        handler_ = rhs.handler_;
    }
    messages_ = rhs.messages_;
}

CglCliqueMerging &CglCliqueMerging::operator=(const CglCliqueMerging &rhs) {
    if (this != &rhs) {
        gutsOfDestructor();
        maxExtensions_ = rhs.maxExtensions_;
        maxItBK_ = rhs.maxItBK_;
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

CglCliqueMerging::~CglCliqueMerging() {
    gutsOfDestructor();
}

void CglCliqueMerging::gutsOfDestructor() {
    nExtended_ = 0;
    nDominated_ = 0;

    if (defaultHandler_) {
        delete handler_;
        handler_ = NULL;
    }
}

// Pass in Message handler (not deleted at end)
void CglCliqueMerging::passInMessageHandler(CoinMessageHandler *handler) {
    if (defaultHandler_)
        delete handler_;
    defaultHandler_ = false;
    handler_ = handler;
}

// Set language
void CglCliqueMerging::newLanguage(CoinMessages::Language language) {
    messages_ = CglMessage(language);
}

/* tries to extend every clique in model using
 * its conflict graph. Dominated cliques are removed */
void CglCliqueMerging::mergeCliques(OsiSolverInterface &model) {
    nExtended_ = nDominated_ = 0;

    const CoinPackedMatrix *cpmRow = model.getMatrixByRow();
    const CoinBigIndex *starts = cpmRow->getVectorStarts();
    const double *Arhs = model.getRightHandSide();
	const CGraph *cgraph = model.getCGraph();
    clock_t startExtend;

    double *rc = new double[model.getNumCols() * 2]();
    size_t *idx = new size_t[model.getNumCols()];

    // list of rows which are cliques
    std::vector<size_t> cliques;

    // if it is a clique or not and if it is dominated
    std::vector<enum CliqueType> cliqueState(model.getNumRows(), NotAClique);

    CliqueSet *newCliques = clq_set_create(); // where new cliques will be stored

    size_t *nColClqs = NULL; // number of cliques that a column is involved
    size_t *allCollClqs = NULL; // contiguous vector to store column cliques
    size_t **colClqs = NULL; // cliques that a column is involved

    /* here cliques are stored only with column indexes
       * complimentary variables are stores with indexes
       * numCols + idxVar, so that coefficients are not necessary
       * to store */
    // clique elements per row
    SpsMtx origClqs(model.getNumRows(), model.getNumElements());

    bool *iv = new bool[model.getNumCols() * 2](); // incidence vector for columns (including compliment)
    bool *ivr = new bool[model.getNumRows()](); // incidence vector for rows

    size_t *rtc = new size_t[model.getNumRows()]; // rows to check per thread

    std::vector<CliqueSize> clqsSize;

    OsiSolverInterface::OsiNameVec clqNames;

    /* new rows */
    std::vector<int> nrStart; nrStart.push_back(0);
    std::vector<int> nrIdx;
    std::vector<double> nrCoef;
    std::vector<double> nrLB;
    std::vector<double> nrUB;

    int *clqIdxs = new int[model.getNumCols()];
    double *clqCoefs = new double[model.getNumCols()];

    /* row names */
    OsiSolverInterface::OsiNameVec rowNames;

    const size_t nCliques = detect_cliques(model, cliqueState, cliques, origClqs);

    if (nCliques == 0) {
        goto TERMINATE;
    }

    /* filling cliques per col */
    {
        size_t totnz = 0;
        nColClqs = new size_t[model.getNumCols()]();
        for (size_t i = 0; i < nCliques; i++) {
            size_t ir = cliques[i];
            size_t nzr = origClqs.nz(ir);
            const size_t *row = origClqs.row(ir);
#ifdef DEBUG
            assert(nzr >= 2);
#endif

            for (size_t j = 0; (j < nzr); ++j) {
                size_t col = row[j];
                if (col >= model.getNumCols()) {
                    col -= model.getNumCols();
                }
                ++nColClqs[col];
                ++totnz;
            }
        }
        allCollClqs = new size_t[totnz];
        colClqs = new size_t*[model.getNumCols()];
        // start of each column
        colClqs[0] = allCollClqs;
        for (size_t i = 1; i < model.getNumCols(); i++) {
            colClqs[i] = colClqs[i - 1] + nColClqs[i - 1];
        }

        std::fill(nColClqs, nColClqs + model.getNumCols(), 0);

        // filling cliques of each col
        for (size_t i = 0; i < nCliques; i++) {
            size_t ir = cliques[i];
            size_t nzr = origClqs.nz(ir);
            const size_t *row = origClqs.row(ir);
#ifdef DEBUG
            assert(nzr >= 2);
#endif

            for (size_t j = 0; j < nzr; j++) {
                size_t col = row[j];

                if (col >= model.getNumCols()) {
                    col -= model.getNumCols();
                }

                colClqs[col][nColClqs[col]++] = ir;
            }
        }
    }

    startExtend = clock();
    for (size_t iclq = 0; iclq < nCliques; iclq++) {
        size_t row = cliques[iclq];
        size_t nz = cpmRow->getVectorLengths()[row];

        const double *ccoef = cpmRow->getElements() + starts[row];
        const int *cidx = cpmRow->getIndices() + starts[row];

        std::copy(cidx, cidx + nz, idx);

        if (cliqueState[row] != Dominated) {
            for (size_t ii = 0; ii < nz; ii++) {
                if (ccoef[ii] <= -1e-5) {
                    idx[ii] += model.getNumCols();
                }
            }

            /* extending */
            {
                CliqueExtender *clqe = clqe_create(cgraph);
                //                printf("aaa cgraph: %p\n", cgraph ); fflush(stdout); fflush(stderr);
                //                printf("rc %p cgraph size %d\n", (void*) rc, cgraph_size(cgraph) ); fflush(stdout); fflush(stderr);
                clqe_set_max_it_bk(clqe, maxItBK_);
                clqe_set_costs(clqe, rc, cgraph_size(cgraph));

                clock_t startext = clock();
                size_t status = clqe_extend(clqe, idx, nz, model.getNumCols(), CLQEM_EXACT);
                //clqMergeSecsExtend += ((double)clock()-startext)/(double)CLOCKS_PER_SEC;
                if (status > 0) {
                    if (clqMergeVerbose >= 2) {
                        printf("-> constraint %s extended: \n", model.getRowName(row).c_str());
                    }

                    //++clqMergeNExtended;
                    ++nExtended_;

                    // to sort cliques found per size
                    clqsSize.clear();

                    const CliqueSet *clqs = clqe_get_cliques(clqe);
#ifdef DEBUG
                    assert(clq_set_number_of_cliques(clqs) >= 1);
#endif

                    size_t nCliquesToInsert = std::min(maxExtensions_, clq_set_number_of_cliques(clqs));
                    if (nCliquesToInsert) {
                        cliqueState[row] = Dominated;
                        if (nCliquesToInsert == 1) {
#if __cplusplus >= 201103L
                            clqsSize.emplace_back(0, clq_set_clique_size(clqs, 0));
#else
                            clqsSize.push_back(CliqueSize(0, clq_set_clique_size(clqs, 0)));
#endif
                        } else {
                            for (size_t ic = 0; ic < clq_set_number_of_cliques(clqs); ic++) {
#if __cplusplus >= 201103L
                                clqsSize.emplace_back(ic, clq_set_clique_size(clqs, ic));
#else
                                clqsSize.push_back(CliqueSize(ic, clq_set_clique_size(clqs, ic)));
#endif
                            }

                            sort(clqsSize.begin(), clqsSize.end());
#ifdef DEBUG
                            assert(clqsSize[0].size >= clqsSize[1].size);
#endif
                        }
                    }

                    for (size_t ic = 0; ic < nCliquesToInsert; ic++) {
                        size_t idxclique = clqsSize[ic].clqIdx;
                        size_t size = clq_set_clique_size(clqs, idxclique);
                        //printf("   adding %d\n", size );
                        const size_t *el = clq_set_clique_elements(clqs, idxclique);

                        add_clique(model, newCliques, size, el, cliqueState, origClqs, iv, ivr, rtc,
                                (const size_t**)colClqs, nColClqs);
                        char origrname[256] = "";
                        strncpy(origrname, model.getRowName(row).c_str(), 256);
                        char extn[64] = "";
                        if (nCliquesToInsert > 1) {
                            sprintf(extn, "xt(%ld)", ic);
                            strcat(origrname, extn);
                        }
                        if (clqMergeVerbose >= 2) {
                            printf("\t%s : ", origrname);
                        }
                        clqNames.push_back(origrname);

                        if (clqMergeVerbose >= 2) {
                            for (size_t k = 0; k < size; k++) {
                                char strneg[3] = "";
                                size_t icol = el[k];
#ifdef DEBUG
                                assert(icol < std::numeric_limits<int>::max());
#endif

                                if (icol >= model.getNumCols()) {
                                    strcpy(strneg, "!");
                                    icol -= model.getNumCols();
                                }
                                printf("%s%s ", strneg, model.getColName(icol).c_str());
                            }
                            printf("\n");
                        }
                    }
                }
                clqe_free(&clqe);
            }
        } // non dominated clique constraints
    } // all clique constraints

    //clqMergeSecsExtendAndDominate = ((double)clock()-(double)startExtend)/((double)CLOCKS_PER_SEC);

    {
        clock_t startRemAdd = clock();
        /* removing dominated cliques */
        {
            int nToRemove = 0;
            int *toRemove = new int[model.getNumRows()];

            for (size_t i = 0; i < model.getNumRows(); i++) {
                if (cliqueState[i] == Dominated) {
                    if (model.getRowSense()[i] == 'E') {
                        //++clqMergeNDominatedEqPart;
                        // adding >= part, <= part will be added separately
                        int nz = cpmRow->getVectorLengths()[i];
                        double rhs = Arhs[i];
                        char rname[256];
                        strncpy(rname, model.getRowName(i).c_str(), 256);
                        char nrname[512];
                        sprintf(nrname, "%sEp", rname);

                        const double *crcoef = cpmRow->getElements() + starts[i];
                        const int *cridx = cpmRow->getIndices() + starts[i];

                        nrIdx.insert(nrIdx.end(), cridx, cridx + nz);
                        nrCoef.insert(nrCoef.end(), crcoef, crcoef + nz);
                        nrLB.push_back(rhs);
                        nrUB.push_back(DBL_MAX);
#ifdef DEBUG
                        assert(*(nrStart.rbegin()) + nz < std::numeric_limits<int>::max());
#endif
                        nrStart.push_back(*(nrStart.rbegin()) + nz);
                        rowNames.push_back(nrname);


                    } else {
                        //++clqMergeNDominatedFull;
                        ++nDominated_;
                    }
                    toRemove[nToRemove++] = i;
                }
            }

            if (nToRemove) {
                model.deleteRows(nToRemove, toRemove);
            }

            delete[] toRemove;
        }

        /* adding stronger cliques */
        if (clq_set_number_of_cliques(newCliques) > 0) {
            size_t nNewCliques = clq_set_number_of_cliques(newCliques);

            for (size_t ic = 0; ic < nNewCliques; ic++) {
                size_t size = clq_set_clique_size(newCliques, ic);
                const size_t *el = clq_set_clique_elements(newCliques, ic);
                double rhs = 1.0;

                for (size_t i = 0; i < size; i++) {
                    if (el[i] >= model.getNumCols()) {
#ifdef DEBUG
                        assert(el[i] - model.getNumCols() < std::numeric_limits<int>::max());
#endif
                        clqIdxs[i] = el[i] - model.getNumCols();
                        clqCoefs[i] = -1.0;
                        rhs -= 1.0;
                    } else {
#ifdef DEBUG
                        assert(el[i] < std::numeric_limits<int>::max());
#endif
                        clqIdxs[i] = el[i];
                        clqCoefs[i] = 1.0;
                    }
                }

                // adding
                nrIdx.insert(nrIdx.end(), clqIdxs, clqIdxs + size);
                nrCoef.insert(nrCoef.end(), clqCoefs, clqCoefs + size);
                nrLB.push_back(-DBL_MAX);
                nrUB.push_back(rhs);
#ifdef DEBUG
                assert(*(nrStart.rbegin()) + size < std::numeric_limits<int>::max());
#endif
                nrStart.push_back(*(nrStart.rbegin()) + size);
                rowNames.push_back(clqNames[ic]);
            }
        }

        //clqMergeSecsAddAndRemove = ((double)clock()-(double)startRemAdd)/((double)CLOCKS_PER_SEC);
    }

    /* flushing rows */
    {
        int crIdx = model.getNumRows();
#ifdef DEBUG
        assert(nrStart.size() - 1 < std::numeric_limits<int>::max());
#endif
        model.addRows(nrStart.size() - 1, &nrStart[0], &nrIdx[0], &nrCoef[0], &nrLB[0], &nrUB[0]);
        for (int i = 0; i < (int) rowNames.size(); ++i) {
            model.setRowName(crIdx + i, rowNames[i]);
        }
    }

    /*
      if (clqMergeVerbose)
      {
          printf("%d extended, %d dom full, %d dom eq.\n", clqMergeNExtended, clqMergeNDominatedFull, clqMergeNDominatedEqPart );
          fflush( stdout );
      }*/

    TERMINATE:
    clq_set_free(&newCliques);

    delete[] iv;
    delete[] ivr;
    delete[] rtc;
    delete[] clqIdxs;
    delete[] clqCoefs;
    delete[] idx;
    delete[] rc;

    if (nColClqs) {
        delete[] nColClqs;
    }
    if (allCollClqs) {
        delete[] allCollClqs;
    }
    if (colClqs) {
        delete[] colClqs;
    }

    //remover
    handler_->message(CGL_PROCESS_CLQMRG, messages_) << nExtended_ << nDominated_ << CoinMessageEol;
}
