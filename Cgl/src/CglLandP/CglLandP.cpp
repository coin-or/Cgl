// Copyright (C) 2005-2008, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     07/21/05
//---------------------------------------------------------------------------
#include "CglLandP.hpp"
#include "CglLandPSimplex.hpp"
#ifdef COIN_USE_XPR
#include "OsiXprSimplexInterface.hpp"
#endif
#define INT_INFEAS(value) fabs(value - floor(value+0.5))

#include "OsiClpSolverInterface.hpp"

#define CLONE_SI //Solver is cloned between two cuts

#include "CoinTime.hpp"
#include "CglGomory.hpp"
#include "CoinFactorization.hpp"
#include <fstream>
namespace LAP
{
//Setup output messages
LapMessages::LapMessages( )
        :CoinMessages(LAP_MESSAGES_DUMMY_END)
{
    strcpy(source_,"Lap");
    addMessage(BEGIN_ROUND,CoinOneMessage( 1, 2,"Starting %s round %d variable considered for separation."));
    addMessage(END_ROUND,CoinOneMessage(2, 2,"End ouf %s round %d cut generated in %g seconds."));
    addMessage(DURING_SEP,CoinOneMessage(3,1,"After %g seconds, separated %d cuts."));
    addMessage(CUT_REJECTED, CoinOneMessage(4,1,"Cut rejected for %s."));
    addMessage(CUT_FAILED,CoinOneMessage(5,1,"Generation failed."));
    addMessage(LAP_CUT_FAILED_DO_MIG, CoinOneMessage(3006,1,"Failed to generate a cut generate a Gomory cut instead"));
}
}
using namespace LAP;
CglLandP::Parameters::Parameters():
        CglParam(),
        pivotLimit(20),
        pivotLimitInTree(10),
        maxCutPerRound(50),
        failedPivotLimit(1),
        degeneratePivotLimit(0),
        extraCutsLimit(5),
        pivotTol(1e-4),
        away(5e-4),
        timeLimit(DBL_MAX),
        singleCutTimeLimit(DBL_MAX),
        rhsWeight(1.),
        useTableauRow(true),
        modularize(false),
        strengthen(true),
        countMistakenRc(false),
        sepSpace(Fractional),
        perturb(true),
        normalization(Unweighted),
        rhsWeightType(Fixed),
        lhs_norm(L1),
        generateExtraCuts(none),
        pivotSelection(mostNegativeRc)
{
    EPS = 1e-08;
}

CglLandP::Parameters::Parameters(const Parameters &other):
        CglParam(other),
        pivotLimit(other.pivotLimit),
        pivotLimitInTree(other.pivotLimitInTree),
        maxCutPerRound(other.maxCutPerRound),
        failedPivotLimit(other.failedPivotLimit),
        degeneratePivotLimit(other.degeneratePivotLimit),
        extraCutsLimit(other.extraCutsLimit),
        pivotTol(other.pivotTol),
        away(other.away),
        timeLimit(other.timeLimit),
        singleCutTimeLimit(other.singleCutTimeLimit),
        rhsWeight(other.rhsWeight),
        useTableauRow(other.useTableauRow),
        modularize(other.modularize),
        strengthen(other.strengthen),
        countMistakenRc(other.countMistakenRc),
        sepSpace(other.sepSpace),
        perturb(other.perturb),
        normalization(other.normalization),
        rhsWeightType(other.rhsWeightType),
        lhs_norm(other.lhs_norm),
        generateExtraCuts(other.generateExtraCuts),
        pivotSelection(other.pivotSelection)
{}

CglLandP::Parameters & CglLandP::Parameters::operator=(const Parameters &other)
{
    if (this != &other) {
        CglParam::operator=(other);
        pivotLimit = other.pivotLimit;
        pivotLimitInTree = other.pivotLimitInTree;
        maxCutPerRound = other.maxCutPerRound;
        failedPivotLimit = other.failedPivotLimit;
        degeneratePivotLimit = other.failedPivotLimit;
        extraCutsLimit = other.extraCutsLimit;
        pivotTol = other.pivotTol;
        away = other.away;
        timeLimit = other.timeLimit;
        singleCutTimeLimit = other.singleCutTimeLimit;
        rhsWeight = other.rhsWeight;
        useTableauRow = other.useTableauRow;
        modularize = other.modularize;
        strengthen = other.strengthen;
        countMistakenRc = other.countMistakenRc;
        sepSpace = other.sepSpace;
        perturb = other.perturb;
        normalization = other.normalization;
        rhsWeightType = other.rhsWeightType;
        lhs_norm = other.lhs_norm;
        generateExtraCuts = other.generateExtraCuts;
        pivotSelection = other.pivotSelection;
    }
    return *this;
}

CglLandP::CachedData::CachedData(int nBasics, int nNonBasics):
        basics_(NULL), nonBasics_(NULL), nBasics_(nBasics),
        nNonBasics_(nNonBasics), basis_(NULL), colsol_(NULL),
        slacks_(NULL), integers_(NULL)
{
    if (nBasics_>0) {
        basics_ = new int[nBasics_];
        integers_ = new bool [nNonBasics_ + nBasics_];
    }
    if (nNonBasics_>0)
        nonBasics_ = new int[nNonBasics_];
    if (nBasics_ + nNonBasics_ > 0) {
        colsol_ = new double[nBasics_ + nNonBasics_];
        slacks_ = &colsol_[nNonBasics_];
    }
}

CglLandP::CachedData::CachedData(const CachedData &source):
        basics_(NULL), nonBasics_(NULL), nBasics_(source.nBasics_),
        nNonBasics_(source.nNonBasics_), basis_(NULL),
        colsol_(NULL), slacks_(NULL), integers_(NULL)
{
    if (nBasics_>0) {
        basics_ = new int[nBasics_];
        CoinCopyN(source.basics_, nBasics_, basics_);
        integers_ = new bool [nNonBasics_ + nBasics_];
        CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);
    }
    if (nNonBasics_>0) {
        nonBasics_ = new int[nNonBasics_];
        CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
    }
    if (nBasics_ + nNonBasics_ > 0) {
        colsol_ = new double[nBasics_ + nNonBasics_];
        slacks_ = &colsol_[nNonBasics_];
        CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
    }
    if (source.basis_!=NULL)
        basis_ = new CoinWarmStartBasis(*source.basis_);
}

CglLandP::CachedData& CglLandP::CachedData::operator=(const CachedData &source)
{
    if (this != &source) {
        nBasics_ = source.nBasics_;
        nNonBasics_ = source.nNonBasics_;
        if (basics_ == NULL) delete [] basics_;
        basics_ = NULL;
        if (nonBasics_ == NULL) delete [] nonBasics_;
        nonBasics_ = NULL;
        if (basis_ == NULL) delete [] basis_;
        basis_ = NULL;
        if (colsol_ == NULL) delete [] colsol_;
        colsol_ = NULL;
        if (slacks_ == NULL) delete [] slacks_;
        slacks_ = NULL;
        if (integers_ == NULL) delete [] integers_;
        integers_ = NULL;
        if (nBasics_>0) {
            basics_ = new int[nBasics_];
            CoinCopyN(source.basics_, nBasics_, basics_);
            integers_ = new bool [nBasics_ + nNonBasics_];
            CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);
        }
        if (nNonBasics_>0) {
            nonBasics_ = new int[nNonBasics_];
            CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
        }
        if (nBasics_ + nNonBasics_ > 0) {
            colsol_ = new double[nBasics_ + nNonBasics_];
            slacks_ = &colsol_[nNonBasics_];
            CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
        }
        if (source.basis_!=NULL)
            basis_ = new CoinWarmStartBasis(*source.basis_);
    }
    return *this;
}

void
CglLandP::CachedData::getData(const OsiSolverInterface &si)
{
    int nBasics = si.getNumRows();
    int nNonBasics = si.getNumCols();
    if (basis_ != NULL)
        delete basis_;
    basis_ = dynamic_cast<CoinWarmStartBasis *> (si.getWarmStart());
    if (!basis_)
        throw NoBasisError();

    if (nBasics_ > 0 || nBasics != nBasics_) {
        delete [] basics_;
        basics_ = NULL;
    }
    if (basics_ == NULL) {
        basics_ = new int[nBasics];
        nBasics_ = nBasics;
    }

    if (nNonBasics_ > 0 || nNonBasics != nNonBasics_) {
        delete [] nonBasics_;
        nonBasics_ = NULL;
    }
    if (nonBasics_ == NULL) {
        nonBasics_ = new int[nNonBasics];
        nNonBasics_ = nNonBasics;
    }
    int n = nBasics + nNonBasics;
    if ( nBasics_ + nNonBasics_ > 0 || nBasics_ + nNonBasics_ != n) {
        delete [] colsol_;
        delete [] integers_;
        integers_ = NULL;
        colsol_ = NULL;
        slacks_ = NULL;
    }
    if (colsol_ == NULL) {
        colsol_ = new double[n];
        slacks_ = &colsol_[nNonBasics];
    }

    if (integers_ == NULL) {
        integers_ = new bool[n];
    }

    const double * rowLower = si.getRowLower();
    const double * rowUpper = si.getRowUpper();
    //determine which slacks are integer
    const CoinPackedMatrix * m = si.getMatrixByCol();
    const double * elems = m->getElements();
    const int * inds = m->getIndices();
    const CoinBigIndex * starts = m->getVectorStarts();
    const int * lenghts = m->getVectorLengths();
    //    int numElems = m->getNumElements();
    int numCols = m->getNumCols();
    assert(numCols == nNonBasics_);
    //   int numRows = m->getNumRows();
    CoinFillN(integers_ ,n, true);
    for (int i = 0 ;  i < numCols ; i++) {
        if (si.isContinuous(i))
            integers_[i] = false;
    }
    bool * integerSlacks = integers_ + numCols;
    for (int i = 0 ; i < nBasics ; i++) {
        if (rowLower[i] > -1e50 && INT_INFEAS(rowLower[i]) > 1e-15)
            integerSlacks[i] = false;
        if (rowUpper[i] < 1e50 && INT_INFEAS(rowUpper[i]) > 1e-15)
            integerSlacks[i] = false;
    }
    for (int i = 0 ;  i < numCols ; i++) {
        CoinBigIndex end = starts[i] + lenghts[i];
        if (integers_[i]) {
            for (CoinBigIndex k=starts[i] ; k < end; k++) {
                if (integerSlacks[inds[k]] && INT_INFEAS(elems[k])>1e-15 )
                    integerSlacks[inds[k]] = false;
            }
        } else {
            for (CoinBigIndex k=starts[i] ; k < end; k++) {
                if (integerSlacks[inds[k]])
                    integerSlacks[inds[k]] = false;
            }
        }
    }

    CoinCopyN(si.getColSolution(), si.getNumCols(), colsol_);
    CoinCopyN(si.getRowActivity(), si.getNumRows(), slacks_);
    for (int i = 0 ; i < si.getNumRows() ; i++) {
        slacks_[i]*=-1;
        if (rowLower[i]>-1e50) {
            slacks_[i] += rowLower[i];
        } else {
            slacks_[i] += rowUpper[i];
        }
    }
    //Now get the fill the arrays;
    nNonBasics = 0;
    nBasics = 0;



    //For having the index variables correctly ordered we need to access to OsiSimplexInterface
    {
        OsiSolverInterface * ncSi = (const_cast<OsiSolverInterface *>(&si));
        ncSi->enableSimplexInterface(0);
        ncSi->getBasics(basics_);
        ncSi->disableSimplexInterface();
    }

    int numStructural = basis_->getNumStructural();
    for (int i = 0 ; i < numStructural ; i++) {
        if (basis_->getStructStatus(i)== CoinWarmStartBasis::basic) {
            nBasics++;
            //Basically do nothing
#ifdef LANDP_DEBUG
            if (nBasics>nBasics_) {
                std::cerr<<"Error in number of basic variables"<<std::endl;
                throw CoinError("Unexpected number of basic variables (what is going on)","CglLandP::CachedData","GetData");
            }
#endif
        } else {
            nonBasics_[nNonBasics++] = i;
#ifdef LANDP_DEBUG
            if (nNonBasics>nNonBasics_) {
                std::cerr<<"Error in number of non-basic variables"<<std::endl;
                throw CoinError("Unexpected number of non-basic variables (what is going on)","CglLandP::CachedData","GetData");
            }
#endif
        }
    }

    int numArtificial = basis_->getNumArtificial();
    for (int i = 0 ; i < numArtificial ; i++) {
        if (basis_->getArtifStatus(i)== CoinWarmStartBasis::basic) {
            //Just check number of basics
            nBasics++;

#ifdef LANDP_DEBUG
            if (nBasics>nBasics_) {
                std::cerr<<"Error in number of basic variables"<<std::endl;
                throw CoinError("Unexpected number of basic variables (what is going on)","CglLandP::CachedData","GetData");
            }
#endif
        } else {
            nonBasics_[nNonBasics++] = i + basis_->getNumStructural();
#ifdef LANDP_DEBUG
            if (nNonBasics>nNonBasics_) {
                std::cerr<<"Error in number of non-basic variables"<<std::endl;
                throw CoinError("Unexpected number of non-basic variables (what is going on)","CglLandP::CachedData","GetData");
            }
#endif
        }
    }
#ifdef LANDP_DEBUG
    //Check that the expected number of basics and non-basics is found
    if (nBasics!=nBasics_) {
        std::cerr<<"Warning Number of basics variable is not the one expected"<<std::endl;
        while (nBasics<nBasics_)
            basics_[nBasics++] = -1;
    }
    if (nNonBasics!=nNonBasics_) {
        std::cerr<<"Warning Number of basics variable is not the one expected"<<std::endl;
        while (nNonBasics<nNonBasics_)
            nonBasics_[nNonBasics++] = -1;
    }
#endif
}

CglLandP::CachedData::~CachedData()
{
    if (basics_!=NULL)
        delete [] basics_;
    if (nonBasics_!=NULL)
        delete [] nonBasics_;
    if (colsol_ != NULL)
        delete [] colsol_;
    delete basis_;
    if (integers_)
        delete [] integers_;
}

CglLandP::CglLandP(const CglLandP::Parameters &params,
                   const LAP::Validator &validator):
        params_(params), cached_(), validator_(validator), numcols_(-1),
        originalColLower_(NULL), originalColUpper_(NULL),
        canLift_(false),
        extraCuts_()
#ifdef DO_STAT
        , roundsStats_()
#endif
{
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(0);
    messages_ = LapMessages();
}


CglLandP::~CglLandP()
{
    delete handler_;
    if (originalColLower_ != NULL)
        delete [] originalColLower_;
    if (originalColUpper_ != NULL)
        delete [] originalColUpper_;
}

CglLandP::CglLandP(const CglLandP & source):
        params_(source.params_), cached_(source.cached_),
        validator_(source.validator_), numcols_(source.numcols_),
        originalColLower_(NULL), originalColUpper_(NULL),
        canLift_(source.canLift_),
        extraCuts_(source.extraCuts_)
#ifdef DO_STAT
        ,
        roundsStats_(source.roundsStats_)
#endif
{
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(source.handler_->logLevel());
    messages_ = LapMessages();
    if (numcols_ != -1) {
        assert(numcols_ > 0);
        assert(originalColLower_!=NULL);
        assert(originalColUpper_!=NULL);
        originalColLower_ = new double[numcols_];
        originalColUpper_ = new double[numcols_];
        CoinCopyN(source.originalColLower_,numcols_,originalColLower_);
        CoinCopyN(source.originalColUpper_,numcols_,originalColUpper_);
    }
}

/** Assignment operator */
CglLandP& CglLandP::operator=(const CglLandP &rhs)
{
    if (this != &rhs) {
        params_ = rhs.params_;
        cached_ = rhs.cached_;
        validator_ = rhs.validator_;
        extraCuts_ = rhs.extraCuts_;
    }
    return *this;
}


CglCutGenerator *
CglLandP::clone() const
{
    return new CglLandP(*this);
}

extern double restaurationTime;


//Quick if vector are sorted otherwise falls back on CoinPackedVector procedure
bool equals(const CoinPackedVector& x, const CoinPackedVector& y, const CoinRelFltEq &eq
            , double yScaling)
{
    const double * xVal = x.getElements();
    const double * yVal = y.getElements();

    const int * xInd = x.getIndices();
    const int * yInd = y.getIndices();

    const int& xNelem = x.getNumElements();
    const int& yNelem = y.getNumElements();
    //check that two vetors are ordered
    for (int i = 1 ; i < xNelem ; i++) {
        if (xInd[i] <= xInd[i - 1]) {
            return x.isEquivalent(y,eq);
        }
    }
    for (int i = 1 ; i < yNelem ; i++) {
        if (yInd[i] <= yInd[i - 1]) {
            return x.isEquivalent(y,eq);
        }
    }
    int i =0, j = 0;
    for ( ; i < xNelem && j < yNelem; i++,j++) {
        if (yInd[j]*yScaling < xInd[i]) {
            if (!eq(yVal[j],0.))
                return false;
            i--;
        } else if (yInd[j]*yScaling > xInd[i]) {
            if (!eq(xVal[i],0.))
                return false;
            j--;
        } else if (!eq(xVal[i],yVal[j]*yScaling)) {
            return false;
        }
    }
    for (; i < xNelem ; i++) {
        if (!eq(xVal[i],0.))
            return false;
    }
    for (;j < yNelem ; j++) {
        if (!eq(yVal[j]*yScaling,0.))
            return false;
    }
    return true;
}

//Quick if vector are sorted otherwise falls back on CoinPackedVector procedure
bool similarInfNorm(const CoinPackedVector& x, const CoinPackedVector& y,
                    double limit, double yScaling)
{
    const double * xVal = x.getElements();
    const double * yVal = y.getElements();

    const int * xInd = x.getIndices();
    const int * yInd = y.getIndices();

    const int& xNelem = x.getNumElements();
    const int& yNelem = y.getNumElements();
    //check that two vetors are ordered
    for (int i = 1 ; i < xNelem ; i++) {
        if (xInd[i] <= xInd[i - 1]) {
            CoinPackedVector diff = x - yScaling * y;
            return diff.infNorm() < limit;
        }
    }
    for (int i = 1 ; i < yNelem ; i++) {
        if (yInd[i] <= yInd[i - 1]) {
            CoinPackedVector diff = x - yScaling * y;
            return diff.infNorm() < limit;
        }
    }
    int i =0, j = 0;
    for ( ; i < xNelem && j < yNelem; i++,j++) {
        if (yInd[j] < xInd[i]) {
            if (fabs(yVal[j]*yScaling) > limit) return false;
            i--;
        } else if (yInd[j] > xInd[i]) {
            if (fabs(xVal[i]) > limit) return false;
            j--;
        } else {
            if (fabs(xVal[i] - yVal[j]*yScaling) > limit) return false;
        }
    }
    for (; i < xNelem ; i++) {
        if (fabs(xVal[i]) > limit) return false;
    }
    for (;j < yNelem ; j++) {
        if (fabs(yVal[j]*yScaling) > limit) return false;
    }
    return true;
}

//Quick if vector are sorted otherwise falls back on CoinPackedVector procedure
double similarTwoNorm(const CoinPackedVector& x, const CoinPackedVector& y,
                      double limit, double yScaling)
{
    const double * xVal = x.getElements();
    const double * yVal = y.getElements();

    const int * xInd = x.getIndices();
    const int * yInd = y.getIndices();

    const int& xNelem = x.getNumElements();
    const int& yNelem = y.getNumElements();
    //check that two vetors are ordered
    for (int i = 1 ; i < xNelem ; i++) {
        if (xInd[i] <= xInd[i - 1]) {
            CoinPackedVector diff = x - y;
            std::cout<<"Vector are not sorted"<<std::endl;
            return diff.twoNorm() < limit;
        }
    }
    for (int i = 1 ; i < yNelem ; i++) {
        if (yInd[i] <= yInd[i - 1]) {
            CoinPackedVector diff = x - y;
            std::cout<<"Vector are not sorted"<<std::endl;
            return diff.twoNorm() < limit;
        }
    }
    limit *= limit;
    int i =0, j = 0;
    double v = 0;
    double yScalingSquare = yScaling * yScaling;
    for ( ; i < xNelem && j < yNelem; i++,j++) {
        if (yInd[j] < xInd[i]) {
            v += yVal[j] * yVal[j] * yScalingSquare;
            if (v>limit) {
                return false;
            }
            i--;
        } else if (yInd[j] > xInd[i]) {
            v += xVal[i] * xVal[i];
            if (v>limit) {
                return false;
            }
            j--;
        } else {
            double diff = xVal[i] - yVal[j] * yScaling;
            v += diff*diff;
            if (v>limit) {
                return false;
            }
        }
    }
    for (; i < xNelem ; i++) {
        v += xVal[i] * xVal[i];
        if (v>limit) {
            return false;
        }
    }
    for (;j < yNelem ; j++) {
        v += yVal[j] * yVal[j] * yScalingSquare;
        if (v>limit) {
            return false;
        }
    }
    return true;
}

struct cutsCos {
    int i;
    int j;
    double angle;
    cutsCos(int  i_, int j_ , double angle_):i(i_), j(j_), angle(angle_) {
    }
    bool operator<(const cutsCos&other)const {
        return angle > other.angle;
    }
};


void
CglLandP::scanExtraCuts(OsiCuts& cs, const double * colsol) const
{
    int numAdded = 0;
    for (int i = extraCuts_.sizeRowCuts() - 1; i > -1 ; i--) {
        double violation = extraCuts_.rowCut(i).violated(colsol);
        if (violation > 0.) {
            cs.insert(extraCuts_.rowCut(i));
            numAdded++;
//      std::cout<<"A cut computed in a previous iteration is violated by "<<violation<<"."<<std::endl;
            //extraCuts_.eraseRowCut(i);
        }
    }
//  std::cout<<"Added "<<numAdded<<" previously generated cuts."<<std::endl;
}

void
CglLandP::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
                       const CglTreeInfo info ) const
{
    if (info.pass == 0 && !info.inTree) {
        numrows_ = si.getNumRows();
    }
    scanExtraCuts(cs, si.getColSolution());
    Parameters params = params_;
    params.rhsWeight = numrows_ + 2;


#if DO_STAT
    printf("CUTGAP After %i pass objective is %g\n",info.pass, si.getObjValue());
    if (info.pass == 0 && !info.inTree) { //lookup miplib problem
        std::string name;
        si.getStrParam(OsiProbName,name);
        roundsStats_.lookupProblem(name.c_str());
    }
#endif

    if (info.inTree) { //put lower pivot limit
        params.pivotLimit = std::min(params.pivotLimit, params.pivotLimitInTree);
        params.countMistakenRc = true;
    }
    if (params.timeLimit < 0) {
        params.pivotLimit = 0;
    }
#ifdef LANDP_DEBUG
    if (!si.basisIsAvailable()) {
        std::cerr<<"No basis!!!"<<std::endl;
        throw -1;
    }
#endif

    cached_.getData(si);
    CglLandPSimplex landpSi(si,cached_, params, validator_);
    if (params.generateExtraCuts == CglLandP::AllViolatedMigs) {
        landpSi.genThisBasisMigs(cached_, params);
    }
    landpSi.setLogLevel(handler_->logLevel());
    int nCut = 0;

#ifdef DO_STAT
    roundsStats_.addStatistic();
    roundStatistic & stat = roundsStats_.rounds.back();
    stat.doubleStat[roundStatistic::Time] = -CoinCpuTime();
#endif
    //get integer fractionals and sort
#if 0
    CoinPackedVector xFrac;
    int nFrac = getSortedFractionals(xFrac, cached_, params);
    const int * indices = xFrac.getIndices();
#endif

    std::vector<int> indices;
    getSortedFractionalIndices(indices,cached_, params);

    int numrows = si.getNumRows();


#ifdef DO_STAT
    //Get informations on current optimum
    {
        OsiSolverInterface * gapTester = si.clone();
        gapTester->resolve();

        roundsStats_.analyseOptimalBasis(gapTester,info.pass, numrows_);
        delete gapTester;
    }
#endif

    params_.timeLimit += CoinCpuTime();
    CoinRelFltEq eq(1e-04);

    for (unsigned int i = 0; i < indices.size() && nCut < params.maxCutPerRound &&
            nCut < cached_.nBasics_ ; i++) {

        //Check for time limit
        int iRow = indices[i];
        assert(iRow < numrows);
        OsiRowCut cut;
        int code=1;
        OsiSolverInterface * ncSi = NULL;
#ifndef DO_STAT
        if (params.pivotLimit != 0)
#endif
        {
            ncSi = si.clone();
            landpSi.setSi(ncSi);
            ncSi->setDblParam(OsiDualObjectiveLimit, DBL_MAX);
            ncSi->messageHandler()->setLogLevel(0);
        }

        int generated = 0;
#ifndef DO_STAT
        if (params.pivotLimit == 0) {
            generated = landpSi.generateMig(iRow, cut, cached_, params);
        } else
#endif
        {
            generated = landpSi.optimize(iRow, cut, cached_, params);
            if (params.generateExtraCuts == CglLandP::AllViolatedMigs) {
                landpSi.genThisBasisMigs(cached_, params);
            }
            landpSi.resetSolver(cached_.basis_);
        }
        code = 0;
        if (generated)
            code = validator_(cut, cached_.colsol_, si, params, originalColLower_, originalColUpper_);
        if (!generated || code) {
#ifdef DO_STAT
            if (!code)
                stat.intStat[roundStatistic::NumFailures]++;
                std::cout<<"Failed to generate a cut"<<std::endl;
            else {
                stat.intStat[roundStatistic::NumRejected]++;
                handler_->message(CUT_REJECTED, messages_)<<validator_.failureString(code)<<CoinMessageEol;
            }
#endif
            if (params.pivotLimit !=0) {
                handler_->message(LAP_CUT_FAILED_DO_MIG, messages_)<<validator_.failureString(code)<<CoinMessageEol;
                landpSi.freeSi();
                OsiSolverInterface * ncSi = si.clone();
                landpSi.setSi(ncSi);
                params.pivotLimit = 0;
                if (landpSi.optimize(iRow, cut, cached_, params)) {
                    code = validator_(cut, cached_.colsol_, si, params, originalColLower_, originalColUpper_);
                }
                params.pivotLimit = params_.pivotLimit;
            }
        }

#ifndef DO_STAT
        if (params.pivotLimit != 0)
#endif
        {
            landpSi.freeSi();
        }
        if (code) {
            handler_->message(CUT_REJECTED, messages_)<<
                validator_.failureString(code)<<CoinMessageEol;
#ifdef DO_STAT
            stat.intStat[roundStatistic::NumRejected]++;
#endif
        } else {
            if (canLift_) {
                cut.setGloballyValid(true);
            }
            cs.insertIfNotDuplicate(cut, eq);
            //cs.insert(cut);
          {
                //std::cout<<"Violation "<<cut.violated(cached_.colsol_)<<std::endl;
#ifdef DO_STAT
                roundStatistic oneCut;
#if LandP_DEBUG > 1
                oneCut.meanNnz = landpSi.extra.Nnz;
                oneCut.meanNneg = landpSi.extra.Nneg;
#endif
                oneCut.doubleStat[roundStatistic::NumPivots] = landpSi.extra.numPivots;;
                oneCut.intStat[roundStatistic::NumCuts] = 1;
                oneCut.doubleStat[roundStatistic::CutViolation] = landpSi.extra.violation;
                oneCut.doubleStat[roundStatistic::CutDepth] = cut.violated(si.getColSolution())/cut.row().twoNorm();
                oneCut.doubleStat[roundStatistic::AngleToObj] = acos(cos(cut.row(),
                        si.getObjCoefficients(), si.getNumCols()));
                if (landpSi.extra.numPivots) {
                    oneCut.doubleStat[roundStatistic::NegativeRcRows] = landpSi.extra.nNegativeRcRows;
                    oneCut.doubleStat[roundStatistic::BestRow] = landpSi.extra.bestRow;
                    oneCut.doubleStat[roundStatistic::MaxBestRow] = landpSi.extra.bestRow;
                    oneCut.doubleStat[roundStatistic::BestRc] = landpSi.extra.bestRc;
                    oneCut.doubleStat[roundStatistic::MaxRc] = landpSi.extra.maxRc;
                }

                roundsStats_.addCut(oneCut);
#endif
                nCut++;
            }
        }
    }

    Cuts& extra = landpSi.extraCuts();
    for (int i = 0 ; i < cached_.nNonBasics_; i++) {
        OsiRowCut * cut = extra.rowCut(i);
        if (cut == NULL) continue;
        int code = validator_(*cut, cached_.colsol_, si, params,
                              originalColLower_, originalColUpper_);
        if (code) {
            handler_->message(LAP_CUT_FAILED_DO_MIG, messages_)<<validator_.failureString(code)<<CoinMessageEol;
#ifdef DO_STAT
            stat.intStat[roundStatistic::NumRejected]++;
#endif
        } else {
            cs.insertIfNotDuplicate(*cut, eq); {
                nCut++;
#ifdef DO_STAT
                stat.intStat[roundStatistic::NumExtra]++;
#endif
            }
        }
        delete cut;
    }


    params_.timeLimit -= CoinCpuTime();

#ifdef DO_STAT
    if (stat.doubleStat[roundStatistic::NumPivots] > 0.) {
        double numPivots = (double) stat.doubleStat[roundStatistic::NumPivots];
        stat.doubleStat[roundStatistic::BestRow] /= numPivots;
        stat.doubleStat[roundStatistic::NegativeRcRows] /= numPivots;
        stat.doubleStat[roundStatistic::BestRc] /= numPivots;
    }
    if (nCut > 0) {
        stat.intStat[roundStatistic::NumCuts] = nCut;
        stat.doubleStat[roundStatistic::Time] += CoinCpuTime();

        OsiSolverInterface * gapTester = si.clone();
        gapTester->applyCuts(cs);
        gapTester->resolve();



        stat.doubleStat[roundStatistic::Bound] = gapTester->getObjValue();
        delete gapTester;
        if (!info.inTree) {
            roundsStats_.displayRound();
        }
    }
#endif

}


template < class S, class T, class U >
class StableCompare {
public:
    inline bool operator()(const CoinTriple<S,T,U>& t1,
                           const CoinTriple<S,T,U>& t2) const
    {
        return (t1.third < t2.third) ||
               ((t1.third == t2.third) && (t1.second < t2.second));
    }

};

template <class T1,class T2>
struct StableExternalComp{
    const std::vector<T1> &vec_1_;
    const std::vector<T2> &vec_2_;
    StableExternalComp(const std::vector<T1> &vec_1,
                       const std::vector<T2> &vec_2):
            vec_1_(vec_1),
            vec_2_(vec_2){
    }
    CoinRelFltEq eq;
    bool operator()(int i, int j){
        bool result = (vec_1_[i] < vec_1_[j]) ||
                      ( ((vec_1_[i]== vec_1_[j]))
                        && (vec_2_[i] < vec_2_[j]));
        return result;
    }

};
void
CglLandP::getSortedFractionalIndices(std::vector<int> &frac_indices,
                                     const CachedData &data,
                                     const CglLandP::Parameters & params) const{
    std::vector<int> colIndices;
    std::vector<double> values;
    std::vector<int> indices;
    for (int i = 0 ; i < data.nBasics_ ; i++) {
        const int& iCol = data.basics_[i];
        if (iCol >= data.nNonBasics_ ||
                !data.integers_[iCol] ||
                INT_INFEAS(data.colsol_[iCol]) <= params.away)
            continue;
        const double value = INT_INFEAS(data.colsol_[iCol]);

        frac_indices.push_back(i);
        indices.push_back(values.size());
        values.push_back(- value);
        colIndices.push_back(iCol);
    }
    std::sort(indices.begin(), indices.end(),StableExternalComp<double, int>(values,colIndices));
    colIndices = frac_indices;
    for (unsigned int i = 0; i < indices.size() ; i++){
        frac_indices[i] = colIndices[indices[i]];
    }

}


