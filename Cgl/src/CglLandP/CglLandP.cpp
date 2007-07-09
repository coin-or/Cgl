// Copyright (C) 2005, Pierre Bonami and others.  All Rights Reserved.
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
    addMessage(LAP_CUT_FAILED_DO_MIG, CoinOneMessage(3006,2,"Failed to generate a cut generate a Gomory cut instead"));
  }
}
using namespace LAP;
  CglLandP::Parameters::Parameters():
    pivotLimit(20),
    pivotLimitInTree(10),
    maxCutPerRound(50),    
    failedPivotLimit(1),
    degeneratePivotLimit(0),
    pivotTol(1e-4),
    away(5e-4),
    timeLimit(DBL_MAX),
    singleCutTimeLimit(DBL_MAX),
    useTableauRow(true),
    modularize(false),
    strengthen(true),
    reducedSpace(true),
    perturb(true),
    scaleCuts(false),
  pivotSelection(mostNegativeRc)
  {
  }
  
  CglLandP::Parameters::Parameters(const Parameters &other):
    pivotLimit(other.pivotLimit),
    pivotLimitInTree(other.pivotLimitInTree),
    maxCutPerRound(other.maxCutPerRound),
    failedPivotLimit(other.failedPivotLimit),
    degeneratePivotLimit(other.degeneratePivotLimit),
    pivotTol(other.pivotTol),
    away(other.away),
    timeLimit(other.timeLimit),
    singleCutTimeLimit(other.singleCutTimeLimit),
    useTableauRow(other.useTableauRow),
    modularize(other.modularize),
    strengthen(other.strengthen),
    reducedSpace(other.reducedSpace),
    perturb(other.perturb),
    scaleCuts(other.scaleCuts),
  pivotSelection(other.pivotSelection)
{}
 
CglLandP::Parameters & CglLandP::Parameters::operator=(const Parameters &other){
  if(this != &other)
  {
    pivotLimit = other.pivotLimit;
    pivotLimitInTree = other.pivotLimitInTree;
    maxCutPerRound = other.maxCutPerRound;
    failedPivotLimit = other.failedPivotLimit;
    degeneratePivotLimit = other.failedPivotLimit;
    pivotTol = other.pivotTol;
    away = other.away;
    timeLimit = other.timeLimit;
    singleCutTimeLimit = other.singleCutTimeLimit;
    useTableauRow = other.useTableauRow;
    modularize = other.modularize;
    strengthen = other.strengthen;
    reducedSpace = other.reducedSpace;
    perturb = other.perturb;
    scaleCuts = other.scaleCuts;
    pivotSelection = other.pivotSelection;
  }
  return *this;
}

  CglLandP::CachedData::CachedData(int nBasics, int nNonBasics):
    basics_(NULL), nonBasics_(NULL), nBasics_(nBasics),
    nNonBasics_(nNonBasics), basis_(NULL), colsol_(NULL),
    slacks_(NULL), integers_(NULL)
  {
    if(nBasics_>0)
      {
	basics_ = new int[nBasics_];
	integers_ = new bool [nNonBasics_ + nBasics_];
      }
    if(nNonBasics_>0)
      nonBasics_ = new int[nNonBasics_];
    if(nBasics_ + nNonBasics_ > 0)
      {
	colsol_ = new double[nBasics_ + nNonBasics_];
	slacks_ = &colsol_[nNonBasics_];
      }
  }
  
  CglLandP::CachedData::CachedData(const CachedData &source):
    basics_(NULL), nonBasics_(NULL), nBasics_(source.nBasics_),
    nNonBasics_(source.nNonBasics_), basis_(NULL),
    colsol_(NULL), slacks_(NULL), integers_(NULL)
  {
    if(nBasics_>0)
      {
	basics_ = new int[nBasics_];
	CoinCopyN(source.basics_, nBasics_, basics_);
	integers_ = new bool [nNonBasics_ + nBasics_];
	CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);
      }
    if(nNonBasics_>0)
      {
	nonBasics_ = new int[nNonBasics_];
	CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
      }
    if(nBasics_ + nNonBasics_ > 0)
      {
	colsol_ = new double[nBasics_ + nNonBasics_];
	slacks_ = &colsol_[nNonBasics_];
	CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
      }
    if(source.basis_!=NULL)
      basis_ = new CoinWarmStartBasis(*source.basis_);
  }

CglLandP::CachedData& CglLandP::CachedData::operator=(const CachedData &source){
  if(this != &source)
  {
    nBasics_ = source.nBasics_;
    nNonBasics_ = source.nNonBasics_; 
    if(basics_ == NULL) delete [] basics_; basics_ = NULL;
    if(nonBasics_ == NULL) delete [] nonBasics_; nonBasics_ = NULL;
    if(basis_ == NULL) delete [] basis_; basis_ = NULL;
    if(colsol_ == NULL) delete [] colsol_; colsol_ = NULL;
    if(slacks_ == NULL) delete [] slacks_; slacks_ = NULL;
    if(integers_ == NULL) delete [] integers_; integers_ = NULL;
    if(nBasics_>0)
    {
      basics_ = new int[nBasics_];
      CoinCopyN(source.basics_, nBasics_, basics_);
      integers_ = new bool [nBasics_ + nNonBasics_];
      CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);     
    }
    if(nNonBasics_>0)
    {
      nonBasics_ = new int[nNonBasics_];
      CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
    }
    if(nBasics_ + nNonBasics_ > 0)
    {
      colsol_ = new double[nBasics_ + nNonBasics_];
      slacks_ = &colsol_[nNonBasics_];
      CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
  }
    if(source.basis_!=NULL)
      basis_ = new CoinWarmStartBasis(*source.basis_);
  }
  return *this;
}

  void 
  CglLandP::CachedData::getData(const OsiSolverInterface &si)
{
    int nBasics = si.getNumRows();
    int nNonBasics = si.getNumCols();
    if(basis_ != NULL)
      delete basis_;
    basis_ = dynamic_cast<CoinWarmStartBasis *> (si.getWarmStart());
    if(!basis_)
      throw NoBasisError();
    
    if(nBasics_ > 0 || nBasics != nBasics_)
    {
      delete [] basics_;
      basics_ = NULL;
    }
    if(basics_ == NULL)
    {
      basics_ = new int[nBasics];
      nBasics_ = nBasics;
    }
    
    if(nNonBasics_ > 0 || nNonBasics != nNonBasics_)
    {
      delete [] nonBasics_;
      nonBasics_ = NULL;
    }
    if(nonBasics_ == NULL)
    {
      nonBasics_ = new int[nNonBasics];
      nNonBasics_ = nNonBasics;
    }
    int n = nBasics + nNonBasics;
    if( nBasics_ + nNonBasics_ > 0 || nBasics_ + nNonBasics_ != n)
    {
      delete [] colsol_;
      delete [] integers_;
      integers_ = NULL;
      colsol_ = NULL;
      slacks_ = NULL;
    }
    if(colsol_ == NULL)
    {
      colsol_ = new double[n];
      slacks_ = &colsol_[nNonBasics];
    }
    
    if(integers_ == NULL){
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
    for(int i = 0 ;  i < numCols ; i++)
    {
      if(si.isContinuous(i))
        integers_[i] = false;
    }
    bool * integerSlacks = integers_ + numCols;
    for(int i = 0 ; i < nBasics ; i++)
      {
	if(rowLower[i] > -1e50 && INT_INFEAS(rowLower[i]) > 1e-15)
	  integerSlacks[i] = false;
	if(rowUpper[i] < 1e50 && INT_INFEAS(rowUpper[i]) > 1e-15)
	  integerSlacks[i] = false;
      }
    for(int i = 0 ;  i < numCols ; i++)
    {
      CoinBigIndex end = starts[i] + lenghts[i];
      if(integers_[i])
      {
        for(CoinBigIndex k=starts[i] ; k < end; k++)
	      {
          if(integerSlacks[inds[k]] && INT_INFEAS(elems[k])>1e-15 )
            integerSlacks[inds[k]] = false;
	      }
      }
      else
      {
        for(CoinBigIndex k=starts[i] ; k < end; k++)
	      {
		if(integerSlacks[inds[k]])
		  integerSlacks[inds[k]] = false;
	      }
      }
    }
    
    CoinCopyN(si.getColSolution(), si.getNumCols(), colsol_);
    CoinCopyN(si.getRowActivity(), si.getNumRows(), slacks_);
    for(int i = 0 ; i < si.getNumRows() ; i++)
    {
      slacks_[i]*=-1;
      if(rowLower[i]>-1e50)
      {
        slacks_[i] += rowLower[i];
      }
      else
      {
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
    
    for(int i = 0 ; i < basis_->getNumStructural() ; i++)
    {
      if(basis_->getStructStatus(i)== CoinWarmStartBasis::basic)
      {
        nBasics++;
        //Basically do nothing
#ifdef LANDP_DEBUG
        if(nBasics>nBasics_)
	      {
          std::cerr<<"Error in number of basic variables"<<std::endl;
          throw CoinError("Unexpected number of basic variables (what is going on)","CglLandP::CachedData","GetData");
	      }
#endif
      }
      else
      {
        nonBasics_[nNonBasics++] = i;
#ifdef LANDP_DEBUG
        if(nNonBasics>nNonBasics_)
	      {
          std::cerr<<"Error in number of non-basic variables"<<std::endl;
          throw CoinError("Unexpected number of non-basic variables (what is going on)","CglLandP::CachedData","GetData");
	      }
#endif
      }
    }
    
    for(int i = 0 ; i < basis_->getNumArtificial() ; i++)
    {
      if(basis_->getArtifStatus(i)== CoinWarmStartBasis::basic)
      {
        //Just check number of basics
        nBasics++;
        
#ifdef LANDP_DEBUG
        if(nBasics>nBasics_)
	      {
          std::cerr<<"Error in number of basic variables"<<std::endl;
          throw CoinError("Unexpected number of basic variables (what is going on)","CglLandP::CachedData","GetData");
	      }
#endif
      }
      else
      {
        nonBasics_[nNonBasics++] = i + basis_->getNumStructural();
#ifdef LANDP_DEBUG
        if(nNonBasics>nNonBasics_)
	      {
          std::cerr<<"Error in number of non-basic variables"<<std::endl;
          throw CoinError("Unexpected number of non-basic variables (what is going on)","CglLandP::CachedData","GetData");
	      }
#endif
      }
    }
#ifdef LANDP_DEBUG
    //Check that the expected number of basics and non-basics is found
    if(nBasics!=nBasics_)
    {
      std::cerr<<"Warning Number of basics variable is not the one expected"<<std::endl;
      while(nBasics<nBasics_)
        basics_[nBasics++] = -1;
    }
    if(nNonBasics!=nNonBasics_)
    {
      std::cerr<<"Warning Number of basics variable is not the one expected"<<std::endl;
      while(nNonBasics<nNonBasics_)
        nonBasics_[nNonBasics++] = -1;
    }
#endif
}
  
  CglLandP::CachedData::~CachedData()
  {
    if(basics_!=NULL)
      delete [] basics_;
    if (nonBasics_!=NULL)
      delete [] nonBasics_;
    if(colsol_ != NULL)
      delete [] colsol_;
    delete basis_;
    if(integers_)
      delete [] integers_;
  }
  
  CglLandP::CglLandP(const CglLandP::Parameters &params,
		     const CglValidator &validator):
    params_(params), cached_(), validator_(validator)
#ifdef DO_STAT
    ,
    roundsStats_(),
    used(false)
#endif
  {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(0);
    messages_ = LapMessages();
#ifdef DO_STAT
    logStream = &std::cout;
    logAssignCount = NULL;
#endif
  }
  
  
  CglLandP::~CglLandP(){
#if DO_STAT
    releaseLogStream();
    if(used) displayStats();
#endif
    delete handler_;
  }
  
  CglLandP::CglLandP(const CglLandP & source):
    params_(source.params_), cached_(source.cached_), 
    validator_(source.validator_)
#ifdef DO_STAT
    ,
    roundsStats_(source.roundsStats_),
    logAssignCount(source.logAssignCount),
    used(false)
#endif
  {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(source.handler_->logLevel());
    messages_ = LapMessages();  
#ifdef DO_STAT
    logStream = source.logStream;
    if(logAssignCount && *logAssignCount) (*logAssignCount)++;
#endif
  }

/** Assignment operator */
CglLandP& CglLandP::operator=(const CglLandP &rhs)
{
  if(this != &rhs)
  {
    params_ = rhs.params_;
    cached_ = rhs.cached_;
    validator_ = rhs.validator_;
#ifdef DO_STAT
    releaseLogStream();
    logAssignCount = rhs.logAssignCount;
    logStream = rhs.logStream;
    if(logAssignCount && *logAssignCount) (*logAssignCount)++;
#endif
  }
  return *this;
}


  CglCutGenerator * 
  CglLandP::clone() const
  { 
    return new CglLandP(*this);
  }
  
  extern double restaurationTime;

  void
  CglLandP::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			 const CglTreeInfo info ) const
{
#if DO_STAT
  used = true;
#endif
  Parameters params = params_;
#if DO_STAT
  if(info.pass == 0 && !info.inTree)//lookup miplib problem
  {
    std::string name;
    si.getStrParam(OsiProbName,name);
    roundsStats_.lookupProblem(name.c_str(), si.getObjValue());
  }
#endif
  if(info.inTree)//put lower pivot limit
  {
    params.pivotLimit = min(params.pivotLimit, params.pivotLimitInTree);
  }
  if(params.timeLimit < 0){
    params.pivotLimit = 0;
  }
#ifdef LANDP_DEBUG
  if(!si.basisIsAvailable())
  {
    std::cerr<<"No basis!!!"<<std::endl;
    throw -1;
  }
#endif
 
  cached_.getData(si);
  CglLandPSimplex landpSi(si,cached_, params.reducedSpace, params.pivotLimit);
  landpSi.setLogLevel(handler_->logLevel());
  int nCut = 0;
  
#ifdef DO_STAT 
  roundStatistic stat;
  stat.time = -CoinCpuTime();
#endif
  //get integer fractionals and sort
  CoinPackedVector xFrac;
  int numcols = si.getNumCols();
  for(int i = 0 ; i < numcols ; i++)
  {
    const double &value = si.getColSolution()[i];
    if(si.isInteger(i) &&
       cached_.basis_->getStructStatus(i) == CoinWarmStartBasis::basic &&
       INT_INFEAS(value) > params_.away)
    {
      xFrac.insert(i, - INT_INFEAS(value));
    }
  }
  xFrac.sortIncrElement();
  int nFrac =xFrac.getNumElements();
  const int * indices = xFrac.getIndices();
  int numrows = si.getNumRows();
  int i = 0;
 
  params_.timeLimit += CoinCpuTime();
  
  for(; i < nFrac && nCut < params.maxCutPerRound && 
      nCut < cached_.nBasics_ ; i++)
  {
    
    //Check for time limit
    int iRow = 0;     
    for(; iRow < numrows ; iRow++)
    {
      if(cached_.basics_[iRow]== indices[i])
        break;
    }
    assert(iRow < numrows);
    OsiRowCut cut;
    int code=1;
    OsiSolverInterface * ncSi = NULL;
#ifndef DO_STAT
    if(params.pivotLimit != 0)
#endif
   {
      ncSi = si.clone();
      landpSi.setSi(ncSi);
      ncSi->setDblParam(OsiDualObjectiveLimit, DBL_MAX);
      ncSi->messageHandler()->setLogLevel(0);
    }
    int generated = 0;
#ifndef DO_STAT
    if(params.pivotLimit == 0)
    {
      generated = landpSi.generateMig(iRow, cut, cached_, params);
    }
    else
#endif
    {
      generated = landpSi.findBestCut(iRow, cut, cached_, params);
    }
    if(generated)
      code = validator_(cut, cached_.colsol_, si);
    if(!generated || code) 
    {
#ifdef DO_STAT
      if(!code)
      stat.numberFailures++;
      else stat.numberRejected++;
#endif
      if(params.pivotLimit !=0){
        handler_->message(CUT_REJECTED, messages_)<<CoinMessageEol;
        
        landpSi.freeSi();
        OsiSolverInterface * ncSi = si.clone();
        landpSi.setSi(ncSi);
        params.pivotLimit = 0;
        if(landpSi.findBestCut(iRow, cut, cached_, params))
        {
          code = validator_(cut, cached_.colsol_, si);
        }
        params.pivotLimit = params_.pivotLimit;
      }
    }
    
#ifndef DO_STAT 
    if(params.pivotLimit != 0)
#endif
    {    
      landpSi.freeSi();
    }
    if(code)
    {
	      handler_->message(LAP_CUT_FAILED_DO_MIG, messages_)<<validator_.failureString(code)<<CoinMessageEol;
#ifdef DO_STAT
	      stat.numberRejected++;
#endif	      
	    }
      else
	    {
              CoinRelFltEq eq(1e-04);
	      cs.insertIfNotDuplicate(cut, eq);
#ifdef DO_STAT
#if LandP_DEBUG > 1
	      stat.meanNnz += landpSi.extra.Nnz;
	      stat.meanNneg += landpSi.extra.Nneg;
#endif
	      stat.meanAngleToObj += landpSi.extra.AngleToObj;
	      stat.meanNumPivots += landpSi.extra.numPivots;
	      stat.averageCutViolation += landpSi.extra.depth;
	      if(landpSi.extra.numPivots)
        {
          stat.meanNegativeRcRows += ((double) landpSi.extra.nNegativeRcRows);
          stat.meanBestRow += ((double) landpSi.extra.bestRow);
          stat.maxBestRow = max((double) landpSi.extra.bestRow, stat.maxBestRow);
          stat.meanBestRc += landpSi.extra.bestRc;
          stat.maxRc = max((double) stat.maxRc, landpSi.extra.maxRc);
        }
#endif
	      nCut++;
	    }
  }
  params_.timeLimit -= CoinCpuTime();
  
#ifdef DO_STAT
  if(stat.meanNumPivots > 0.)
  {
    stat.meanBestRow /= (double) stat.meanNumPivots;
    stat.meanNegativeRcRows /= (double) stat.meanNumPivots;
    stat.meanBestRc /= (double) stat.meanNumPivots;
  }
  if(nCut > 0)
  {
    stat.numberCuts = nCut;
    stat.meanNnz /=(double) nCut;
    stat.meanNneg /= (double) nCut;
    stat.meanAngleToObj /= (double) nCut;
    stat.meanNumPivots /= (double) nCut;
    stat.averageCutViolation /= (double) nCut;
    stat.time += CoinCpuTime();
    
    
    OsiSolverInterface * gapTester = si.clone();
    gapTester->applyCuts(cs);
    gapTester->resolve();
    stat.bound = gapTester->getObjValue();
    delete gapTester;
    if(!info.inTree){ 
      roundsStats_.addStatistic(stat); 
      roundsStats_.displayRound(*logStream);
    }
  }
#endif

#ifdef LAP_DEBUG
  std::cout<<"Total restauration time :"<<restaurationTime<<std::endl;
#endif
}

#ifdef DO_STAT
void 
CglLandP::lookupProblem(const std::string &s){
  OsiClpSolverInterface si;
  si.readMps(s.c_str());
  si.initialSolve();
  roundsStats_.lookupProblem(s.c_str(), si.getObjValue());
}
#endif

