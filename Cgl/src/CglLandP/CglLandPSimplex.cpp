// Copyright (C) 2005, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     21/07/05
//---------------------------------------------------------------------------
#include "CglLandPSimplex.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinIndexedVector.hpp"
#define THROW_ON_BAD_REDUCED_COST  //throw -1

#define REMOVE_LOG 0

#define RED_COST_CHECK 1e-6

#define OLD_COMPUTATION
#include <algorithm>
//#define TEST_M3
namespace LAP
{
  
  LandPMessages::LandPMessages()
  :
  CoinMessages(DUMMY_END)
{
    strcpy(source_,"Lap");
    addMessage(Separating, CoinOneMessage(1,3+REMOVE_LOG,"Starting separation on variable %d, initial depth of cut %f"));
    addMessage(FoundImprovingRow, CoinOneMessage(2,4,"Found improving row (leaving variable). Row %d (basic var %d), leaving status %d, sign of gamma %d"));
    addMessage(FoundBestImprovingCol, CoinOneMessage(3,4," Found best improvement (entering variable). Var %d, value of gamma %f, expected depth of next cut %f"));
    addMessage(WarnFailedBestImprovingCol, CoinOneMessage(6001,3,"Failed to find an improving entering variable while reduced cost was %f, depth of current cut %f, best cut depth with pivot %f"));
    //Log line is cut number time pivot number,  cut depth, leaving, incoming gamma degenerate
    addMessage(LogHead, CoinOneMessage(5,3+REMOVE_LOG, "Pivot no \t cut depth \t leaving var \t incoming var \t direction \t gamma \t degenerate"));
    addMessage(PivotLog, CoinOneMessage(6,3+REMOVE_LOG,"%8d\t %9f\t %11d \t %11d \t %11d \t %8f \t %12d \t %.5g \t %11d"));
    addMessage(FinishedOptimal, CoinOneMessage(7,2,"Found optimal lift-and-project cut, depth %f number of pivots performed %d"));
    addMessage(HitLimit, CoinOneMessage(8,2,"Stopping lift-and-project optimization hit %s limit. Number of pivots %d"));
    addMessage(WarnBadSigmaComputation, CoinOneMessage(6002,1,"Cut depth after pivot is not what was expected by computations before, difference %.15f"));
    addMessage(WarnBadRowComputation, CoinOneMessage(6003, 1, "Row obtained after pivot is not what was expected (distance between the two %f in norm inf)."));
    addMessage(WarnGiveUpRow,CoinOneMessage(6004,1,"Limit of %d negative reduced costs with no strict improvement"));
    addMessage(PivotFailedSigmaUnchanged,CoinOneMessage(6005,1,"A pivot failed to be performed (probably refactorization was performed) but sigma is unchanged continue..."));
    addMessage(PivotFailedSigmaIncreased,CoinOneMessage(6006,1,"A pivot failed to be performed, and sigma has changed exit without generating cut"));
    addMessage(FailedSigmaIncreased,CoinOneMessage(6006,1,"Cut violation has increased in last pivot"));
    addMessage
      (WarnBadRhsComputation, 
       CoinOneMessage(6007,1,
                      "rhs obtained  after pivot is not what was expected (distance between the two %f)."));
    addMessage(WarnFailedPivotTol, 
               CoinOneMessage(6008, 2,"All pivots are below tolerance")); 
    addMessage(WarnFailedPivotIIf, 
               CoinOneMessage(6009, 2,"There is no possible pivot within tolerance (every pivot make rhs for current row %f too close to integer feasibility")); 
}

double restaurationTime;



void CglLandPSimplex::printTableau(std::ostream & os)
{
  int width = 9;
  os<<"Tableau at current basis"<<std::endl;
  os<<"    ";
  //Head with non basics indices
  for(int i = 0 ; i < numcols_ ; i++)
  {
    os.width(width);
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    std::cout<<nonBasics_[i]<<" ";
  }
  
  os.width(width);
  os.setf(std::ios_base::right, std::ios_base::adjustfield);
  std::cout<<'b';
  
  os<<std::endl;
  
  //print row by row
  for(int i = 0 ; i < numrows_ ; i++)
  {
    //    int ind = basics_[i];
    row_i_.num = i;
    pullTableauRow(row_i_);
    row_i_.print(os, width, nonBasics_, numcols_);
    
  }
  
}



CglLandPSimplex::CglLandPSimplex(const OsiSolverInterface &si,
                                 const CglLandP::CachedData &cached,
                                 bool reducedSpace, int pivotLimit):
gammas_(false),
rWk1_(NULL),
rWk2_(NULL),
rWk3_(NULL),
rWk4_(NULL),
rIntWork_(NULL),
rowFlags_(NULL),
colHasToComputeContrib_(NULL),
colCandidateToLeave_(NULL),
basics_(NULL), nonBasics_(NULL),
inM1_(NULL), inM2_(NULL), inM3_(NULL),
tau_(0), sigma_(0), colsolToCut_(NULL),
colsol_(NULL),
numcols_(0),numrows_(0),
loBounds_(NULL),
upBounds_(NULL),
inDegenerateSequence_(false),
chosenReducedCostVal_(1e100),
si_(NULL),
nNegativeRc_(0) 
#ifdef LandP_DEBUG
, 
debug_(si.getNumCols(),si.getNumRows())
#endif
{
  restaurationTime = 0.;
  numcols_ = si.getNumCols();
  numrows_ = si.getNumRows();
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  messages_ = LandPMessages();
  
  si_ = const_cast<OsiSolverInterface *>(&si);
  
  OsiClpSolverInterface * clpSi = dynamic_cast<OsiClpSolverInterface *>(si_);
  if(clpSi) solver_ = clp;
#ifdef COIN_USE_XPR
  else
  {
    OsiXprSolverInterface * xprSi = dynamic_cast<OsiXprSolverInterface *>(si_);
    if(xprSi) solver_ = xpress;
#endif
#ifdef COIN_USE_CPX
    else
    {
      OsiXprSolverInterface * cpxSi = dynamic_cast<OsiXprSolverInterface *>(si_);
      if(cpxSi) solver_ = cplex;
    }
#ifdef COIN_USE_XPR
  }
#endif
#endif
  int rowsize = numcols_ + numrows_ + 1;
  row_k_.size(rowsize);
  loBounds_ = new double[numcols_ + numrows_];
  upBounds_ = new double[numcols_ + numrows_];
  
  CoinCopyN(si.getColLower(),numcols_, loBounds_);
  CoinCopyN(si.getColUpper(),numcols_,upBounds_);
  const double * rowUpper = si.getRowUpper();
  const double * rowLower = si.getRowLower();
  double infty = si.getInfinity();
  int i=numcols_;
  for(int iRow = 0; iRow < numrows_ ;iRow++, i++)
  {
    if(rowUpper[iRow] < infty)
      loBounds_[i]=0.;
    else loBounds_[i]= - infty;
    if(rowLower[iRow] <= - infty)
      upBounds_[i] = infty;
    else if(rowUpper[iRow] < infty)
      upBounds_[i] = rowUpper[iRow] - rowLower[iRow];
    else
      upBounds_[i] = 0.;
  }
#ifndef DO_STAT
  if(pivotLimit != 0)
#endif
  {
    own_ = true;
    rWk1_ = new double [numrows_];
    rWk2_ = new double [numrows_];
    rWk3_ = new double [numrows_];
    rWk4_ = new double [numrows_];
    rIntWork_ = new int[numrows_];

    row_i_.size(rowsize);
    newRow_.size(rowsize);
    rowFlags_ = new bool[numrows_];
    colHasToComputeContrib_ = new bool[numcols_ + numrows_];
    colCandidateToLeave_ = new bool[numcols_];
    basics_ = new int[numrows_];
    nonBasics_ = new int[numcols_];
    inM1_ = new int[numcols_];
    inM2_ = new int[numcols_];
    inM3_ = new int[numcols_];
    colsolToCut_ = new double[numcols_ + numrows_];
    colsol_ = new double[numcols_ + numrows_];
  }
#ifndef DO_STAT
  else{
    own_ = false;
    si_->enableSimplexInterface(0);
    basis_ = *cached.basis_;
  }
#endif
  cacheUpdate(cached,reducedSpace);
}

CglLandPSimplex::~CglLandPSimplex()
{
  delete handler_;
  delete [] loBounds_;
  delete [] upBounds_;
  if(own_){
  delete [] rowFlags_;
  delete [] colHasToComputeContrib_;
  delete [] colCandidateToLeave_;
  delete [] basics_;
  delete [] nonBasics_;
  delete [] colsolToCut_;
  delete [] colsol_;
  delete [] inM1_;
  delete [] inM2_;
  delete [] inM3_;
  delete []rWk1_;
  delete []rWk2_;
  delete []rWk3_;
  delete []rWk4_;
  delete []rIntWork_;
  }
  else{
    si_->disableSimplexInterface();
  }
}

void
CglLandPSimplex::cacheUpdate(const CglLandP::CachedData &cached, bool reducedSpace)
{
  integers_ = cached.integers_;
  if(own_){
    CoinCopyN(cached.basics_, numrows_, basics_);
    CoinCopyN(cached.nonBasics_, numcols_, nonBasics_);
    CoinCopyN(cached.colsol_, numrows_ + numcols_, colsol_);
    for(int i = 0 ; i < numcols_ ; i++)
    {
      colsol_[nonBasics_[i]] = 0;   
    }
    CoinCopyN(cached.colsol_, numrows_ + numcols_, colsolToCut_);
#if LandP_DEBUG > 1
    CoinCopyN(cached.colsol_, numrows_ + numcols_, debug_.trueInitialSol_);
#endif
    //Zero all non basics in colsol setup the reduced space
    CoinFillN(colHasToComputeContrib_,numcols_+numrows_,true);
    for(int i = 0 ; i < numcols_ ; i++)
    {
      colsolToCut_[nonBasics_[i]] = colsol_[nonBasics_[i]] = 0;
    }
    /** Mark the variables at zero in solution to cut so that we know that their contribution to reduced cost has to be computed*/
    if(reducedSpace)
    {
      for(int ii = 0; ii < numcols_ ; ii++)
      {
        if(colsolToCut_[ii] - upBounds_[ii] > 1e-08 || colsolToCut_[ii] - loBounds_[ii] < 1e-08)
        {
          colHasToComputeContrib_[ii]=false;
        }
      }
    }
  }
  else
  {
    basics_ = cached.basics_;
    nonBasics_ = cached.nonBasics_;
  }
#if LandP_DEBUG > 1
  si_->enableSimplexInterface(0);
  basis_ = *cached.basis_;
  debug_.getCurrentTableau(*si_,*this);
  si_->disableSimplexInterface();
#endif
}

bool CglLandPSimplex::resetSolver(const CoinWarmStartBasis * basis)
{
  double time = - CoinCpuTime();
  si_->disableSimplexInterface();
  //  if(!si_->setWarmStart(basis))
  //    {std::cerr<<"Error in basis restoration"<<std::endl;
  //    throw ;}
  
  //si.solveFromHotStart();
  //si.unmarkHotStart();
  //si_->resolve();
  restaurationTime += CoinCpuTime() + time;
  return 0;
}


//Erase last k cuts
void eraseLastCuts(OsiCuts & cuts, int k = 2)
{
  
  int numCuts = cuts.sizeRowCuts();
  int begin = numCuts - 1;
  int end = max(numCuts - 3,0);
  for(int ii = begin ; ii > end ; ii--)
  {
    cuts.eraseRowCut(ii);
  }
}

bool
CglLandPSimplex::generateMig(int row, OsiRowCut & cut,const CglLandP::CachedData &cached,const CglLandP::Parameters & params) const{
  row_k_.num = row;
  pullTableauRow(row_k_);
  row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);
  if(params.strengthen || params.modularize)
    createMIG(row_k_, cut);
  else
    createIntersectionCut(row_k_, cut);  
#ifdef DO_STAT
    double sigma = computeCglpObjective(row_k_);
    extra.depth = -sigma;
#endif

#ifdef LANDP_DEBUG
  CglLandPSimplex debug(*si_,cached, params.reducedSpace, 1);
  OsiSolverInterface * ncSi = si_->clone();
  debug.setSi(ncSi);
 // ncSi->disableSimplexInterface();
  ncSi->setDblParam(OsiDualObjectiveLimit, DBL_MAX);
  OsiRowCut cut2;
  debug.findBestCut(row, cut2, cached, params);
  if(cut!= cut2){
  cut.print();
  cut2.print();
  }
  debug.freeSi();
#endif
  return 1;//At this point nothing failed, always generate a cut
}

bool 
CglLandPSimplex::findBestCut
(int row, OsiRowCut & cut,const CglLandP::CachedData &cached,const CglLandP::Parameters & params)
{
  bool optimal = false;
  int nRowFailed = 0;
  
  double timeLimit = CoinMin(params.timeLimit, params.singleCutTimeLimit);
  timeLimit += CoinCpuTime();
  // double timeBegin = CoinCpuTime();
  OsiClpSolverInterface * lclClp = 
    dynamic_cast<OsiClpSolverInterface *>(si_);
  /** Copy the cached information */
  CoinCopyN(cached.basics_, numrows_, basics_);
  CoinCopyN(cached.nonBasics_, numcols_, nonBasics_);
  CoinCopyN(cached.colsol_, numrows_ + numcols_, colsol_);
  CoinCopyN(cached.colsol_, numrows_ + numcols_, colsolToCut_);


#ifdef ROW_PETURBATION
  if(params.perturb && perturbed_==NULL)
    perturbed_ = new int[numcols_ + numrows_];
  if(params.perturb)
    CoinFillN(perturbed_,numcols_ + numrows_, 0);
  else if(perturbed_!=NULL)
  {
    delete perturbed_;
    perturbed_ =NULL;
  }
#endif

  basis_ = *cached.basis_;
  
  for(int i = 0 ; i < numcols_ ; i++)
  {
    colsolToCut_[nonBasics_[i]] = colsol_[nonBasics_[i]] = 0;
  }
  
  //  if(params.pivotLimit)
  si_->enableSimplexInterface(0);

#ifdef LandP_DEBUG 
  TOTO
  //Check that basics_ is correct
  int * basic2 = new int [numrows_];
  si_->getBasics(basic2);
  for(int i = 0; i < numrows_ ; i++)
    assert(basics_[i]==basic2[i]);
  delete [] basic2;
#endif

  row_k_.num = row;
  pullTableauRow(row_k_);
  row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);
  if(params.modularize)
    modularizeRow(row_k_);
#ifdef ROW_PETURBATION
  if(params.perturb)
    perturbRow(row_k_,params.perturbationEpsilon,perturbed_);
#endif
  
  //  row_k_.print(numcols_ + numrows_);
  double fullSpaceSigma = computeCglpObjective(row_k_);
  
  if(fullSpaceSigma > 0)
  {
    std::cerr<<"Non separating Gomory cut?????"<<  sigma_<<std::endl;
    
    resetSolver(cached.basis_);
    return 0;
  }
  
  updateM1_M2_M3(row_k_, 0./*params.pivotTol*/, (params.reducedSpace!=0), params.perturb);
  sigma_ = computeCglpObjective(row_k_);
  
  handler_->message(Separating,messages_)<<basics_[row_k_.num]<<sigma_<<CoinMessageEol<<CoinMessageEol;
  handler_->message(LogHead, messages_)<<CoinMessageEol<<CoinMessageEol;
  
  
  extra.nNegativeRcRows = 0;
  extra.bestRow = 0;
  extra.maxBestRow = 0;
  extra.bestRc = 0;
  extra.maxRc = -DBL_MAX;
  
  //Save the variable basic in this row
  //  int var_k = cached.basics_[k_];
  
  //Put a flag on each row to say if we want to continue trying to use it
  CoinFillN(rowFlags_,numrows_,true);
  
  int numberConsecutiveDegenerate = 0;
  bool allowDegeneratePivot = numberConsecutiveDegenerate < params.degeneratePivotLimit;
  bool beObstinate = 0;
  int numPivots = 0;
  int numFailedPivots = 0;
  bool hasFlagedRow = false;
  int maxTryRow = 5;
  while (  !optimal && numPivots < params.pivotLimit)
  {
    if(timeLimit - CoinCpuTime() < 0.) break;
#ifdef GENERATE_ALL //Code to generate all the cuts found during the procedure
    {
      OsiRowCut cut;
      TabRow aRow(row_k_);
      //	pullTableauRow(aRow);
      //aRow.rhs = aRow.rhs - floor(aRow.rhs);
      if(params.strengthen || params.modularize)
        createMIG(aRow, cut);
      else
        createIntersectionCut(aRow, cut);
      cuts.insert(cut);
    }
#endif
    
#if LandP_DEBUG > 6
    if(handler_->logLevel()>=4)//Output current cut 
    {
      TabRow temp(row_k_);
      OsiRowCut cut;
      double normalization = normCoef(row_k_);
      if(params.strengthen || params.modularize)
        createMIG(temp, cut);
      else
        createIntersectionCut(temp, cut);
      if(params.scaleCuts)
        scale(cut, normalization);
      put_in_non_basic_init_space(cut);
      if(params.strengthen || params.modularize)
        std::cout<<"MIG violation :"<<- cut.violated(debug_.initialColsol_)
          <<std::endl;
      else
        std::cout<<"Intersection cut violation :"<< - cut.violated(debug_.initialColsol_)
          <<std::endl;
      
      if(handler_->logLevel()>=5)
      {
        cut.print();	  
        std::cout<<std::endl<<std::endl;
      }
      
    }
#endif
    
    
    updateM1_M2_M3(row_k_, params.pivotTol, (params.reducedSpace!=0), params.perturb);
    sigma_ = computeCglpObjective(row_k_);
    int direction = 0;
    int gammaSign = 0;
    int leaving = -1;
    int incoming = -1;
    leaving = fastFindCutImprovingPivotRow(direction, gammaSign, params.pivotTol);
    //    handler_->message(FoundImprovingRow, messages_)
    //<<row_i_.num<<basics_[row_i_.num]<<direction<<gammaSign<<CoinMessageEol;
    double bestSigma;
    if(leaving >= 0)
    {
      if(params.pivotSelection == CglLandP::mostNegativeRc)
      {
          incoming = fastFindBestPivotColumn(direction, gammaSign, 
                                             params.pivotTol, params.away, 
                                             (params.reducedSpace!=0),
                                             allowDegeneratePivot,
                                             bestSigma
                                             );
//          incoming = findBestPivotColumn(direction, params.pivotTol,
//                                         params.reducedSpace,allowDegeneratePivot, params.modularize);
        while(incoming == -1 && !optimal &&
              nRowFailed++ < maxTryRow)// if no improving was found rescan the tables of reduced cost to find a good one
	      {
          //	      std::cerr<<"Although reduced cost was <=0, no strictly improving pivot found"<<std::endl;
          rowFlags_[leaving] = false;
          if(!hasFlagedRow)
            hasFlagedRow = true;
          leaving = rescanReducedCosts(direction, gammaSign, params.pivotTol);
          if(leaving >= 0)
          {
              incoming = fastFindBestPivotColumn(direction, gammaSign, 
                                                 params.pivotTol, 
                                                 params.away,
                                                 (params.reducedSpace!=0),
                                                 allowDegeneratePivot,
                                                 bestSigma
                                                 );
//              incoming = findBestPivotColumn(direction,params.pivotTol,
//                                             params.reducedSpace,allowDegeneratePivot, params.modularize); 
          }
          else optimal = true;
        }
      }
      else
      {
        incoming = findBestPivot(leaving, direction, params);
      }
      if(incoming >= 0 && !optimal)
      {
        if(inDegenerateSequence_)//flag leaving row
        {
          numberConsecutiveDegenerate++;
          allowDegeneratePivot = numberConsecutiveDegenerate < params.degeneratePivotLimit;
          rowFlags_[leaving] = false;
        }
        else
        {
          beObstinate = 0;
          numberConsecutiveDegenerate = 0;
          allowDegeneratePivot = numberConsecutiveDegenerate < params.degeneratePivotLimit;
        }
        double gamma = - row_k_.row[nonBasics_[incoming]] / row_i_.row[nonBasics_[incoming]];
        if(fabs(row_i_.row[nonBasics_[incoming]])<=1e-12)
          std::cerr<<"Too small a pivot"<<std::endl;
        int savedStatus = -1;
        if((nonBasics_[incoming] < numcols_ && basis_.getStructStatus(nonBasics_[incoming])==CoinWarmStartBasis::atUpperBound) ||
           (nonBasics_[incoming] >= numcols_ && basis_.getArtifStatus(nonBasics_[incoming] - numcols_)==CoinWarmStartBasis::atUpperBound))
          savedStatus = 1;
        bool pivoted = changeBasis(incoming,leaving,direction);
        if(pivoted)
        {
          numPivots++;
#ifdef LandP_DEBUG
#if LandP_DEBUG > 10
           double lastSigma = sigma_;
          sigma_ = computeCglpObjective(row_k_);
          if(fabs(sigma_ - debug_.newSigma_) > RED_COST_CHECK)
          {
            handler_->message(WarnBadSigmaComputation, messages_)
            <<fabs(sigma_ - debug_.newSigma_)<<CoinMessageEol<<CoinMessageEol;
            throw -1;
          }
#endif
#endif
          if(params.modularize)
            modularizeRow(row_k_);
          sigma_ = computeCglpObjective(row_k_);

#ifdef ROW_PETURBATION
          if(perturbed_)
            perturbRow(row_k_,params.perturbationEpsilon,perturbed_);
#endif

#ifdef LandP_DEBUG          
          if(sigma_-lastSigma>1e-3*(-sigma_))
          {
            std::cerr<<"sigma has increased!!! : "<<sigma_-lastSigma<<", direction: "<<direction<<std::endl;
            changeBasis(incoming,leaving,savedStatus);
            if(params.modularize)
            {
              row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);
              for(int i = 0 ; i < numcols_ ; i++)
              {
                if(nonBasics_[i] >= numcols_) continue;
                int &ni = nonBasics_[i];
                if(si_->isInteger(ni))
                  row_k_.row[ni] = modularizedCoef(row_k_.row[ni],row_k_.rhs);
              }
            }
            std::cout<<"Reverted old, failed and new sigma "<<lastSigma<<", "<<sigma_<<computeCglpObjective(row_k_)<<std::endl;
            assert(fabs(lastSigma - computeCglpObjective(row_k_)) < 1e-08);
            sigma_= lastSigma;
            handler_->message(FailedSigmaIncreased,messages_)<<CoinMessageEol<<CoinMessageEol;  
            //		  resetSolver(cached.basis_);
            //	  eraseLastCuts(cuts);
            break;
            return 0;          
          }
#endif
          handler_->message(PivotLog,messages_)<<numPivots<<sigma_<<
            nonBasics_[incoming]<<basics_[leaving]<<direction<<gamma<<inDegenerateSequence_<<CoinMessageEol<<CoinMessageEol;
        }
        else//pivot failed
        {
          numFailedPivots++;
          //check wether sigma has changed if it has exit cut generation if it has not continue
          double lastSigma = sigma_;
          sigma_ = computeCglpObjective(row_k_);
          if( sigma_-lastSigma>1e-8)
          {
            handler_->message(PivotFailedSigmaIncreased,messages_)<<CoinMessageEol<<CoinMessageEol;
            //break;
            resetSolver(cached.basis_);
            //	    eraseLastCuts(cuts);
            return 0;
          }
          handler_->message(PivotFailedSigmaUnchanged,messages_)<<CoinMessageEol<<CoinMessageEol;
          numFailedPivots = params.failedPivotLimit + 1;
          resetSolver(cached.basis_);
          return 0;
          if(numFailedPivots > params.failedPivotLimit)
            break;
          
        }
	    }
      else //attained max number of leaving vars tries with no improvement
      {
        handler_->message(WarnGiveUpRow,messages_)<<nRowFailed<<CoinMessageEol<<CoinMessageEol;
        break;    
      }
    }
    else
    {
      if(hasFlagedRow && beObstinate)
      {
        //Reset row flags
        CoinFillN(rowFlags_,numrows_,true);
        hasFlagedRow = false;
        if(inDegenerateSequence_)
        {
          allowDegeneratePivot = false;
          beObstinate = false;
        }
      }
      else
      {
        //could perturb but Ionut skipped that will see later
        optimal = true;
        handler_->message(FinishedOptimal, messages_)<<sigma_<<numPivots<<CoinMessageEol<<CoinMessageEol;
      }
    }
  }
  
  if(!optimal && numPivots >= params.pivotLimit)
  {
    std::string limit="pivots";
    handler_->message(HitLimit, messages_)<<limit<<numPivots<<CoinMessageEol<<CoinMessageEol;
  }
  if(!optimal && numFailedPivots >= params.failedPivotLimit)
  {
    std::string limit="failed pivots";
    handler_->message(HitLimit, messages_)<<limit<<numPivots<<CoinMessageEol<<CoinMessageEol;
  }
  //Create the cut
  
  //Retake the row from the factorization
  if(0 && lclClp)//Refactorize check that everything is ok
  {
    int status = lclClp->getModelPtr()->factorize(); //statusOfProblem(0);
    if(status)
	  {
	    std::cerr<<"Problem refactorizing "<<status<<std::endl;
	    throw -1;
	  }
    status = lclClp->getModelPtr()->getSolution();
    if(status)
	  {
	    std::cerr<<"Problem refactorizing "<<status<<std::endl;
	    throw -1;
	  }
  }
  pullTableauRow(row_k_);
  row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);
  
  //  double normalization = 100*normCoef(row_k_);
#ifdef GENERATE_ALL  // Code to be used when genrating all cuts found during procedure
  if(!optimal || !params.generateAll)
#endif
  {
    if(params.strengthen || params.modularize)
      createMIG(row_k_, cut);
    else
      createIntersectionCut(row_k_, cut);
  }  
  
  extra.numPivots = numPivots;
  extra.depth = - sigma_;
  
#ifdef LandP_DEBUG 
  double cosToObj = 0;
  //(cut.cosToVec(si_->getObjCoefficients(), si_->getNumCols()));
  extra.AngleToObj = acos(cosToObj);
#if LandP_DEBUG > 1
  if(handler_->logLevel()>=3)//Output optimal
  {
    OsiRowCut cut2(cut);
    put_in_non_basic_init_space(cut2);
    if(params.strengthen || params.modularize)
      std::cout<<"MIG violation :"<<- cut.violated(debug_.trueInitialSol_)
        <<std::endl;
    else
      std::cout<<"Intersection cut violation :"<<- cut.violated(debug_.trueInitialSol_)
        <<std::endl;
    
    CoinPackedVector row = cut2.row();
    //    const int * inds = row.getIndices();
    const double * elems = row.getElements();
    int size = row.getNumElements();
    int &nNegative = extra.Nneg = 0;
    for(int k = 0 ; k < size ; k++)
    {
      if(elems[k] < -1e-10) nNegative++;
    }	
    extra.Nnz = size;
    std::cout<<"Number of non zero coef in optimal cut "<<size
      <<" Number of negative coefs"<<nNegative<<std::endl;
    std::cout<<"Cos to obj " <<cosToObj
      <<", Angle to objective "<<extra.AngleToObj<<std::endl;
    std::cout<<"Number of pivots "<<extra.numPivots<<std::endl;
    if(handler_->logLevel()>=5)
      cut2.print();
  }
#endif
#endif
  
  resetSolver(cached.basis_);
  return 1;//At this point nothing failed, always generate a cut
}


bool
CglLandPSimplex::changeBasis(int incoming, int leaving, int leavingStatus)
{
  double infty = si_->getInfinity();
  int clpLeavingStatus = leavingStatus;
  
#ifndef OLD_COMPUTATION
  adjustTableauRow(basics_[row_i_.num], row_i_, leavingStatus);
  double gamma = - row_k_.row[nonBasics_[incoming]] / row_i_.row[nonBasics_[incoming]];
  //Update row k
  for(int i = 0 ; i < numcols_; i++)
  {
    row_k_.row[nonBasics_[i]] += gamma * row_i_.row[nonBasics_[i]];
  }
  row_k_.row[basics_[leaving]] += gamma;
  row_k_.rhs += gamma * row_i_.rhs;
  resetOriginalTableauRow(basics_[row_i_.num], row_i_, leavingStatus);
#endif
  
  if(solver_ == clp)
  {
    if(basics_[leaving] >= numcols_)
      clpLeavingStatus = - leavingStatus;
  }
  int code = 0;
  try
  {
    code = si_->pivot(nonBasics_[incoming],basics_[leaving], clpLeavingStatus);
  }
  catch(int i)
  {
    if(i==-1)
      pullTableauRow(row_i_);
    throw i;
    
  }
  if (code) 
  {
    pullTableauRow(row_k_);
    row_k_.rhs = row_k_.rhs - floor(row_k_.rhs);
    return 0;
  }
  
  
  
  
  //Update basics_, non-basics_ and basis
  
  //swap bounds
	if(!code)
  {
    
    int & indexLeaving = basics_[leaving];
    colsolToCut_[indexLeaving] = leavingStatus==1?  getUpBound(indexLeaving) - colsolToCut_[indexLeaving]:colsolToCut_[indexLeaving] - getLoBound(indexLeaving);
    
    if(basics_[leaving] < numcols_)
    {
      basis_.setStructStatus(indexLeaving, leavingStatus==1 ? CoinWarmStartBasis::atUpperBound : CoinWarmStartBasis::atLowerBound);
    }
    else
    {
      int iRow = basics_[leaving] - numcols_;
      basis_.setArtifStatus(iRow,  leavingStatus==1 ? CoinWarmStartBasis::atUpperBound : CoinWarmStartBasis::atLowerBound);
      //    assert(leavingStatus==-1 || (rowLower_[iRow]>-1e50 && rowUpper_[iRow] < 1e50));
    }
    
    if(nonBasics_[incoming] < numcols_)
    {
      int & indexIncoming = nonBasics_[incoming];
      CoinWarmStartBasis::Status status = basis_.getStructStatus(indexIncoming);
      colsolToCut_[indexIncoming] = status==CoinWarmStartBasis::atUpperBound?  getUpBound(indexIncoming)
        - colsolToCut_[indexIncoming]:colsolToCut_[indexIncoming] + getLoBound(indexIncoming);
      basis_.setStructStatus(indexIncoming, CoinWarmStartBasis::basic);
    }
    else
    {
      int iRow = nonBasics_[incoming] - numcols_;
      int & indexIncoming = nonBasics_[incoming];
      colsolToCut_[indexIncoming] = basis_.getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound? getUpBound(indexIncoming)
        - colsolToCut_[indexIncoming]: colsolToCut_[indexIncoming] + getLoBound(indexIncoming);
      basis_.setArtifStatus(iRow,  CoinWarmStartBasis::basic);
    }
    
    int swap = basics_[leaving];
    basics_[leaving] = nonBasics_[incoming];
    nonBasics_[incoming] = swap;
    //update solution of leaving variable
    colsol_[nonBasics_[incoming]] = 0;
  }
  //update solution for basics
  const double * lpSol = si_->getColSolution();
  const double * rowAct = si_->getRowActivity();
  const double * rowLower = si_->getRowLower();
  const double * rowUpper = si_->getRowUpper();
  
  for(int i = 0 ; i < numrows_ ; i++)
  {
    int& iCol = basics_[i];
    if(iCol<numcols_)
      colsol_[iCol] = lpSol[iCol];
    else // Osi does not give direct acces to the value of slacks
    {
      int iRow = iCol - numcols_;
      colsol_[iCol] = - rowAct[iRow];
      if(rowLower[iRow]> -infty)
	    {
	      colsol_[iCol] += rowLower[iRow];
	    }
      else
	    {
	      colsol_[iCol] += rowUpper[iRow];
	    }
    }
  }
  // basics_ may unfortunately change reload
  int k = basics_[row_k_.num];
  si_->getBasics(basics_);
  
  
  
  
  if(basics_[row_k_.num] != k)//scan basics to find k again
  {
    for(int ii = 0 ; ii < numrows_ ; ii++)
    {
      if(basics_[ii]==k)
	    {
	      row_k_.num=ii;
	      break;
	    }
    }
  }
  
#ifdef OLD_COMPUTATION
  pullTableauRow(row_k_);
  row_k_.rhs =  row_k_.rhs - floor(row_k_.rhs);
  
#else
#ifdef LandP_DEBUG
  {
    //Check that row_k_ is the same has what it should be
    CglLandPSimplex::TabRow myrow;
    int rowsize = numcols_ + numrows_;
    myrow.size(rowsize + 1);
    
    myrow.num = row_k_.num;
    pullTableauRow(myrow);
    myrow.rhs = myrow.rhs - floor(myrow.rhs);
    CoinRelFltEq eq;
    for(int ii = 0 ; ii < rowsize ; ii++)
    {
      if(!perturbed_ || !perturbed_[ii])
        assert(eq(row_k_.row[ii],myrow.row[ii]));
      else//coefficient has been perturbed get unperturbed
      {
        double coef= 0; //it should be 0
        assert(eq(coef,myrow.row[ii]));
      }
    }
    assert(eq(row_k_.rhs,myrow.rhs));
  }
#endif
#endif
#ifdef LandP_DEBUG
#if LandP_DEBUG > 2
	if(!code)
  {
    //Check that row_k_ is what we had computed in newRow_
    CoinRelFltEq eq;
    
    if(!eq(row_k_.rhs,debug_.bestNewRhs_))
    {
      double error = fabs(row_k_.rhs - debug_.bestNewRhs_);
      handler_->message(WarnBadRhsComputation,  messages_)<<error<<CoinMessageEol<<CoinMessageEol;
      if(error > 1e-06)
        throw -1;
    }
    double norm_inf = 0.;
    for(int i = 0 ; i < numcols_ + numrows_ ; i++)
    {
      if(!eq(row_k_.row[i],debug_.bestNewRow_[i]))
      {
        norm_inf = max(fabs(row_k_.row[i]-debug_.bestNewRow_[i]),norm_inf);
        std::cerr<<"Variable" <<i<<", coefficient : "<<row_k_.row[i]<<" should have been "<<debug_.bestNewRow_[i]<<std::endl;
        std::cerr<<"lo bound "<<getLoBound(i)<<" up bound "<< getUpBound(i)<<std::endl;
      }
    }
    
    if(norm_inf > 1e-06)
    {
      handler_->message(WarnBadRowComputation,  messages_)<<norm_inf<<CoinMessageEol<<CoinMessageEol;
      throw -1;
    }
  }
#endif
#endif
  if(!code)
    return true;
  else return false;
}

/** Find a row which can be used to perform an improving pivot return index of the cut or -1 if none exists
* (i.e., find the leaving variable).*/
int 
CglLandPSimplex::findCutImprovingPivotRow( int &direction, int &gammaSign, double tolerance)
{
  bool bestRed = 0;
  tolerance = -10*tolerance;
  int bestRow = -1;
  int bestDirection = 0;
  int bestGamma = 0;
  double infty = si_->getInfinity();
  for(row_i_.num = 0 ; row_i_.num < numrows_ ; row_i_.num++)
  {
    if(row_i_.num != row_k_.num//obviously not necessary to combine row k with itself
       && rowFlags_[row_i_.num] //row has not been flaged
                                //   && fabs(getUpBound(basics_[row_i_.num]) - getLoBound(basics_[row_i_.num]))>1e-09 //variable is not fixed
       ) 
    {
      // 	  int outVar = 227;
      // 	  for(int i = 0 ; i < numrows_ ; i++)
      // 	    {
      // 	      if(outVar == basics_[i])
      // 		row_i_.num = i;
      // 	    }
      pullTableauRow(row_i_);
      computeRedCostConstantsInRow();
      
      if(getLoBound(basics_[row_i_.num]) > -infty)
        // variable can leave at its lower bound	    
        //Compute reduced cost with basics_[i] value decreasing
	    {
	      direction = -1;
        
	      gammaSign = -1;
	      double redCost = computeCglpRedCost(direction, gammaSign);
	      if(redCost<tolerance)
        {
          if(bestRed)
          {
            tolerance = redCost;
            bestRow = row_i_.num;
            bestDirection = direction;
            bestGamma = gammaSign;
          }
          else return row_i_.num;
        }
	      gammaSign = 1;
	      redCost = computeCglpRedCost(direction, gammaSign);
	      if(redCost<tolerance)
        {
          if(bestRed)
          {
            tolerance = redCost;
            bestRow = row_i_.num;
            bestDirection = direction;
            bestGamma = gammaSign;
          }
          else return row_i_.num;
        }
	    }
      if( getUpBound(basics_[row_i_.num])<infty) // variable can leave at its upper bound	    
                                                 //Compute reduced cost with basics_[i] value decreasing
	    {
	      direction = 1;
	      // 	      adjustTableauRow(i_, row_i_, rhs_i_, direction);
	      gammaSign = -1;	
	      double redCost = computeCglpRedCost(direction, gammaSign);
	      if(redCost<tolerance)
        {
          if(bestRed)
          {
            tolerance = redCost;
            bestRow = row_i_.num;
            bestDirection = direction;
            bestGamma = gammaSign;
          }
          else return row_i_.num;
        }
	      gammaSign = 1;
	      redCost = computeCglpRedCost(direction, gammaSign);
		      if(redCost<tolerance)
          {
            if(bestRed)
            {
              tolerance = redCost;
              bestRow = row_i_.num;
              bestDirection = direction;
              bestGamma = gammaSign;
            }
            else return row_i_.num;
          }
	    }
      rowFlags_[row_i_.num]=false;
    }
  }
  direction = bestDirection;
  gammaSign = bestGamma;
  row_i_.num=bestRow;
  if(row_i_.num>=0 && row_i_.num!= bestRow)
  {
    row_i_.num=bestRow;
    pullTableauRow(row_i_);
  }
  return bestRow;
}


/** Find a row which can be used to perform an improving pivot the fast way
* (i.e., find the leaving variable).
\return index of the cut or -1 if none exists. */
int 
CglLandPSimplex::fastFindCutImprovingPivotRow( int &direction, int &gammaSign, double tolerance)
{
  bool modularize = false;
  //Fill vector to compute contribution to reduced cost of variables in M1 and M2 (nz non-basic vars in row_k_.row).
  // 1. Put the values
  // 2. Post multiply by basis inverse
  //   for(int i = 0 ; i < numrows_ ; i++)
  //     rWk1_[i]=0;
  double * rWk1bis_ =NULL;
  CoinFillN(rWk1_,numrows_,(double) 0);
  if(modularize)
    CoinFillN(rWk1bis_, numrows_, (double) 0);
  int capacity = 0;
  const CoinPackedMatrix* mat = si_->getMatrixByCol();
  
  const CoinBigIndex* starts = mat->getVectorStarts();
  const CoinBigIndex * lengths = mat->getVectorLengths();
  const int * indices = mat->getIndices();
  const double * elements = mat->getElements();
  
  for(int i = 0 ;  i < numcols_ && inM1_[i]!= -1; i++)
  {
    if(inM1_[i]<numcols_)
    {
      const CoinBigIndex& begin = starts[inM1_[i]];
      const CoinBigIndex end = begin + lengths[inM1_[i]];
      bool swap = false;
      if(basis_.getStructStatus(inM1_[i])==CoinWarmStartBasis::atUpperBound) swap = true;
      for(CoinBigIndex k = begin ; k < end ;  k++)
      {
        if(swap)
        {
          rWk1_[indices[k]] += elements[k] * sigma_;
          if(modularize)
            rWk1bis_[indices[k]] += elements[k] * (colsolToCut_[inM1_[i]] - sigma_);
        }
        else
        {
          rWk1_[indices[k]] -= elements[k] * sigma_;
          if(modularize)
            rWk1bis_[indices[k]] -= elements[k] * (colsolToCut_[inM1_[i]] - sigma_);
        }
      }	
    }
    else
    {
      bool swap = false;
      if(basis_.getArtifStatus(inM1_[i] - numcols_)==CoinWarmStartBasis::atUpperBound) swap = true;
      if(swap)
      {
        rWk1_[inM1_[i] - numcols_] += sigma_;
        if(modularize)
          rWk1bis_[inM2_[i] - numcols_] += (colsolToCut_[inM1_[i]] - sigma_);
      }
      else
      {
        rWk1_[inM1_[i] - numcols_] -= sigma_;
        if(modularize)
          rWk1bis_[inM2_[i] - numcols_] -= (colsolToCut_[inM1_[i]] - sigma_);
      }
    }
  }
  for(int i = 0 ;  i < numcols_ && inM2_[i]!= -1; i++)
  {
    if(inM2_[i]<numcols_)
    {
      const CoinBigIndex& begin = starts[inM2_[i]];
      const CoinBigIndex end = begin + lengths[inM2_[i]];
      bool swap = false;
      if(basis_.getStructStatus(inM2_[i])==CoinWarmStartBasis::atUpperBound) swap = true;
      for(CoinBigIndex k = begin ; k < end ;  k++)
      {
        if(swap)
        {
          rWk1_[indices[k]] += elements[k] * (colsolToCut_[inM2_[i]] - sigma_);
          if(modularize)
            rWk1bis_[indices[k]] += elements[k] * sigma_;
        }
        else
        {
          rWk1_[indices[k]] -= elements[k] * (colsolToCut_[inM2_[i]] - sigma_);
          if(modularize)
            rWk1bis_[indices[k]] -= elements[k] * sigma_;
        }
      }
    }
    else
    {
      bool swap = false;
      if(basis_.getArtifStatus(inM2_[i] - numcols_)==CoinWarmStartBasis::atUpperBound) swap = true;
      if(swap)
      {
        rWk1_[inM2_[i] - numcols_] += (colsolToCut_[inM2_[i]] - sigma_);
        if(modularize)
          rWk1bis_[inM2_[i] - numcols_] += sigma_;
      }
      else
      {
        rWk1_[inM2_[i] - numcols_] -= (colsolToCut_[inM2_[i]] - sigma_);
        if(modularize)
          rWk1bis_[inM2_[i] - numcols_] -= sigma_;
      }
    }
  }
  
  for(int i = 0 ; i < numrows_ ; i++)
  {
    if(rWk1_[i])
      rIntWork_[capacity++] = i;
  }
  //  CoinFillN(&rIntWork_[capacity],numrows_-capacity,0);
  CoinIndexedVector indexed;
  indexed.borrowVector(numrows_, capacity, rIntWork_, rWk1_);

  OsiClpSolverInterface * clpSi = dynamic_cast<OsiClpSolverInterface *>(si_);
  if(clpSi)
    clpSi->getBInvACol(&indexed);
  else
    throw CoinError("Function not implemented in this OsiSolverInterface",
                    "getBInvACol","CglLandpSimplex");
  indexed.returnVector();
  if(modularize)
  {
    capacity = 0;
    for(int i = 0 ; i < numrows_ ; i++)
    {
      if(rWk1bis_[i])
        rIntWork_[capacity++] = i;
    }
    //  CoinFillN(&rIntWork_[capacity],numrows_-capacity,0);
    CoinIndexedVector indexed;
    indexed.borrowVector(numrows_, capacity, rIntWork_, rWk1bis_);
    OsiClpSolverInterface * clpSi = dynamic_cast<OsiClpSolverInterface *>(si_);
    if(clpSi)
      clpSi->getBInvACol(&indexed);
    else
      indexed.returnVector();
  }
  //Now compute the contribution of the variables in M3_
  //Need to get the column of the tableau in rW3_ for each of these and 
  // add up with correctly in storage for multiplier for negative gamma (named rW3_) and 
  // for positive gamma (which is named rW4_)
  if(inM3_[0]!=-1)
  {
    CoinFillN(rWk3_,numrows_,0.);
    CoinFillN(rWk4_,numrows_,0.);
    if(modularize)
    {
      double * rWk3bis_ = NULL;
      double * rWk4bis_ = NULL;
      CoinFillN(rWk3bis_,numrows_,0.);    
      CoinFillN(rWk4bis_,numrows_,0.);
    }
  }
  for(int i = 0 ; i < numcols_ && inM3_[i]!=-1 ;i++)
  {
    si_->getBInvACol(inM3_[i], rWk2_);
    bool swap = false;
    if(inM3_[i] < numcols_ && basis_.getStructStatus(inM3_[i])==CoinWarmStartBasis::atUpperBound) swap = true;
    if(inM3_[i]>= numcols_ && basis_.getArtifStatus(inM3_[i] - numcols_)==CoinWarmStartBasis::atUpperBound) swap = true;
    
    for(int j = 0 ; j < numrows_ ; j++)
    {
#if LandP_DEBUG >1
      if(fabs(rWk2_[j])>1e20)
        throw -1;
#endif
      if(swap)
        rWk2_[j] = - rWk2_[j];
      if(rWk2_[j] > 0.)
      {
        //is in M1 for multiplier with negative gamma
        rWk3_[j] -= sigma_*rWk2_[j];
        //is in M2 for multiplier with positive gamma
        rWk4_[j] -= (colsolToCut_[inM3_[i]] - sigma_)*rWk2_[j];
        
#if LandP_DEBUG >1
        if(fabs(rWk3_[j])>1e20)
          throw -1;
        if(fabs(rWk4_[j])>1e20)
          throw -1;
#endif      
        
      }
      else if (rWk2_[j] < 0.)
      {
        //is in M2 for multiplier with negative gamma
        rWk3_[j] -=(colsolToCut_[inM3_[i]] - sigma_)*rWk2_[j];
        
        //is in M1 for multiplier with positive gamma
        rWk4_[j] -= sigma_*rWk2_[j];
        
#if LandP_DEBUG >1
        if(fabs(rWk3_[j])>1e20)
          throw -1;
        if(fabs(rWk4_[j])>1e20)
          throw -1;
#endif
      }
    }
    
  }
  //Now, we just need to add up everything correctly for each of the reduced 
  //cost. Compute the Tau in rWk2_ which is not used anymore then compute the reduced cost in rWk1_ for u^l_j 
  // rwk2_ in u^u_j rWk3_ for v^u_j rWk4_ for v^u_j 
  // Let's rename not to get too much confused
  double * ul_i = rWk1_;
  double * uu_i = rWk2_;
  double * vl_i = rWk3_;
  double * vu_i = rWk4_;
  int bestRow = -1;
  int bestDirection;
  int bestGammaSign;
  // double infty = si_->getInfinity();
  
  
  nNegativeRc_ = 0;//counter
    nNegativeRcRows_ = 0;//counter
      
      double fzero = colsolToCut_[basics_[row_k_.num]] - floor(colsolToCut_[basics_[row_k_.num]]);
      fzero = row_k_.rhs;
      for(int i = 0 ; i < numcols_ ; i++)
      {
        fzero -= colsolToCut_[nonBasics_[i]] * row_k_.row[nonBasics_[i]];
      }
      
      double bestReducedCost = -tolerance;
      for(int i = 0 ; i < numrows_ ; i++)
      {
        if(i == row_k_.num//obviously not necessary to combine row k with itself
                          //   && fabs(getUpBound(basics_[row_i_.num]) - getLoBound(basics_[row_i_.num]))>1e-09 //variable is not fixed
           )
        {
          ul_i[i]=uu_i[i]=vl_i[i]=vu_i[i]=10.;
          continue;
        }
#if LandP_DEBUG > 3 //compute r.c. with old function and double check (VERY SLOW!!!)
        row_i_.num = i;
        pullTableauRow(row_i_);
        computeRedCostConstantsInRow();
        if(fabs(rWk1_[i]-tau_)/( (fabs(rWk1_[i])>1) ? fabs(rWk1_[i]) : 1.)> RED_COST_CHECK)
        {
          std::cout<<"constant : "<<rWk1_[i]<<", "<<tau_<<std::endl;
          THROW_ON_BAD_REDUCED_COST;
        }
#endif
        
        double tau1 = rWk1_[i];
        double tau2 = rWk1_[i];
        double tau3 = rWk1_[i];
        double tau4 = rWk1_[i];
        if(inM3_[0]!=-1)
        {
        tau1 += rWk3_[i];
        tau2 += rWk4_[i];
        tau3 += rWk4_[i];
        tau4 += rWk3_[i];
        }
        if(modularize)
        {
          tau1 = rWk1_[i] + rWk3_[i];
          tau2 = rWk1bis_[i] + rWk3_[i];
          tau3 = - rWk1_[i] - rWk4_[i];
          tau4 = - rWk1bis_[i] - rWk3_[i];
        }
        
        
        double redCost;
        bool hasNegativeRc = false;
        
        double loBound = getLoBound(basics_[i]);
        
        if(loBound > -1e50)
        {      
          redCost = -sigma_ + (tau1)
          + (1 - fzero) * ( colsol_[basics_[i]]   
                            - loBound);
          
          
#if LandP_DEBUG > 3
          double alternate = computeCglpRedCost(-1, -1);
          if(fabs(alternate-redCost)/((fabs(alternate) > 1) ? fabs(alternate) : 1.)> RED_COST_CHECK)
          {
            std::cout<<"ul_i Old function : "<<alternate
            <<" new one : "<<redCost<<std::endl;
            alternate = computeCglpRedCost(-1, -1);
            THROW_ON_BAD_REDUCED_COST;
          }
#endif    
          
          if(redCost < -tolerance)
          {
            ul_i[i] = redCost;
            nNegativeRc_++;
            hasNegativeRc = true;
          }
          else
            ul_i[i] = 10.;
          if(redCost < bestReducedCost 
             && rowFlags_[i] )//row has not been flaged
          {
            bestDirection = -1;
            bestGammaSign = -1;
            bestReducedCost = redCost;
            bestRow = i;
          }
          
          
          redCost = -sigma_ - (tau2)
            - (1 - fzero) * ( colsol_[basics_[i]]   
                              - loBound)
            - loBound + colsolToCut_[basics_[i]];
          
#if LandP_DEBUG > 3
          alternate = computeCglpRedCost(-1, 1);
          if(fabs(alternate-redCost)/((fabs(alternate) > 1) ? fabs(alternate) : 1.)> RED_COST_CHECK)
          {
            std::cout<<"vl_i Old function : "<<alternate
            <<" new one : "<<redCost<<std::endl;
            alternate = computeCglpRedCost(-1, 1);
            THROW_ON_BAD_REDUCED_COST;
          }
#endif
          if(redCost < -tolerance)
          {
            vl_i[i] = redCost;
            nNegativeRc_++;
            hasNegativeRc = true;
          }
          else
            vl_i[i]=10.;
          
          if(redCost < bestReducedCost
             && rowFlags_[i]) //row has not been flaged
          {
            bestDirection = -1;
            bestGammaSign = 1;
            bestReducedCost = redCost;
            bestRow = i;
          }
          
          
          
        }
        else
        {
          ul_i[i] = 10.;
          vl_i[i] = 10.;
        }
        double upBound = getUpBound(basics_[i]);
        if(getUpBound(basics_[i]) < 1e50)
        {
          double redCost = - sigma_ - (tau3)
          + (1 - fzero) * ( - colsol_[basics_[i]]   
                            + upBound);
          
#if LandP_DEBUG > 3
          double alternate = computeCglpRedCost(1, -1);
          if(fabs(alternate-redCost)/((fabs(alternate) > 1) ? fabs(alternate) : 1.)> RED_COST_CHECK)
          {
            std::cout<<"uu_i Old function : "<<alternate
            <<" new one : "<<redCost<<std::endl;
            THROW_ON_BAD_REDUCED_COST;
          }
#endif      
          if(redCost < -tolerance)
          {
            nNegativeRc_++;
            uu_i[i] = redCost;
            hasNegativeRc = true;
          }
          else
            uu_i[i] = 10.;
          
          if(redCost < bestReducedCost
             && rowFlags_[i]) //row has not been flaged
          {
            bestDirection = 1;
            bestGammaSign = -1;
            bestReducedCost = redCost;
            bestRow = i;
          }
          
          redCost = -sigma_ + (tau4)
            - (1 - fzero) * ( - colsol_[basics_[i]]
                              + upBound)
            + upBound - colsolToCut_[basics_[i]];
          
#if LandP_DEBUG > 3
          alternate = computeCglpRedCost(1, 1);
          if(fabs(alternate-redCost)/((fabs(alternate) > 1) ? fabs(alternate) : 1.)>1e-7)
          {
            std::cout<<"vu_i Old function : "<<alternate
            <<" new one : "<<redCost<<std::endl;
            THROW_ON_BAD_REDUCED_COST;
          }
#endif
          if(redCost < -tolerance)
          {
            vu_i[i] = redCost;
            nNegativeRc_++;
            hasNegativeRc = true;
          }
          
          
          else
            vu_i[i] = 10.;
          
          if(redCost < bestReducedCost
             && rowFlags_[i]) //row has not been flaged
          {
            bestDirection = 1;
            bestGammaSign = 1;
            bestReducedCost = redCost;
            bestRow = i;
          }     
        }
        else
        {
          uu_i[i] = 10.;
          vu_i[i] = 10.;
        }
        if(hasNegativeRc) nNegativeRcRows_ ++;
      }
      //  throw -1;
      direction = bestDirection;
      gammaSign = bestGammaSign;
      //  row_i_.num=bestRow;
      if(bestRow != -1)
      {
        chosenReducedCostVal_ = bestReducedCost;
        row_i_.num=bestRow;
        pullTableauRow(row_i_);
      }
      return bestRow;
}


/** Find a row which can be used to perform an improving pivot tables are already filled.
\return index of the cut or -1 if none exists. */
int 
CglLandPSimplex::rescanReducedCosts( int &direction, int &gammaSign, double tolerance)
{
  // The reduced cost are already here in rWk1_ is u^l_j 
  // rwk2_ is u^u_j rWk3_ is v^u_j rWk4_ is v^u_j 
  // Let's rename not to get too much confused
  double * ul_i = rWk1_;
  double * uu_i = rWk2_;
  double * vl_i = rWk3_;
  double * vu_i = rWk4_;
  int bestRow = -1;
  int bestDirection;
  int bestGammaSign;
  // double infty = si_->getInfinity();
  double bestReducedCost = -tolerance;
  for(int i = 0 ; i < numrows_ ; i++)
  {
    if(i == row_k_.num//obviously not necessary to combine row k with itself
       || !rowFlags_[i] //row has not been flaged
                        //   && fabs(getUpBound(basics_[row_i_.num]) - getLoBound(basics_[row_i_.num]))>1e-09 //variable is not fixed
       )
      continue;
    
    if(ul_i[i] < bestReducedCost
       && rowFlags_[i]) //row has not been flaged
    {
      bestDirection = -1;
      bestGammaSign = -1;
      bestReducedCost = ul_i[i];
      bestRow = i;
    }
    
    if(vl_i[i] < bestReducedCost
       && rowFlags_[i]) //row has not been flaged
    {
      bestDirection = -1;
      bestGammaSign = 1;
      bestReducedCost = vl_i[i];
      bestRow = i;
    }
    
    
    if(uu_i[i] < bestReducedCost
       && rowFlags_[i]) //row has not been flaged
    {
      bestDirection = 1;
      bestGammaSign = -1;
      bestReducedCost = uu_i[i];
      bestRow = i;
    }
    if(vu_i[i] < bestReducedCost
       && rowFlags_[i]) //row has not been flaged
    {
      bestDirection = 1;
      bestGammaSign = 1;
      bestReducedCost = vu_i[i];
      bestRow = i;
    }
    
  }
  //  throw -1;
  direction = bestDirection;
  gammaSign = bestGammaSign;
  //  row_i_.num=bestRow;
  if(bestRow != -1)
  {
    chosenReducedCostVal_ = bestReducedCost;
    row_i_.num=bestRow;
    pullTableauRow(row_i_);
  }
  return bestRow;
}

/** Find the column which leads to the best cut (i.e., find incoming variable).*/
int 
CglLandPSimplex::fastFindBestPivotColumn(int direction, int gammaSign,
                                         double pivotTol, double rhsTol, bool reducedSpace, bool allowDegenerate, double & bestSigma)
{
  gammas_.clear();
  
#if 0 //when have to do a deep debuging might want to stop for a specific row
  if(basics_[row_i_.num] == 343)
    std::cout<<"There"<<std::endl;
#endif
  adjustTableauRow(basics_[row_i_.num], row_i_, direction);
  double fzero = colsolToCut_[basics_[row_k_.num]] - floor(colsolToCut_[basics_[row_k_.num]]);
  fzero = row_k_.rhs;
  for(int i = 0 ; i < numcols_ ; i++)
  {
    fzero -= colsolToCut_[nonBasics_[i]] * row_k_.row[nonBasics_[i]];
  }
  double p = -row_k_.rhs * (1 - fzero);
  double q = row_i_.rhs * fzero;
  
  if(gammaSign < 0)
    q -= row_i_.rhs;
  double r = 1.;
  double s = (double) gammaSign;
  
#ifdef LandP_DEBUG
  // One can recompute the reduced cost based on the f+ computation
  // It is usefull for double checking correctness of computations
  // but, since in this procedure we only compute for pivot > gammaTolerance (non degenerate pivots)
  // the test can fail.
  
  // Those are used to recompute the true reduced cost.*/
  double pTrue = p;
  double qTrue = q;
  double rTrue = r;
  double sTrue = s;
#endif
  double diff = 0;
  
  bool haveSmallGammaPivot = false;
  double gammaTolerance = 0;
  if(allowDegenerate)
    gammaTolerance = 0;
  //fill the array with the gammas of correct sign
  for(int i = 0 ; i < numcols_ ; i++)
  {
    if(colCandidateToLeave_[i]==false) continue;
    double gamma = 1;
    if(fabs(row_i_.row[nonBasics_[i]]) > pivotTol)
    {
      gamma = - row_k_.row[nonBasics_[i]]/row_i_.row[nonBasics_[i]];
      if(gamma * gammaSign >= gammaTolerance)
	    {
	      gammas_.insert(i,gamma*gammaSign);
	    }
      
    }
    
    gamma = fabs(gamma); //  we already know the sign of gamma, its absolute value is more usefull
    if(row_k_.row[nonBasics_[i]]>0. && gamma >= gammaTolerance )
    {
      if(gammaSign > 0)
	    {
	      p += row_k_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	      
	    }
      else
	    {
        //	      if(fabs(colsolToCut_[nonBasics_[i]]) > 0)
      {
        p += row_k_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
        q += row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
      }
	    }
      r += row_k_.row[nonBasics_[i]];
      s += row_i_.row[nonBasics_[i]];
    }
    else if (row_k_.row[nonBasics_[i]]<0 && gamma >= gammaTolerance)
    {
      if(gammaSign > 0)
	    {
	      q -= row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	    }
      r -= row_k_.row[nonBasics_[i]];
      s -= row_i_.row[nonBasics_[i]];
    }
    else
    {
      haveSmallGammaPivot |= true;
      if(gammaSign > 0 && row_i_.row[nonBasics_[i]] < 0)
	    {
	      q -= row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	    }
      else if (gammaSign < 0 && row_i_.row[nonBasics_[i]] < 0)
	    {
	      q += row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	    }
      
      s += gammaSign*fabs(row_i_.row[nonBasics_[i]]);     
      diff += gammaSign*fabs(row_i_.row[nonBasics_[i]]);
    }
#ifdef LandP_DEBUG
    //fill the true p0 q0 r0 s0
    rTrue += fabs(row_k_.row[nonBasics_[i]]);
    if(row_k_.row[nonBasics_[i]]>0.)
    {
      if(gammaSign > 0)
	    {
	      pTrue += row_k_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	    }
      else 
	    {
	      pTrue += row_k_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	      qTrue += row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	    }
      sTrue += row_i_.row[nonBasics_[i]];
    }
    else if (row_k_.row[nonBasics_[i]]<0)
    {
      if(gammaSign > 0)
	    {
	      qTrue -= row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	    }
      sTrue -= row_i_.row[nonBasics_[i]];
    }
    else
    {
      if(gammaSign > 0 && row_i_.row[nonBasics_[i]] < 0)
	    {
	      qTrue -= row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
	    }
      else if (gammaSign < 0)
	    {
	      if (row_i_.row[nonBasics_[i]] < 0)
        {
          qTrue += row_i_.row[nonBasics_[i]] * colsolToCut_[nonBasics_[i]];
        }
	      
	    }
      sTrue += gammaSign*fabs(row_i_.row[nonBasics_[i]]);
    }
#endif
  }
  
#ifdef LandP_DEBUG
  if(!debug_.eq(pTrue/rTrue, sigma_) && !debug_.req(pTrue/rTrue, sigma_))
  {
    std::cout<<" Row k :";
    row_k_.print(std::cout, 9 , nonBasics_,numcols_);
    double sigma2 = computeCglpObjective(row_k_);
    std::cerr<<"Error in initial sigma HERE"<<pTrue/rTrue<<"   "<<sigma_<<"    "<<sigma2<<"          "<<p/r - sigma_<<std::endl;
    std::cout<<"f0 is "<<fzero<<", colsol : "<<colsolToCut_[basics_[row_k_.num]]<<std::endl;
    std::cerr<<"p :"<<p<<"  r :"<<r<<std::endl;
    throw -1;
  }
  
  double reducedCostTrue = gammaSign * ( qTrue -  (pTrue*sTrue)/rTrue);
  double reducedCost = gammaSign * ( q -  (p*s)/r);
  if(!debug_.eq(reducedCostTrue, chosenReducedCostVal_) && !debug_.req(reducedCostTrue, chosenReducedCostVal_))
  {
    std::cerr<<"Error in initial values "<<reducedCost <<"   "<<reducedCostTrue <<"   "
    <<chosenReducedCostVal_<<"  "<<reducedCost - chosenReducedCostVal_<<std::endl;
    
    std::cout<<" Row k :";
    row_k_.print(std::cout, 9 , nonBasics_,numcols_);
    
    std::cout<<" Row i :";
    row_i_.print(std::cout, 9, nonBasics_, numcols_);
    //#ifdef THINGS_ARE_REALLY_WRONG  
    fastFindCutImprovingPivotRow(direction, gammaSign, pivotTol);
    resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);	
    //#endif
    throw -1;
  }
  
#if LandP_DEBUG > 2
         newRow_.row[basics_[row_k_.num]] = 1.;
         newRow_.rhs = row_k_.rhs;
     double sigma2 = computeCglpObjective(0, false);
     if(!debug_.eq(p/r, sigma2))
       {
         std::cerr<<"Error in initial values"<<std::endl;
         throw -1;
       }
#endif
#endif
  
  int n = gammas_.getNumElements();
  gammas_.sortIncrElement();
  const int* inds = gammas_.getIndices();
  const double * elements = gammas_.getElements();
  int bestColumn = -1;
  double newSigma = 1e100;
  bestSigma = sigma_ = p/r;
  int lastValid = -1;
  bool sigmaWorst = false;
  for(int i = 0 ; i < n ; i++)
  {
    if(elements[i]< gammaTolerance)
    {
      std::cout<<"Something is very wrong in new procedure, in the loop there is a gamma below the tolerance"<<std::endl;
    }
    if( (row_k_.rhs + gammaSign * elements[i] * row_i_.rhs> (1 - rhsTol *gammaTolerance) ||
         row_k_.rhs + gammaSign * elements[i] * row_i_.rhs< (rhsTol)*gammaTolerance)) //getting to close to integer feasibility (is that really bad????)
    {
#ifdef LandP_DEBUG 
      std::cerr<<"Getting to close to integer feasibility get last valid entering variable"<<std::endl;
#endif
      //break;
    }
    if(fabs(row_i_.row[nonBasics_[inds[i]]])< pivotTol)
    {
      std::cout<<"Pivot below tolerance, should we break??"<<std::endl;
      //   break;
    }
    else
    {
      newSigma = (p + gammaSign * elements[i] * q)/(r + gammaSign*elements[i] * s);
      if(newSigma > bestSigma + pivotTol)
      {
        break;
      }
      else if(newSigma <= bestSigma && colCandidateToLeave_[inds[i]])
      {
        bestColumn = inds[i];
        bestSigma = newSigma;
        lastValid = i;
      }
      else
      {
        sigmaWorst=true;
      }
#if LandP_DEBUG > 2
      newRow_.row[basics_[row_k_.num]] = 1.;
      newRow_.rhs = row_k_.rhs + gammaSign * elements[i] * row_i_.rhs;
      debug_.newSigma_ = computeCglpObjective(gammaSign * elements[i], false);
             if(!debug_.eq(debug_.newSigma_, newSigma) && !debug_.req(debug_.newSigma_, newSigma_))//Ok we really disagree
             {
               std::cout<<"Disagreement in leaving column computation : "<<std::endl
               <<"Old method says new objective will be "<<debug_.newSigma_<<std::endl
               <<"New method says new objective will be "<<newSigma<<std::endl
       		 <<"Delta : "<<debug_.newSigma_-newSigma<<std::endl;
               throw -1;
              }
#endif
    }
    if(gammaSign*(q * r - p * s) >= 0)/* function is starting to increase stop here*/
    {
      break;
    }
    
    int col = nonBasics_[inds[i]];
    if(row_i_.row[col] > 0)
    {
      if(gammaSign > 0)
      {
        p += row_k_.row[col] * colsolToCut_[col];
        q += row_i_.row[col] * colsolToCut_[col];
        r += row_k_.row[col]*2;
        s += row_i_.row[col]*2;
      }
      else
      {
        p -= row_k_.row[col] * colsolToCut_[col];
        q -= row_i_.row[col] * colsolToCut_[col];
        r -= row_k_.row[col]*2;
        s -= row_i_.row[col]*2;
      }
    }
    else if(row_i_.row[col]<0)
    {
      if(gammaSign > 0)
      {
        p -= row_k_.row[col] * colsolToCut_[col];
        q -= row_i_.row[col] * colsolToCut_[col];
        r -= row_k_.row[col]*2;
        s -= row_i_.row[col]*2;
      }
      else
      {
        p += row_k_.row[col] * colsolToCut_[col];
        q += row_i_.row[col] * colsolToCut_[col];
        r += row_k_.row[col]*2;
        s += row_i_.row[col]*2;
      }
    }
  }
#if LandP_DEBUG > 2
  //Do we do slow full check
#if LandP_DEBUG > 10 //Desactivated better to check all the computations as they are done
  
  //Check with older and slower method that we still get the same result.
  int oldEntering = findBestPivotColumn(direction, gammaTolerance, pivotTol, reducedSpace, allowDegenerate);
  if(oldEntering != bestColumn)//oh oh disagreement!
  {
    //There may be several pivots with best depth
    if(fabs(debug_.newSigma_ - bestSigma)>1e-08)//Ok we really disagree
    {
      std::cout<<"Disagreement in leaving column computation : "<<std::endl
      <<"Old method, entering variable : "<<oldEntering<<" for a new objective of "<<debug_.newSigma_<<std::endl
      <<"New method, entering variable : "<<bestColumn<<" for a new objective of "<<bestSigma<<std::endl;
      //      throw -1;
    }
  }
#else
  //Just do a fast check that the computed sigma is correct and store the next row for later pivot validity check
  if(lastValid!=-1)
  {
    CoinFillN(newRow_.row,numcols_ + numrows_ + 1, 0.);
    newRow_.row[basics_[row_k_.num]] = 1.;
    newRow_.rhs = row_k_.rhs + gammaSign * elements[lastValid] * row_i_.rhs;
    debug_.newSigma_ = computeCglpObjective(gammaSign * elements[lastValid], false);
          if(!debug_.eq(debug_.newSigma_ , bestSigma))//Ok we really disagree
     	{
     	  std::cout<<"Disagreement in leaving column computation : "<<std::endl
     		   <<"Old method says new objective will be "<<debug_.newSigma_<<std::endl
     		   <<"New method says new objective will be "<<bestSigma<<std::endl;
     	  //	  	  throw -1;
     	}
  }
  
  debug_.bestNewRhs_ = newRow_.rhs;
  CoinCopyN(newRow_.row, numcols_ + numrows_, debug_.bestNewRow_);
#endif
#endif
  resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
  //Get the results did we find a valid pivot? is it degenerate?
  if(bestColumn == -1)// Apparently no pivot is within the tolerances
  {
    handler_->message(WarnFailedPivotTol, messages_)<<CoinMessageEol<<CoinMessageEol;
    return -1;
  }
  if(bestSigma < sigma_ - 1e-07)//everything has gone ok
  {
    handler_->message(FoundBestImprovingCol, messages_)<<nonBasics_[bestColumn]<<gammaSign * elements[lastValid]<<bestSigma<<CoinMessageEol<<CoinMessageEol;
    inDegenerateSequence_ = false;
    return bestColumn;
  }
  if(bestSigma > sigma_)//Oups pivot makes the cut worse
  {
    handler_->message(WarnFailedBestImprovingCol, messages_)<<chosenReducedCostVal_<<sigma_<<bestSigma<<CoinMessageEol<<CoinMessageEol;
    throw -1;
    return -1;
  }
  else if(allowDegenerate)//Pivot is degenerate and we allow
  {
    inDegenerateSequence_ = true;
    return bestColumn;
  }
  else//we don't accept a degenerate pivot
  {    
    handler_->message(WarnFailedBestImprovingCol, messages_)<<sigma_<<bestSigma<<CoinMessageEol<<CoinMessageEol;
    return -1;
  }
}


struct point{
double x;
double y;
point(double a, double b):x(a), y(b)
{}
bool operator<(const point & other){return x < other.x;}
};



int 
CglLandPSimplex::findBestPivotColumn(int direction, 
                                     double pivotTol, bool reducedSpace, bool allowDegenerate, bool modularize)
{
  CoinFillN(newRow_.row,numcols_ + numrows_ + 1, 0.);
  int varOut=-1;
  
  adjustTableauRow(basics_[row_i_.num], row_i_, direction);
  
  double m = si_->getInfinity();
  
  int j = 0;
  double gamma = 0.;
  
  bool pivotTolFail = false;//Failed to find a pivot because gamma tolerance
    bool rhsIIFail = true; // Failed to find a pivot because of integer infeasibility
    for(;j< numcols_ ; j++)
    {
      if(reducedSpace && 
         !colCandidateToLeave_[j]
         )
        continue;
      if(fabs(row_i_.row[nonBasics_[j]])< pivotTol)
      {
        continue;
      }
      gamma = - row_k_.row[nonBasics_[j]]/row_i_.row[nonBasics_[j]];

      newRow_.row[basics_[row_k_.num]] = 1.;
      newRow_.rhs = row_k_.rhs + gamma * row_i_.rhs;
      if(newRow_.rhs > 1e-5  && newRow_.rhs < 1 - 1e-5 )
      {
        rhsIIFail = false;
        double m_j = computeCglpObjective(gamma, modularize);
        if(m_j < m)
        {
          varOut = j;
          m = m_j;
#ifdef LandP_DEBUG
          debug_.bestNewRhs_ = newRow_.rhs;
          CoinCopyN(newRow_.row, numcols_ + numrows_, debug_.bestNewRow_);
          debug_.newSigma_ = m_j;
#endif
        }
      }
    }
    resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
    
    if(m < sigma_ )
    {
      handler_->message(FoundBestImprovingCol, messages_)<<nonBasics_[varOut]<<gamma<<m<<CoinMessageEol<<CoinMessageEol;
      inDegenerateSequence_ = false;
      return varOut;
    }
    else if(allowDegenerate && m<=sigma_)
    {
      inDegenerateSequence_ = true;
    }
    else
    {
      //      std::cout<<" leaving candidate : "<<row_i_.num<<", "<<basics_[row_i_.num]<<std::endl;
      if(pivotTolFail)
        handler_->message(WarnFailedPivotTol, messages_)<<CoinMessageEol<<CoinMessageEol;
      else if(rhsIIFail)
        handler_->message(WarnFailedPivotIIf, messages_)<<row_k_.rhs<<CoinMessageEol<<CoinMessageEol;
      else
        handler_->message(WarnFailedBestImprovingCol, messages_)<<chosenReducedCostVal_<<sigma_<<m<<CoinMessageEol<<CoinMessageEol;
      return -1;
    }
    return -1;
}

struct reducedCost
{
  /** To avoid computing two times the same row direction will have strange
  values, direction is -1 or 1 if for only one of the two direction
  rc is <0 and is -2 or 2 if the two direction have one <0 rc with the sign
  indicating which one of the two directions is the best.<br>
  Note that by theory only one reduced cost (for u_i, or v_i) 
  maybe negative for each direction.
  */
  int direction;
  /** gammSign is the sign of gamma (corresponding to taking rc for u_i or v_i)
  for the best of the two rc for this row.*/
  int gammaSign;
  /** gammaSign2 is the sign of gamma for the worst of the two rc for this row.*/
  int gammaSign2;
  /** if both reduced costs are <0 value is the smallest of the two.*/
  double value;
  /** greatest of the two reduced costs */
  double value2;
  /** index of the row.*/
  int row;
  bool operator<(const reducedCost & other)
  { return (value>other.value);}
};

/** Find incoming and leaving variables which lead to the most violated
adjacent normalized lift-and-project cut. 
\remark At this point reduced costs should be already computed.
\return incoming variable variable,
\param leaving variable
\param direction leaving direction
\param numTryRows number rows tried
\param pivotTol pivot tolerance
\param reducedSpace separaration space (reduced or full)
\param allowNonStrictlyImproving wether or not to allow non stricly improving pivots.      
*/
int CglLandPSimplex::findBestPivot(int &leaving, int & direction, 
                                   const CglLandP::Parameters & params){
  // 1. Sort <0 reduced costs in increasing order
  // 2. for numTryRows reduced costs call findBestPivotColumn
  // if better record
  
  // The reduced cost are already here in rWk1_ is u^l_j 
  // rwk2_ is u^u_j rWk3_ is v^u_j rWk4_ is v^u_j 
  // Let's rename not to get too much confused
  double * ul_i = rWk1_;
  double * uu_i = rWk2_;
  double * vl_i = rWk3_;
  double * vu_i = rWk4_;
  
  reducedCost * rc = new reducedCost[nNegativeRcRows_];
  int k = 0;
  rc[k].direction = 0;//initialize first rc
    int k2 = 0;
    for(int i = 0 ; i < numrows_ ; i++)
    {
      if(ul_i[i] < -params.pivotTol)
        //       && rowFlags_[i]) //row has not been flaged
      {
        rc[k].direction = -1;
        rc[k].gammaSign = -1;
        rc[k].value = ul_i[i];
        rc[k].row = i;
        k2++;
      }
      
      if(vl_i[i] < -params.pivotTol)
        //&& rowFlags_[i]) //row has not been flaged
      {
        rc[k].direction = -1;
        rc[k].gammaSign = 1;
        rc[k].value = vl_i[i];
        rc[k].row = i;
        k2++;
      }
      
      
      if(uu_i[i] < -params.pivotTol)
        //       && rowFlags_[i]) //row has not been flaged
      {
        if(rc[k].direction == 0)
        {
          rc[k].direction = 1;
          rc[k].gammaSign = -1;         
          rc[k].value = uu_i[i];
          rc[k].row = i;
        }
        else
        {
          if(uu_i[i] <rc[k].value)//this one is better
          {
            rc[k].direction = 2;
            rc[k].gammaSign2 = rc[k].gammaSign;
            rc[k].gammaSign = -1;
            rc[k].value2 = rc[k].value;
            rc[k].value = uu_i[i];
          }
          else
          {
            rc[k].direction = -2;
            rc[k].gammaSign2 = -1;
            rc[k].value2 = uu_i[i];
          }
        }
        k2++;
      }
      if(vu_i[i] < -params.pivotTol)
        //&& rowFlags_[i]) //row has not been flaged
      {
        if(rc[k].direction==0)
        {
          rc[k].direction = 1;
          rc[k].gammaSign = 1;
          rc[k].value = vu_i[i];
          rc[k].row = i;
        }
        else
        {
          if(vu_i[i] < rc[k].value)//this one is better
          {
            rc[k].direction = 2;
            rc[k].gammaSign2 = rc[k].gammaSign;
            rc[k].gammaSign = 1;
            rc[k].value2 = rc[k].value;
            rc[k].value = vu_i[i];
          }
          else //the other one is better
          {
            rc[k].direction = -2;
            rc[k].gammaSign2 = 1;
            rc[k].value2 = vu_i[i];
          }
        }
        k2++;
      }
      if(rc[k].direction!=0) //We have added a row with < 0 rc during
                             //last iteration
      {
        k++;
        if(k<nNegativeRcRows_)
          rc[k].direction = 0;
        else
          break;
      }
      
  }
    
    assert(k2==nNegativeRc_);
    assert(k==nNegativeRcRows_);
    
    //now make a heap
    std::make_heap(rc, rc + k);
    //  assert(rc[0].value==chosenReducedCostVal_);  
    int bestLeaving = -1;
    int bestIncoming = -1;
    int bestDirection = 0;
    int bestGammaSign = 0;
    
    double bestSigma = DBL_MAX;
    double bestRc;
    //now scan the heap
    int best_l = 0;
    int notImproved = 0;
    for(int l = 0; l < k && l < 10	; l++, notImproved++)
    {
      if(!rowFlags_[rc[l].row]) continue;//this row has been marked to be skipped
                                         //     if(bestLeaving != -1 && rc[l].value > -1e-02) break;
      if(rc[l].value > -1e-02) break;
      row_i_.num=rc[l].row;
      pullTableauRow(row_i_);//Get the tableau row
        
        //compute f+ or f- for the best negative rc corresponding to this row
        chosenReducedCostVal_ = rc[l].value;
        double sigma;
        int incoming = 
          fastFindBestPivotColumn
          (rc[l].direction, rc[l].gammaSign, 
           params.pivotTol, params.away, 
           (params.reducedSpace!=0), 
           0, 
           sigma);
        if(incoming!=-1 && bestSigma > sigma)
        {
//          std::cout<<"I found a better pivot "<<sigma - sigma_<< " for indice number "<<l<<std::endl;
          best_l = l;
          bestSigma = sigma;
          bestIncoming = incoming;
          bestLeaving = rc[l].row;
          bestDirection = rc[l].direction > 0 ? 1 : -1;
          bestGammaSign = rc[l].gammaSign;
          bestRc = rc[l].value;
          notImproved = 0;
        }
        
        //Now evenutally compute f+ or f- for the other negative rc (if if exists)
        if(rc[l].direction == 2 || rc[l].direction == -2)
        {
          rc[l].direction/= -2;//Reverse the direction
          chosenReducedCostVal_ = rc[l].value2;//need to set this for debug double
                                               //checks
            int incoming = fastFindBestPivotColumn
              (rc[l].direction, rc[l].gammaSign2, 
               params.pivotTol, params.away, 
               (params.reducedSpace!=0), 
               0,
               sigma);
            if(incoming!=-1 && bestSigma > sigma)
            {
              // std::cout<<"I found a better pivot "<<sigma - sigma_<<std::endl;
              best_l = l;
              bestSigma = sigma;
              bestIncoming = incoming;
              bestLeaving = rc[l].row;
              bestDirection = rc[l].direction;
              bestGammaSign = rc[l].gammaSign2;
              bestRc = rc[l].value2;
              notImproved = 0;
            }
        }
    }
    
    row_i_.num = leaving = bestLeaving;
    chosenReducedCostVal_ = bestRc;
    assert(best_l <= nNegativeRcRows_);
    if(bestLeaving!=-1)
    {
      //      std::cout<<"Best pivot pivot "<<best_l<<std::endl;
      pullTableauRow(row_i_);
      extra.nNegativeRcRows += nNegativeRcRows_;
#if LandP_DEBUG > 2
      //Recall findBestPivotColumn to have good debug information
      double sigma=DBL_MAX;
      assert(bestIncoming == fastFindBestPivotColumn
             (bestDirection, bestGammaSign, 
              params.gammaTol,
              params.pivotTol, params.rhsTol, 
              params.reducedSpace, 
              params.allowNonStrictlyImproving, 
              sigma));
      assert(debug_.eq(sigma, bestSigma) || debug_.req(sigma, bestSigma));
#endif
      
      extra.bestRow += (best_l+1);
      extra.maxBestRow = max (extra.maxBestRow, (best_l+1));
      extra.bestRc += chosenReducedCostVal_;
      extra.maxRc = max(extra.maxRc, chosenReducedCostVal_);
    }
    direction = bestDirection;
    delete [] rc;
    return bestIncoming;
    
}
double 
CglLandPSimplex::computeCglpObjective(TabRow &row) const
{
  double numerator = -row.rhs * (1 - row.rhs);
  double denominator = 1;
  
  for(int j = 0 ; j < numcols_ ; j++)
  {
    denominator += fabs(row.row[ nonBasics_[j]]);
    numerator += (row.row[ nonBasics_[j]] > 0. ? 
                  row.row[ nonBasics_[j]] *(1- row.rhs):
                  - row.row[ nonBasics_[j]] * row.rhs)* colsolToCut_[nonBasics_[j]];
  }
  return numerator/denominator;
}


double 
CglLandPSimplex::computeCglpObjective(double gamma, bool strengthen)
{
  double numerator = -newRow_.rhs * (1 - newRow_.rhs);
  double denominator = 1;
  
  newRow_.row[ basics_[row_i_.num]] = newRowCoefficient(basics_[row_i_.num], gamma);
  if(strengthen && row_i_.num < numcols_ && si_->isInteger(row_i_.num))
    newRow_.row[ basics_[row_i_.num]] = modularizedCoef(newRow_.row[ basics_[row_i_.num]], newRow_.rhs);
  denominator += fabs(newRow_.row[ basics_[row_i_.num]]);
  numerator += (newRow_.row[ basics_[row_i_.num]] > 0 ? 
                newRow_.row[ basics_[row_i_.num]] *(1- newRow_.rhs):
                - newRow_.row[ basics_[row_i_.num]] * newRow_.rhs)* 
    colsolToCut_[basics_[row_i_.num]];
  
  for(int j = 0 ; j < numcols_ ; j++)
  {
    
    newRow_.row[nonBasics_[j]] = newRowCoefficient(nonBasics_[j], gamma);
    if(strengthen && nonBasics_[j] < numcols_ &&  si_->isInteger(j))
      newRow_.row[ nonBasics_[j]] = modularizedCoef(newRow_.row[ nonBasics_[j]], newRow_.rhs);
//    if(colCandidateToLeave_[nonBasics_[j]] == false) continue;
    denominator += fabs(newRow_.row[ nonBasics_[j]]);
    numerator += (newRow_.row[ nonBasics_[j]] > 0 ? 
                  newRow_.row[ nonBasics_[j]] *(1- newRow_.rhs):
                  - newRow_.row[ nonBasics_[j]] * newRow_.rhs)* 
      colsolToCut_[nonBasics_[j]];
  }
  //  assert (fabs(numerator/denominator - computeCglpObjective(newRow_))<1e-04);
  return numerator/denominator;
}

/** Modularize row.*/
void 
CglLandPSimplex::modularizeRow(TabRow & row)
{
  for(int i = 0 ; i < numcols_ ; i++)
  {
    int &ni = nonBasics_[i];
    if(ni >= numcols_) continue;
    if(si_->isInteger(ni))
      row.row[ni] = modularizedCoef(row.row[ni],row.rhs);
  }
}


/** Compute the reduced cost of Cglp */
double 
CglLandPSimplex::computeCglpRedCost(int direction, int gammaSign)
{
  double toBound;
  toBound = direction == -1 ? getLoBound(basics_[row_i_.num])  : getUpBound(basics_[row_i_.num]);
  
  double value =0;
  int sign = gammaSign * direction;
  double tau1 = 0;
  double tau2 = 0;
  for(int i = 0 ; i < numcols_ && inM3_[i]!= -1 ; i++)
  {
    tau1 += fabs(row_i_ .row[inM3_[i]]);
    if (sign == 1 && row_i_.row[inM3_[i]] < 0)
    {
      tau2 += row_i_.row[inM3_[i]] * colsolToCut_ [inM3_[i]];
    }
    else if(sign == -1 && row_i_.row[inM3_[i]] > 0)
    {
      tau2 += row_i_.row[inM3_[i]] * colsolToCut_ [inM3_[i]];
    }
  }
  double Tau = - sign * (tau_ + tau2) - tau1 * sigma_;
  value = - sigma_ + Tau
    + (1 - colsolToCut_[basics_[row_k_.num]]) * sign * (row_i_.rhs   
                                                        -  toBound)
    + (gammaSign == 1)*direction*(toBound - colsolToCut_[basics_[row_i_.num]]);
  
  return value;
}


/** Compute the value of sigma and thau (which are constants for a row i as defined in Mike Perregaard thesis */
void 
CglLandPSimplex::computeRedCostConstantsInRow()
{
  double tau1 = 0; //the part which will be multiplied by sigma
  double tau2 = 0;//the rest
    
    for(int i = 0 ; i < numcols_ && inM1_[i] != -1 ; i++)
    {
      tau1 += row_i_.row[inM1_[i]];      
    }
    for(int i = 0 ; i < numcols_ && inM2_[i] != -1 ; i++)
    {      
      tau1 -= row_i_.row[inM2_[i]];
      tau2 += row_i_.row[inM2_[i]] * colsolToCut_[inM2_[i]];
    }
    tau_ = sigma_ * tau1 + tau2;
    
}

void 
CglLandPSimplex::updateM1_M2_M3(TabRow & row, double tolerance, bool reducedSpace, bool perturb)
{
  tolerance = 0;
  int nM1 = 0;
  int nM2 = 0;
  int nM3 = 0;
  for(int i = 0; i<numcols_ ; i++)
  {
    if(row.row[nonBasics_[i]]< -tolerance) 
    {
      if(colHasToComputeContrib_[nonBasics_[i]])
      {
        inM1_[nM1++] = nonBasics_[i];
        colCandidateToLeave_[i]=1;
      }
      else
        colCandidateToLeave_[i]=0;
    }
    else if(row.row[nonBasics_[i]]>tolerance)
    {
      if(colHasToComputeContrib_[nonBasics_[i]])
      {
        inM2_[nM2++] = nonBasics_[i];
        colCandidateToLeave_[i]=1;
      }
      else
        colCandidateToLeave_[i]=0;
    }
    else
    {
      if(colHasToComputeContrib_[nonBasics_[i]])
      {
        if(perturb)//assign to M1 or M2 at random
        {
          int sign = CoinDrand48() > 0.5 ? 1 : -1; // NOT on a thread by thread basis
          if(sign == -1)//put into M1
          {
            inM1_[nM1++] = nonBasics_[i];
            colCandidateToLeave_[i]=1;
          }
          else//put into M2
          {
            inM2_[nM2++] = nonBasics_[i];
            colCandidateToLeave_[i]=1;
          }
        }
        else 
        {
          inM3_[nM3++] = nonBasics_[i];
          colCandidateToLeave_[i] = 1;
        }
      }
      else
        colCandidateToLeave_[i] = 0;
    }
   }
  if(nM1<numcols_) inM1_[nM1]=-1;
  if(nM2<numcols_) inM2_[nM2]=-1;
  if(nM3<numcols_) inM3_[nM3]=-1;
}

/** Create the intersection cut of row k*/
void 
CglLandPSimplex::createIntersectionCut(const TabRow & row, OsiRowCut &cut) const
{
  const double * colLower = si_->getColLower();
  const double * rowLower = si_->getRowLower();
  const double * colUpper = si_->getColUpper();
  const double * rowUpper = si_->getRowUpper();
  // double f_0 = row.rhs;
  //put the row back into original form
  for(int j = 0; j < numcols_ ; j++)
  {
    if((nonBasics_[j] < numcols_))
    {
      CoinWarmStartBasis::Status status =
	    (nonBasics_[j] < numcols_) ? 
	    basis_.getStructStatus(nonBasics_[j]) :
	    basis_.getArtifStatus(nonBasics_[j] - numcols_);
      
      if(status==CoinWarmStartBasis::atLowerBound)
	    {
	      //	      row.rhs += getLoBound(nonBasics_[j]) * row.row[nonBasics_[j]];
	    }
      else if(status==CoinWarmStartBasis::atUpperBound)
	    {
	      row.row[nonBasics_[j]] = - row.row[nonBasics_[j]];
	      //	      row.rhs += getUpBound(nonBasics_[j]) * row.row[nonBasics_[j]];
	    }
      else
	    {
	      throw;
	    }
    }
  }
  
  
  
  //  return ;
  
  cut.setUb(DBL_MAX);
  double * vec = new double[numcols_+ numrows_ ];
  CoinFillN(vec, numcols_ + numrows_, 0.);
  double infty = si_->getInfinity();
  double cutRhs = row.rhs;
  cutRhs = cutRhs ;//* (1 - cutRhs);
    for(int j = 0; j < numcols_ ; j++)
    {
      if(fabs(row.row[nonBasics_[j]])>1e-10)
      {
        double value = intersectionCutCoef(row.row[nonBasics_[j]], row.rhs);
        
        if(nonBasics_[j]<numcols_)
        {
          CoinWarmStartBasis::Status status = //CoinWarmStartBasis::basic; 
          basis_.getStructStatus(nonBasics_[j]);
	      	if(status==CoinWarmStartBasis::atUpperBound)
	      	  {
            value = - intersectionCutCoef(- row.row[nonBasics_[j]], row.rhs) ;
	      	    cutRhs += value * colUpper[nonBasics_[j]];
	      	  }
          else
            cutRhs += value * colLower[nonBasics_[j]];
          vec[nonBasics_[j]] += value;
        }
        else if(nonBasics_[j]>=numcols_)
        {
          int iRow = nonBasics_[j] - numcols_;
          
          if(rowLower[iRow] > -infty)
          {
            value = -value;
            cutRhs -= value*rowLower[iRow];	
            assert(basis_.getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound ||
                   (fabs(rowLower[iRow] - rowUpper[iRow]) < 1e-08));
          }
          else
          {
            cutRhs -= value*rowUpper[iRow];
            assert(basis_.getArtifStatus(iRow)==CoinWarmStartBasis::atLowerBound);
          }
          vec[nonBasics_[j]] = value;
          assert(fabs(cutRhs)<1e100);
        }
      }
    }

    const CoinPackedMatrix * mat = si_->getMatrixByCol();
    const CoinBigIndex * starts = mat->getVectorStarts();
    const int * lengths = mat->getVectorLengths();
    const double * values = mat->getElements();
    const CoinBigIndex * indices = mat->getIndices();
    for(int j = 0 ; j < numcols_ ; j++)
    {
      const int& start = starts[j];
      int end = start + lengths[j];
          for(int k = start ; k < end ; k++)
          {
            vec[j] -= vec[numcols_ + indices[k]] * values[k];	      
          }
      
    }
    
    //Pack vec into the cut
    int * inds = new int [numcols_];
    int nelem = 0;
    for(int i = 0 ; i < numcols_ ; i++)
    {
     if(fabs(vec[i]) > COIN_INDEXED_TINY_ELEMENT)
     {
       vec[nelem] = vec[i];
       inds[nelem++] = i;
     }
    }
      
    cut.setLb(cutRhs);
    cut.setRow(nelem, inds, vec, false);
    delete [] vec;
    
}

/** Create MIG cut from row k*/
void 
CglLandPSimplex::createMIG( TabRow &row, OsiRowCut &cut) const
{
  const double * colLower = si_->getColLower();
  const double * rowLower = si_->getRowLower();
  const double * colUpper = si_->getColUpper();
  const double * rowUpper = si_->getRowUpper();
  
  if(1)
  {
    
    //Clean up row
    //       for(int i = 0 ; i < numcols_ + numrows_ ; i++)
    // 	{
    // 	  if(fabs(row.row[i])<1e-6)
    // 	    row.row[i] = 0.;
    // 	}
    //       row.print(numcols_ + numrows_);
    double f_0 = row.rhs - floor(row.rhs);
    //put the row back into original form
    for(int j = 0; j < numcols_ ; j++)
    {
      if((nonBasics_[j] < numcols_))
	    {
	      CoinWarmStartBasis::Status status =
        (nonBasics_[j] < numcols_) ? 
        basis_.getStructStatus(nonBasics_[j]) :
        basis_.getArtifStatus(nonBasics_[j] - numcols_);
	      
	      if(status==CoinWarmStartBasis::atLowerBound)
        {
          //	      row.rhs += getLoBound(nonBasics_[j]) * row.row[nonBasics_[j]];
        }
	      else if(status==CoinWarmStartBasis::atUpperBound)
        {
          row.row[nonBasics_[j]] = - row.row[nonBasics_[j]];
          //	      row.rhs += getUpBound(nonBasics_[j]) * row.row[nonBasics_[j]];
        }
	      else
        {
          throw;
        }
	    }
    }
    
    row.rhs = f_0;
    
    cut.setUb(DBL_MAX);
    double * vec = new double[numcols_ + numrows_];
    CoinFillN(vec, numcols_ + numrows_, 0.);
    //f_0 = row.rhs - floor(row.rhs);
    double infty = si_->getInfinity();
    double cutRhs = row.rhs - floor(row.rhs);
    cutRhs = cutRhs * (1 - cutRhs);
    assert(fabs(cutRhs)<1e100);
    for(int j = 0; j < numcols_ ; j++)
    {
      int iRow = nonBasics_[j] - numcols_;
      if(fabs(row.row[nonBasics_[j]])>0.)
        //	 &&! (iRow >= 0 && (rowLower[iRow] - rowUpper[iRow]) > -1e-08))
      {
        double value = strengthenedIntersectionCutCoef(nonBasics_[j], row.row[nonBasics_[j]], row.rhs);
        if(fabs(value) > 1e9)
        {
          std::cout<<"Big coefficient:" <<value<<std::endl;
        }
        //	  if(fabs(value)>0.)
        {
          if(nonBasics_[j]<numcols_)
          {
            CoinWarmStartBasis::Status status = //CoinWarmStartBasis::basic; 
            basis_.getStructStatus(nonBasics_[j]);
            if(status==CoinWarmStartBasis::atUpperBound)
	      	  {
              value = - strengthenedIntersectionCutCoef(nonBasics_[j], - row.row[nonBasics_[j]], row.rhs) ;
	      	    cutRhs += value * colUpper[nonBasics_[j]];
	      	  }
            else
              cutRhs += value * colLower[nonBasics_[j]];
            assert(fabs(cutRhs)<1e100);
            vec[nonBasics_[j]] += value;
          }
          else if(nonBasics_[j]>=numcols_)
          {
            
            if(rowLower[iRow] > -infty)
            {
              value = -value;
              cutRhs -= value*rowLower[iRow];	
              assert(basis_.getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound ||
                     (rowUpper[iRow] < infty));
            }
            else
            {
              cutRhs -= value*rowUpper[iRow];
              //		  assert(basis_.getArtifStatus(iRow)==CoinWarmStartBasis::atLowerBound);
            }
            vec[nonBasics_[j]] = value;
            assert(fabs(cutRhs)<1e100);
          }
        }
      }
      }
      
    //Eliminate slacks
      const CoinPackedMatrix * mat = si_->getMatrixByCol();
      const CoinBigIndex * starts = mat->getVectorStarts();
      const int * lengths = mat->getVectorLengths();
      const double * values = mat->getElements();
      const CoinBigIndex * indices = mat->getIndices();
      const double * vecSlacks = vec + numcols_;
      for(int j = 0 ; j < numcols_ ; j++)
      {
        const CoinBigIndex& start = starts[j];
        CoinBigIndex end = start + lengths[j];
        double & val = vec[j];
        for(CoinBigIndex k = start ; k < end ; k++)
        {
          val -= vecSlacks[indices[k]] * values[k];	      
        }
      }
    
    //Pack vec into the cut
    int * inds = new int [numcols_];
    int nelem = 0;
    for(int i = 0 ; i < numcols_ ; i++)
    {
     if(fabs(vec[i]) > COIN_INDEXED_TINY_ELEMENT)
     {
       vec[nelem] = vec[i];
       inds[nelem++] = i;
     }
    }
      
    cut.setLb(cutRhs);
    cut.setRow(nelem, inds, vec, false);
    delete [] vec;
    delete [] inds;
    
    }
  
}

double 
CglLandPSimplex::normCoef(TabRow &row)
{
  double res = 1;
  for(int i = 0 ; i < numcols_ ; i++)
    res += fabs(row_k_.row[nonBasics_[i]]);
  return res/(1-row.rhs);
}


/** scale the cut passed as argument using provided normalization factor*/
void 
CglLandPSimplex::scale(OsiRowCut &cut, double norma)
{
  assert(norma >0.);
  CoinPackedVector row;
  row.reserve(cut.row().getNumElements());
  for(int i = 0 ; i < cut.row().getNumElements() ; i++)
  {
    row.insert(cut.row().getIndices()[i], cut.row().getElements()[i]/norma);
  }
  cut.setLb(cut.lb()/norma);
  cut.setRow(row);
}

/** scale the cut passed as argument*/
void 
CglLandPSimplex::scale(OsiRowCut &cut)
{
  double rhs = fabs(cut.lb());
  CoinPackedVector row;
  row.reserve(cut.row().getNumElements());
  for(int i = 0 ; i < cut.row().getNumElements() ; i++)
  {
    row.insert(cut.row().getIndices()[i], cut.row().getElements()[i]/rhs);
  }
  cut.setLb(cut.lb()/rhs);
  cut.setRow(row);
}

/** Get the row i of the tableau */
void 
CglLandPSimplex::pullTableauRow(TabRow &row) const
{
  const double * rowLower = si_->getRowLower();
  const double * rowUpper = si_->getRowUpper();
  
  CoinFillN(row.row, numcols_ + numrows_ ,0.);

  double infty = si_->getInfinity();
  /* Get the row */
  si_->getBInvARow(row.num,row.row,&row.row[numcols_]);
//  row.row[basics_[row.num]]=1;
  /* get the rhs */
  if(solver_ == xpress)//rhs was computed in BInvARow
    row.rhs = row.row[numcols_ + numrows_];
  else
  {
    int iCol = basics_[row.num];
    if(iCol<numcols_)
      row.rhs = si_->getColSolution()[iCol];
    else // Osi does not give direct acces to the value of slacks
    {
      iCol -= numcols_;
      row.rhs = - si_->getRowActivity()[iCol];
      if(rowLower[iCol]> -infty)
	    {
	      row.rhs += rowLower[iCol];
	    }
      else
	    {
	      row.rhs+= rowUpper[iCol];
	    }
    }
  }
  //Now adjust the row of the tableau to reflect non-basic variables activity
  for(int j = 0; j < numcols_ ; j++)
  {
    if(nonBasics_[j]<numcols_)
    {
      if(basis_.getStructStatus(nonBasics_[j])==CoinWarmStartBasis::atLowerBound)
	    {
	      if( solver_ == xpress && fabs(loBounds_[nonBasics_[j]])>0)
        {
          if(loBounds_[nonBasics_[j]] > -infty)
            row.rhs -= row.row[nonBasics_[j]] * loBounds_[nonBasics_[j]];
          else
          {
            throw CoinError("Structural at lower bound while there is no upper bound","CglLandPSimplex","pullTableauRow");
          }
        }
	    }
      else if(basis_.getStructStatus(nonBasics_[j])==CoinWarmStartBasis::atUpperBound)
	    {
	      if(solver_ == xpress && fabs(upBounds_[nonBasics_[j]])>0)
        {
          if(upBounds_[nonBasics_[j]]<infty)
            row.rhs -= row.row[nonBasics_[j]] * upBounds_[nonBasics_[j]];
          else
          {
            throw CoinError("Structural at upper bound while there is no upper bound","CglLandPSimplex","pullTableauRow");
          }
        }
	      row.row[nonBasics_[j]] = -row.row[nonBasics_[j]];
	    }
      else
	    {
	      throw CoinError("Invalid basis","CglLandPSimplex","pullTableauRow");
	    }
    }
    else
    {
      int iRow = nonBasics_[j] - numcols_;

	    if(basis_.getArtifStatus(iRow)==CoinWarmStartBasis::atUpperBound)
	    {
	      row.row[nonBasics_[j]] = -row.row[nonBasics_[j]];
	    }
    }
  }
}
/** Adjust the row of the tableau to reflect leaving variable direction */
void 
CglLandPSimplex::adjustTableauRow(int var, TabRow & row, int direction)
{
  double bound;
  if(direction > 0)
  {
//    row.row[basics_[row.num]] = - row.row[basics_[row.num]];

    for(int j = 0 ; j < numcols_ ; j++)
      row.row[nonBasics_[j]] = - row.row[nonBasics_[j]];

    row.rhs = -row.rhs;
    bound = getUpBound(var);
    colsolToCut_[var] = - colsolToCut_[var];
  }
  else if(direction < 0)
  {
    bound = - getLoBound(var);
  }
  row.rhs += bound;
  //  assert(fabs(row.rhs)<1e100);
  colsolToCut_[var] += bound;
}


/** reset the tableau row after a call to adjustTableauRow */
void 
CglLandPSimplex::resetOriginalTableauRow(int var, TabRow & row, int direction)
{
  if (direction > 0){
    adjustTableauRow(var, row, direction);
  }
  else
  {
    double bound = getLoBound(var);
    row.rhs += bound;
    colsolToCut_[var] += bound;
  }
}



#if LandP_DEBUG > 1
/** Create MIG cut from row k and express it in the original non-basic space*/
void 
CglLandPSimplex::put_in_non_basic_init_space(OsiRowCut &cut)
{
  const CoinPackedVector & row = cut.row();
  double lb = cut.lb();
  cut.setUb(DBL_MAX);
  double * vec = new double[numcols_ + numrows_];
  CoinFillN(vec, numcols_ + numrows_, 0.);
  //f_0 = row.rhs - floor(row.rhs);
  
  const int * ind = row.getIndices();
  const double * vals = row.getElements();
  int n = row.getNumElements();
  for(int j = 0; j < n ; j++)
  {
    //Now re-express coefficient int terms of initial non-basics
    int iRow = -1;
    CoinWarmStartBasis::Status status =
	  (ind[j] < numcols_) ? 
	  debug_.initialBasis_->getStructStatus(ind[j]) :
	  debug_.initialBasis_->getArtifStatus(ind[j] - numcols_);
    if(status == CoinWarmStartBasis::atUpperBound)
	  {
	    //swap the coefficient and reflect in rhs
	    vec[ind[j]] -= vals[j];
	    lb -= vals[j]*getUpBound(ind[j]);
	  }
    else if (status == CoinWarmStartBasis::atLowerBound)
	  { 
	    vec[ind[j]] += vals[j];
	    lb +=  vals[j]*getLoBound(ind[j]);
	  }
    else
	  {
	    // nonBasics_[j] was basic in initial tableau
	    //find the row in which it was
	    for (int i = 0 ; i < numrows_ ; i++)
	      if(debug_.initialBasics_[i] == ind[j])
        {
          iRow = i;
          break;
        }
          if(iRow==-1)//could not find row something is wrong 
          {
            std::cerr<<"Something is wrong in basis"<<std::endl;
            throw -1;
          }
          const CoinPackedMatrix &mat = debug_.initialTableau_;
	    lb -= vals[j] * debug_.initialColsol_[ind[j]];
	    int start = mat.getVectorStarts()[iRow];
	    int end = start + mat.getVectorLengths()[iRow];
	    for(int k = start ; k < end ; k++)
      {
        vec[debug_.initialNonBasics_[mat.getIndices()[k]]] -= vals[j] * mat.getElements()[k];	      
      }
	  }      
  }
  
  
  //Pack vec into the cut
  CoinPackedVector cutRow;
  for(int i = 0 ; i < numcols_ +numrows_; i++)
  {
    if(fabs(vec[i]) > 1e-10)//COIN_INDEXED_TINY_ELEMENT)
      cutRow.insert(i, vec[i]);
    }
  delete [] vec;
  cut.setLb(lb);
  cut.setRow(cutRow);
  
  
}


/** Get the tableau corresponding to the current basis.*/
void 
CglLandPSimplex::DebugData::getCurrentTableau(OsiSolverInterface &si, CglLandPSimplex &lap)
{
  si.getBasics(initialBasics_);
  initialBasis_=dynamic_cast<CoinWarmStartBasis *>(si.getWarmStart());
  CoinCopyN(lap.basics_, lap.numrows_, initialBasics_);
  CoinCopyN(lap.nonBasics_, lap.numcols_, initialNonBasics_);
  CoinCopyN(lap.colsol_, lap.numcols_ + lap.numrows_, initialColsol_);
  //Now get the tableau
  TabRow row;
  row.size(lap.numcols_ + lap.numrows_);
  CoinBigIndex maxNnz = si.getMatrixByCol()->getNumElements() * 2 + 100;
  int * indices = new int[maxNnz];
  double * values = new double[maxNnz];
  CoinBigIndex * starts = new CoinBigIndex[lap.numrows_ + 1];
  CoinBigIndex * lengths = new CoinBigIndex[lap.numrows_];
  CoinBigIndex nnz = 0;
  for(row.num = 0 ; row.num < lap.numrows_ ; row.num++)
  {
    starts[row.num] = nnz;
    lengths[row.num] = 0;
    lap.pullTableauRow(row);
    for(int i = 0 ; i < lap.numcols_ ; i++)
    {
      if(fabs(row.row[lap.nonBasics_[i]])>1e-10)
      {
        values[nnz] = row.row[lap.nonBasics_[i]];
        indices[nnz++] = i;
        lengths[row.num]++;
        if(nnz>=maxNnz-1)
        {
          maxNnz = 2*maxNnz + 100;
          int * indices2 = new int[maxNnz];
          double * values2 = new double[maxNnz];
          CoinCopyN(indices, nnz, indices2);
          delete [] indices;
          indices = indices2;
          
          CoinCopyN(values, nnz, values2);
          delete [] values;
          values = values2;
        }
      }
    }
  }
  starts[lap.numrows_] = nnz;
  initialTableau_.assignMatrix(0, lap.numcols_, lap.numrows_, nnz, values, indices, starts, lengths);
}
#endif



#if 0
int 
CglLandPSimplex::plotCGLPobj(int direction, double gammaTolerance,
                                     double pivotTol, bool reducedSpace, bool allowDegenerate, bool modularize)
{
  row_i_.num=2;
  pullTableauRow(row_i_);
  std::cout<<"Row : "<<row_i_.num<<std::endl;
  CoinFillN(newRow_.row,numcols_ + numrows_ + 1, 0.);
  int varOut=-1;
  std::list< point > points; 
  adjustTableauRow(basics_[row_i_.num], row_i_, -1);
  
  double m = si_->getInfinity();
  
  int j = 0;
  double gamma = 0.;
  
  bool pivotTolFail = false;//Failed to find a pivot because gamma tolerance
    bool rhsIIFail = true; // Failed to find a pivot because of integer infeasibility
    for(;j< numcols_ ; j++)
    {
      if(reducedSpace && 
         !colCandidateToLeave_[j]
         )
        continue;
      if(fabs(row_i_.row[nonBasics_[j]])< pivotTol)
      {
        continue;
      }
      gamma = - row_k_.row[nonBasics_[j]]/row_i_.row[nonBasics_[j]];
      if(0 && fabs(gamma)< gammaTolerance) 
      {
        continue;
        pivotTolFail = true;
      }
      newRow_.row[basics_[row_k_.num]] = 1.;
      newRow_.rhs = row_k_.rhs + gamma * row_i_.rhs;
      if(newRow_.rhs > 1e-5 * gammaTolerance && newRow_.rhs < 1 - 1e-5 * gammaTolerance)
      {
        rhsIIFail = false;
        double m_j = computeCglpObjective(gamma, modularize);
        points.push_back(point(gamma, m_j));

      }
    }
    resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
   points.sort();
   {
   std::ofstream os("cglpObj_low");
   for(std::list< point >::iterator i = points.begin();
       i != points.end() ; i++)
   {
      os<<i->x<<"    "<<i->y<<std::endl;
   }
   }
   CoinFillN(newRow_.row,numcols_ + numrows_ + 1, 0.);
  varOut=-1;
  points.clear(); 
  adjustTableauRow(basics_[row_i_.num], row_i_, 1);
  
   m = si_->getInfinity();
  
  j = 0;
  gamma = 0.;
  
  pivotTolFail = false;//Failed to find a pivot because gamma tolerance
    rhsIIFail = true; // Failed to find a pivot because of integer infeasibility
    for(;j< numcols_ ; j++)
    {
      if(reducedSpace && 
         !colCandidateToLeave_[j]
         )
        continue;
      if(fabs(row_i_.row[nonBasics_[j]])< pivotTol)
      {
        continue;
      }
      gamma = - row_k_.row[nonBasics_[j]]/row_i_.row[nonBasics_[j]];
      if(0 && fabs(gamma)< gammaTolerance) 
      {
        continue;
        pivotTolFail = true;
      }
      newRow_.row[basics_[row_k_.num]] = 1.;
      newRow_.rhs = row_k_.rhs + gamma * row_i_.rhs;
      if(newRow_.rhs > 1e-5 * gammaTolerance && newRow_.rhs < 1 - 1e-5 * gammaTolerance)
      {
        rhsIIFail = false;
        double m_j = computeCglpObjective(gamma, modularize);
        points.push_back(point(gamma, m_j));

      }
    }
    resetOriginalTableauRow(basics_[row_i_.num], row_i_, direction);
{
   std::ofstream os("cglpObj_up");
   points.sort();
   for(std::list< point >::iterator i = points.begin();
       i != points.end() ; i++)
   {
      os.precision(10);
      os<<i->x<<"    "<<i->y<<std::endl;
   }
   }
}
#endif
}
