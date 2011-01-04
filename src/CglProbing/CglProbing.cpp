
/*
  Copyright (C) 2002, International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#define PROBING100 0

// #define PRINT_DEBUG
// #define CGL_DEBUG 1
// #undef NDEBUG

#include "CoinPackedMatrix.hpp"
#include "CglProbing.hpp"
#include "CglProbingRowCut.hpp"
#include "CglProbingDebug.hpp"

/*
  PROBING_EXTRA_STUFF has been false since June 2006, so anything guarded by
  this symbol is likely to be extremely dusty.
*/
//#define PROBING_EXTRA_STUFF true
#define PROBING_EXTRA_STUFF false
// Seems to be unused  -- lh, 101202 --
// #define FIXED_ALLOWANCE 10
#define SIZE_ROW_MULT 4
#define SIZE_ROW_ADD 2000

typedef struct {double infeasibility;int sequence;} double_int_pair;
class double_int_pair_compare {
public:
  bool operator() (double_int_pair x , double_int_pair y) const
  {
    return ( x.infeasibility < y.infeasibility);
  }
};


#if CGL_DEBUG > 0
static int nPath=0;
#endif

//-------------------------------------------------------------------
// Generate disaggregation cuts
//------------------------------------------------------------------- 
void CglProbing::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info2) const
{

# if CGL_DEBUG > 0
  std::cout << "Entering CglProbing::generateCuts." << std::endl ;
/*
  Note that the debugger must be activated before generateCuts is called.
  Even then, getRowCutDebugger will return a pointer only if we're on the
  optimal path.
*/
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  if (debugger) {
    std::cout << "On optimal path " << nPath << std::endl ;
    nPath++ ;
  }
# if CGL_DEBUG > 3
  dumpSoln(si) ;
# endif
# endif
/*
  Enforce the rule that a negative value prohibits row cut generation in the
  search tree. (Code 0x4 is `column cuts only'.)
*/
  int saveRowCuts = rowCuts_ ;
  if (rowCuts_ < 0) {
    if (info2.inTree)
      rowCuts_ = 0x4 ;
    else
      rowCuts_ = -rowCuts_ ;
  }
/*
  Setup. Create arrays that will hold row and column bounds while we're
  working.

  TODO: !rowCopy_ is used lots of places as surrogate for `we have no
	snapshot'. Snapshot creation is user-initiated; it won't just
	happen.   -- lh, 100924 --

  TODO: The distinction between colLower, colUpper, vs. colLower_, colUpper_,
	is that the latter can be accessed by the client after we return.
	See generateCutsAndModify.
*/
  int nRows = si.getNumRows() ; 
  double *rowLower = new double[nRows+1] ;
  double *rowUpper = new double[nRows+1] ;

  int nCols = si.getNumCols() ;
  if (!rowCopy_) {
    numberRows_ = nRows ;
    numberColumns_ = nCols ;
  }
  double *colLower = new double[nCols] ;
  double *colUpper = new double[nCols] ;
/*
  We need a modifiable copy of the CglTreeInfo object.
*/
  CglTreeInfo info = info2 ;
/*
  Do the work.
*/
  int ninfeas =
	gutsOfGenerateCuts(si,cs,rowLower,rowUpper,colLower,colUpper,&info) ;
/*
  Infeasibility is indicated by a stylized infeasible column cut.

  If we're debugging, do a sanity check that this cut really is infeasible.
*/
  if (ninfeas) {
    OsiRowCut rc ;
    rc.setLb(DBL_MAX) ;
    rc.setUb(0.0) ;   
    cs.insert(rc) ;
#   if CGL_DEBUG > 0
    const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
    if (debugger)
      assert(!debugger->invalidCut(rc)) ; 
#   endif
  }
/*
  Delete the working arrays and restore the row cut mode.
*/
  delete [] rowLower ;
  delete [] rowUpper ;
  delete [] colLower ;
  delete [] colUpper ;
  delete [] colLower_ ;
  delete [] colUpper_ ;
  colLower_ = NULL ;
  colUpper_ = NULL ;
  rowCuts_ = saveRowCuts ;

# if CGL_DEBUG > 0
  std::cout << "Leaving CglProbing::generateCuts." << std::endl ;
# endif

  return ;
}


/*
  TODO: Note that this is not a const method (contrast with generateCuts).
	-- lh, 100924
  
  TODO: Note that the CglTreeInfo object is not const here. In generateCuts,
	it's passed in as a const and a nonconst clone is created. I need to
	look into that.  -- lh, 100924 --

  Generally, check generateCuts for additional comment. I went through it first
  and didn't repeat comments on common code blocks here.
*/
int CglProbing::generateCutsAndModify (const OsiSolverInterface &si, 
				       OsiCuts &cs,
				       CglTreeInfo *info) 
{
# if CGL_DEBUG > 0
  std::cout << "Entering CglProbing::generateCutsAndModify." << std::endl ;
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  if (debugger) {
    std::cout << "On optimal path " << nPath << std::endl ;
    nPath++ ;
  }
# if CGL_DEBUG > 3
  dumpSoln(si) ;
# endif
# endif

/*
  Enforce the rule that a negative value prohibits row cut generation in the
  search tree. (Code 0x4 is `column cuts only'.)
*/
  int saveRowCuts = rowCuts_ ;
  if (rowCuts_ < 0) {
    if (info->inTree)
      rowCuts_ = 4 ;
    else
      rowCuts_ = -rowCuts_ ;
  }
/*
  Tweak the cut generation mode. Boiled down: If the user sets mode to 0, then
  at the root we'll ignore it for three passes, using mode 1 instead. Then
  we'll do one pass with mode 0, but permanently convert the mode to 1. If
  we're not at the root, mode 0 is automatically converted to mode 1 for cut
  generation, but we don't change the visible setting.

  TODO: Ok, specifying mode 0 is a fiction, except for pass 4 at the root.
	Why pass 4? And why are we restricted to just once?  -- lh, 100924 --
*/
  int saveMode = mode_ ;
  bool rowCliques = false ;
  if (!mode_) {
    if (info->pass != 4 || info->inTree) {
      mode_ = 1 ;
    } else {
      saveMode = 1 ;
      rowCliques = true ;
    }
  }
/*
  Create working arrays for row and column bounds.
*/
  int nRows = si.getNumRows() ; 
  double * rowLower = new double[nRows+1] ;
  double * rowUpper = new double[nRows+1] ;

  int nCols = si.getNumCols() ; 
  double * colLower = new double[nCols] ;
  double * colUpper = new double[nCols] ;
/*
  Do the work.
*/
  int ninfeas = gutsOfGenerateCuts(si,cs,
				   rowLower,rowUpper,colLower,colUpper,info) ;
/*
  Infeasibility is indicated by a stylized infeasible column cut.

  If we're debugging, check whether this is a plausible outcome.

  TODO: Remind me --- Why can't I use the debugger created at the top of the
	method. And don't I need to dispose of both of these?
*/
  if (ninfeas) {
    OsiRowCut rc ;
    rc.setLb(DBL_MAX) ;
    rc.setUb(0.0) ;   
    cs.insert(rc) ;
#   if CGL_DEBUG > 0
    const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
    if (debugger)
      assert(!debugger->invalidCut(rc)) ; 
#   endif
  }
/*
  Restore modes.
*/
  rowCuts_ = saveRowCuts ;
  mode_ = saveMode ;
/*
  jjf: move bounds so can be used by user

  Mode 3 is just like mode 2, except that the tightened row bounds will be
  returned to the client.
*/
  if (mode_ == 3) {
    delete [] rowLower_ ;
    delete [] rowUpper_ ;
    rowLower_ = rowLower ;
    rowUpper_ = rowUpper ;
  } else {
    delete [] rowLower ;
    delete [] rowUpper ;
  }
/*
  But given that this is generateCutsAndModify, we certainly want to return the
  tightened column bounds.
*/
  delete [] colLower_ ;
  delete [] colUpper_ ;
  colLower_ = colLower ;
  colUpper_ = colUpper ;
/*
  jjf: Setup information 

  Boil down detailed clique information into a form that's easily useable when
  calculating row bounds.
*/
  if (rowCliques && numberRows_ && numberColumns_)
    setupRowCliqueInformation(si) ;

# if CGL_DEBUG > 0
  std::cout
    << "Leaving CglProbing::generateCutsAndModify, return "
    << ninfeas << std::endl ;
# endif

  return ninfeas ;
}

bool analyze(const OsiSolverInterface * solverX, char * intVar,
             double * lower, double * upper)
{
  OsiSolverInterface * solver = solverX->clone();
  const double *objective = solver->getObjCoefficients() ;
  int numberColumns = solver->getNumCols() ;
  int numberRows = solver->getNumRows();
  double direction = solver->getObjSense();
  int iRow,iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix  matrixByCol(*solver->getMatrixByCol());
  const double * element = matrixByCol.getElements();
  const int * row = matrixByCol.getIndices();
  const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
  const int * columnLength = matrixByCol.getVectorLengths();

  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();

  char * ignore = new char [numberRows];
  int * which = new int[numberRows];
  double * changeRhs = new double[numberRows];
  memset(changeRhs,0,numberRows*sizeof(double));
  memset(ignore,0,numberRows);
  int numberChanged=0;
  bool finished=false;
  while (!finished) {
    int saveNumberChanged = numberChanged;
    for (iRow=0;iRow<numberRows;iRow++) {
      int numberContinuous=0;
      double value1=0.0,value2=0.0;
      bool allIntegerCoeff=true;
      double sumFixed=0.0;
      int jColumn1=-1,jColumn2=-1;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (upper[jColumn] > lower[jColumn]+1.0e-8) {
          if (!intVar[jColumn]) {
            if (numberContinuous==0) {
              jColumn1=jColumn;
              value1=value;
            } else {
              jColumn2=jColumn;
              value2=value;
            }
            numberContinuous++;
          } else {
            if (fabs(value-floor(value+0.5))>1.0e-12)
              allIntegerCoeff=false;
          }
        } else {
          sumFixed += lower[jColumn]*value;
        }
      }
      double low = rowLower[iRow];
      if (low>-1.0e20) {
        low -= sumFixed;
        if (fabs(low-floor(low+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      double up = rowUpper[iRow];
      if (up<1.0e20) {
        up -= sumFixed;
        if (fabs(up-floor(up+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      if (!allIntegerCoeff)
        continue; // can't do
      if (numberContinuous==1) {
        // see if really integer
        // This does not allow for complicated cases
        if (low==up) {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // integer
              numberChanged++;
              intVar[jColumn1]=77;
            }
          }
        } else {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // This constraint will not stop it being integer
              ignore[iRow]=1;
            }
          }
        }
      } else if (numberContinuous==2) {
        if (low==up) {
          /* need general theory - for now just look at 2 cases -
             1 - +- 1 one in column and just costs i.e. matching objective
             2 - +- 1 two in column but feeds into G/L row which will try and minimize
             (take out 2 for now - until fixed)
          */
          if (fabs(value1)==1.0&&value1*value2==-1.0&&!lower[jColumn1]
              &&!lower[jColumn2]&&columnLength[jColumn1]==1&&columnLength[jColumn2]==1) {
            int n=0;
            int i;
            double objChange=direction*(objective[jColumn1]+objective[jColumn2]);
            double bound = CoinMin(upper[jColumn1],upper[jColumn2]);
            bound = CoinMin(bound,1.0e20);
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                which[n++]=jRow;
                changeRhs[jRow]=value;
              }
            }
            for ( i=columnStart[jColumn2];i<columnStart[jColumn2]+columnLength[jColumn2];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                if (!changeRhs[jRow]) {
                  which[n++]=jRow;
                  changeRhs[jRow]=value;
                } else {
                  changeRhs[jRow]+=value;
                }
              }
            }
            if (objChange>=0.0) {
              // see if all rows OK
              bool good=true;
              for (i=0;i<n;i++) {
                int jRow = which[i];
                double value = changeRhs[jRow];
                if (value) {
                  value *= bound;
                  if (rowLength[jRow]==1) {
                    if (value>0.0) {
                      double rhs = rowLower[jRow];
                      if (rhs>0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    } else {
                      double rhs = rowUpper[jRow];
                      if (rhs<0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    }
                  } else if (rowLength[jRow]==2) {
                    if (value>0.0) {
                      if (rowLower[jRow]>-1.0e20)
                        good=false;
                    } else {
                      if (rowUpper[jRow]<1.0e20)
                        good=false;
                    }
                  } else {
                    good=false;
                  }
                }
              }
              if (good) {
                // both can be integer
                numberChanged++;
                intVar[jColumn1]=77;
                numberChanged++;
                intVar[jColumn2]=77;
              }
            }
            // clear
            for (i=0;i<n;i++) {
              changeRhs[which[i]]=0.0;
            }
          }
        }
      }
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8&&!intVar[iColumn]) {
        double value;
        value = upper[iColumn];
        if (value<1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        value = lower[iColumn];
        if (value>-1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        bool integer=true;
        for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          if (!ignore[iRow]) {
            integer=false;
            break;
          }
        }
        if (integer) {
          // integer
          numberChanged++;
          intVar[iColumn]=77;
        }
      }
    }
    finished = numberChanged==saveNumberChanged;
  }
  bool feasible=true;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (intVar[iColumn]==77) {
      if (upper[iColumn]>1.0e20) {
        upper[iColumn] = 1.0e20;
      } else {
        upper[iColumn] = floor(upper[iColumn]+1.0e-5);
      }
      if (lower[iColumn]<-1.0e20) {
        lower[iColumn] = -1.0e20;
      } else {
        lower[iColumn] = ceil(lower[iColumn]-1.0e-5);
        if (lower[iColumn]>upper[iColumn])
          feasible=false;
      }
      if (lower[iColumn]==0.0&&upper[iColumn]==1.0)
        intVar[iColumn]=1;
      else if (lower[iColumn]==upper[iColumn])
        intVar[iColumn]=0;
      else
        intVar[iColumn]=2;
    }
  }
  delete [] which;
  delete [] changeRhs;
  delete [] ignore;
  //if (numberChanged)
  //printf("%d variables could be made integer\n",numberChanged);
  delete solver;
  return feasible;
}


/*
  NOTE: It's assumed that the arrays passed in for rowLower and rowUpper
	have an extra slot that can be used for bounds on the objective,
	in the event that gutsOfGenerateCuts decides to put it in the matrix.
*/
int CglProbing::gutsOfGenerateCuts (const OsiSolverInterface &si, 
				    OsiCuts &cs ,
				    double *rowLower, double *rowUpper,
                                    double *colLower, double *colUpper,
                                    CglTreeInfo *info) const
{

# if CGL_DEBUG > 0
  std::cout
    << "Entering CglProbing::gutsOfGenerateCuts, mode "
    << mode_ << std::endl ;
# endif

  int nRows ;
  
  CoinPackedMatrix *rowCopy = NULL ;
  int numberRowCutsBefore = cs.sizeRowCuts() ;

/*
  Establish a cutoff -- nontrivial, if possible. The client can forbid use of
  the objective by setting usingObjective_ to -1.
*/
  double cutoff ;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff) ;
  if (!cutoff_available || usingObjective_ < 0) {
    cutoff = si.getInfinity() ;
  }
  cutoff *= si.getObjSense() ;
/*
  This looks like a sure failure if we're maximising and grabbed the solver's
  infinity.  In any event, why is large and positive ok but large and
  negative is not? --lh, 100917 --
*/
  if (fabs(cutoff) > 1.0e30)
    assert (cutoff > 1.0e30) ;

  int mode = mode_ ;
  int nCols = si.getNumCols(); 
/*
  Get a modifiable array of variable types. Force reevaluation of the
  integer classification (binary vs. general integer) based on current
  bounds.
*/
  const char *intVarOriginal = si.getColType(true) ;
  char *intVar = CoinCopyOfArray(intVarOriginal,nCols) ;
/*
  Get modifiable bound arrays and the current solution.
*/
  CoinMemcpyN(si.getColLower(),nCols,colLower) ;
  CoinMemcpyN(si.getColUpper(),nCols,colUpper) ;
  const double * colsol = si.getColSolution() ;
/*
  Stage 1: Preliminary processing.

  During Stage 1 we'll set `reasonable' bounds on variables and constraints,
  get rid of constraints that are judged unsuitable for further probing (too
  dense, or not amenable to propagation or implication), and do one bit of
  really simple coefficient tightening. At the end of Stage 1, we'll have a
  row- and column-major copy of the reduced and strengthened constraint
  system.
*/
/*
  Put reasonable bounds on general integer variables.
  TODO: Here's another source of bogus bounds, distinct from 1e20!
	-- lh, 101009 --
*/
  int numberIntegers = 0 ;
  for (int i = 0 ; i < nCols ; i++) {
    if (intVar[i]) {
      numberIntegers++ ;
      if (intVar[i] == 2) {
	if (colsol[i] < 1.0e10 && colUpper[i] > 1.0e12) 
	  colUpper[i] = CGL_REASONABLE_INTEGER_BOUND ;
	if (colsol[i] > -1.0e10 && colLower[i] < -1.0e12) 
	  colLower[i] = -CGL_REASONABLE_INTEGER_BOUND ;
      }
    }
  }
  bool feasible = true ;
/*
  If we're at the root of the search tree, scan the constraint system to see
  if there are naturally integer variables that are not declared as integer.
  Convert them to integer type.

  TODO: If we lose feasibility here, we should clean up and bail. This is not
	handled well in subsequent blocks of code.  -- lh, 101021 --
*/
  if (!info->inTree && !info->pass) {
    feasible = analyze(&si,intVar,colLower,colUpper) ;
  }
/*
  Attempt to tighten bounds based on reduced costs. Consider movement by +/-1
  and use that as a filter. For integer variables, this will fix the variable
  at the current bound. For continuous variables, we can restrict the gap
  between lower and upper bound to less than 1.0.

  jjf: Should be in CbcCutGenerator and check if basic.

  lh: (100917) Not sure why that's relevant, except that reduced cost is
      zero by definition for basic variables. Why not just test for d<j> = 0?
  
  lh: (100917) Seems like a prime candidate for a subroutine (inlined, perhaps).
      And that's a really good explanation for the redundant code. This
      probably *was* a subroutine at one time.

  This block disabled as of 100917 (PROBING_EXTRA_STUFF == false). Probably
  disabled further back than that.

  lh: (101021) And I'm asking myself `Why why?' in the print statements below.
      Apparently John was of the opinion that anything that could be fixed here
      should have been fixed elsewhere. Did this functionality actually move
      elsewhere, so that those print statements are just telltales in case
      something slipped by?
*/
  if (feasible && PROBING_EXTRA_STUFF) {
    const double * djs = si.getReducedCost() ;

    // Begin redundant code block -- we did this just above!
    const double * colsol = si.getColSolution() ;
    double direction = si.getObjSense() ;
    double cutoff ;
    bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff) ;
    if (!cutoff_available || usingObjective_ < 0) {
      cutoff = si.getInfinity() ;
    }
    cutoff *= direction ;
    if (fabs(cutoff) > 1.0e30)
      assert (cutoff > 1.0e30) ;
    // End redundant code block

    double current = si.getObjValue() ;
    current *= direction ;
    double gap = CoinMax(cutoff-current,1.0e-1) ;
    for (int i = 0; i < nCols; ++i) {
      double djValue = djs[i]*direction ;
      // Consider only variables that are not fixed.
      if (colUpper[i]-colLower[i] > 1.0e-8) {
	// If we're currently at lower bound ...
        if (colsol[i] < colLower[i]+primalTolerance_) {
	  // ... and +1 will drive us over the cutoff ...
          if (djValue > gap) {
	    // ... and we have to move at least +1 ...
            if (si.isInteger(i)) {
	      // ... fix at lower bound.
              printf("why int %d not fixed at lb\n",i) ;
              colUpper[i] = colLower[i] ;
            } else {
	      // ... if a fractional move is possible, calculate the max
	      // movement and crank down the upper bound.
              double newUpper = colLower[i] + gap/djValue ;
              if (newUpper < colUpper[i]) {
                //printf("%d ub from %g to %g\n",i,colUpper[i],newUpper) ;
                colUpper[i] = CoinMax(newUpper,colLower[i]+1.0e-5) ;
              }
            }
          }
        } else if (colsol[i] > colUpper[i]-primalTolerance_) {
          if (-djValue > gap) {
            if (si.isInteger(i)) {
              printf("why int %d not fixed at ub\n",i) ;
              colLower[i] = colUpper[i] ;
            } else {
              double newLower = colUpper[i] + gap/djValue ;
              if (newLower>colLower[i]) {
                //printf("%d lb from %g to %g\n",i,colLower[i],newLower) ;
                colLower[i]= CoinMin(newLower,colUpper[i]-1.0e-5) ;
              }
            }
          }
        }
      }
    }
  }  // End cost bound tightening using reduced costs.

  int ninfeas = 0 ;
  // Set up maxes
  int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_ ;
  int maxElements = info->inTree ? maxElements_ : maxElementsRoot_ ;

  // Forcing probing to consider really dense constraints wasn't good,
  // even at the root.
  // if (!info->inTree && !info->pass) maxElements = nCols ;

  // Get objective offset
  double offset ;
  si.getDblParam(OsiObjOffset,offset) ;
#ifdef COIN_DEVELOP
  // print it just once
  if (offset && !info->inTree && !info->pass) 
    printf("CglProbing obj offset %g\n",offset) ;
#endif
/*
  If we have a snapshot handy, most of the setup is done. Otherwise, we need to
  create row- and column-major copies of the constraint matrix. If the client
  has specified mode 0 (work from snapshot), it is confused, but we silently
  correct it.

  TODO: Ah ... maybe a `You're confused' message and an error return would be
	more appropriate? Seems like a serious algorithm error at a higher
	level.  -- lh, 100917 --

	Or maybe the client is within the documentation. After all, John's
	comment for mode 0 says `only available with snapshot; otherwise as
	mode 1'.  -- lh, 101021 --

  TODO: Note that the column-major copy is not created until much further
	down, after we've settled on the set of rows to process and applied
	various techniques for row strengthening.   -- lh, 101021 --
*/
  if (!rowCopy_) {
    nRows = si.getNumRows(); 
    if (mode == 0)
      mode = 1 ;
/*
  Make the row copy and create modifiable row bound arrays. If we have a
  nontrivial cutoff and are not forbidden from using it, add in the objective
  coefficients.

  TODO: For the objective, we make an entry in rowLower and rowUpper and
	install our bounds. That implies that the arrays passed in must have
	sufficient capacity.  -- lh, 100917 --

  TODO: Really, there should be a boolean for using the objective, so we don't
	keep reevaluating this condition. And why do we generate a boolean,
	maximise, instead of using sense*objcoeff, as is common elsewhere?
	-- lh, 100917 --

  TODO: No column copy?  -- lh, 100917 --
	It's made far down below, at the end of the row strengthening block,
	after we've decided on the rows to keep.
*/
    if (cutoff < 1.0e30 && usingObjective_ > 0) {
      rowCopy = new CoinPackedMatrix(*si.getMatrixByRow(),1,nCols,false) ;
    } else {
      rowCopy = new CoinPackedMatrix(*si.getMatrixByRow()) ;
    }
    if (cutoff < 1.0e30 && usingObjective_ > 0) {
      int *columns = new int[nCols] ;
      double *elements = new double[nCols] ;
      int n = 0 ;
      const double *objective = si.getObjCoefficients() ;
      bool maximize = (si.getObjSense() == -1) ;
      for (int i = 0 ; i < nCols ; i++) {
	if (objective[i]) {
	  elements[n] = (maximize) ? -objective[i] : objective[i] ;
	  columns[n++] = i ;
	}
      }
      rowCopy->appendRow(n,columns,elements) ;
      delete [] columns ;
      delete [] elements ;
      CoinMemcpyN(si.getRowLower(),nRows,rowLower) ;
      CoinMemcpyN(si.getRowUpper(),nRows,rowUpper) ;
      rowLower[nRows] = -DBL_MAX ;
      rowUpper[nRows] = cutoff+offset ;
      nRows++ ;
    } else {
      CoinMemcpyN(si.getRowLower(),nRows,rowLower) ;
      CoinMemcpyN(si.getRowUpper(),nRows,rowUpper) ;
    }
  }
/*
  If we're starting from a snapshot, make a copy for further processing. Make
  sure we adjust rowLower and rowUpper to accommodate the objective (the
  coefficients are already in).

  TODO: Another good reason to hive off snapshot as a separate class.
	-- lh, 100917 --

  TODO: Ummmm ... what if the objective isn't already in? I see, I'm checking
	usingObjective_, which (one hopes) is correctly set from before. But
	it's publically accessible. This should be an internal state thing.
	-- lh, 100917 --
  
  TODO: Note that we don't make a working copy of the column copy held by the
	snapshot.  -- lh, 101021 --

  TODO: Why do we need to avoid modifying the snapshot? Could we get by without
	making yet another copy?  -- lh, 101125 --
*/
  else {
    nRows = numberRows_ ;
    // Sanity -- our cached copy must have the same number of columns!
    assert(nCols == numberColumns_) ;
    rowCopy = new CoinPackedMatrix(*rowCopy_) ;
    // More sanity.
    assert (rowCopy_->getNumRows() == numberRows_) ;
    rowLower = new double[nRows] ;
    rowUpper = new double[nRows] ;
    CoinMemcpyN(rowLower_,nRows,rowLower) ;
    CoinMemcpyN(rowUpper_,nRows,rowUpper) ;
    if (usingObjective_ > 0) {
      rowLower[nRows-1] = -DBL_MAX ;
      rowUpper[nRows-1] = cutoff+offset ;
    }
  }
/*
  Remove rows with too many elements and rows that are deemed useless for
  further probing.  Try to strengthen rows in place, if allowed. Deal with
  fixed variables. We're going to do this in place by editing the internal
  data structures of rowCopy.

  TODO: At first glance, it looks like nRealRows is just so we can avoid
	tossing the objective (that was John's comment). But then why not
	just set the upper bound to nRows-1 or nRows-2, depending on
	usingObjective_?  -- lh, 100917 --

  TODO: This block of code looks like an excellent candidate for a subroutine.
	Maybe it once was, given that it's wrapped in its own code block.
	-- lh, 100923 --

  TODO: OUTRUBBISH is disabled, so either it didn't improve performance or
	it's buggy. -- lh, 100923
*/
  // int i ;
  CoinBigIndex * rowStartPos = NULL ;
  int * realRows = NULL ;
  // Begin row strengthen block
  {
    //#define OUTRUBBISH

    int *rowLength = rowCopy->getMutableVectorLengths(); 
    double * elements = rowCopy->getMutableElements() ;
    int * column = rowCopy->getMutableIndices() ;
    CoinBigIndex *rowStart = rowCopy->getMutableVectorStarts() ;
#ifdef OUTRUBBISH
    double large = 1.0e3 ;
#endif
    int nDelete = 0 ;
    int nKeep = 0 ;
    int *which = new int[nRows] ;
    int nElements = rowCopy->getNumElements() ;
    int nTotalOut = 0 ;
    int nRealRows = si.getNumRows() ;
/*
  Preliminary assessment for dense rows and free rows.  Count the number of
  coefficients that we'll remove. We're looking for rows that are too long or
  have row bounds so weak as to be ineffective. Don't toss the objective
  (index larger than nRealRows.
*/
    for (int i = 0 ; i < nRows ; i++) {
      if (rowLength[i] > maxElements ||
	  (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20)) {
	if (i < nRealRows) 
	  nTotalOut += rowLength[i] ;
      }
    }
/*
  If the total coefficients in dense rows amounts to less than 1/10 of the
  matrix, just deal with it (maxElements set to nCols won't exclude anything).
*/
    if (nTotalOut*10 < nElements)
      maxElements = nCols ;
#ifdef OUTRUBBISH
    int nExtraDel = 0 ;
#endif
/*
  Main processing loop(1).
  
  Tossing the junk will significantly reduce the constraint system, so let's
  get to it.

  Examine each row for length and effectiveness.  We can summarily toss free
  rows and rows that are too dense.

  TODO: We carefully exclude the objective from deletion for density, but note
	that 1e30 is the break for adding it and 1e20 is the break for removing
	it. Admittedly, between 1e20 and 1e30 there's not likely to be much
	effect, but it'd be nice to be consistent.  -- lh, 100923 --
*/
    for (int i = 0 ; i < nRows ; i++) {
      if ((rowLength[i] > maxElements && i < nRealRows) ||
	  (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20)) {
	which[nDelete++] = i ;
      } else {
/*
  This row passes the coefficient limit and has at least one finite bound.
  But is it useful for bound propagation and probing?

  TODO: This is lhs bound evaluation. Should be done up front.
	-- lh, 100917 --

  TODO: The OUTRUBBISH and !OUTRUBBISH blocks are in no way mutually exclusive,
	even though that's the configuration here. It's easy to conceive of a
	blend in which the OUTRUBBISH block runs to remove `weak' constraints
	at some threshold. Constraints that pass are then considered for
	strengthening.  -- lh, 100923 --
*/
#ifdef OUTRUBBISH
/* 
  jjf: out all rows with infinite plus and minus

  It's actually a bit more extreme than that; large is set to 1.0e3 above,
  which is a fair bit short of infinity. Very aggressive removal of constraints
  perceived to be weak. -- lh, 100923 --
*/
	int nPlus = (rowUpper[i] > -large) ? 0 : 1 ;
	int nMinus = (rowLower[i] < large) ? 0 : 1 ;
	CoinBigIndex start = rowStart[i] ;
	CoinBigIndex end = start + rowLength[i] ;
	for (CoinBigIndex j = start ; j < end ; j++) { 
	  int iColumn = column[j] ;
	  if (colUpper[iColumn] > large) {
	    if (elements[j] > 0.0)
	      nPlus++ ;
	    else
	      nMinus++ ;
	  }
	  if (colLower[iColumn] < -large) {
	    if (elements[j] < 0.0)
	      nPlus++ ;
	    else
	      nMinus++ ;
	  }
	}
	if (!nPlus||!nMinus) {
	  rowLower[nKeep] = rowLower[i] ;
	  rowUpper[nKeep] = rowUpper[i] ;
	  nKeep++ ;
	} else {
	  nExtraDel++ ;
	  which[nDelete++] = i ;
	}
#else  // OUTRUBBISH
/*
  The (way) less aggressive approach. We won't actually delete anything
  here.  We're only looking at rows that can be strengthened in place.
  Call this Strengthen(1). See notes for details.

  Process this row only if info says we can strengthen rows in place; only on
  the first pass for this node; only if the row has at least one infinite
  bound (lower or upper)
*/
	if (info->strengthenRow && !info->pass &&
	    (rowLower[i] < -1.0e20 || rowUpper[i] > 1.0e20)) {
/*
  Assess the row. If we find anything other than binary variables, we're
  immediately done.

  TODO: Really, the loop should abort as soon as nPlus > 1 & nMinus > 1.
	Presumably it doesn't because the test would take more time than
	just finishing the processing?
	-- lh, 100923 --

  TODO: coefficients of zero are lumped into nMinus, but arguably a coeff
	of zero shouldn't exist in a sparse matrix structure.
	-- lh, 100923 --

  TODO: Is the code not consistent about labelling binary variables? This is
	precisely the condition used to set the code for binary in intVar.
	-- lh, 100923 --

  TODO: Doesn't allow for nPlus == 1 and nMinus == 1. Maybe that's a special
	case that gets caught elsewhere? And why increment effectiveness
	rather than simply assign to it?   -- lh, 100923 --
*/
	  int nPlus = 0 ;
	  int nMinus = 0 ;
	  for (CoinBigIndex j = rowStart[i] ; j < rowStart[i+1] ; j++) {
	    int jColumn = column[j] ;
	    if (intVar[jColumn] &&
		colLower[jColumn] == 0.0 && colUpper[jColumn] == 1.0) {
	      double value = elements[j] ;
	      if (value > 0.0) {
		nPlus++ ;
	      } else {
		nMinus++ ;
	      }
	    } else {
	      nPlus = 2 ;
	      nMinus = 2 ;
	      break ;
	    }
	  }
/*
  Look for the single positive coefficient in the row and reduce it. (Similarly
  for negative coefficient in the next block.

  TODO: So why don't we just break from the loop after we've tweaked the
	coefficient?  There's no real loop here to accumulate effectiveness.
	-- lh, 100923 --
*/
	  double effectiveness = 0.0 ;
	  if (nPlus == 1 && rowUpper[i] > 0.0 && rowUpper[i] < 1.0e10) {
	    // can make element smaller
	    for (CoinBigIndex j = rowStart[i] ; j < rowStart[i+1] ; j++) {
	      double value = elements[j] ;
	      if (value > 0.0) {
		elements[j] -= rowUpper[i] ;
		//printf("pass %d row %d plus el from %g to %g\n",info->pass,
		//     i,elements[j]+rowUpper[i],elements[j]) ;
	      }
	      effectiveness += fabs(elements[j]) ;
	    }
	    rowUpper[i] = 0.0 ;
	  } else
	  if (nMinus == 1 && rowLower[i] < 0.0 && rowLower[i] > -1.0e10) {
	    // can make element smaller in magnitude
	    for (CoinBigIndex j = rowStart[i] ; j < rowStart[i+1] ; j++) {
	      double value = elements[j] ;
	      if (value < 0.0) {
		elements[j] -= rowLower[i] ;
		//printf("pass %d row %d minus el from %g to %g\n",info->pass,
		//     i,elements[j]+rowLower[i],elements[j]) ;
	      }
	      effectiveness += fabs(elements[j]) ;
	    }
	    rowLower[i]=0.0 ;
	  }
/*
  If the coefficient is nonzero, remember the result. We need to remember the
  result because (as with everything done here) we're working with a copy.
  This business of strengthening the coefficient must be propagated back out
  to the client, who gets to make the decision about whether to use it.

  TODO: Why do we create a block-local object, then create a clone for later
	use? Why not just use new?  -- lh, 100923 --
*/
	  if (effectiveness) {
	    OsiRowCut rc ;
	    rc.setLb(rowLower[i]) ;
	    rc.setUb(rowUpper[i]) ;
	    int start = rowStart[i] ;
	    rc.setRow(rowLength[i],column+start,elements+start,false) ;
	    rc.setEffectiveness(effectiveness) ;
	    assert (!info->strengthenRow[i]) ;
	    info->strengthenRow[i] = rc.clone() ;
	  }
	}
	rowLower[nKeep] = rowLower[i] ;
	rowUpper[nKeep] = rowUpper[i] ;
	nKeep++ ;
#endif   // OUTRUBBISH
      }
    }  // End main processing loop(1)
/*
  Did we delete anything in processing loop(1)? If so, deal with it. We need
  to delete the rows, of course. There may also be bookkeeping.  In
  particular, strengthened constraints are held in strengthenRow by the
  original row index. We need a cross-reference, realRows, to keep track.

  TODO: Huh? What's the deal with the objective here? Why does it need
	special treatment in realRows?  -- lh, 100923 -- 

  TODO: If we have no constraints remaining (nKeep == 0), then surely we could
	skip this bit of work.  -- lh, 101021 --

  TODO: Because realRows is not allocated unless we delete rows from the
	system, we will not be able to strengthen in place. See the matching
	TODO in probe().  -- lh, 101209 --
*/
    if (nDelete) {
#ifdef OUTRUBBISH
      if (nExtraDel) {
	printf("%d rows deleted (extra %d)\n",nDelete,nExtraDel) ;
      }
#else
#endif
      if (info->strengthenRow) {
	realRows = new int [nRows] ;
	CoinZeroN(realRows,nRows) ;
	for (int i = 0 ; i < nDelete ; i++)
	  realRows[which[i]] = -1 ;
	int k = 0 ;
	for (int i = 0 ; i < nRows ; i++) {
	  if (!realRows[i]) {
	    if (i < nRealRows)
	      realRows[k++] = i ; // keep
	    else
	      realRows[k++] = -1 ; // objective - discard
	  }
	}
      }
      rowCopy->deleteRows(nDelete,which) ;
      nRows = nKeep ;
    }
    delete [] which ;
/*
  Do we have any constraints left to work with? If not, clean up and go home.
*/
    if (!nRows) {
#ifdef COIN_DEVELOP
      printf("All rows too long for probing\n") ;
#endif
      delete rowCopy ;
      // Read `if we're working from a snapshot'
      if (rowCopy_) {
	delete [] rowLower ;
	delete [] rowUpper ;
      }
      delete [] intVar ;
/*
  jjf: and put back unreasonable bounds on integer variables

  TODO: Do we need to be so careful here? Why not just do a block copy?
	-- lh, 100923 --
*/
      const double * trueLower = si.getColLower() ;
      const double * trueUpper = si.getColUpper() ;
      for (int i = 0 ; i < nCols ; i++) {
	if (intVarOriginal[i] == 2) {
	  if (colUpper[i] == CGL_REASONABLE_INTEGER_BOUND) 
	    colUpper[i] = trueUpper[i] ;
	  if (colLower[i] == -CGL_REASONABLE_INTEGER_BOUND) 
	    colLower[i] = trueLower[i] ;
	}
      }
      delete [] realRows ;
      return 0 ;
    }
/*
  We're still inside the row strengthen block ...

  jjf: Out elements for fixed columns and sort

  In a few more words: Process the rows, adjusting rowLower and rowUpper for
  fixed variables; the columns are physically removed. As we rewrite the
  coefficients, sort them so that negative coefficients precede positive
  coefficients.

  We need to refresh the mutable pointers; deleteRows above may have rendered
  the previous pointers invalid.

  TODO: This code duplicates a block in snapshot(), except that the code in
	snapshot doesn't substitute for the values of fixed variables. The
	two blocks should probably be pulled out into a method.
*/
    elements = rowCopy->getMutableElements() ;
    column = rowCopy->getMutableIndices() ;
    rowStart = rowCopy->getMutableVectorStarts() ;
    rowLength = rowCopy->getMutableVectorLengths(); 

#if 0
    // For debugging?
    int nFixed = 0 ;
    for (i = 0 ; i < nCols ; i++) {
      if (colUpper[i] == colLower[i])
	nFixed++ ;
    }
    printf("%d columns fixed\n",nFixed) ;
#endif

    CoinBigIndex newSize = 0 ;
    int * column2 = new int[nCols] ;
    double * elements2 = new double[nCols] ;
    rowStartPos = new CoinBigIndex [nRows] ;
    for (int i = 0 ; i < nRows ; i++) {
      double offset = 0.0 ;
      CoinBigIndex start = rowStart[i] ;
      CoinBigIndex end = start+rowLength[i] ;
      rowStart[i] = newSize ;
      CoinBigIndex save = newSize ;
      int nOther = 0 ;
      for (CoinBigIndex j = start ; j < end ; j++) { 
	int iColumn = column[j] ;
	if (colUpper[iColumn] > colLower[iColumn]) {
	  double value = elements[j] ;
	  if (value < 0.0) {
	    elements[newSize] = value ;
	    column[newSize++] = iColumn ;
	  } else {
	    elements2[nOther] = value ;
	    column2[nOther++] = iColumn ;
	  }
	} else {
	  offset += colUpper[iColumn]*elements[j] ;
	}
      }
      rowStartPos[i] = newSize ;
      for (int k = 0 ; k < nOther ; k++) {
	elements[newSize] = elements2[k] ;
	column[newSize++] = column2[k] ;
      }
      rowLength[i] = newSize-save ;
      if (offset) {
	if (rowLower[i] > -1.0e20)
	  rowLower[i] -= offset ;
	if (rowUpper[i] < 1.0e20)
	  rowUpper[i] -= offset ;
      }
    }
    delete [] column2 ;
    delete [] elements2 ;
    rowStart[nRows] = newSize ;
    rowCopy->setNumElements(newSize) ;
  } 
  // End row strengthen block

/*
  Here, finally, we're making a column-major copy, after we've winnowed the
  rows for the set we want to work with.

  TODO: Still, I'm asking myself, ``Why does the snapshot keep a colum-major
	copy?  -- lh, 101021 --
*/
  CoinPackedMatrix * columnCopy = new CoinPackedMatrix(*rowCopy,0,0,true) ;

/*
  End of Stage 1.
*/

#if CGL_DEBUG > 0
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger() ;
  if (debugger && !debugger->onOptimalPath(si))
    debugger = NULL ;
#else
  const OsiRowCutDebugger * debugger = NULL ;
#endif

/*
  Begin Stage 2
*/

  const int *column = rowCopy->getIndices() ;
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *rowLength = rowCopy->getVectorLengths(); 
  const double *rowElements = rowCopy->getElements() ;
/*
  jjf: Arrays so user can find out what happened

  TODO: Misleading? lookedAt_ seems to be used to hold the set of variables
	selected for probing. The question is whether the selection activity
	is duplicated in probe() and/or probeCliques().
*/
  if (!lookedAt_) {
    lookedAt_ = new int[nCols] ;
  }
  numberThisTime_ = 0 ;
/*
  jjf: Let us never add more than twice the number of rows worth of row cuts
  jjf: Keep cuts out of cs until end so we can find duplicates quickly

  TODO: Well, ok, but none of the values for nRowsFake agree with that comment
	(assuming nRowsFake is the limit on the number of row cuts). And the
	row_cut constructor also adjusts the value.
	-- lh, 100923 --

  TODO: There's a similar calculation in probe() to size the local container
	used there. It allocates much more space.
*/
  int nRowsSafe = CoinMin(nRows,si.getNumRows()) ;
  int nRowsFake = info->inTree ? nRowsSafe/3 : nRowsSafe ;
  if (!info->inTree && !info->pass) 
    nRowsFake *= 5 ;
  CglProbingRowCut rowCut(nRowsFake,!info->inTree) ;

  int * markR = new int [nRows] ;
  double * minR = new double [nRows] ;
  double * maxR = new double [nRows] ;
/*
  In modes 1 and 2, we're allowed to tighten bounds. Call tighten() to do bound
  propagation (two pass limit), then create `column cuts' to record the
  tightened bounds.

  TODO: Surely we could skip this if we've already lost feasibility! Or does
	tighten() bail immediately? Not immediately obvious. This method
	really needs a cleanupAndGoHome() method that can be called at the
	point where we discover infeasibility.
	-- lh, 100923 --

  NOTE: Split tighten() into tighten() and tightenClique() -- lh, 101124 --
*/
  if (mode) {
    if (cliqueRowStart_ && numberRows_ && cliqueRowStart_[numberRows_]) {
      ninfeas = tightenClique(colLower,colUpper,column,rowElements,
			      rowStart,rowStartPos,rowLength,
			      rowLower,rowUpper,nRows,nCols,
			      intVar,2,primalTolerance_);
    } else {
      ninfeas = tighten(colLower,colUpper,column,rowElements,
			rowStart,rowStartPos,rowLength,
			rowLower,rowUpper,nRows,nCols,
			intVar,2,primalTolerance_);
    }
    if (!feasible)
      ninfeas = 1 ;
    if (!ninfeas) {   // block extends to end of if (mode)
/*
  jjf: create column cuts where integer bounds have changed

  TODO: Looks like another good candidate for a method. Already a block.
	-- lh, 100923 --
  
  TODO: ifCut could just as well be a boolean.  -- lh, 100923 --
*/
      {
	const double *lower = si.getColLower() ;
	const double *upper = si.getColUpper() ;
	const double *colsol = si.getColSolution() ;
	int numberChanged = 0 ;
	int ifCut = 0 ;
	CoinPackedVector lbs ;
	CoinPackedVector ubs ;
	for (int i = 0 ; i < nCols ; ++i) {
	  if (intVar[i]) {
	    colUpper[i] = CoinMin(upper[i],floor(colUpper[i]+1.0e-4)) ;
	    if (colUpper[i] < upper[i]-1.0e-8) {
	      if (colUpper[i] < colsol[i]-1.0e-8)
		ifCut = 1 ;
	      ubs.insert(i,colUpper[i]) ;
	      numberChanged++ ;
	    }
	    colLower[i] = CoinMax(lower[i],ceil(colLower[i]-1.0e-4)) ;
	    if (colLower[i] > lower[i]+1.0e-8) {
	      if (colLower[i] > colsol[i]+1.0e-8)
		ifCut = 1 ;
	      lbs.insert(i,colLower[i]) ;
	      numberChanged++ ;
	    }
	  }
	}
	if (numberChanged) {
	  OsiColCut cc ;
	  cc.setUbs(ubs) ;
	  cc.setLbs(lbs) ;
	  if (ifCut) {
	    cc.setEffectiveness(100.0) ;
	  } else {
	    cc.setEffectiveness(1.0e-5) ;
	  }
#if CGL_DEBUG > 0
	  checkBounds(si,cc) ;
#endif
	  cs.insert(cc) ;
	}
      }  // end create column cuts block

      // This next block extends to end of mode 1,2 block
      if (maxProbe > 0) {
        numberThisTime_ = 0 ;
/*
  jjf: get min max etc for rows

  Calculate L<i> and U<i> for each row. If one or both are sufficiently small
  to process further, markR[i] is set to -1. If both are too large to warrant
  further processing, markR[i] is set to -2.

  TODO: Why do we need a second run for this? Why not do it as part of
	tighten()? My old notes indicate that tighten calculates the row
	bounds (and indeed it must to propagate column bounds). On a purely
	mechanical level, tighten does not fill in minR, maxR, markR.
	-- lh, 100923 --
*/
        tighten2(colLower,colUpper,column,rowElements,
                 rowStart,rowLength,rowLower,rowUpper,
                 minR,maxR,markR,nRows) ;
/*
  Decide what to look at. Mode 1 is unsatisfied variables only; mode 2 is all
  variables.

  At the root, mode 1 will still look at all integer variables, but the order
  is sorted by infeasibility and flipped on alternate passes. In mode 2,
  there's no sort.

  array holds pairs (index,away*multiplier). As explained below, away is
  small when infeasibility is large, so large negative values are variables
  that are close to integrality.
*/
        if (mode == 1) {
          const double *colsol = si.getColSolution() ;
          double_int_pair *array = new double_int_pair [nCols] ;
#	  ifdef ZEROFAULT
	  std::memset(array,0,sizeof(double_int_pair)*nCols) ;
#	  endif
	  double multiplier = -1.0 ;
/*
  Examine unfixed integer variables and record the ones that are unsatisfied
  (i.e., not at an integral value).

  TODO: Why do we change the sign of the infeasibility between even and odd
	passes at the root? Best guess is that it reverses the order in which
	integer variables are processed. -- lh, 100923 --

  TODO: the constant .49999 is actually a hardwired integer feasibility
	tolerance of 1e-5. This should be a parameter; better, an attribute
	of the CglProbing object. Presumably the game with .49999-(*) and
	multiplier is a convenience for sorting (most infeasible are close
	to zero, nearly integral are large and negative). But the comparison
	method (double_int_pair_compare) isn't used elsewhere, so maybe this
	can be simplified.  -- lh, 100924 --

  TODO: Note that if we're at the root, we look at all variables.
	-- lh, 100924 --

  TODO: At one time, numberThisTime_ was capped at maxProbe, but that's
	disabled. Enforced elsewhere?  -- lh, 100924 --
*/
	  if (info->inTree || (info->pass&1) != 0)
	    multiplier = 1.0 ;
	  // const int * columnLength = si.getMatrixByCol()->getVectorLengths() ;
          for (int i = 0 ; i < nCols ; i++) {
            if (intVar[i] && (colUpper[i]-colLower[i]) > 1.0e-8) {
              double away = fabs(0.5-(colsol[i]-floor(colsol[i]))) ;
              if (away < 0.49999 || !info->inTree) {
                //array[numberThisTime_].infeasibility = away ;
                array[numberThisTime_].infeasibility = away*multiplier ;
                //array[numberThisTime_].infeasibility = -columnLength[i] ;
                array[numberThisTime_++].sequence = i ;
              }
            }
          }
	  // printf("maxP %d num %d\n",maxProbe,numberThisTime_) ;
          std::sort(array,array+numberThisTime_,double_int_pair_compare()) ;
          // numberThisTime_ = CoinMin(numberThisTime_,maxProbe) ;
          for (int i = 0 ; i < numberThisTime_ ; i++) {
            lookedAt_[i] = array[i].sequence ;
          }
          delete [] array ;
        }
/*
  For mode 2 or 3, probe all integer variables that are not fixed.
*/
	else {
          for (int i = 0 ; i < nCols ; i++) {
            if (intVar[i] && (colUpper[i]-colLower[i]) > 1.0e-8) {
              lookedAt_[numberThisTime_++] = i ;
            }
          }
        }
/*
  lookedAt_ now contains a list of indices of variables to probe. Rows worth
  processing are tagged with -1 in markR 
*/
#if 0
/*
  TODO: This bit of code seems to have migrated up to the head of Stage 1,
	main processing loop (1), where it's used to limit the constraints
	installed in rowCopy. Presumably this is the origin of the abort().
	We should never see a dense row in the working matrix.
*/
        // Only look at short rows
        for (i=0;i<nRows;i++) {
          if (rowLength[i]>maxElements)
            abort(); //markR[i]=-2 ;
        }
#endif
/*
  TODO: Eh? Sort by index for repeatability? Seems pointless, and it is
	commented out. -- lh, 100924 --
*/
	// sort to be clean
	// std::sort(lookedAt_,lookedAt_+numberThisTime_) ;
/*
  And do the probing, with or without cliques.
*/
        if (!numberCliques_) {
          ninfeas = probe(si,cs,colLower,colUpper,rowCopy,columnCopy,
                          rowStartPos,realRows,rowLower,rowUpper,
                          intVar,minR,maxR,markR,info) ;
        } else {
          ninfeas = probeCliques(si,debugger,cs,
				 colLower,colUpper,rowCopy,columnCopy,
                                 realRows,rowLower,rowUpper,
                                 intVar,minR,maxR,markR,info) ;
        }
      } // end maxProbe > 0
    } // end !nInfeas
  } // end mode != 0
/*
  End of block for mode = 1 or 2; begin mode = 0. This works with the snapshot
  and does not do bound propagation.

  TODO: In fact, it's unclear to me just what we're doing here. The activity
	is completely unsymmetric with the previous case. -- lh, 100924 --
*/
  else
  if (maxProbe > 0) {
/*
  jjf: global cuts from previous calculations
       could check more thoroughly that integers are correct (sanity)
       make up list of new variables to look at 

  Initial hypothesis: We're working from a snapshot, so we've done at least
  one round of probing and cut generation. In particular, we should have
  disaggregation cuts, stashed in cutVector_. So we're only going to probe
  variables where we have no information. We're looking at binary
  variables because that's the natural index for a collection of
  disaggregation cuts generated by splitting a constraint where a single
  binary variable controls whether a bunch of other variables can be
  nonzero.

  The comment in CglProbing for mode 0: `if no information exists for
  that variable then probing will be done' makes more sense now. And if you
  take that together with the comment `only with snapshot; otherwise as mode
  1' then the tweaking of mode at the head of generateCutsAndModify almost
  makes sense.

  -- lh, 101021 --

  TODO: Check out how we can get away with reordering the major operations
	here. In the mode != 0 block, it's tighten2; select & sort; probe.
	Here, it's select and sort; tighten2; probe.  -- lh, 101021 --

	Now I think about it, there's no problem. Tighten2 flags rows worth
	processing, which is independent of columns worth processing. Order
	makes no difference.  -- lh, 101125 --

*/
    assert(numberIntegers == numberIntegers_) ;
    numberThisTime_ = 0 ;
    const double *colsol = si.getColSolution() ;
    double_int_pair *array = new double_int_pair [nCols] ;
#   ifdef ZEROFAULT
    std::memset(array,0,sizeof(double_int_pair)*nCols) ;
#   endif
    for (int i = 0 ; i < number01Integers_ ; i++) {
      int j = cutVector_[i].sequence ;
      if (!cutVector_[i].index && (colUpper[j]-colLower[j]) > 1.0e-8) {
	double away = fabs(0.5-(colsol[j]-floor(colsol[j]))) ;
        array[numberThisTime_].infeasibility = away ;
        array[numberThisTime_++].sequence = i ;
      }
    }
    std::sort(array,array+numberThisTime_,double_int_pair_compare()) ;
    numberThisTime_ = CoinMin(numberThisTime_,maxProbe) ;
    for (int i = 0 ; i < numberThisTime_ ; i++) {
      lookedAt_[i] = array[i].sequence ;
    }
    // sort to be clean
    // std::sort(lookedAt_,lookedAt_+numberThisTime_) ;
    delete [] array ;
    // get min max etc for rows
    tighten2(colLower, colUpper, column, rowElements,
	     rowStart, rowLength, rowLower, rowUpper,
	     minR , maxR , markR, nRows) ;
/*
  jjf: don't do cuts at all if 0 (i.e. we are just checking bounds)

  TODO: The comment is wrong, at least for generateCuts. What happens there is
	that any negative value of rowCuts_ is converted to 0x4, which is
	supposedly column cuts only. But it's definitely not 0, so this block
	will still execute.  -- lh, 100924 --

  TODO: Why don't we get a choice between probe and probeCliques? And why did
	we create a new cut collection (csNew) instead of using parameter cs?
	-- lh, 100924 --
*/
    OsiCuts csNew ;
    if (rowCuts_) {
#if 0
/*
  As for mode = 1 or 2; this bit of code migrated up to the head of the Stage 1
  main processing loop.
*/
      // Only look at short rows
      for (i=0;i<nRows;i++) {
        if (rowLength[i]>maxElements)
          abort(); //markR[i]=-2 ;
      }
#endif
      ninfeas = probeCliques(si,debugger,csNew,colLower,colUpper,
			     rowCopy,columnCopy,
			     realRows,rowLower,rowUpper,
			     intVar,minR,maxR,markR,info) ;
    }
/*
  TODO: And I have no idea what we're getting into here --- completely
	asymmetric to the mode 1 and 2 case. I'm going to bail from here for a
	bit and look over the startup code and probing methods, then come back
	to this.  -- lh, 100924 --

	Made it back this far with an understanding of cliques, but have been
	spending far too much time at home with Windows. As a working theory,
	it's this next bit that does disaggregation cuts, or strengthens
	disaggregation cuts with implication cuts from probing.

	If I interpret the slides from the 3/97 presentation correctly,
	probing should do implication cuts, and disaggregation would be a
	separate activity.  So I'm off to look at the probing code and
	understand what it does. There's 500 lines of uncommented code coming
	up, followed by another 200 lines with the comment `see if any
	disaggregation cuts are violated.' I'll be back. -- lh, 101021 --
*/
    if (!ninfeas) {
      // go through row cuts
      int nCuts = csNew.sizeRowCuts() ;
      int iCut ;
      // need space for backward lookup
      // just for ones being looked at
      int *backward = new int [2*nCols] ;
      int *onList = backward+nCols ;
      for (int i = 0 ; i < nCols ; i++) {
        backward[i] = -1 ;
        onList[i] = 0 ;
      }
      for (int i = 0 ; i < number01Integers_ ; i++) {
	int j = cutVector_[i].sequence ;
	backward[j] = i ;
	onList[j] = 1 ;
      }
      // first do counts
      // we know initialized to zero
      for (iCut=0;iCut<nCuts;iCut++) {
	OsiRowCut rcut ;
	CoinPackedVector rpv ;
	rcut = csNew.rowCut(iCut) ;
	rpv = rcut.row() ;
	assert(rpv.getNumElements()==2) ;
	const int * indices = rpv.getIndices() ;
	double* elements = rpv.getElements() ;
	double lb=rcut.lb() ;
	// find out which integer
        int which=0 ;
        int i=backward[indices[0]] ;
        if (i<0||!onList[indices[0]]) {
          which=1 ;
          i=backward[indices[1]] ;
          // Just possible variable was general integer but now 0-1
          if (!onList[indices[which]])
            continue ;
        }
        int other = indices[1-which] ;
	if (lb==-DBL_MAX) {
          if (!rcut.ub()) {
            // UB
            if (elements[which]<0.0) {
              //assert (elements[1-which]>0.0) ;
              // delta to 0 => x to 0.0
              cutVector_[i].length++ ;
            } else {
              if (elements[1-which]<0.0 && fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++ ;
              } else {
                if (onList[other]) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other] ;
                    cutVector_[i].length++ ;
                    cutVector_[j].length++ ;
                  } else {
                    continue ;
                  }
                }
              }
            }
          } else {
            if (onList[other]) {
              if (elements[0]==1.0 && elements[1]==1.0&&rcut.ub()==1.0) {
                // can do something ?
                int j=backward[other] ;
                cutVector_[i].length++ ;
                cutVector_[j].length++ ;
              } else {
                continue ;
              }
            }
          }
	} else {
          assert(rcut.ub()==DBL_MAX) ;
          if (!lb) {
            // LB
            if (elements[which]>0.0) {
              //assert (elements[1-which]<0.0) ;
              // delta to 0 => x to 0.0
              // flip so same as UB
              cutVector_[i].length++; 
            } else {
              if (elements[1-which]<0.0 && fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++ ;
              } else {
                if (onList[other]) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other] ;
                    cutVector_[i].length++ ;
                    cutVector_[j].length++ ;
                  } else {
                    continue ;
                  }
                }
              }
            }
          }
	} 
      }
      // allocate space
      for (int i=0;i<number01Integers_;i++) {
	int j=cutVector_[i].sequence ;
	if (onList[j] && !cutVector_[i].index) {
	  disaggregation thisOne=cutVector_[i] ;
	  cutVector_[i].index=new disaggregationAction [thisOne.length] ;
          cutVector_[i].length=0 ;
	}
      }
      // now put in
      for (iCut=0;iCut<nCuts;iCut++) {
	OsiRowCut rcut ;
	CoinPackedVector rpv ;
	int iput ;
	rcut = csNew.rowCut(iCut) ;
	rpv = rcut.row() ;
	assert(rpv.getNumElements()==2) ;
	const int * indices = rpv.getIndices() ;
	double* elements = rpv.getElements() ;
	double lb=rcut.lb() ;
	// find out which integer
	// find out which integer
        int which=0 ;
        int i=backward[indices[0]] ;
        if (i<0||!onList[indices[0]]) {
          which=1 ;
          i=backward[indices[1]] ;
          // Just possible variable was general integer but now 0-1
          if (!onList[indices[which]])
            continue ;
        }
        int other = indices[1-which] ;
        int j = other ? backward[other] : -1 ;
	if (lb==-DBL_MAX) {
          if (!rcut.ub()) {
            // UB
            if (elements[which]<0.0) {
              iput=cutVector_[i].length ;
              if (j>=0)
                setAffectedInDisaggregation(cutVector_[i].index[iput],j) ;
              else
                setAffectedInDisaggregation(cutVector_[i].index[iput],other) ;
              setWhenAtUBInDisaggregation(cutVector_[i].index[iput],false) ;
              setAffectedToUBInDisaggregation(cutVector_[i].index[iput],false) ;
	      setZeroOneInDisaggregation(cutVector_[i].index[iput],onList[other]!=0) ;
              cutVector_[i].length++ ;
            } else { 
              if (elements[1-which]<0.0 && fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                // delta to 1 => x to upper bound
                iput=cutVector_[i].length ;
                if (j>=0)
                  setAffectedInDisaggregation(cutVector_[i].index[iput],j) ;
                else
                  setAffectedInDisaggregation(cutVector_[i].index[iput],other) ;
                setWhenAtUBInDisaggregation(cutVector_[i].index[iput],true) ;
                setAffectedToUBInDisaggregation(cutVector_[i].index[iput],true) ;
                setZeroOneInDisaggregation(cutVector_[i].index[iput],onList[other]!=0) ;
                cutVector_[i].length++ ;
              } else {
                if (onList[other]) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other] ;
                    assert (j>=0) ;
                    // flip so value0 1.0
                    if (value1==1.0) {
                      j=i ;
                      i=backward[other] ;
                      value1=value0 ;
                      value0=1.0 ;
                    }
                    assert (value0==1.0) ;
                    assert (value1==-1.0) ;
                    iput=cutVector_[i].length ;
                    setAffectedInDisaggregation(cutVector_[i].index[iput],j) ;
                    setWhenAtUBInDisaggregation(cutVector_[i].index[iput],true) ;
                    setAffectedToUBInDisaggregation(cutVector_[i].index[iput],true) ;
                    setZeroOneInDisaggregation(cutVector_[i].index[iput],true) ;
                    cutVector_[i].length++ ;
                    iput=cutVector_[j].length ;
                    setAffectedInDisaggregation(cutVector_[j].index[iput],i) ;
                    setWhenAtUBInDisaggregation(cutVector_[j].index[iput],false) ;
                    setAffectedToUBInDisaggregation(cutVector_[j].index[iput],false) ;
                    setZeroOneInDisaggregation(cutVector_[j].index[iput],true) ;
                    cutVector_[j].length++ ;
                  }
                }
              }
            }
          } else {
            if (onList[other]) {
              if (elements[0]==1.0 && elements[1]==1.0&&rcut.ub()==1.0) {
                // can do something ?
                int j=backward[other] ;
                assert (j>=0) ;
                iput=cutVector_[i].length ;
                setAffectedInDisaggregation(cutVector_[i].index[iput],j) ;
                setWhenAtUBInDisaggregation(cutVector_[i].index[iput],true) ;
                setAffectedToUBInDisaggregation(cutVector_[i].index[iput],false) ;
                setZeroOneInDisaggregation(cutVector_[i].index[iput],true) ;
                cutVector_[i].length++ ;
                iput=cutVector_[j].length ;
                setAffectedInDisaggregation(cutVector_[j].index[iput],i) ;
                setWhenAtUBInDisaggregation(cutVector_[j].index[iput],true) ;
                setAffectedToUBInDisaggregation(cutVector_[j].index[iput],false) ;
                setZeroOneInDisaggregation(cutVector_[j].index[iput],true) ;
                cutVector_[j].length++ ;
              } else {
#ifdef COIN_DEVELOP
                abort() ;
#endif
		continue ;
              }
            }
          }
	} else {
          assert(rcut.ub()==DBL_MAX) ;
          if (!lb) {
            // LB
            if (elements[which]>0.0) {
              iput=cutVector_[i].length ;
              if (j>=0)
                setAffectedInDisaggregation(cutVector_[i].index[iput],j) ;
              else
                setAffectedInDisaggregation(cutVector_[i].index[iput],other) ;
              setWhenAtUBInDisaggregation(cutVector_[i].index[iput],false) ;
              setAffectedToUBInDisaggregation(cutVector_[i].index[iput],false) ;
              setZeroOneInDisaggregation(cutVector_[i].index[iput],onList[other]!=0) ;
              cutVector_[i].length++ ;
            } else { 
              if (elements[1-which]<0.0 && fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                iput=cutVector_[i].length ;
                if (j>=0)
                  setAffectedInDisaggregation(cutVector_[i].index[iput],j) ;
                else
                  setAffectedInDisaggregation(cutVector_[i].index[iput],other) ;
                setWhenAtUBInDisaggregation(cutVector_[i].index[iput],true) ;
                setAffectedToUBInDisaggregation(cutVector_[i].index[iput],true) ;
                setZeroOneInDisaggregation(cutVector_[i].index[iput],onList[other]!=0) ;
                cutVector_[i].length++ ;
              } else {
                if (onList[other]) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other] ;
                    assert (j>=0) ;
                    // flip so value0 -1.0
                    if (value1==-1.0) {
                      j=i ;
                      i=backward[other] ;
                      value1=value0 ;
                      value0=-1.0 ;
                    }
                    assert (value0==-1.0) ;
                    assert (value1==1.0) ;
                    iput=cutVector_[i].length ;
                    setAffectedInDisaggregation(cutVector_[i].index[iput],j) ;
                    setWhenAtUBInDisaggregation(cutVector_[i].index[iput],true) ;
                    setAffectedToUBInDisaggregation(cutVector_[i].index[iput],true) ;
                    setZeroOneInDisaggregation(cutVector_[i].index[iput],true) ;
                    cutVector_[i].length++ ;
                    iput=cutVector_[j].length ;
                    setAffectedInDisaggregation(cutVector_[j].index[iput],i) ;
                    setWhenAtUBInDisaggregation(cutVector_[j].index[iput],false) ;
                    setAffectedToUBInDisaggregation(cutVector_[j].index[iput],false) ;
                    setZeroOneInDisaggregation(cutVector_[j].index[iput],true) ;
                    cutVector_[j].length++ ;
                  }
                }
              }
            }
          }
	} 
      }
      delete [] backward ;
      // Now sort and get rid of duplicates
      // could also see if any are cliques
      int longest=0 ;
      for (int i=0;i<number01Integers_;i++) 
        longest = CoinMax(longest, cutVector_[i].length) ;
      unsigned int * sortit = new unsigned int[longest] ;
      for (int i=0;i<number01Integers_;i++) {
        disaggregation & thisOne=cutVector_[i] ;
        int k ;
        int number = thisOne.length ;
        for (k=0;k<number;k++) {
          int affected = affectedInDisaggregation(thisOne.index[k]) ;
          int zeroOne = zeroOneInDisaggregation(thisOne.index[k]) ? 1 : 0 ;
          int whenAtUB = whenAtUBInDisaggregation(thisOne.index[k]) ? 1 : 0 ;
          int affectedToUB = affectedToUBInDisaggregation(thisOne.index[k]) ? 1: 0 ;
          sortit[k]=(affected<<3)|(zeroOne<<2)|(whenAtUB<<1)|affectedToUB ;
        }
        std::sort(sortit,sortit+number) ;
        int affectedLast = 0xffffffff ;
        int zeroOneLast = 0 ;
        int whenAtUBLast = 0 ;
        int affectedToUBLast = 0; 
        int put=0 ;
        for (k=0;k<number;k++) {
          int affected = sortit[k]>>3 ;
          int zeroOne = (sortit[k]&4)>>2 ;
          int whenAtUB = (sortit[k]&2)>>1 ;
          int affectedToUB = sortit[k]&1 ;
          disaggregationAction action ;
	  action.affected=0 ;
          setAffectedInDisaggregation(action,affected) ;
          setZeroOneInDisaggregation(action,zeroOne!=0) ;
          setWhenAtUBInDisaggregation(action,whenAtUB!=0) ;
          setAffectedToUBInDisaggregation(action,affectedToUB!=0) ;
          if (affected!=affectedLast||zeroOne!=zeroOneLast) {
            // new variable
            thisOne.index[put++]=action ;
          } else if (whenAtUB!=whenAtUBLast||affectedToUB!=affectedToUBLast) {
            // new action - what can we discover
            thisOne.index[put++]=action ;
            int j=cutVector_[i].sequence ;
            int k=affected ;
            if (zeroOne) {
              k=cutVector_[k].sequence ;
              if (logLevel_>1)
                printf("For %d %d 0-1 pair",j,k) ;
            } else {
              if (logLevel_>1)
                printf("For %d %d pair",j,k) ;
            }
            if (logLevel_>1)
              printf(" old whenAtUB, affectedToUB %d %d, new whenAtUB, affectedToUB %d %d\n",
                     whenAtUBLast, affectedToUBLast,whenAtUB, affectedToUB) ;
          }
          affectedLast=affected ;
          zeroOneLast=zeroOne ;
          whenAtUBLast=whenAtUB ;
          affectedToUBLast=affectedToUB ;
        }
        if (put<number) {
          //printf("%d reduced from %d to %d\n",i,number,put) ;
          thisOne.length=put ;
        }
      }
      // And look at all where two 0-1 variables involved
      for (int i=0;i<number01Integers_;i++) {
        disaggregation & thisOne=cutVector_[i] ;
        int k ;
        int number = thisOne.length ;
        for (k=0;k<number;k++) {
          int affected = affectedInDisaggregation(thisOne.index[k]) ;
          bool zeroOne = zeroOneInDisaggregation(thisOne.index[k]) ;
          if (zeroOne && static_cast<int>(affected)>i) {
            bool whenAtUB = whenAtUBInDisaggregation(thisOne.index[k]) ;
            bool affectedToUB = affectedToUBInDisaggregation(thisOne.index[k]) ;
            disaggregation otherOne=cutVector_[affected] ;
            int numberOther = otherOne.length ;
            // Could do binary search if a lot
            int lastAction=-1 ;
            for (int j=0;j<numberOther;j++) {
              if (affectedInDisaggregation(otherOne.index[j])==i) {
                bool whenAtUBOther = whenAtUBInDisaggregation(otherOne.index[j]) ;
                bool affectedToUBOther = affectedToUBInDisaggregation(otherOne.index[j]) ;
                /* action -
                   0 -> x + y <=1 (1,1 impossible)
                   1 -> x - y <=0 (1,0 impossible)
                   2 -> -x + y <=0 (0,1 impossible)
                   3 -> -x -y <= -1 (0,0 impossible)
                  10 -> x == y
                  11 -> x + y == 1
                  20 -> x == 0
                  21 -> x == 1
                  22 -> y == 0
                  23 -> y == 1
                */
                int action=-1 ;
                if (whenAtUB) {
                  if (affectedToUB) {
                    // x -> 1 => y -> 1
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=10; // x,y must be same
                      } else {
                        // y -> 1 => x -> 0 
                        action=20; // If x is 1 then contradiction
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=23; // if y is 0 then contradiction
                      } else {
                        // y -> 0 => x -> 0
                        action=1; // x,y 1,0 impossible 
                      }
                    }
                  } else {
                    // x -> 1 => y -> 0
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=22; // If y is 1 then contradiction
                      } else {
                        // y -> 1 => x -> 0 
                        action=0; 
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=11; // x,y with same values impossible
                      } else {
                        // y -> 0 => x -> 0
                        action=20; // If x is 1 then contradiction
                      }
                    }
                  }
                } else {
                  if (affectedToUB) {
                    // x -> 0 => y -> 1
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=21; // If x is 0 then contradiction
                      } else {
                        // y -> 1 => x -> 0 
                        action=11; // x,y must be different
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=3; // one of x,y must be 1
                      } else {
                        // y -> 0 => x -> 0
                        action=23; // if y is 0 then contradiction
                      }
                    }
                  } else {
                    // x -> 0 => y -> 0
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=2; // x,y 0,1 impossible
                      } else {
                        // y -> 1 => x -> 0 
                        action=22; // If y is 1 then contradiction
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=21; // if x is 0 then contradiction
                      } else {
                        // y -> 0 => x -> 0
                        action=10; // x,y must be same
                      }
                    }
                  }
                }
                assert (action>=0) ;
                if (action<4) {
                  // clique - see if there
                  if (oneFixStart_) {
                    switch (action) {
                    case 0:
                      break ;
                    case 1:
                      break ;
                    case 2:
                      break ;
                    case 3:
                      break ;
                    }
                    // If not can we add or strengthen
                  }
                  // check last action
                  if (lastAction>=0) {
                    if (logLevel_>1)
                      printf("XX lastAction %d, this %d\n",lastAction,action) ;
                  }
                } else if (action<12) {
                  if (logLevel_>1)
                    printf("XX Could eliminate one of %d %d 0-1 variables %c\n",i,affected,
                           (lastAction>=0) ? '*' : ' ') ;
                  if (info->strengthenRow) {
                    OsiRowCut rc ;
                    int index[2] ;
                    double element[2] ;
                    index[0]=cutVector_[i].sequence ;
                    element[0]=1.0 ;
                    index[1]=cutVector_[affected].sequence ;
                    if (action==10) {
                      // 10 -> x == y
                      rc.setLb(0.0) ;
                      rc.setUb(0.0);   
                      element[1]= -1.0 ;
                    } else {
                      // 11 -> x + y == 1
                      rc.setLb(1.0) ;
                      rc.setUb(1.0);   
                      element[1]= 1.0 ;
                    }
                    rc.setRow(2,index,element,false) ;
                    cs.insert(rc) ;
                  }
                } else {
                  if (action<22) {
                    if (logLevel_>1)
                      printf("XX Could fix a 0-1 variable %d\n",i) ;
                  } else {
                    if (logLevel_>1)
                      printf("XX Could fix a 0-1 variable %d\n",affected) ;
		  }
                }
                //printf("%d when %d forces %d to %d , %d when %d forces %d to %d\n",
                //     i,whenAtUB,affected,affectedToUB, 
                //     affected, whenAtUBOther,i, affectedToUBOther) ;
              }
            }
          }
        }
      }
      delete [] sortit ;
    }
    if (cutVector_) {
      // now see if any disaggregation cuts are violated
      for (int i=0;i<number01Integers_;i++) {
	int j=cutVector_[i].sequence ;
	double solInt=colsol[j] ;
	double  upper, solValue ;
	int icol ;
	int index[2] ;
	double element[2] ;
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  double away = fabs(0.5-(solInt-floor(solInt))) ;
	  if (away<0.4999999) {
	    disaggregation thisOne=cutVector_[i] ;
	    int k ;
	    OsiRowCut rc ;
	    for (k=0;k<thisOne.length;k++) {
	      icol = affectedInDisaggregation(thisOne.index[k]) ;
              if (zeroOneInDisaggregation(thisOne.index[k]))
                icol = cutVector_[icol].sequence ;
	      solValue=colsol[icol] ;
	      upper=colUpper_[icol] ;
              double infeasibility=0.0 ;
              if (!whenAtUBInDisaggregation(thisOne.index[k])) {
                if (!affectedToUBInDisaggregation(thisOne.index[k])) {
                  // delta -> 0 => x to lb (at present just 0)
                  infeasibility = solValue - upper * solInt ;
                  if (infeasibility > 1.0e-3) {
                    rc.setLb(-DBL_MAX) ;
                    rc.setUb(0.0) ;
                    index[0]=icol ;
                    element[0]=1.0 ;
                    index[1]=j ;
                    element[1]= -upper ;
                  } else {
                    infeasibility=0.0 ;
                  }
                } else {
                  // delta -> 0 => x to ub
                  abort() ;
                }
              } else {
                if (affectedToUBInDisaggregation(thisOne.index[k])) {
                  // delta -> 1 => x to ub (?)
                  icol = affectedInDisaggregation(thisOne.index[k]) ;
                  if (zeroOneInDisaggregation(thisOne.index[k]))
                    icol = cutVector_[icol].sequence ;
                  solValue=colsol[icol] ;
                  upper=colUpper_[icol] ;
                  if (!colLower[icol]) {
                    infeasibility = upper * solInt - solValue ;
                    if (infeasibility > 1.0e-3) {
                      rc.setLb(-DBL_MAX) ;
                      rc.setUb(0.0) ;
                      index[0]=icol ;
                      element[0]=-1.0 ;
                      index[1]=j ;
                      element[1]= upper ;
                    } else {
                      infeasibility=0.0 ;
                    }
                  } else {
                    assert (upper==colLower[icol]) ;
                    infeasibility=0.0 ;
                  }
                } else {
                  // delta + delta2 <= 1
                  assert (zeroOneInDisaggregation(thisOne.index[k])) ;
                  // delta -> 1 => delta2 -> 0
                  icol = affectedInDisaggregation(thisOne.index[k]) ;
                  icol = cutVector_[icol].sequence ;
                  // only do if icol > j
                  if (icol >j && colUpper[icol] ) {
                    solValue=colsol[icol] ;
                    if (!colLower[icol]) {
                      infeasibility = solInt + solValue - 1.0 ;
                      if (infeasibility > 1.0e-3) {
                        rc.setLb(-DBL_MAX) ;
                        rc.setUb(1.0) ;
                        index[0]=icol ;
                        element[0]=1.0 ;
                        index[1]=j ;
                        element[1]= 1.0 ;
                      } else {
                        infeasibility=0.0 ;
                      }
                    } else {
                      assert (upper==colLower[icol]) ;
                      infeasibility=0.0 ;
                    }
                  }
                }
              }
              if (infeasibility) {
                rc.setEffectiveness(infeasibility) ;
                rc.setRow(2,index,element,false) ;
                if (logLevel_>1)
                  printf("%g <= %g * x%d + %g * x%d <= %g\n",
                         rc.lb(),element[0],index[0],element[1],index[1],rc.ub()) ;
#if CGL_DEBUG > 0
                if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                rowCut.addCutIfNotDuplicate(rc) ;
              }
	    }
	  }
	}
      }
    }
  }
/*
  End probing, mode 0. Time to clean up and go home.
*/
  delete [] markR ;
  delete [] minR ;
  delete [] maxR ;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info->strengthenRow,0) ;
  }
  // delete stuff
  delete rowCopy ;
  delete columnCopy ;
  if (rowCopy_) {
    delete [] rowLower ;
    delete [] rowUpper ;
  }
  delete [] intVar ;
  delete [] rowStartPos ;
  delete [] realRows ;
  // and put back unreasonable bounds on integer variables
  const double * trueLower = si.getColLower() ;
  const double * trueUpper = si.getColUpper() ;
  if (!ninfeas) {
    for (int i=0;i<nCols;i++) {
      if (intVarOriginal[i]==2) {
	if (colUpper[i] == CGL_REASONABLE_INTEGER_BOUND) 
	  colUpper[i] = trueUpper[i] ;
	if (colLower[i] == -CGL_REASONABLE_INTEGER_BOUND) 
	  colLower[i] = trueLower[i] ;
      }
    }
  } else {
    memcpy(colLower,trueLower,nCols*sizeof(double)) ;
    memcpy(colUpper,trueUpper,nCols*sizeof(double)) ;
  }
  if (!info->inTree && ((info->options&4)==4||((info->options&8)&&!info->pass))) {
    int numberRowCutsAfter = cs.sizeRowCuts() ;
    for (int i=numberRowCutsBefore;i<numberRowCutsAfter;i++)
      cs.rowCutPtr(i)->setGloballyValid() ;
  }

# if CGL_DEBUG > 0
  std::cout
    << "Leaving CglProbing::gutsOfGenerateCuts, return "
    << ninfeas << std::endl ;
# endif

  return ninfeas ;
}



// Does probing and adding cuts
int CglProbing::probeCliques( const OsiSolverInterface & si, 
                              const OsiRowCutDebugger *
#if CGL_DEBUG > 0
			 debugger
#endif
                              ,OsiCuts & cs, 
                              double * colLower, double * colUpper, 
		       CoinPackedMatrix *rowCopy,
			      CoinPackedMatrix *columnCopy, const int * realRows,
		       double * rowLower, double * rowUpper,
		       char * intVar, double * minR, double * maxR, 
		       int * markR, 
                       CglTreeInfo * info) const
{
  // Set up maxes
  int maxStack = info->inTree ? maxStack_ : maxStackRoot_;
  int nRows=rowCopy->getNumRows();
  int nCols=rowCopy->getNumCols();
  double * colsol = new double[nCols];
  double * djs = new double[nCols];
  const double * currentColLower = si.getColLower();
  const double * currentColUpper = si.getColUpper();
  double * tempL = new double [nCols];
  double * tempU = new double [nCols];
  int * markC = new int [nCols];
  int * stackC = new int [2*nCols];
  int * stackR = new int [nRows];
  double * saveL = new double [2*nCols];
  double * saveU = new double [2*nCols];
  double * saveMin = new double [nRows];
  double * saveMax = new double [nRows];
  double * element = new double[nCols];
  int * index = new int[nCols];
  // For trying to extend cliques
  int * cliqueStack=NULL;
  int * cliqueCount=NULL;
  int * to_01=NULL;
  if (!mode_) {
    to_01 = new int[nCols];
    cliqueStack = new int[numberCliques_];
    cliqueCount = new int[numberCliques_];
    int i;
    for (i=0;i<numberCliques_;i++) {
      cliqueCount[i]=cliqueStart_[i+1]-cliqueStart_[i];
    }
    for (i=0;i<nCols;i++) 
      to_01[i]=-1;
    for (i=0;i<number01Integers_;i++) {
      int j=cutVector_[i].sequence;
      to_01[j]=i;
    }
  }
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info->inTree ? nRows/3 : nRows;
  CglProbingRowCut rowCut(nRowsFake, !info->inTree);
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements();
  double movement;
  int i, j, k,kk,jj;
  int kcol,krow;
  bool anyColumnCuts=false;
  double dbound, value, value2;
  int ninfeas=0;
  int rowCuts;
  double disaggEffectiveness;
  if (mode_) {
    /* clean up djs and solution */
    CoinMemcpyN(si.getReducedCost(),nCols,djs);
    CoinMemcpyN( si.getColSolution(),nCols,colsol);
    disaggEffectiveness=1.0e-3;
    rowCuts=rowCuts_;
  } else {
    // need to go from a neutral place
    memset(djs,0,nCols*sizeof(double));
    CoinMemcpyN( si.getColSolution(),nCols,colsol);
    disaggEffectiveness=-1.0e10;
    if (rowCuts_!=4)
      rowCuts=1;
    else
      rowCuts=4;
  }
  for (i = 0; i < nCols; ++i) {
    /* was if (intVar[i]) */
    if (1) {
      if (colUpper[i]-colLower[i]>1.0e-8) {
	if (colsol[i]<colLower[i]+primalTolerance_) {
	  colsol[i]=colLower[i];
	  djs[i] = CoinMax(0.0,djs[i]);
	} else if (colsol[i]>colUpper[i]-primalTolerance_) {
	  colsol[i]=colUpper[i];
	  djs[i] = CoinMin(0.0,djs[i]);
	} else {
	  djs[i]=0.0;
	}
	/*if (fabs(djs[i])<1.0e-5) 
	  djs[i]=0.0;*/
      }
    }
  }

  int ipass=0,nfixed=-1;

  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff);
  if (!cutoff_available||usingObjective_<0) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  cutoff *= si.getObjSense();
  if (fabs(cutoff)>1.0e30)
    assert (cutoff>1.0e30);
  double current = si.getObjValue();
  // make irrelevant if mode is 0
  if (!mode_)
    cutoff=DBL_MAX;
  /* for both way coding */
  int nstackC0=-1;
  int * stackC0 = new int[maxStack];
  double * lo0 = new double[maxStack];
  double * up0 = new double[maxStack];
  int nstackR,nstackC;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3;
    } else {
      markC[i]=0;
    }
  }
  double tolerance = 1.0e1*primalTolerance_;
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_;
  // If we are going to replace coefficient then we don't need to be effective
  double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3;
  while (ipass<maxPass&&nfixed) {
    int iLook;
    ipass++;
    nfixed=0;
    for (iLook=0;iLook<numberThisTime_;iLook++) {
      double solval;
      double down;
      double up;
      int awayFromBound=1;
      j=lookedAt_[iLook];
      solval=colsol[j];
      down = floor(solval+tolerance);
      up = ceil(solval-tolerance);
      if(colUpper[j]-colLower[j]<1.0e-8) markC[j]=3;
      if (markC[j]||!intVar[j]) continue;
      double saveSolval = solval;
      if (solval>=colUpper[j]-tolerance||solval<=colLower[j]+tolerance||up==down) {
	awayFromBound=0;
	if (solval<=colLower[j]+2.0*tolerance) {
	  solval = colLower[j]+1.0e-1;
	  down=colLower[j];
	  up=down+1.0;
	} else if (solval>=colUpper[j]-2.0*tolerance) {
	  solval = colUpper[j]-1.0e-1;
	  up=colUpper[j];
	  down=up-1;
	} else {
          // odd
          up=down+1.0;
          solval = down+1.0e-1;
        }
      }
      assert (up<=colUpper[j]);
      assert (down>=colLower[j]);
      assert (up>down);
      if ((solval-down>1.0e-6&&up-solval>1.0e-6)||mode_!=1) {
	int istackC,iway, istackR;
	int way[]={1,2,1};
	int feas[]={1,2,4};
	int feasible=0;
	int notFeasible;
	for (iway=0;iway<3;iway ++) {
	  int fixThis=0;
	  double objVal=current;
	  int goingToTrueBound=0;
	  stackC[0]=j;
	  markC[j]=way[iway];
          double solMovement;
	  if (way[iway]==1) {
	    movement=down-colUpper[j];
            solMovement = down-colsol[j];
	    assert(movement<-0.99999);
	    if (fabs(down-colLower[j])<1.0e-7) {
	      goingToTrueBound=2;
	      down=colLower[j];
	    }
	  } else {
	    movement=up-colLower[j];
            solMovement = up-colsol[j];
	    assert(movement>0.99999);
	    if (fabs(up-colUpper[j])<1.0e-7) {
	      goingToTrueBound=2;
	      up=colUpper[j];
	    }
	  }
	  if (goingToTrueBound&&(colUpper[j]-colLower[j]>1.5||colLower[j]))
	    goingToTrueBound=1;
	  // switch off disaggregation if not wanted
	  if ((rowCuts&1)==0)
	    goingToTrueBound=0;
#ifdef PRINT_DEBUG
	  if (fabs(movement)>1.01) {
	    printf("big %d %g %g %g\n",j,colLower[j],solval,colUpper[j]);
	  }
#endif
	  if (solMovement*djs[j]>0.0)
	    objVal += solMovement*djs[j];
	  nstackC=1;
	  nstackR=0;
	  saveL[0]=colLower[j];
	  saveU[0]=colUpper[j];
          assert (saveU[0]>saveL[0]);
	  notFeasible=0;
	  if (movement<0.0) {
	    colUpper[j] += movement;
	    colUpper[j] = floor(colUpper[j]+0.5);
#ifdef PRINT_DEBUG
	    printf("** Trying %d down to 0\n",j);
#endif
	  } else {
	    colLower[j] += movement;
	    colLower[j] = floor(colLower[j]+0.5);
#ifdef PRINT_DEBUG
	    printf("** Trying %d up to 1\n",j);
#endif
	  }
	  if (fabs(colUpper[j]-colLower[j])<1.0e-6) 
	    markC[j]=3; // say fixed
	  istackC=0;
	  /* update immediately */
	  for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
	    int irow = row[k];
	    value = columnElements[k];
	    assert (markR[irow]!=-2);
	    if (markR[irow]==-1) {
	      stackR[nstackR]=irow;
	      markR[irow]=nstackR;
	      saveMin[nstackR]=minR[irow];
	      saveMax[nstackR]=maxR[irow];
	      nstackR++;
#if 0
	    } else if (markR[irow]==-2) {
	      continue;
#endif
	    }
	    /* could check immediately if violation */
	    if (movement>0.0) {
	      /* up */
	      if (value>0.0) {
		/* up does not change - down does */
                if (minR[irow]>-1.0e10)
                  minR[irow] += value;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      } else {
		/* down does not change - up does */
                if (maxR[irow]<1.0e10)
                  maxR[irow] += value;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      }
	    } else {
	      /* down */
	      if (value<0.0) {
		/* up does not change - down does */
                if (minR[irow]>-1.0e10)
                  minR[irow] -= value;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      } else {
		/* down does not change - up does */
                if (maxR[irow]<1.0e10)
                  maxR[irow] -= value;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      }
	    }
	  }
	  while (istackC<nstackC&&nstackC<maxStack) {
	    int jway;
	    int jcol =stackC[istackC];
	    jway=markC[jcol];
	    // If not first and fixed then skip
	    if (jway==3&&istackC) {
	      //istackC++;
	      //continue;
              //printf("fixed %d on stack\n",jcol);
	    }
	    // Do cliques
	    if (oneFixStart_&&oneFixStart_[jcol]>=0) {
	      int start;
	      int end;
	      if (colLower[jcol]>saveL[istackC]) {
		// going up
		start = oneFixStart_[jcol];
		end = zeroFixStart_[jcol];
	      } else {
		assert (colUpper[jcol]<saveU[istackC]);
		// going down
		start = zeroFixStart_[jcol];
		end = endFixStart_[jcol];
	      }
	      for (int i=start;i<end;i++) {
		int iClique = whichClique_[i];
		for (int k=cliqueStart_[iClique];k<cliqueStart_[iClique+1];k++) {
		  int kcol = sequenceInCliqueEntry(cliqueEntry_[k]);
                  if (jcol==kcol)
                    continue;
		  int kway = oneFixesInCliqueEntry(cliqueEntry_[k]);
                  if (kcol!=jcol) {
                    if (!markC[kcol]) {
                      // not on list yet
                      if (nstackC<2*maxStack) {
                        markC[kcol] = 3; // say fixed
                        fixThis++;
                        stackC[nstackC]=kcol;
                        saveL[nstackC]=colLower[kcol];
                        saveU[nstackC]=colUpper[kcol];
                        assert (saveU[nstackC]>saveL[nstackC]);
                        nstackC++;
                        if (!kway) {
                          // going up
                          double solMovement=1.0-colsol[kcol];
                          if (solMovement>0.0001) {
                            assert (djs[kcol]>=0.0);
                            objVal += djs[kcol]*solMovement;
                          }
                          colLower[kcol]=1.0;
                          /* update immediately */
                          for (int jj =columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                            krow = row[jj];
                            value = columnElements[jj];
			    assert (markR[krow]!=-2);
                            if (markR[krow]==-1) {
                              stackR[nstackR]=krow;
                              markR[krow]=nstackR;
                              saveMin[nstackR]=minR[krow];
                              saveMax[nstackR]=maxR[krow];
                              nstackR++;
#if 0
                            } else if (markR[krow]==-2) {
                              continue;
#endif
                            }
                            /* could check immediately if violation */
                            /* up */
                            if (value>0.0) {
                              /* up does not change - down does */
                              if (minR[krow]>-1.0e10)
                                minR[krow] += value;
                              if (minR[krow]>rowUpper[krow]+1.0e-5) {
                                colUpper[kcol]=-1.0e50; /* force infeasible */
                                break;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow]<1.0e10)
                                maxR[krow] += value;
                              if (maxR[krow]<rowLower[krow]-1.0e-5) {
                                notFeasible=1;
                                break;
                              }
                            }
                          }
                        } else {
                          // going down
                          double solMovement=0.0-colsol[kcol];
                          if (solMovement<-0.0001) {
                            assert (djs[kcol]<=0.0);
                            objVal += djs[kcol]*solMovement;
                          }
                          colUpper[kcol]=0.0;
                          /* update immediately */
                          for (int jj =columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                            krow = row[jj];
                            value = columnElements[jj];
			    assert (markR[krow]!=-2);
                            if (markR[krow]==-1) {
                              stackR[nstackR]=krow;
                              markR[krow]=nstackR;
                              saveMin[nstackR]=minR[krow];
                              saveMax[nstackR]=maxR[krow];
                              nstackR++;
#if 0
                            } else if (markR[krow]==-2) {
                              continue;
#endif
                            }
                            /* could check immediately if violation */
                            /* down */
                            if (value<0.0) {
                              /* up does not change - down does */
                              if (minR[krow]>-1.0e10)
                                minR[krow] -= value;
                              if (minR[krow]>rowUpper[krow]+1.0e-5) {
                                notFeasible=1;
                                break;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow]<1.0e10)
                                maxR[krow] -= value;
                              if (maxR[krow]<rowLower[krow]-1.0e-5) {
                                notFeasible=1;
                                break;
                              }
                            }
                          }
                        }
                      }
                    } else if (markC[kcol]==1) {
                      // marked as going to 0
                      assert (!colUpper[kcol]);
                      if (!kway) {
                        // contradiction
                        notFeasible=1;
                        break;
                      }
                    } else if (markC[kcol]==2) {
                      // marked as going to 1
                      assert (colLower[kcol]);
                      if (kway) {
                        // contradiction
                        notFeasible=1;
                        break;
                      }
                    } else {
                      // marked as fixed
                      assert (markC[kcol]==3);
                      int jkway;
                      if (colLower[kcol])
                        jkway=1;
                      else
                        jkway=0;
                      if (kway==jkway) {
                        // contradiction
                        notFeasible=1;
                        break;
                      }
                    }
                  }
		}
		if (notFeasible)
		  break;
	      }
	      if (notFeasible)
		istackC=nstackC+1;
	    }
	    for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
	      // break if found not feasible
	      if (notFeasible)
		break;
	      int irow = row[k];
	      /*value = columnElements[k];*/
	      assert (markR[irow]!=-2);
#if 0
	      if (markR[irow]!=-2) {
#endif
		/* see if anything forced */
		for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
		  double moveUp=0.0;
		  double moveDown=0.0;
		  double newUpper=-1.0,newLower=1.0;
		  kcol=column[kk];
		  bool onList = (markC[kcol]!=0);
		  if (markC[kcol]!=3) {
		    value2=rowElements[kk];
                    int markIt=markC[kcol];
		    if (value2 < 0.0) {
		      if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colUpper[kcol]+
			  (rowUpper[irow]-minR[irow])/value2;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
                            markIt |= 2;
			    newLower = ceil(dbound-primalTolerance_);
			  } else {
			    newLower=dbound;
			    if (newLower+primalTolerance_>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol];
                              markIt |= 2;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
                            }
			  }
			  moveUp = newLower-colLower[kcol];
			}
		      }
		      if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colLower[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 1;
			    newUpper = floor(dbound+primalTolerance_);
			  } else {
			    newUpper=dbound;
			    if (newUpper-primalTolerance_<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol];
                              markIt |= 1;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
                            }
			  }
			  moveDown = newUpper-colUpper[kcol];
			}
		      }
		    } else {
		      /* positive element */
		      if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colUpper[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 2;
			    newLower = ceil(dbound-primalTolerance_);
			  } else {
			    newLower=dbound;
			    if (newLower+primalTolerance_>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol];
                              markIt |= 2;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
			    }
			  }
			  moveUp = newLower-colLower[kcol];
			}
		      }
		      if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colLower[kcol] + 
			  (rowUpper[irow]-minR[irow])/value2;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 1;
			    newUpper = floor(dbound+primalTolerance_);
			  } else {
			    newUpper=dbound;
			    if (newUpper-primalTolerance_<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol];
                              markIt |= 1;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
			    }
			  }
			  moveDown = newUpper-colUpper[kcol];
			}
		      }
		    }
		    if (nstackC<2*maxStack) {
                      markC[kcol] = markIt;
		    }
		    if (moveUp&&nstackC<2*maxStack) {
		      fixThis++;
#ifdef PRINT_DEBUG
		      printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]);
		      value=0.0;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]);
			} else {
			  value += rowElements[jj]*colLower[ii];
			}
		      }
		      printf(" - fixed %g\n",value);
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
			}
		      }
		      printf("\n");
#endif
		      if (!onList) {
			stackC[nstackC]=kcol;
			saveL[nstackC]=colLower[kcol];
			saveU[nstackC]=colUpper[kcol];
                        assert (saveU[nstackC]>saveL[nstackC]);
			nstackC++;
			onList=true;
		      }
		      if (newLower>colsol[kcol]) {
			if (djs[kcol]<0.0) {
			  /* should be infeasible */
			  assert (newLower>colUpper[kcol]+primalTolerance_);
			} else {
			  objVal += moveUp*djs[kcol];
			}
		      }
		      if (intVar[kcol])
			newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4));
		      colLower[kcol]=newLower;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]=3; // say fixed
		      }
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj];
			value = columnElements[jj];
			assert (markR[krow]!=-2);
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow;
			  markR[krow]=nstackR;
			  saveMin[nstackR]=minR[krow];
			  saveMax[nstackR]=maxR[krow];
			  nstackR++;
#if 0
			} else if (markR[krow]==-2) {
			  continue;
#endif
			}
			/* could check immediately if violation */
			/* up */
			if (value>0.0) {
			  /* up does not change - down does */
                          if (minR[krow]>-1.0e10)
                            minR[krow] += value*moveUp;
			  if (minR[krow]>rowUpper[krow]+1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			} else {
			  /* down does not change - up does */
                          if (maxR[krow]<1.0e10)
                            maxR[krow] += value*moveUp;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			}
		      }
		    }
		    if (moveDown&&nstackC<2*maxStack) {
		      fixThis++;
#ifdef PRINT_DEBUG
		      printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]);
		      value=0.0;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]);
			} else {
			  value += rowElements[jj]*colLower[ii];
			}
		      }
		      printf(" - fixed %g\n",value);
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
			}
		      }
		      printf("\n");
#endif
		      if (!onList) {
			stackC[nstackC]=kcol;
			saveL[nstackC]=colLower[kcol];
			saveU[nstackC]=colUpper[kcol];
                        assert (saveU[nstackC]>saveL[nstackC]);
			nstackC++;
			onList=true;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_);
			} else {
			  objVal += moveDown*djs[kcol];
			}
		      }
		      if (intVar[kcol])
			newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4));
		      colUpper[kcol]=newUpper;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]=3; // say fixed
		      }
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj];
			value = columnElements[jj];
			assert (markR[krow]!=-2);
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow;
			  markR[krow]=nstackR;
			  saveMin[nstackR]=minR[krow];
			  saveMax[nstackR]=maxR[krow];
			  nstackR++;
#if 0
			} else if (markR[krow]==-2) {
#endif
			  continue;
			}
			/* could check immediately if violation */
			/* down */
			if (value<0.0) {
			  /* up does not change - down does */
                          if (minR[krow]>-1.0e10)
                            minR[krow] += value*moveDown;
			  if (minR[krow]>rowUpper[krow]+1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			} else {
			  /* down does not change - up does */
                          if (maxR[krow]<1.0e10)
                            maxR[krow] += value*moveDown;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			}
		      }
		    }
		    if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
		      notFeasible=1;;
		      k=columnStart[jcol]+columnLength[jcol];
		      istackC=nstackC+1;
#ifdef PRINT_DEBUG
		      printf("** not feasible this way\n");
#endif
		      break;
		    }
		  }
		}
#if 0
	      }
#endif
	    }
	    istackC++;
	  }
	  if (!notFeasible) {
	    if (objVal<=cutoff) {
	      feasible |= feas[iway];
	    } else {
#ifdef PRINT_DEBUG
	      printf("not feasible on dj\n");
#endif
	      notFeasible=1;
	      if (iway==1&&feasible==0) {
		/* not feasible at all */
		ninfeas=1;
		j=nCols-1;
		break;
	      }
	    }
	  } else if (iway==1&&feasible==0) {
	    /* not feasible at all */
	    ninfeas=1;
	    j=nCols-1;
            iLook=numberThisTime_;
	    ipass=maxPass;
	    break;
	  }
	  if (notFeasible)
	    goingToTrueBound=0;
	  if (iway==2||(iway==1&&feasible==2)) {
	    /* keep */
	    iway=3;
	    nfixed++;
            if (mode_) {
	      OsiColCut cc;
	      int nTot=0,nFix=0,nInt=0;
	      bool ifCut=false;
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC];
		if (intVar[icol]) {
		  if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
		    element[nFix]=colUpper[icol];
		    index[nFix++]=icol;
		    nInt++;
		    if (colsol[icol]>colUpper[icol]+primalTolerance_) {
		      ifCut=true;
		      anyColumnCuts=true;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot=nFix;
		cc.setUbs(nFix,index,element);
		nFix=0;
	      }
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC];
		if (intVar[icol]) {
		  if (colLower[icol]>currentColLower[icol]+1.0e-4) {
		    element[nFix]=colLower[icol];
		    index[nFix++]=icol;
		    nInt++;
		    if (colsol[icol]<colLower[icol]-primalTolerance_) {
		      ifCut=true;
		      anyColumnCuts=true;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot+=nFix;
		cc.setLbs(nFix,index,element);
	      }
	      // could tighten continuous as well
	      if (nInt) {
		if (ifCut) {
		  cc.setEffectiveness(100.0);
		} else {
		  cc.setEffectiveness(1.0e-5);
		}
#if CGL_DEBUG > 0
		checkBounds(si,cc);
#endif
		cs.insert(cc);
	      }
	    }
	    for (istackC=0;istackC<nstackC;istackC++) {
	      int icol=stackC[istackC];
	      if (colUpper[icol]-colLower[icol]>primalTolerance_) {
		markC[icol]=0;
	      } else {
		markC[icol]=3;
	      }
	    }
	    for (istackR=0;istackR<nstackR;istackR++) {
	      int irow=stackR[istackR];
	      markR[irow]=-1;
	    }
	  } else {
	    /* is it worth seeing if can increase coefficients
	       or maybe better see if it is a cut */
	    if (iway==0) {
	      nstackC0=CoinMin(nstackC,maxStack);
	      double solMove = saveSolval-down;
	      double boundChange;
	      if (notFeasible) {
		nstackC0=0;
	      } else {
		for (istackC=0;istackC<nstackC0;istackC++) {
		  int icol=stackC[istackC];
		  stackC0[istackC]=icol;
		  lo0[istackC]=colLower[icol];
		  up0[istackC]=colUpper[icol];
		}
	      }
	      /* restore all */
              int nCliquesAffected=0;
              assert (iway==0);
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC];
		double oldU=saveU[istackC];
		double oldL=saveL[istackC];
		if(goingToTrueBound==2&&istackC) {
                  // Work for extending cliques
                  if (!mode_&&numberCliques_) {
                    int i_01 = to_01[icol];
                    if (i_01>=0) {
                      int start;
                      int end;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol];
                        end = endFixStart_[icol];
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol];
                        end = zeroFixStart_[icol];
                      }
                      //if (end>start)
                      //printf("j %d, other %d is in %d cliques\n",
                      //     j,i_01,end-start);
                      for (int i=start;i<end;i++) {
                        int iClique = whichClique_[i];
                        int size = cliqueStart_[iClique+1]-cliqueStart_[iClique];
                        if (cliqueCount[iClique]==size) {
                          // first time
                          cliqueStack[nCliquesAffected++]=iClique;
                        }
                        // decrement counts
                        cliqueCount[iClique]--;
                      }
                    }
                  }
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (jval-down)
		  boundChange = oldU-colUpper[icol];
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(-DBL_MAX);
		    rc.setUb(colUpper[icol]-down*boundChange);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= - boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange;
		    if (mode_) 
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false);
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (jval-down)
		  boundChange = oldL-colLower[icol];
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(colLower[icol]-down*boundChange);
		    rc.setUb(DBL_MAX);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]=- boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false);
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
#if 0
		      printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,oldL,oldU,colLower[icol],colsol[icol], j, colsol[j]);
		      printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
			     icol,boundChange,j,rc.ub());
#endif
		    }
		  }
		}
		colUpper[icol]=oldU;
		colLower[icol]=oldL;
		markC[icol]=0;
	      }
              if (nCliquesAffected) {
                for (int i=0;i<nCliquesAffected;i++) {
                  int iClique = cliqueStack[i];
                  int size = cliqueCount[iClique];
                  // restore
                  cliqueCount[iClique]= cliqueStart_[iClique+1]-cliqueStart_[iClique];
                  if (!size) {
                    if (logLevel_>1)
                      printf("** could extend clique by adding j!\n");
                  }
                }
              }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR];
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts;
		  double gap = rowUpper[irow]-maxR[irow];
		  double sum=0.0;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]];
		    }
		    if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_||(info->strengthenRow&&rowLower[irow]<-1.0e20)) {
		      // can be a cut
		      // subtract gap from upper and integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL);
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        // If strengthenRow point to row
                        //if(info->strengthenRow)
                        //printf("a point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      if (rowLower[irow]<-1.0e20) {
			printf("5Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
		      int realRow = (rowLower[irow]<-1.0e20) ? irow : -1;
		      if (realRows&&realRow>0)
			realRow=realRows[realRow];
		      rowCut.addCutIfNotDuplicate(rc,realRow);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum+gap*colsol[j]<minR[irow]+primalTolerance_||(info->strengthenRow&&rowUpper[irow]>1.0e20)) {
		      // can be a cut
		      // add gap to lower and integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL);
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("b point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      if (rowUpper[irow]>1.0e20) {
			printf("6Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
		      int realRow = (rowUpper[irow]>1.0e20) ? irow : -1;
		      if (realRows&&realRow>0)
			realRow=realRows[realRow];
		      rowCut.addCutIfNotDuplicate(rc,realRow);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    } else {
	      if (iway==1&&feasible==3) {
		iway=3;
		/* point back to stack */
		for (istackC=nstackC-1;istackC>=0;istackC--) {
		  int icol=stackC[istackC];
		  markC[icol]=istackC+1000;
		}
		if (mode_) {
		  OsiColCut cc;
		  int nTot=0,nFix=0,nInt=0;
		  bool ifCut=false;
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC];
		    int istackC1=markC[icol]-1000;
		    if (istackC1>=0) {
		      if (CoinMin(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
			saveL[istackC1]=CoinMin(lo0[istackC],colLower[icol]);
			if (intVar[icol]) {
			  element[nFix]=saveL[istackC1];
			  index[nFix++]=icol;
			  nInt++;
			  if (colsol[icol]<saveL[istackC1]-primalTolerance_)
			    ifCut=true;
			}
			nfixed++;
		      }
		    }
		  }
		  if (nFix) {
		    nTot=nFix;
		    cc.setLbs(nFix,index,element);
		    nFix=0;
		  }
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC];
		    int istackC1=markC[icol]-1000;
		    if (istackC1>=0) {
		      if (CoinMax(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
			saveU[istackC1]=CoinMax(up0[istackC],colUpper[icol]);
			if (intVar[icol]) {
			  element[nFix]=saveU[istackC1];
			  index[nFix++]=icol;
			  nInt++;
			  if (colsol[icol]>saveU[istackC1]+primalTolerance_)
			    ifCut=true;
			}
			nfixed++;
		      }
		    }
		  }
		  if (nFix) {
		    nTot+=nFix;
		    cc.setUbs(nFix,index,element);
		  }
		  // could tighten continuous as well
		  if (nInt) {
		    if (ifCut) {
		      cc.setEffectiveness(100.0);
		    } else {
		      cc.setEffectiveness(1.0e-5);
		    }
#if CGL_DEBUG > 0
		    checkBounds(si,cc);
#endif
		    cs.insert(cc);
		  }
		}
	      } else {
		goingToTrueBound=0;
	      }
	      double solMove = up-saveSolval;
	      double boundChange;
	      /* restore all */
              int nCliquesAffected=0;
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC];
		double oldU=saveU[istackC];
		double oldL=saveL[istackC];
		if(goingToTrueBound==2&&istackC) {
                  // Work for extending cliques
                  if (!mode_&&numberCliques_&&iway==3) {
                    int i_01 = to_01[icol];
                    if (i_01>=0) {
                      int start;
                      int end;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol];
                        end = endFixStart_[icol];
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol];
                        end = zeroFixStart_[icol];
                      }
                      //if (end>start)
                      //printf("up j %d, other %d is in %d cliques\n",
                      //     j,i_01,end-start);
                      for (int i=start;i<end;i++) {
                        int iClique = whichClique_[i];
                        int size = cliqueStart_[iClique+1]-cliqueStart_[iClique];
                        if (cliqueCount[iClique]==size) {
                          // first time
                          cliqueStack[nCliquesAffected++]=iClique;
                        }
                        // decrement counts
                        cliqueCount[iClique]--;
                      }
                    }
                  }
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (up-jval)
		  boundChange = oldU-colUpper[icol];
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(-DBL_MAX);
		    rc.setUb(colUpper[icol]+up*boundChange);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= + boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false);
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (up-jval)
		  boundChange = oldL-colLower[icol];
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(colLower[icol]+up*boundChange);
		    rc.setUb(DBL_MAX);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= + boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false);
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
		    }
		  }
		}
		colUpper[icol]=oldU;
		colLower[icol]=oldL;
                if (oldU>oldL+1.0e-4)
                  markC[icol]=0;
                else
                  markC[icol]=3;
	      }
              if (nCliquesAffected) {
                for (int i=0;i<nCliquesAffected;i++) {
                  int iClique = cliqueStack[i];
                  int size = cliqueCount[iClique];
                  // restore
                  cliqueCount[iClique]= cliqueStart_[iClique+1]-cliqueStart_[iClique];
                  if (!size) {
                    if (logLevel_>1)
                      printf("** could extend clique by adding j!\n");
                  }
                }
              }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR];
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts;
		  double gap = rowUpper[irow]-maxR[irow];
		  double sum=0.0;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]];
		    }
		    if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_||(info->strengthenRow&&rowLower[irow]<-1.0e20)) {
		      // can be a cut
		      // add gap to integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL);
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("c point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      if (rowLower[irow]<-1.0e20) {
			printf("7Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
		      int realRow = (rowLower[irow]<-1.0e20) ? irow : -1;
		      if (realRows&&realRow>0)
			realRow=realRows[realRow];
			rowCut.addCutIfNotDuplicate(rc,realRow);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum-gap*colsol[j]<rowLower[irow]+primalTolerance_||(info->strengthenRow&&rowUpper[irow]>1.0e20)) {
		      // can be a cut
		      // subtract gap from integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL);
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]-gap*(colUpper[j]-1));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("d point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      if (rowUpper[irow]>1.0e20) {
			printf("8Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
		      int realRow = (rowUpper[irow]>1.0e20) ? irow : -1;
		      if (realRows&&realRow>0)
			realRow=realRows[realRow];
			rowCut.addCutIfNotDuplicate(rc,realRow);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    }
	  }
	}
      }
    }
  }
  delete [] cliqueStack;
  delete [] cliqueCount;
  delete [] to_01;
  delete [] stackC0;
  delete [] lo0;
  delete [] up0;
  delete [] tempL;
  delete [] tempU;
  delete [] markC;
  delete [] stackC;
  delete [] stackR;
  delete [] saveL;
  delete [] saveU;
  delete [] saveMin;
  delete [] saveMax;
  delete [] index;
  delete [] element;
  delete [] djs;
  delete [] colsol;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info->strengthenRow,0);
  }
  return (ninfeas);
}


// Does probing and adding cuts for clique slacks
int 
CglProbing::probeSlacks( const OsiSolverInterface & si, 
                          const OsiRowCutDebugger * 
#if CGL_DEBUG > 0
			 debugger
#endif
			 ,OsiCuts & cs, 
                          double * colLower, double * colUpper, CoinPackedMatrix *rowCopy,
			 CoinPackedMatrix *columnCopy,
                          double * rowLower, double * rowUpper,
                          char * intVar, double * minR, double * maxR,int * markR,
                          CglTreeInfo * info) const
{
  if (!numberCliques_)
    return 0;
  // Set up maxes
  int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_;
  int maxStack = info->inTree ? maxStack_ : maxStackRoot_;
  int nRows=rowCopy->getNumRows();
  int nCols=rowCopy->getNumCols();
  double * colsol = new double[nCols];
  CoinMemcpyN( si.getColSolution(),nCols,colsol);
  int rowCuts=rowCuts_;
  double_int_pair * array = new double_int_pair [numberCliques_];
  // look at <= cliques
  int iClique;
  int nLook=0;
  for (iClique=0;iClique<numberCliques_;iClique++) {
    if (!cliqueType_[iClique].equality) {
      double sum=0.0;
      for (int j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        double value = colsol[iColumn];
        if (oneFixesInCliqueEntry(cliqueEntry_[j]))
          sum += value;
        else
          sum -= value;
      }
      double away = fabs(0.5-(sum-floor(sum)));
      if (away<0.49999) {
        array[nLook].infeasibility=away;
        array[nLook++].sequence=iClique;
      }
    }
  }
  std::sort(array,array+nLook,double_int_pair_compare());
  nLook=CoinMin(nLook,maxProbe);
  const double * currentColLower = si.getColLower();
  const double * currentColUpper = si.getColUpper();
  double * tempL = new double [nCols];
  double * tempU = new double [nCols];
  int * markC = new int [nCols];
  int * stackC = new int [2*nCols];
  int * stackR = new int [nRows];
  double * saveL = new double [2*nCols];
  double * saveU = new double [2*nCols];
  double * saveMin = new double [nRows];
  double * saveMax = new double [nRows];
  double * element = new double[nCols];
  int * index = new int[nCols];
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info->inTree ? nRows/3 : nRows;
  CglProbingRowCut rowCut(nRowsFake, !info->inTree);
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements();
  double movement;
  int i, j, k,kk,jj;
  int kcol,irow,krow;
  bool anyColumnCuts=false;
  double dbound, value, value2;
  int ninfeas=0;
  for (i = 0; i < nCols; ++i) {
    if (colUpper[i]-colLower[i]>1.0e-8) {
      if (colsol[i]<colLower[i]+primalTolerance_) {
        colsol[i]=colLower[i];
      } else if (colsol[i]>colUpper[i]-primalTolerance_) {
        colsol[i]=colUpper[i];
      }
    }
  }

  int ipass=0,nfixed=-1;

  /* for both way coding */
  int nstackC0=-1;
  int * stackC0 = new int[maxStack];
  double * lo0 = new double[maxStack];
  double * up0 = new double[maxStack];
  int nstackR,nstackC;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3;
    } else {
      markC[i]=0;
    }
  }
  double tolerance = 1.0e1*primalTolerance_;
  // If we are going to replace coefficient then we don't need to be effective
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_;
  double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3;
  while (ipass<maxPass&&nfixed) {
    int iLook;
    ipass++;
    nfixed=0;
    for (iLook=0;iLook<nLook;iLook++) {
      double solval;
      double down;
      double up;
      int iClique=array[iLook].sequence;
      solval=0.0;
      j=0;
      for (j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        double value = colsol[iColumn];
        if (oneFixesInCliqueEntry(cliqueEntry_[j]))
          solval += value;
        else
          solval -= value;
      }
      down = floor(solval+tolerance);
      up = ceil(solval-tolerance);
      int istackC,iway, istackR;
      int way[]={1,2,1};
      int feas[]={1,2,4};
      int feasible=0;
      int notFeasible;
      for (iway=0;iway<3;iway ++) {
        int fixThis=0;
        stackC[0]=j;
        markC[j]=way[iway];
        if (way[iway]==1) {
          movement=down-colUpper[j];
          assert(movement<-0.99999);
          down=colLower[j];
        } else {
          movement=up-colLower[j];
          assert(movement>0.99999);
          up=colUpper[j];
        }
        nstackC=1;
        nstackR=0;
        saveL[0]=colLower[j];
        saveU[0]=colUpper[j];
        assert (saveU[0]>saveL[0]);
        notFeasible=0;
        if (movement<0.0) {
          colUpper[j] += movement;
          colUpper[j] = floor(colUpper[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d down to 0\n",j);
#endif
        } else {
          colLower[j] += movement;
          colLower[j] = floor(colLower[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d up to 1\n",j);
#endif
        }
        if (fabs(colUpper[j]-colLower[j])<1.0e-6) 
          markC[j]=3; // say fixed
        istackC=0;
        /* update immediately */
        for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
          int irow = row[k];
          value = columnElements[k];
          if (markR[irow]==-1) {
            stackR[nstackR]=irow;
            markR[irow]=nstackR;
            saveMin[nstackR]=minR[irow];
            saveMax[nstackR]=maxR[irow];
            nstackR++;
          } else if (markR[irow]==-2) {
            continue;
          }
          /* could check immediately if violation */
          if (movement>0.0) {
            /* up */
            if (value>0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] += value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] += value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          } else {
            /* down */
            if (value<0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] -= value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] -= value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          }
        }
        while (istackC<nstackC&&nstackC<maxStack) {
          int jway;
          int jcol =stackC[istackC];
          jway=markC[jcol];
          // If not first and fixed then skip
          if (jway==3&&istackC) {
            //istackC++;
            //continue;
            //printf("fixed %d on stack\n",jcol);
          }
          // Do cliques
          if (oneFixStart_&&oneFixStart_[jcol]>=0) {
            int start;
            int end;
            if (colLower[jcol]>saveL[istackC]) {
              // going up
              start = oneFixStart_[jcol];
              end = zeroFixStart_[jcol];
            } else {
              assert (colUpper[jcol]<saveU[istackC]);
              // going down
              start = zeroFixStart_[jcol];
              end = endFixStart_[jcol];
            }
            for (int i=start;i<end;i++) {
              int iClique = whichClique_[i];
              for (int k=cliqueStart_[iClique];k<cliqueStart_[iClique+1];k++) {
                int kcol = sequenceInCliqueEntry(cliqueEntry_[k]);
                if (jcol==kcol)
                  continue;
                int kway = oneFixesInCliqueEntry(cliqueEntry_[k]);
                if (kcol!=jcol) {
                  if (!markC[kcol]) {
                    // not on list yet
                    if (nstackC<2*maxStack) {
                      markC[kcol] = 3; // say fixed
                      fixThis++;
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      if (!kway) {
                        // going up
                        colLower[kcol]=1.0;
                        /* update immediately */
                        for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                          krow = row[jj];
                          value = columnElements[jj];
                          if (markR[krow]==-1) {
                            stackR[nstackR]=krow;
                            markR[krow]=nstackR;
                            saveMin[nstackR]=minR[krow];
                            saveMax[nstackR]=maxR[krow];
                            nstackR++;
                          } else if (markR[krow]==-2) {
                            continue;
                          }
                          /* could check immediately if violation */
                          /* up */
                          if (value>0.0) {
                            /* up does not change - down does */
                            if (minR[krow]>-1.0e10)
                              minR[krow] += value;
                            if (minR[krow]>rowUpper[krow]+1.0e-5) {
                              colUpper[kcol]=-1.0e50; /* force infeasible */
                              break;
                            }
                          } else {
                            /* down does not change - up does */
                            if (maxR[krow]<1.0e10)
                              maxR[krow] += value;
                            if (maxR[krow]<rowLower[krow]-1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          }
                        }
                      } else {
                        // going down
                        colUpper[kcol]=0.0;
                        /* update immediately */
                        for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                          krow = row[jj];
                          value = columnElements[jj];
                          if (markR[krow]==-1) {
                            stackR[nstackR]=krow;
                            markR[krow]=nstackR;
                            saveMin[nstackR]=minR[krow];
                            saveMax[nstackR]=maxR[krow];
                            nstackR++;
                          } else if (markR[krow]==-2) {
                            continue;
                          }
                          /* could check immediately if violation */
                          /* down */
                          if (value<0.0) {
                            /* up does not change - down does */
                            if (minR[krow]>-1.0e10)
                              minR[krow] -= value;
                            if (minR[krow]>rowUpper[krow]+1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          } else {
                            /* down does not change - up does */
                            if (maxR[krow]<1.0e10)
                              maxR[krow] -= value;
                            if (maxR[krow]<rowLower[krow]-1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          }
                        }
                      }
                    }
                  } else if (markC[kcol]==1) {
                    // marked as going to 0
                    assert (!colUpper[kcol]);
                    if (!kway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  } else if (markC[kcol]==2) {
                    // marked as going to 1
                    assert (colLower[kcol]);
                    if (kway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  } else {
                    // marked as fixed
                    assert (markC[kcol]==3);
                    int jkway;
                    if (colLower[kcol])
                      jkway=1;
                    else
                      jkway=0;
                    if (kway==jkway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  }
                }
              }
              if (notFeasible)
                break;
            }
            if (notFeasible)
              istackC=nstackC+1;
          }
          for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
            // break if found not feasible
            if (notFeasible)
              break;
            irow = row[k];
            /*value = columnElements[k];*/
            if (markR[irow]!=-2) {
              /* see if anything forced */
              for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
                double moveUp=0.0;
                double moveDown=0.0;
                double newUpper=-1.0,newLower=1.0;
                kcol=column[kk];
                bool onList = (markC[kcol]!=0);
                if (markC[kcol]!=3) {
                  value2=rowElements[kk];
                  int markIt=markC[kcol];
                  if (value2 < 0.0) {
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colUpper[kcol]+
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  } else {
                    /* positive element */
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colUpper[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  }
                  if (nstackC<2*maxStack) {
                    markC[kcol] = markIt;
		  }
                  if (moveUp&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      onList=true;
                    }
                    if (intVar[kcol])
                      newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4));
                    colLower[kcol]=newLower;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
                      markC[kcol]=3; // say fixed
		    }
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* up */
                      if (value>0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveUp;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveUp;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (moveDown&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      onList=true;
                    }
                    if (intVar[kcol])
                      newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4));
                    colUpper[kcol]=newUpper;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
                      markC[kcol]=3; // say fixed
		    }
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* down */
                      if (value<0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveDown;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveDown;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
                    notFeasible=1;;
                    k=columnStart[jcol]+columnLength[jcol];
                    istackC=nstackC+1;
#ifdef PRINT_DEBUG
                    printf("** not feasible this way\n");
#endif
                    break;
                  }
                }
              }
            }
          }
          istackC++;
        }
        if (!notFeasible) {
          feasible |= feas[iway];
        } else if (iway==1&&feasible==0) {
          /* not feasible at all */
          ninfeas=1;
          j=nCols-1;
          iLook=nLook;
          ipass=maxPass;
          break;
        }
        if (iway==2||(iway==1&&feasible==2)) {
          /* keep */
          iway=3;
          nfixed++;
          if (mode_) {
            OsiColCut cc;
            int nTot=0,nFix=0,nInt=0;
            bool ifCut=false;
            for (istackC=0;istackC<nstackC;istackC++) {
              int icol=stackC[istackC];
              if (intVar[icol]) {
                if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
                  element[nFix]=colUpper[icol];
                  index[nFix++]=icol;
                  nInt++;
                  if (colsol[icol]>colUpper[icol]+primalTolerance_) {
                    ifCut=true;
                    anyColumnCuts=true;
                  }
                }
              }
            }
            if (nFix) {
              nTot=nFix;
              cc.setUbs(nFix,index,element);
              nFix=0;
            }
            for (istackC=0;istackC<nstackC;istackC++) {
              int icol=stackC[istackC];
              if (intVar[icol]) {
                if (colLower[icol]>currentColLower[icol]+1.0e-4) {
                  element[nFix]=colLower[icol];
                  index[nFix++]=icol;
                  nInt++;
                  if (colsol[icol]<colLower[icol]-primalTolerance_) {
                    ifCut=true;
                    anyColumnCuts=true;
                  }
                }
              }
            }
            if (nFix) {
              nTot+=nFix;
              cc.setLbs(nFix,index,element);
            }
            // could tighten continuous as well
            if (nInt) {
              if (ifCut) {
                cc.setEffectiveness(100.0);
              } else {
                cc.setEffectiveness(1.0e-5);
              }
#if CGL_DEBUG > 0
              checkBounds(si,cc);
#endif
              cs.insert(cc);
            }
          }
          for (istackC=0;istackC<nstackC;istackC++) {
            int icol=stackC[istackC];
            if (colUpper[icol]-colLower[icol]>primalTolerance_) {
              markC[icol]=0;
            } else {
              markC[icol]=3;
            }
          }
          for (istackR=0;istackR<nstackR;istackR++) {
            int irow=stackR[istackR];
            markR[irow]=-1;
          }
        } else {
          /* is it worth seeing if can increase coefficients
             or maybe better see if it is a cut */
          if (iway==0) {
            nstackC0=CoinMin(nstackC,maxStack);
            if (notFeasible) {
              nstackC0=0;
            } else {
              for (istackC=0;istackC<nstackC0;istackC++) {
                int icol=stackC[istackC];
                stackC0[istackC]=icol;
                lo0[istackC]=colLower[icol];
                up0[istackC]=colUpper[icol];
              }
            }
            /* restore all */
            assert (iway==0);
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              markC[icol]=0;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0) {
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    sum += rowElements[kk]*colsol[column[kk]];
                  }
                  if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_||info->strengthenRow) {
                    // can be a cut
                    // subtract gap from upper and integer coefficient
                    // saveU and saveL spare
                    int * index = reinterpret_cast<int *>(saveL);
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]-gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=-gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(-DBL_MAX);
                    rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0));
                    // effectiveness is how far j moves
                    rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      // If strengthenRow point to row
                      //if(info->strengthenRow)
                      //printf("a point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("9Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
                      rowCut.addCutIfNotDuplicate(rc,irow);
                    }
                  }
                }
                gap = minR[irow]-rowLower[irow];
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  if (!sum) {
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      sum += rowElements[kk]*colsol[column[kk]];
                    }
                  }
                  if (sum+gap*colsol[j]<minR[irow]+primalTolerance_||info->strengthenRow) {
                    // can be a cut
                    // add gap to lower and integer coefficient
                    // saveU and saveL spare
                    int * index = reinterpret_cast<int *>(saveL);
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]+gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0));
                    rc.setUb(DBL_MAX);
                    // effectiveness is how far j moves
                    rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      //if(info->strengthenRow)
                      //printf("b point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("10Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
                      rowCut.addCutIfNotDuplicate(rc,irow);
                    }
                  }
                }
              }
              minR[irow]=saveMin[istackR];
              maxR[irow]=saveMax[istackR];
              markR[irow]=-1;
            }
          } else {
            if (iway==1&&feasible==3) {
              iway=3;
              /* point back to stack */
              for (istackC=nstackC-1;istackC>=0;istackC--) {
                int icol=stackC[istackC];
                markC[icol]=istackC+1000;
              }
              if (mode_) {
                OsiColCut cc;
                int nTot=0,nFix=0,nInt=0;
                bool ifCut=false;
                for (istackC=0;istackC<nstackC0;istackC++) {
                  int icol=stackC0[istackC];
                  int istackC1=markC[icol]-1000;
                  if (istackC1>=0) {
                    if (CoinMin(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
                      saveL[istackC1]=CoinMin(lo0[istackC],colLower[icol]);
                      if (intVar[icol]) {
                        element[nFix]=saveL[istackC1];
                        index[nFix++]=icol;
                        nInt++;
                        if (colsol[icol]<saveL[istackC1]-primalTolerance_)
                          ifCut=true;
                      }
                      nfixed++;
                    }
                  }
                }
                if (nFix) {
                  nTot=nFix;
                  cc.setLbs(nFix,index,element);
                  nFix=0;
                }
                for (istackC=0;istackC<nstackC0;istackC++) {
                  int icol=stackC0[istackC];
                  int istackC1=markC[icol]-1000;
                  if (istackC1>=0) {
                    if (CoinMax(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
                      saveU[istackC1]=CoinMax(up0[istackC],colUpper[icol]);
                      if (intVar[icol]) {
                        element[nFix]=saveU[istackC1];
                        index[nFix++]=icol;
                        nInt++;
                        if (colsol[icol]>saveU[istackC1]+primalTolerance_)
                          ifCut=true;
                      }
                      nfixed++;
                    }
                  }
                }
                if (nFix) {
                  nTot+=nFix;
                  cc.setUbs(nFix,index,element);
                }
                // could tighten continuous as well
                if (nInt) {
                  if (ifCut) {
                    cc.setEffectiveness(100.0);
                  } else {
                    cc.setEffectiveness(1.0e-5);
                  }
#if CGL_DEBUG > 0
                  checkBounds(si,cc);
#endif
                  cs.insert(cc);
                }
              }
            }
            /* restore all */
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              if (oldU>oldL+1.0e-4)
                markC[icol]=0;
              else
                markC[icol]=3;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0) {
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    sum += rowElements[kk]*colsol[column[kk]];
                  }
                  if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_||info->strengthenRow) {
                    // can be a cut
                    // add gap to integer coefficient
                    // saveU and saveL spare
                    int * index = reinterpret_cast<int *>(saveL);
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("c point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("11Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
			rowCut.addCutIfNotDuplicate(rc,irow);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum-gap*colsol[j]<rowLower[irow]+primalTolerance_||info->strengthenRow) {
		      // can be a cut
		      // subtract gap from integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL);
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]-gap*(colUpper[j]-1));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("d point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("12Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
			rowCut.addCutIfNotDuplicate(rc,irow);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    }
	  }
        }
    }
  }
  delete [] stackC0;
  delete [] lo0;
  delete [] up0;
  delete [] tempL;
  delete [] tempU;
  delete [] markC;
  delete [] stackC;
  delete [] stackR;
  delete [] saveL;
  delete [] saveU;
  delete [] saveMin;
  delete [] saveMax;
  delete [] index;
  delete [] element;
  delete [] colsol;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info->strengthenRow,0);
  }
  delete [] array;
  abort();
  return (ninfeas);
}

/* Set mode

  Take care to change only the basic mode, without affecting higher bits.
*/
void CglProbing::setMode(int mode)
{
  if (mode >= 0 && mode < 3) {
    mode_ &= ~15;
    mode_ |= mode;
  }
}

// Return only the basic mode
int CglProbing::getMode() const
{
  return mode_&15;
}


// Set maximum number of passes per node
void CglProbing::setMaxPass(int value)
{
  if (value>0)
    maxPass_=value;
}


// Get maximum number of passes per node
int CglProbing::getMaxPass() const
{
  return maxPass_;
}


// Set log level
void CglProbing::setLogLevel(int value)
{
  if (value>=0)
    logLevel_=value;
}


// Get log level
int CglProbing::getLogLevel() const
{
  return logLevel_;
}


// Set maximum number of unsatisfied variables to look at
void CglProbing::setMaxProbe(int value)
{
  if (value>=0)
    maxProbe_=value;
}


// Get maximum number of unsatisfied variables to look at
int CglProbing::getMaxProbe() const
{
  return maxProbe_;
}


// Set maximum number of variables to look at in one probe
void CglProbing::setMaxLook(int value)
{
  if (value>=0)
    maxStack_=value;
}


// Get maximum number of variables to look at in one probe
int CglProbing::getMaxLook() const
{
  return maxStack_;
}


// Set maximum number of elements in row for scan
void CglProbing::setMaxElements(int value)
{
  if (value>0)
    maxElements_=value;
}


// Get maximum number of elements in row for scan
int CglProbing::getMaxElements() const
{
  return maxElements_;
}


// Set maximum number of passes per node (root node)
void CglProbing::setMaxPassRoot(int value)
{
  if (value>0)
    maxPassRoot_=value;
}


// Get maximum number of passes per node (root node)
int CglProbing::getMaxPassRoot() const
{
  return maxPassRoot_;
}


// Set maximum number of unsatisfied variables to look at (root node)
void CglProbing::setMaxProbeRoot(int value)
{
  if (value>0)
    maxProbeRoot_=value;
}


// Get maximum number of unsatisfied variables to look at (root node)
int CglProbing::getMaxProbeRoot() const
{
  return maxProbeRoot_;
}


// Set maximum number of variables to look at in one probe (root node)
void CglProbing::setMaxLookRoot(int value)
{
  if (value>0)
    maxStackRoot_=value;
}


// Get maximum number of variables to look at in one probe (root node)
int CglProbing::getMaxLookRoot() const
{
  return maxStackRoot_;
}


// Set maximum number of elements in row for scan (root node)
void CglProbing::setMaxElementsRoot(int value)
{
  if (value>0)
    maxElementsRoot_=value;
}


// Get maximum number of elements in row for scan (root node)
int CglProbing::getMaxElementsRoot() const
{
  return maxElementsRoot_;
}


// Set whether to use objective
void CglProbing::setUsingObjective(int yesNo)
{
  usingObjective_=yesNo;
}


// Get whether objective is being used
int CglProbing::getUsingObjective() const
{
  return usingObjective_;
}


// Decide whether to do row cuts
void CglProbing::setRowCuts(int type)
{
  if (type>-5&&type<5)
    rowCuts_=type;
}


// Returns row cuts generation type
int CglProbing::rowCuts() const
{
  return rowCuts_;
}


// Returns tight lower
const double * CglProbing::tightLower() const
{
  return colLower_;
}


// Returns tight upper
const double * CglProbing::tightUpper() const
{
  return colUpper_;
}


// Returns relaxed Row lower
const double * CglProbing::relaxedRowLower() const
{
  return rowLower_;
}


// Returns relaxed Row upper
const double * CglProbing::relaxedRowUpper() const
{
  return rowUpper_;
}


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing ()
:
CglCutGenerator(),
primalTolerance_(1.0e-07),
mode_(1),
rowCuts_(1),
maxPass_(3),
logLevel_(0),
maxProbe_(100),
maxStack_(50),
maxElements_(1000),
maxPassRoot_(3),
maxProbeRoot_(100),
maxStackRoot_(50),
maxElementsRoot_(10000),
usingObjective_(0)
{

  numberRows_=0;
  numberColumns_=0;
  rowCopy_=NULL;
  columnCopy_=NULL;
  rowLower_=NULL;
  rowUpper_=NULL;
  colLower_=NULL;
  colUpper_=NULL;
  numberIntegers_=0;
  number01Integers_=0;
  numberThisTime_=0;
  totalTimesCalled_=0;
  lookedAt_=NULL;
  cutVector_=NULL;
  numberCliques_=0;
  cliqueType_=NULL;
  cliqueStart_=NULL;
  cliqueEntry_=NULL;
  oneFixStart_=NULL;
  zeroFixStart_=NULL;
  endFixStart_=NULL;
  whichClique_=NULL;
  cliqueRow_=NULL;
  cliqueRowStart_=NULL;
  tightenBounds_=NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing (  const CglProbing & rhs)
                                                              :
  CglCutGenerator(rhs),
  primalTolerance_(rhs.primalTolerance_),
  mode_(rhs.mode_),
  rowCuts_(rhs.rowCuts_),
  maxPass_(rhs.maxPass_),
  logLevel_(rhs.logLevel_),
  maxProbe_(rhs.maxProbe_),
  maxStack_(rhs.maxStack_),
  maxElements_(rhs.maxElements_),
  maxPassRoot_(rhs.maxPassRoot_),
  maxProbeRoot_(rhs.maxProbeRoot_),
  maxStackRoot_(rhs.maxStackRoot_),
  maxElementsRoot_(rhs.maxElementsRoot_),
  usingObjective_(rhs.usingObjective_)
{  
  numberRows_=rhs.numberRows_;
  numberColumns_=rhs.numberColumns_;
  numberCliques_=rhs.numberCliques_;
  if (rhs.rowCopy_) {
    rowCopy_= new CoinPackedMatrix(*(rhs.rowCopy_));
    columnCopy_= new CoinPackedMatrix(*(rhs.columnCopy_));
    rowLower_=new double[numberRows_];
    CoinMemcpyN(rhs.rowLower_,numberRows_,rowLower_);
    rowUpper_=new double[numberRows_];
    CoinMemcpyN(rhs.rowUpper_,numberRows_,rowUpper_);
    colLower_=new double[numberColumns_];
    CoinMemcpyN(rhs.colLower_,numberColumns_,colLower_);
    colUpper_=new double[numberColumns_];
    CoinMemcpyN(rhs.colUpper_,numberColumns_,colUpper_);
    int i;
    numberIntegers_=rhs.numberIntegers_;
    number01Integers_=rhs.number01Integers_;
    cutVector_=new disaggregation [number01Integers_];
    CoinMemcpyN(rhs.cutVector_,number01Integers_,cutVector_);
    for (i=0;i<number01Integers_;i++) {
      if (cutVector_[i].index) {
	cutVector_[i].index = CoinCopyOfArray(rhs.cutVector_[i].index,cutVector_[i].length);
      }
    }
  } else {
    rowCopy_=NULL;
    columnCopy_=NULL;
    rowLower_=NULL;
    rowUpper_=NULL;
    colLower_=NULL;
    colUpper_=NULL;
    numberIntegers_=0;
    number01Integers_=0;
    cutVector_=NULL;
  }
  numberThisTime_=rhs.numberThisTime_;
  totalTimesCalled_=rhs.totalTimesCalled_;
  if (numberColumns_)
    lookedAt_=CoinCopyOfArray(rhs.lookedAt_,numberColumns_);
  else
    lookedAt_ = NULL;
  if (numberCliques_) {
    cliqueType_ = new cliqueType [numberCliques_];
    CoinMemcpyN(rhs.cliqueType_,numberCliques_,cliqueType_);
    cliqueStart_ = new int [numberCliques_+1];
    CoinMemcpyN(rhs.cliqueStart_,(numberCliques_+1),cliqueStart_);
    int n = cliqueStart_[numberCliques_];
    cliqueEntry_ = new cliqueEntry [n];
    CoinMemcpyN(rhs.cliqueEntry_,n,cliqueEntry_);
    oneFixStart_ = new int [numberColumns_];
    CoinMemcpyN(rhs.oneFixStart_,numberColumns_,oneFixStart_);
    zeroFixStart_ = new int [numberColumns_];
    CoinMemcpyN(rhs.zeroFixStart_,numberColumns_,zeroFixStart_);
    endFixStart_ = new int [numberColumns_];
    CoinMemcpyN(rhs.endFixStart_,numberColumns_,endFixStart_);
    int n2=-1;
    for (int i=numberColumns_-1;i>=0;i--) {
      if (oneFixStart_[i]>=0) {
	n2=endFixStart_[i];
	break;
      }
    }
    assert (n==n2);
    whichClique_ = new int [n];
    CoinMemcpyN(rhs.whichClique_,n,whichClique_);
    if (rhs.cliqueRowStart_) {
      cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_,numberRows_+1);
      n=cliqueRowStart_[numberRows_];
      cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_,n);
    } else {
      cliqueRow_=NULL;
      cliqueRowStart_=NULL;
    }
  } else {
    cliqueType_=NULL;
    cliqueStart_=NULL;
    cliqueEntry_=NULL;
    oneFixStart_=NULL;
    zeroFixStart_=NULL;
    endFixStart_=NULL;
    cliqueRow_=NULL;
    cliqueRowStart_=NULL;
    whichClique_=NULL;
  }
  if (rhs.tightenBounds_) {
    assert (numberColumns_);
    tightenBounds_=CoinCopyOfArray(rhs.tightenBounds_,numberColumns_);
  } else {
    tightenBounds_=NULL;
  }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglProbing::clone() const
{
  return new CglProbing(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglProbing::~CglProbing ()
{
  // free memory
  delete [] rowLower_;
  delete [] rowUpper_;
  delete [] colLower_;
  delete [] colUpper_;
  delete rowCopy_;
  delete columnCopy_;
  delete [] lookedAt_;
  delete [] cliqueType_;
  delete [] cliqueStart_;
  delete [] cliqueEntry_;
  delete [] oneFixStart_;
  delete [] zeroFixStart_;
  delete [] endFixStart_;
  delete [] whichClique_;
  delete [] cliqueRow_;
  delete [] cliqueRowStart_;
  if (cutVector_) {
    for (int i=0;i<number01Integers_;i++) {
      delete [] cutVector_[i].index;
    }
    delete [] cutVector_;
  }
  delete [] tightenBounds_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglProbing &
CglProbing::operator=(
                                         const CglProbing& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    primalTolerance_=rhs.primalTolerance_;
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    delete [] rowLower_;
    delete [] rowUpper_;
    delete [] colLower_;
    delete [] colUpper_;
    delete rowCopy_;
    delete columnCopy_;
    delete [] lookedAt_;
    delete [] cliqueType_;
    delete [] cliqueStart_;
    delete [] cliqueEntry_;
    delete [] oneFixStart_;
    delete [] zeroFixStart_;
    delete [] endFixStart_;
    delete [] whichClique_;
    delete [] cliqueRow_;
    delete [] cliqueRowStart_;
    delete [] tightenBounds_;
    mode_=rhs.mode_;
    rowCuts_=rhs.rowCuts_;
    maxPass_=rhs.maxPass_;
    logLevel_=rhs.logLevel_;
    maxProbe_=rhs.maxProbe_;
    maxStack_=rhs.maxStack_;
    maxElements_ = rhs.maxElements_;
    maxPassRoot_ = rhs.maxPassRoot_;
    maxProbeRoot_ = rhs.maxProbeRoot_;
    maxStackRoot_ = rhs.maxStackRoot_;
    maxElementsRoot_ = rhs.maxElementsRoot_;
    usingObjective_=rhs.usingObjective_;
    numberCliques_=rhs.numberCliques_;
    if (rhs.rowCopy_) {
      rowCopy_= new CoinPackedMatrix(*(rhs.rowCopy_));
      columnCopy_= new CoinPackedMatrix(*(rhs.columnCopy_));
      rowLower_=new double[numberRows_];
      CoinMemcpyN(rhs.rowLower_,numberRows_,rowLower_);
      rowUpper_=new double[numberRows_];
      CoinMemcpyN(rhs.rowUpper_,numberRows_,rowUpper_);
      colLower_=new double[numberColumns_];
      CoinMemcpyN(rhs.colLower_,numberColumns_,colLower_);
      colUpper_=new double[numberColumns_];
      CoinMemcpyN(rhs.colUpper_,numberColumns_,colUpper_);
      int i;
      numberIntegers_=rhs.numberIntegers_;
      number01Integers_=rhs.number01Integers_;
      for (i=0;i<number01Integers_;i++) {
        delete [] cutVector_[i].index;
      }
      delete [] cutVector_;
      cutVector_=new disaggregation [number01Integers_];
      CoinMemcpyN(rhs.cutVector_,number01Integers_,cutVector_);
      for (i=0;i<number01Integers_;i++) {
        if (cutVector_[i].index) {
          cutVector_[i].index = CoinCopyOfArray(rhs.cutVector_[i].index,cutVector_[i].length);
        }
      }
    } else {
      rowCopy_=NULL;
      columnCopy_=NULL;
      rowLower_=NULL;
      rowUpper_=NULL;
      colLower_=NULL;
      colUpper_=NULL;
      numberIntegers_=0;
      number01Integers_=0;
      cutVector_=NULL;
    }
    numberThisTime_=rhs.numberThisTime_;
    totalTimesCalled_=rhs.totalTimesCalled_;
    if (numberColumns_)
      lookedAt_=CoinCopyOfArray(rhs.lookedAt_,numberColumns_);
    else
      lookedAt_ = NULL;
    if (numberCliques_) {
      cliqueType_ = new cliqueType [numberCliques_];
      CoinMemcpyN(rhs.cliqueType_,numberCliques_,cliqueType_);
      cliqueStart_ = new int [numberCliques_+1];
      CoinMemcpyN(rhs.cliqueStart_,(numberCliques_+1),cliqueStart_);
      int n = cliqueStart_[numberCliques_];
      cliqueEntry_ = new cliqueEntry [n];
      CoinMemcpyN(rhs.cliqueEntry_,n,cliqueEntry_);
      oneFixStart_ = new int [numberColumns_];
      CoinMemcpyN(rhs.oneFixStart_,numberColumns_,oneFixStart_);
      zeroFixStart_ = new int [numberColumns_];
      CoinMemcpyN(rhs.zeroFixStart_,numberColumns_,zeroFixStart_);
      endFixStart_ = new int [numberColumns_];
      CoinMemcpyN(rhs.endFixStart_,numberColumns_,endFixStart_);
      int n2=-1;
      for (int i=numberColumns_-1;i>=0;i--) {
	if (oneFixStart_[i]>=0) {
	  n2=endFixStart_[i];
	  break;
	}
      }
      assert (n==n2);
      whichClique_ = new int [n];
      CoinMemcpyN(rhs.whichClique_,n,whichClique_);
      if (rhs.cliqueRowStart_) {
        cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_,numberRows_+1);
        n=cliqueRowStart_[numberRows_];
        cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_,n);
      } else {
        cliqueRow_=NULL;
        cliqueRowStart_=NULL;
      }
    } else {
      cliqueType_=NULL;
      cliqueStart_=NULL;
      cliqueEntry_=NULL;
      oneFixStart_=NULL;
      zeroFixStart_=NULL;
      endFixStart_=NULL;
      whichClique_=NULL;
      cliqueRow_=NULL;
      cliqueRowStart_=NULL;
    }
    if (rhs.tightenBounds_) {
      assert (numberColumns_);
      tightenBounds_=CoinCopyOfArray(rhs.tightenBounds_,numberColumns_);
    } else {
      tightenBounds_=NULL;
    }
  }
  return *this;
}

/// This can be used to refresh any inforamtion
void 
CglProbing::refreshSolver(OsiSolverInterface * solver)
{
  if (rowCopy_) {
    // snapshot existed - redo
    snapshot(*solver,NULL);
  }
}


/* Creates cliques for use by probing.
   Can also try and extend cliques as a result of probing (root node).
   Returns number of cliques found.
*/
int 
CglProbing::createCliques( OsiSolverInterface & si, 
			  int minimumSize, int maximumSize)
{
  // get rid of what is there
  deleteCliques();
  CoinPackedMatrix matrixByRow(*si.getMatrixByRow());
  int numberRows = si.getNumRows();
  if (!rowCopy_)
    numberRows_=numberRows;
  numberColumns_ = si.getNumCols();

  numberCliques_=0;
  int numberEntries=0;
  int numberIntegers=0;
  int * lookup = new int[numberColumns_];
  int i;
  for (i=0;i<numberColumns_;i++) {
    if (si.isBinary(i))
      lookup[i]=numberIntegers++;
    else
      lookup[i]=-1;
  }

  int * which = new int[numberColumns_];
  int * whichRow = new int[numberRows];
  // Statistics
  int totalP1=0,totalM1=0;
  int numberBig=0,totalBig=0;
  int numberFixed=0;

  // Row copy
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column lengths for slacks
  const int * columnLength = si.getMatrixByCol()->getVectorLengths();

  const double * lower = si.getColLower();
  const double * upper = si.getColUpper();
  const double * rowLower = si.getRowLower();
  const double * rowUpper = si.getRowUpper();
  int iRow;
  for (iRow=0;iRow<numberRows;iRow++) {
    int numberP1=0, numberM1=0;
    int j;
    double upperValue=rowUpper[iRow];
    double lowerValue=rowLower[iRow];
    bool good=true;
    int slack = -1;
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn = column[j];
      int iInteger=lookup[iColumn];
      if (upper[iColumn]-lower[iColumn]<1.0e-8) {
	// fixed
	upperValue -= lower[iColumn]*elementByRow[j];
	lowerValue -= lower[iColumn]*elementByRow[j];
	continue;
      } else if (upper[iColumn]!=1.0||lower[iColumn]!=0.0) {
	good = false;
	break;
      } else if (iInteger<0) {
	good = false;
	break;
      } else {
	if (columnLength[iColumn]==1)
	  slack = iInteger;
      }
      if (fabs(elementByRow[j])!=1.0) {
	good=false;
	break;
      } else if (elementByRow[j]>0.0) {
	which[numberP1++]=iColumn;
      } else {
	numberM1++;
	which[numberIntegers-numberM1]=iColumn;
      }
    }
    int iUpper = static_cast<int> (floor(upperValue+1.0e-5));
    int iLower = static_cast<int> (ceil(lowerValue-1.0e-5));
    int state=0;
    if (upperValue<1.0e6) {
      if (iUpper==1-numberM1)
	state=1;
      else if (iUpper==-numberM1)
	state=2;
      else if (iUpper<-numberM1)
	state=3;
    }
    if (!state&&lowerValue>-1.0e6) {
      if (-iLower==1-numberP1)
	state=-1;
      else if (-iLower==-numberP1)
	state=-2;
      else if (-iLower<-numberP1)
	state=-3;
    }
    if (good&&state) {
      if (abs(state)==3) {
	// infeasible
	numberCliques_ = -99999;
	break;
      } else if (abs(state)==2) {
	// we can fix all
	numberFixed += numberP1+numberM1;
	if (state>0) {
	  // fix all +1 at 0, -1 at 1
	  for (i=0;i<numberP1;i++)
	    si.setColUpper(which[i],0.0);
	  for (i=0;i<numberM1;i++)
	    si.setColLower(which[numberIntegers-i-1],
				 1.0);
	} else {
	  // fix all +1 at 1, -1 at 0
	  for (i=0;i<numberP1;i++)
	    si.setColLower(which[i],1.0);
	  for (i=0;i<numberM1;i++)
	    si.setColUpper(which[numberIntegers-i-1],
				 0.0);
	}
      } else {
	int length = numberP1+numberM1;
        totalP1 += numberP1;
        totalM1 += numberM1;
	if (length >= minimumSize&&length<maximumSize) {
	  whichRow[numberCliques_++]=iRow;
	  numberEntries += length;
	} else if (numberP1+numberM1 >= maximumSize) {
	  // too big
	  numberBig++;
	  totalBig += numberP1+numberM1;
	}
      }
    }
  }
  if (numberCliques_<0) {
    if (logLevel_)
      printf("*** Problem infeasible\n");
  } else {
    if (numberCliques_) {
      if (logLevel_)
        printf("%d cliques of average size %g found, %d P1, %d M1\n",
               numberCliques_,
               (static_cast<double>(totalP1+totalM1))/
	       (static_cast<double> (numberCliques_)),
               totalP1,totalM1);
    } else {
      if (logLevel_>1)
        printf("No cliques found\n");
    }
    if (numberBig) {
      if (logLevel_)
        printf("%d large cliques ( >= %d) found, total %d\n",
	     numberBig,maximumSize,totalBig);
    }
    if (numberFixed) {
      if (logLevel_)
        printf("%d variables fixed\n",numberFixed);
    }
  }
  if (numberCliques_>0) {
    cliqueType_ = new cliqueType [numberCliques_];
    cliqueStart_ = new int [numberCliques_+1];
    cliqueEntry_ = new cliqueEntry [numberEntries];
    oneFixStart_ = new int [numberColumns_];
    zeroFixStart_ = new int [numberColumns_];
    endFixStart_ = new int [numberColumns_];
    whichClique_ = new int [numberEntries];
    numberEntries=0;
    cliqueStart_[0]=0;
    for (i=0;i<numberColumns_;i++) {
      oneFixStart_[i]=-1;
      zeroFixStart_[i]=-1;
      endFixStart_[i]=-1;
    }
    int iClique;
    // Possible some have been fixed
    int numberCliques=numberCliques_;
    numberCliques_=0;
    for (iClique=0;iClique<numberCliques;iClique++) {
      int iRow=whichRow[iClique];
      whichRow[numberCliques_]=iRow;
      int numberP1=0, numberM1=0;
      int j;
      double upperValue=rowUpper[iRow];
      double lowerValue=rowLower[iRow];
      for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	int iColumn = column[j];
	if (upper[iColumn]-lower[iColumn]<1.0e-8) {
	  // fixed
	  upperValue -= lower[iColumn]*elementByRow[j];
	  lowerValue -= lower[iColumn]*elementByRow[j];
	  continue;
	}
	if (elementByRow[j]>0.0) {
	  which[numberP1++]=iColumn;
	} else {
	  numberM1++;
	  which[numberIntegers-numberM1]=iColumn;
	}
      }
      int iUpper = static_cast<int> (floor(upperValue+1.0e-5));
      int iLower = static_cast<int> (ceil(lowerValue-1.0e-5));
      int state=0;
      if (upperValue<1.0e6) {
	if (iUpper==1-numberM1)
	  state=1;
      }
      if (!state&&lowerValue>-1.0e6) {
	state=-1;
      }
      if (abs(state)!=1)
	continue; // must have been fixed
      if (iLower==iUpper) {
	cliqueType_[numberCliques_].equality=1;
      } else {
	cliqueType_[numberCliques_].equality=0;
      }
      if (state>0) {
	for (i=0;i<numberP1;i++) {
	  // 1 is strong branch
	  int iColumn = which[i];
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn);
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],true);
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
	for (i=0;i<numberM1;i++) {
	  // 0 is strong branch
	  int iColumn = which[numberIntegers-i-1];
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn);
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],false);
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
      } else {
	for (i=0;i<numberP1;i++) {
	  // 0 is strong branch
	  int iColumn = which[i];
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn);
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],false);
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
	for (i=0;i<numberM1;i++) {
	  // 1 is strong branch
	  int iColumn = which[numberIntegers-i-1];
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn);
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],true);
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
      }
      numberCliques_++;
      cliqueStart_[numberCliques_]=numberEntries;
    }
    // Now do column lists
    // First do counts
    for (iClique=0;iClique<numberCliques_;iClique++) {
      for (int j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
	int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
	if (oneFixesInCliqueEntry(cliqueEntry_[j]))
	  oneFixStart_[iColumn]++;
	else
	  zeroFixStart_[iColumn]++;
      }
    }
    // now get starts and use which and end as counters
    numberEntries=0;
    for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (oneFixStart_[iColumn]>=0) {
	int n1=oneFixStart_[iColumn];
	int n2=zeroFixStart_[iColumn];
	oneFixStart_[iColumn]=numberEntries;
	which[iColumn]=numberEntries;
	numberEntries += n1;
	zeroFixStart_[iColumn]=numberEntries;
	endFixStart_[iColumn]=numberEntries;
	numberEntries += n2;
      }
    }
    // now put in
    for (iClique=0;iClique<numberCliques_;iClique++) {
      for (int j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
	int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
	if (oneFixesInCliqueEntry(cliqueEntry_[j])) {
	  int put = which[iColumn];
	  which[iColumn]++;
	  whichClique_[put]=iClique;
	} else {
	  int put = endFixStart_[iColumn];
	  endFixStart_[iColumn]++;
	  whichClique_[put]=iClique;
	}
      }
    }
  }
  delete [] which;
  delete [] whichRow;
  delete [] lookup;
  return numberCliques_;
}


// Delete all clique information
void 
CglProbing::deleteCliques()
{
  delete [] cliqueType_;
  delete [] cliqueStart_;
  delete [] cliqueEntry_;
  delete [] oneFixStart_;
  delete [] zeroFixStart_;
  delete [] endFixStart_;
  delete [] whichClique_;
  delete [] cliqueRow_;
  delete [] cliqueRowStart_;
  cliqueType_=NULL;
  cliqueStart_=NULL;
  cliqueEntry_=NULL;
  oneFixStart_=NULL;
  zeroFixStart_=NULL;
  endFixStart_=NULL;
  whichClique_=NULL;
  cliqueRow_=NULL;
  cliqueRowStart_=NULL;
  numberCliques_=0;
}


/*
  Returns true if CglProbing may generate row cuts in the search tree (rather
  than just at the root node).  Provided so the client can know if the matrix
  will change in the tree.  Really meant so column cut generators can still
  be active without worrying code.  Default is true

  The above is enforced in generateCuts and generateCutsAndModify. Note that
  this method may well return an incorrect answer if called during cut
  generation.
*/
bool 
CglProbing::mayGenerateRowCutsInTree() const
{
  return rowCuts_ > 0 ;
}


// Sets up clique information for each row
void 
CglProbing::setupRowCliqueInformation(const OsiSolverInterface & si) 
{
  if (!numberCliques_)
    return;
  CoinPackedMatrix * rowCopy;
  if (!rowCopy_) {
    // create from current
    numberRows_=si.getNumRows(); 
    numberColumns_=si.getNumCols(); 
    rowCopy = new CoinPackedMatrix(*si.getMatrixByRow());
  } else {
    rowCopy = rowCopy_;
    assert(numberRows_<=si.getNumRows()); 
    assert(numberColumns_==si.getNumCols()); 
  }
  assert(numberRows_&&numberColumns_);
  cliqueRowStart_ = new int [numberRows_+1];
  cliqueRowStart_[0]=0;
  // Temporary array while building list
  cliqueEntry ** array = new cliqueEntry * [numberRows_];
  // Which cliques in use
  int * which = new int[numberCliques_];
  int * count = new int[numberCliques_];
  int * back =new int[numberColumns_];
  CoinZeroN(count,numberCliques_);
  CoinFillN(back,numberColumns_,-1);
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * lower = si.getColLower();
  const double * upper = si.getColUpper();
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int j;
    int numberFree=0;
    int numberUsed=0;
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn=column[j];
      if (upper[iColumn]>lower[iColumn]) {
        back[iColumn]=j-rowStart[iRow];
        numberFree++;
        for (int k=oneFixStart_[iColumn];k<endFixStart_[iColumn];k++) {
          int iClique = whichClique_[k];
          if (!count[iClique]) {
            which[numberUsed++]=iClique;
          }
          count[iClique]++;
        }
      }
    }
    // find largest cliques
    bool finished=false;
    int numberInThis=0;
    cliqueEntry * entries = NULL;
    array[iRow]=entries;
    while (!finished) {
      int largest=1;
      int whichClique=-1;
      for (int i=0;i<numberUsed;i++) {
        int iClique = which[i];
        if (count[iClique]>largest) {
          largest=count[iClique];
          whichClique=iClique;
        }
      }
      // Add in if >1 (but not if all as that means clique==row)
      if (whichClique>=0&&largest<numberFree) {
        if (!numberInThis) {
          int length=rowLength[iRow];
          entries = new cliqueEntry [length];
          array[iRow]=entries;
          for (int i=0;i<length;i++) {
            setOneFixesInCliqueEntry(entries[i],false);
            setSequenceInCliqueEntry(entries[i],numberColumns_+1);
          }
        }
        // put in (and take out all counts)
        for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
          int iColumn=column[j];
          if (upper[iColumn]>lower[iColumn]) {
            bool found=false;
            int k;
            for ( k=oneFixStart_[iColumn];k<endFixStart_[iColumn];k++) {
              int iClique = whichClique_[k];
              if (iClique==whichClique) {
                found=true;
                break;
              }
            }
            if (found) {
              for ( k=oneFixStart_[iColumn];k<endFixStart_[iColumn];k++) {
                int iClique = whichClique_[k];
                count[iClique]--;
              }
              for (k=cliqueStart_[whichClique];k<cliqueStart_[whichClique+1];k++) {
                if (sequenceInCliqueEntry(cliqueEntry_[k])==iColumn) {
                  int iback=back[iColumn];
                  setSequenceInCliqueEntry(entries[iback],numberInThis);
                  setOneFixesInCliqueEntry(entries[iback],
					   oneFixesInCliqueEntry(cliqueEntry_[k]));
                  break;
                }
              }
            }
          }
        }
        numberInThis++;
      } else {
        finished=true;
      }
    }
    if (numberInThis) 
      cliqueRowStart_[iRow+1]=cliqueRowStart_[iRow]+rowLength[iRow];
    else
      cliqueRowStart_[iRow+1]=cliqueRowStart_[iRow];
    for (int i=0;i<numberUsed;i++) {
      int iClique = which[i];
      count[iClique]=0;
    }
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn=column[j];
      back[iColumn]=-1;
    }
  }
  delete [] which;
  delete [] count;
  delete [] back;
  // Now put info in one array
  cliqueRow_ = new cliqueEntry [cliqueRowStart_[numberRows_]];
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (array[iRow]) {
      int start = cliqueRowStart_[iRow];
      CoinMemcpyN(array[iRow],rowLength[iRow],cliqueRow_+start);
      delete [] array[iRow];
    }
  }
  delete [] array;
  if (rowCopy!=rowCopy_)
    delete rowCopy;
}


// Mark variables to be tightened
void 
CglProbing::tightenThese(const OsiSolverInterface & solver,int number, const int * which)
{
  delete [] tightenBounds_;
  int numberColumns = solver.getNumCols();
  if (numberColumns_)
    assert (numberColumns_==numberColumns);
  tightenBounds_ = new char [numberColumns];
  memset(tightenBounds_,0,numberColumns);
  for (int i=0;i<number;i++) {
    int k=which[i];
    if (k>=0&&k<numberColumns)
      tightenBounds_[k]=1;
  }
}


// Create C++ lines to get to current state
std::string
CglProbing::generateCpp( FILE * fp) 
{
  CglProbing other;
  fprintf(fp,"0#include \"CglProbing.hpp\"\n");
  fprintf(fp,"3  CglProbing probing;\n");
  if (getMode()!=other.getMode())
    fprintf(fp,"3  probing.setMode(%d);\n",getMode());
  else
    fprintf(fp,"4  probing.setMode(%d);\n",getMode());
  if (getMaxPass()!=other.getMaxPass())
    fprintf(fp,"3  probing.setMaxPass(%d);\n",getMaxPass());
  else
    fprintf(fp,"4  probing.setMaxPass(%d);\n",getMaxPass());
  if (getLogLevel()!=other.getLogLevel())
    fprintf(fp,"3  probing.setLogLevel(%d);\n",getLogLevel());
  else
    fprintf(fp,"4  probing.setLogLevel(%d);\n",getLogLevel());
  if (getMaxProbe()!=other.getMaxProbe())
    fprintf(fp,"3  probing.setMaxProbe(%d);\n",getMaxProbe());
  else
    fprintf(fp,"4  probing.setMaxProbe(%d);\n",getMaxProbe());
  if (getMaxLook()!=other.getMaxLook())
    fprintf(fp,"3  probing.setMaxLook(%d);\n",getMaxLook());
  else
    fprintf(fp,"4  probing.setMaxLook(%d);\n",getMaxLook());
  if (getMaxElements()!=other.getMaxElements())
    fprintf(fp,"3  probing.setMaxElements(%d);\n",getMaxElements());
  else
    fprintf(fp,"4  probing.setMaxElements(%d);\n",getMaxElements());
  if (getMaxPassRoot()!=other.getMaxPassRoot())
    fprintf(fp,"3  probing.setMaxPassRoot(%d);\n",getMaxPassRoot());
  else
    fprintf(fp,"4  probing.setMaxPassRoot(%d);\n",getMaxPassRoot());
  if (getMaxProbeRoot()!=other.getMaxProbeRoot())
    fprintf(fp,"3  probing.setMaxProbeRoot(%d);\n",getMaxProbeRoot());
  else
    fprintf(fp,"4  probing.setMaxProbeRoot(%d);\n",getMaxProbeRoot());
  if (getMaxLookRoot()!=other.getMaxLookRoot())
    fprintf(fp,"3  probing.setMaxLookRoot(%d);\n",getMaxLookRoot());
  else
    fprintf(fp,"4  probing.setMaxLookRoot(%d);\n",getMaxLookRoot());
  if (getMaxElementsRoot()!=other.getMaxElementsRoot())
    fprintf(fp,"3  probing.setMaxElementsRoot(%d);\n",getMaxElementsRoot());
  else
    fprintf(fp,"4  probing.setMaxElementsRoot(%d);\n",getMaxElementsRoot());
  if (rowCuts()!=other.rowCuts())
    fprintf(fp,"3  probing.setRowCuts(%d);\n",rowCuts());
  else
    fprintf(fp,"4  probing.setRowCuts(%d);\n",rowCuts());
  if (getUsingObjective()!=other.getUsingObjective())
    fprintf(fp,"3  probing.setUsingObjective(%d);\n",getUsingObjective());
  else
    fprintf(fp,"4  probing.setUsingObjective(%d);\n",getUsingObjective());
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  probing.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  probing.setAggressiveness(%d);\n",getAggressiveness());
  return "probing";
}

