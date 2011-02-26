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
// #undef NDEBUG

#include "CoinPackedMatrix.hpp"
#include "CglProbing.hpp"
#include "CglProbingRowCut.hpp"
#include "CglProbingDebug.hpp"

/*
  File local namespace for utilities
*/
namespace {

/*
  Helper class for sorting
*/

typedef struct {
  double infeasibility ;
  int sequence ;
} double_int_pair ;

class double_int_pair_compare {
public:
  inline bool operator() (double_int_pair x, double_int_pair y) const
  {
    return ( x.infeasibility < y.infeasibility) ;
  }
} ;

#if CGL_DEBUG > 0
int nPath = 0 ;
#endif

/*
  Compare two sets of upper and lower bound arrays and generate column cuts.
*/

void makeColCuts (int nCols, OsiCuts &cs,
		  const char *const intVar,
		  const double *const colsol,
		  const double *const origlbs,
		  const double *const origubs,
		  const double *const soln,
		  double *const newlbs,
		  double *const newubs
		 )
{
  int numberChanged = 0 ;
  int ifCut = 0 ;
  CoinPackedVector lbs ;
  CoinPackedVector ubs ;
/*
  Scan the bound arrays, recording changed bounds and counting changes and
  actual cuts.
*/
  for (int i = 0 ; i < nCols ; ++i) {
    if (intVar[i]) {
      newubs[i] = CoinMin(origubs[i],floor(newubs[i]+1.0e-4)) ;
      if (newubs[i] < origubs[i]-1.0e-8) {
	if (newubs[i] < colsol[i]-1.0e-8) ifCut++ ;
	ubs.insert(i,newubs[i]) ;
	numberChanged++ ;
      }
      newlbs[i] = CoinMax(origlbs[i],ceil(newlbs[i]-1.0e-4)) ;
      if (newlbs[i] > origlbs[i]+1.0e-8) {
	if (newlbs[i] > colsol[i]+1.0e-8) ifCut++ ;
	lbs.insert(i,newlbs[i]) ;
	numberChanged++ ;
      }
    }
  }
/*
  Stash the changes in a column cut. If the changes actually cut off some part
  of the current solution, boost the effectiveness.
*/
  if (numberChanged) {
    OsiColCut cc ;
    cc.setUbs(ubs) ;
    cc.setLbs(lbs) ;
    if (ifCut) {
      cc.setEffectiveness(100.0) ;
    } else {
      cc.setEffectiveness(1.0e-5) ;
    }
    cs.insert(cc) ;
  }
}



/*
  This method will attempt to identify variables that will naturally take on
  integer values but were not originally identified as such, and mark them as
  integer. Upper and lower bounds are forced to integrality.
  
  Executed once, at first cut generation pass at the root.

  Which begs various questions:
    * Why not execute this independently, before any cut generation occurs?
    * Why not execute repeatedly, as the problem is simplified due to fixed
      variables?
  This would involve calling analyze as part of integer preprocessing.

  The expression abs(x - floor(x+.5)) evaluates to zero when x is integer and
  is nonzero otherwise. Farther down in the code, there's an odd-looking test
  that checks abs(1/aij - floor(1/aij + .5)). What we're looking for is
  coefficients aij = 1/k, k integer. This guarantees that b/aij is integer
  as long as b is integer.

  TODO: Should this routine return a count, instead of a boolean? No
	particular reason other than more information.  -- lh, 110211 --
*/

bool analyze (const OsiSolverInterface *solverX, char *intVar,
	      double *lower, double *upper)
{
# if CGLPROBING_DEBUG > 0
  std::cout << "Entering CglProbing::analyze." << std::endl ;
# endif
  const double e20Inf = 1.0e20 ;
  const int changedToInt = 77 ;

/*
  I don't see why we clone the solver here. There are no modifications to
  anything except the parameters intVar, lower, and upper.
*/
  OsiSolverInterface *solver = solverX->clone() ;
  const double *objective = solver->getObjCoefficients() ;
  int numberColumns = solver->getNumCols() ;
  int numberRows = solver->getNumRows() ;
  double direction = solver->getObjSense() ;
  int iRow,iColumn ;

/*
  Nor do I understand why we make copies of the matrix. Two, in fact! Surely
  we could ask the OSI for these. Aaaaah, but then we might be triggering
  changes inside the OSI? Shouldn't matter, getMatrixBy* are const.
*/
  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow()) ;
  const double *elementByRow = matrixByRow.getElements() ;
  const int *column = matrixByRow.getIndices() ;
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts() ;
  const int *rowLength = matrixByRow.getVectorLengths() ;

  // Column copy
  CoinPackedMatrix matrixByCol(*solver->getMatrixByCol()) ;
  const double *element = matrixByCol.getElements() ;
  const int *row = matrixByCol.getIndices() ;
  const CoinBigIndex *columnStart = matrixByCol.getVectorStarts() ;
  const int *columnLength = matrixByCol.getVectorLengths() ;

  const double * rowLower = solver->getRowLower() ;
  const double * rowUpper = solver->getRowUpper() ;

  char *ignore = new char [numberRows] ;
  int *which = new int[numberRows] ;
  double *changeRhs = new double[numberRows] ;
  memset(changeRhs,0,numberRows*sizeof(double)) ;
  memset(ignore,0,numberRows) ;

  int numberChanged = 0 ;
  bool finished = false ;
/*
  Open main analysis loop. Keep going until nothing changes.
*/
  while (!finished) {
    int saveNumberChanged = numberChanged ;
/*
  Open loop to scan each constraint. With luck, we can conclude that some
  variable must be integer. With less luck, the constraint will not prevent the
  variable from being integer.
*/
    for (iRow = 0 ; iRow < numberRows ; iRow++) {
      int numberContinuous = 0 ;
      double value1 = 0.0,value2 = 0.0 ;
      bool allIntegerCoeff = true ;
      double sumFixed = 0.0 ;
      int jColumn1 = -1,jColumn2 = -1 ;
/*
  Scan the coefficients of the constraint and accumulate some information:
    * Contribution of fixed variables.
    * Count of continuous variables, and column index and coefficients for
      the first and last continuous variable.
    * A boolean, allIntegerCoeff, true if all coefficients of unfixed
      integer variables are integer.
*/
      for (CoinBigIndex j = rowStart[iRow] ;
	   j < rowStart[iRow]+rowLength[iRow] ; j++) {
        int jColumn = column[j] ;
        double value = elementByRow[j] ;
        if (upper[jColumn] > lower[jColumn]+1.0e-8) {
          if (!intVar[jColumn]) {
            if (numberContinuous == 0) {
              jColumn1 = jColumn ;
              value1 = value ;
            } else {
              jColumn2 = jColumn ;
              value2 = value ;
            }
            numberContinuous++ ;
          } else {
            if (fabs(value-floor(value+0.5)) > 1.0e-12)
              allIntegerCoeff = false ;
          }
        } else {
          sumFixed += lower[jColumn]*value ;
        }
      }
/*
  See if the row bounds are integer after adjusting for the contribution of
  fixed variables. Serendipitous cancellation is possible here (if unlikely).
  The fixed contribution could change a non-integer bound to an integer.
*/
      double low = rowLower[iRow] ;
      if (low > -e20Inf) {
        low -= sumFixed ;
        if (fabs(low-floor(low+0.5)) > 1.0e-12)
          allIntegerCoeff = false ;
      }
      double up = rowUpper[iRow] ;
      if (up < e20Inf) {
        up -= sumFixed ;
        if (fabs(up-floor(up+0.5)) > 1.0e-12)
          allIntegerCoeff = false ;
      }
/*
  To make progress, it must be true that the coefficients of all unfixed
  integer variables are integer and the row bounds are integer.
*/
      if (!allIntegerCoeff) continue ;
/*
  Case 1: We have a single continuous variable to consider.

  If a<ij> is of the form 1/k, k integer, then x<j> can be integer. Put another
  way, this constraint will not prevent integrality and we can ignore it in
  further processing.  An equality forces integrality.

  So, for an equality with one continuous variable, shouldn't we mark the
  variable as `cannot be declared integer' so we don't waste time in subsequent
  passes?
*/
      if (numberContinuous == 1) {
        if (low == up) {
          if (fabs(value1) > 1.0e-3) {
            value1 = 1.0/value1 ;
            if (fabs(value1-floor(value1+0.5)) < 1.0e-12) {
              numberChanged++ ;
              intVar[jColumn1] = changedToInt ;
            }
          }
        } else {
          if (fabs(value1) > 1.0e-3) {
            value1 = 1.0/value1 ;
            if (fabs(value1-floor(value1+0.5)) < 1.0e-12) {
              ignore[iRow] = 1 ;
            }
          }
        }
      } else
/*
  Case 2: Two continuous variables. 

  John's comment is `need general theory - for now just look at 2 cases'. I've
  worked out the general theory. See the paper documentation. The code needs
  to be rewritten. What's here is overly restrictive for case 2.1 and wrong for
  case 2.2.

  2.1) Constraint is (integer lhs) + x<1> - x<2> = (integer rhs), x<1> and x<2>
     have lower bounds of zero and do not appear in any other constraint.
     The overall contribution to the objective, (c<1>+c<2>), is unfavourable.
     Consider any given value of (rhs-lhs) = (integer). Some x<i> will head
     towards infinity, the other towards zero, modulo finite bounds and
     feasibility. See the paper documentation.

  2.2) Constraint as above but two coefficients in each column. The second
     coefficients feed into G/L row(s) which will try and minimize. Also
     consider that the second coefficients may be in the same row. Code as
     written is incorrect.

  Case 2.2) has been disabled from the start (with a note that says `take out
  until fixed').
*/
      if (numberContinuous == 2) {
        if (low == up) {
          if (fabs(value1) == 1.0 && value1*value2 == -1.0 &&
	      !lower[jColumn1] && !lower[jColumn2] &&
	      columnLength[jColumn1] == 1 && columnLength[jColumn2] == 1) {
            int n = 0 ;
            int i ;
            double objChange =
		direction*(objective[jColumn1]+objective[jColumn2]) ;
            double bound = CoinMin(upper[jColumn1],upper[jColumn2]) ;
            bound = CoinMin(bound,e20Inf) ;
	    // Since column lengths are 1, there is no second coeff.
            for (i = columnStart[jColumn1] ;
		 i < columnStart[jColumn1]+columnLength[jColumn1] ; i++) {
              int jRow = row[i] ;
              double value = element[i] ;
              if (jRow != iRow) {
                which[n++] = jRow ;
                changeRhs[jRow] = value ;
              }
            }
            for (i = columnStart[jColumn2] ;
		 i < columnStart[jColumn2]+columnLength[jColumn2] ; i++) {
              int jRow = row[i] ;
              double value = element[i] ;
              if (jRow != iRow) {
                if (!changeRhs[jRow]) {
                  which[n++] = jRow ;
                  changeRhs[jRow] = value ;
                } else {
                  changeRhs[jRow] += value ;
                }
              }
            }
            if (objChange >= 0.0) {
              bool good = true ;
	      // Since n = 0, this loop never executes
              for (i = 0 ; i < n ; i++) {
                int jRow = which[i] ;
                double value = changeRhs[jRow] ;
                if (value) {
                  value *= bound ;
                  if (rowLength[jRow] == 1) {
                    if (value>0.0) {
                      double rhs = rowLower[jRow] ;
                      if (rhs>0.0) {
                        double ratio =rhs/value ;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good = false ;
                      }
                    } else {
                      double rhs = rowUpper[jRow] ;
                      if (rhs<0.0) {
                        double ratio =rhs/value ;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good = false ;
                      }
                    }
                  } else if (rowLength[jRow] == 2) {
                    if (value>0.0) {
                      if (rowLower[jRow]>-e20Inf)
                        good = false ;
                    } else {
                      if (rowUpper[jRow]<e20Inf)
                        good = false ;
                    }
                  } else {
                    good = false ;
                  }
                }
              }

              if (good) {
                // both can be integer
                numberChanged++ ;
                intVar[jColumn1] = changedToInt ;
                numberChanged++ ;
                intVar[jColumn2] = changedToInt ;
              }
            }
            // clear (another noop since n is zero)
            for (i = 0 ; i < n ; i++) {
              changeRhs[which[i]] = 0.0 ;
            }
          }
        }
      }
    } // end loop to scan each constraint
/*
  Check for variables x<j> with integer bounds l<j> and u<j> and no
  constraint that prevents x<j> from being integer.
*/
    for (iColumn = 0 ; iColumn < numberColumns ; iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8 && !intVar[iColumn]) {
        double value ;
        value = upper[iColumn] ;
        if (value < e20Inf && fabs(value-floor(value+0.5)) > 1.0e-12) 
          continue ;
        value = lower[iColumn] ;
        if (value > -e20Inf && fabs(value-floor(value+0.5)) > 1.0e-12) 
          continue ;
        bool integer = true ;
        for (CoinBigIndex j = columnStart[iColumn] ;
	     j < columnStart[iColumn]+columnLength[iColumn] ; j++) {
          int iRow = row[j] ;
          if (!ignore[iRow]) {
            integer = false ;
            break ;
          }
        }
        if (integer) {
          numberChanged++ ;
          intVar[iColumn] = changedToInt ;
        }
      }
    }
    finished = (numberChanged == saveNumberChanged) ;
  } // end main loop: while !finished
/*
  Nothing changed on the last pass. Time to clean up the new integer
  variables.  Force bounds to integer values. Feasibility can be lost here.
  The codes assigned here are curious choices. If the variable is fixed at
  zero, declare it binary (1), but if fixed at some other value, declare it
  continuous (0). Otherwise, declare it general integer (2).
*/
  bool feasible = true ;
  for (iColumn = 0 ; iColumn < numberColumns ; iColumn++) {
    if (intVar[iColumn] == changedToInt) {
      if (upper[iColumn] > e20Inf) {
        upper[iColumn] = e20Inf ;
      } else {
        upper[iColumn] = floor(upper[iColumn]+1.0e-5) ;
      }
      if (lower[iColumn] < -e20Inf) {
        lower[iColumn] = -e20Inf ;
      } else {
        lower[iColumn] = ceil(lower[iColumn]-1.0e-5) ;
        if (lower[iColumn] > upper[iColumn])
          feasible = false ;
      }
      if (lower[iColumn] == 0.0 && upper[iColumn] == 1.0)
        intVar[iColumn] = 1 ;
      else if (lower[iColumn] == upper[iColumn])
        intVar[iColumn] = 0 ;
      else
        intVar[iColumn] = 2 ;
    }
  }
/*
  Clean up and return.
*/
  delete [] which ;
  delete [] changeRhs ;
  delete [] ignore ;
  delete solver ;

# if CGLPROBING_DEBUG > 0
  std::cout
    << "CglProbing::analyze: " << numberChanged
    << " variables could be made integer." << std::endl ;
# endif

  return feasible ;
}


/*
  Set up the working row copy.
  
  If we have a snapshot handy, most of the setup is done -- we just make
  a copy of the snapshot.  Otherwise, we need to create a row-major copy
  of the constraint matrix and the rhs and rhslow vectors.  If the client
  has specified mode 0 (work from snapshot), we silently correct it.
*/

CoinPackedMatrix *setupRowCopy (int mode, bool useObj,
				double cutoff, double offset,
				const CoinPackedMatrix *rowCopy_,
				int numberRows_, int numberColumns_,
				const double *const rowLower_,
				const double *const rowUpper_,
				const OsiSolverInterface &si,
				double *rowLower, double *rowUpper)
{
# if CGL_DEBUG > 0
  std::cout << "Entering setupRowCopy, mode " << mode << ", " ;
  if (!rowCopy_) std::cout << "no " ;
  std::cout << "row copy." << std::endl ;
# endif

  int nCols = si.getNumCols() ;
  int nRows = -1 ;
  CoinPackedMatrix *rowCopy = NULL ;
/*
  Make the working copy from the model in the si. If we're using the
  objective, add in the objective coefficients.
*/
  if (!rowCopy_) {

    nRows = si.getNumRows(); 

    CoinMemcpyN(si.getRowLower(),nRows,rowLower) ;
    CoinMemcpyN(si.getRowUpper(),nRows,rowUpper) ;

    if (!useObj) {
      rowCopy = new CoinPackedMatrix(*si.getMatrixByRow()) ;
    } else {
      rowCopy = new CoinPackedMatrix(*si.getMatrixByRow(),1,nCols,false) ;
      const double *objective = si.getObjCoefficients() ;
      const double objmul = si.getObjSense() ;
      int *columns = new int[nCols] ;
      double *elements = new double[nCols] ;
      int kk = 0 ;
      for (int j = 0 ; j < nCols ; j++) {
	if (objective[j]) {
	  elements[kk] = objmul*objective[j] ;
	  columns[kk++] = j ;
	}
      }
      rowCopy->appendRow(kk,columns,elements) ;
      delete[] columns ;
      delete[] elements ;
      rowLower[nRows] = -DBL_MAX ;
      rowUpper[nRows] = cutoff+offset ;
      nRows++ ;
    }
  } else {
/*
  Make the working copy from the snapshot. It must be true that the snapshot
  has the same number of columns as the constraint system in the solver. The
  asserts are sanity checks:
    * the column count could go bad because the user has changed the model
      outside of CglProbing ;
    * rowLower and rowUpper are allocated to the number of rows in the model
      plus 1 (in case we want the objective), so we need to be sure there's
      enough space ;
    * the final check of row counts is simple internal consistency.

  TODO: Another good reason to hive off snapshot as a separate class.
	-- lh, 100917 --
*/
    assert(nCols == numberColumns_) ;
    assert(nRows <= si.getNumRows()+1) ;
    assert(rowCopy_->getNumRows() == numberRows_) ;

    nRows = numberRows_ ;
    rowCopy = new CoinPackedMatrix(*rowCopy_) ;
    CoinMemcpyN(rowLower_,nRows,rowLower) ;
    CoinMemcpyN(rowUpper_,nRows,rowUpper) ;
    if (useObj) {
      rowLower[nRows-1] = -DBL_MAX ;
      rowUpper[nRows-1] = cutoff+offset ;
    }
  }
# if CGL_DEBUG > 0
  if (rowCopy) {
    std::cout << "  verifying matrix." << std::endl ;
    rowCopy->verifyMtx(2) ;
  }
  std::cout << "Leaving setupRowCopy." << std::endl ;
# endif
  return (rowCopy) ;
}


/*
  Remove rows with too many elements and rows that are deemed useless for
  further probing.  Try to strengthen rows in place, if allowed. Deal with
  fixed variables. We're going to do this in place by editing the internal
  data structures of rowCopy. realRows will become a cross-reference from
  row indices in the groomed model back to the original indices.

  Strengthening in place is highly specialised. We're looking for constraints
  that satisfy:
    * All variables binary
    * All coefficients negative (N), except for a single positive (t), and
      b<i> > 0. In this case, we can convert
        a<it>x<t> + SUM{N}a<ik>x<k> <= b<i>
      to
        (a<it>-b<i>)x<t> + SUM{N}a<ik>x<k> <= 0
    * All coefficients positive (P), except for a single negative (t), and
      blow<i> < 0. In this case, we can convert
        a<it>x<t> + SUM{P}a<ik>x<k> >= blow<i>
      to
        (a<it>-blow<i>)x<t> + SUM{P}a<ik>x<k> >= 0
  See the typeset documentation for the full derivation.

  The original code here went to some trouble to avoid deleting the objective
  for density, but it would delete it for an ineffective bound. But we've
  gone to some trouble back in gutsOfGenerateCuts to use the objective only
  if there's a reasonable cutoff. Let's not second guess the decision here.

  While we're here, we will sort the coefficients of each row so that all
  negative coefficients precede all positive coefficients. Record the break
  points in rowStartPos.

  An alternate code block for deleting constraints on bounds would delete
  any constraint where an infinite (1.0e3 :-) bound on a variable would
  lead to an infinite row LB or UB. This is a recipe for deleting the
  entire model, when run prior to an initial round of bound propagation. It's
  been disabled since r675 (2008). I've chopped it entirely. Something useful
  could be done along this line, but only after an initial round of bound
  propagation. Arguably, an initial round of bound propagation would help all
  of CglProbing and should be done in any event.

  Returns true if we're good to proceed, false if the groomed model is not
  worth pursuing.
*/

bool groomModel (
		 bool useObj, int maxRowLen,
		 const OsiSolverInterface &si, const char *const intVar,
		 CoinPackedMatrix *rowCopy,
		 double *const rowLower, double *const rowUpper,
		 const double *const colLower, const double *const colUpper,
		 int *&realRows, CoinBigIndex *&rowStartPos,
		 const CglTreeInfo *info
	        )
{
# if CGL_DEBUG > 0
  std::cout << "Entering groomModel." << std::endl ;
# endif

  realRows = 0 ;
  rowStartPos = 0 ;

  const double weakRhs = 1.0e20 ;
  const double strengthenRhs = 1.0e10 ;
/*
  In the first half (row selection), only elements is modified. In the second
  half (fixed variable removal), all are subject to change.
*/
  CoinBigIndex *rowStart = rowCopy->getMutableVectorStarts() ;
  int *rowLength = rowCopy->getMutableVectorLengths(); 
  double *elements = rowCopy->getMutableElements() ;
  int *column = rowCopy->getMutableIndices() ;
/*
  nRealRows excludes the objective, if we're using it.
*/
  int nRows = rowCopy->getNumRows() ;
  const int nRealRows = (useObj)?(nRows-1):(nRows) ;
  const int nCols = rowCopy->getNumCols() ;
  CoinBigIndex nElements = rowCopy->getNumElements() ;

/*
  Coefficient strengthening is allowed only on the first call of the
  generator, and only if the user has provided an array to record the
  strengthened rows.
*/
  bool allowStrengthen = (info->strengthenRow && !info->pass) ;
  int *which = new int [nRealRows] ;

  int nDense = 0 ;
  int nWeak = 0 ;
/*
  Preliminary assessment for dense rows and free rows.  Count the number of
  coefficients that we'll remove. We're looking for rows that are too long or
  have row bounds so weak as to be ineffective.  If the total coefficients
  in dense rows amounts to less than 1/10 of the matrix, just deal with it
  (maxRowLen set to nCols won't exclude anything).

  At the end of this markup, which[i] holds the row length, or (nCols+1) for
  rows we're going to delete for weak rhs.
*/
  for (int i = 0 ; i < nRealRows ; i++) {
    const int leni = rowLength[i] ;
    if ((rowLower[i] < -weakRhs) && (rowUpper[i] > weakRhs)) {
      nWeak += leni ;
      which[i] = nCols+1 ;
    } else {
      if (leni > maxRowLen) nDense += leni ;
      which[i] = leni ;
    }
  }
  if (nDense*10 < nElements) maxRowLen = nCols ;

/*
  Walk the rows again. Rows that should be deleted will come up positive
  when we subtract maxRowLen and can be summarily dismissed.  Note that we're
  converting which[] as we go, so that at the end it will contain the
  indices of the constraints to be deleted in a block at the beginning.
*/
  int nDelete = 0 ;
  int nKeep = 0 ;
  for (int i = 0 ; i < nRealRows ; i++) {
    if (which[i]-maxRowLen > 0) {
      which[nDelete++] = i ;
      continue ;
    }
/*
  Process the row further only if we can strengthen rows in place. Avoid
  range constraints to avoid complications. See the comments at the head
  for details.
 
  Shift the rhs values here --- it has to be done for any constraint that we
  keep.
*/
    const double bi = rowUpper[i] ;
    const double blowi = rowLower[i] ;
    rowLower[nKeep] = blowi ;
    rowUpper[nKeep] = bi ;
    nKeep++ ;

    if (!allowStrengthen || (blowi > -weakRhs && bi < weakRhs)) continue ;
/*
  Assess the row. If we find anything other than binary variables, we're
  immediately done.

  We're lumping a<ij> = 0 in with a<ij> < 0, but arguably we shouldn't see
  coefficients of zero. Introducing a whole lot of checks into this loop is
  counterproductive; that's why we don't abort as soon as nPlus and nMinus
  both exceed 1.
*/
    int nPlus = 0 ;
    int nMinus = 0 ;
    const CoinBigIndex rstarti = rowStart[i] ;
    const int leni = rowLength[i] ;
    const CoinBigIndex rendi = rstarti+leni ;
    for (CoinBigIndex jj = rstarti ; jj < rendi ; jj++) {
      int j = column[jj] ;
      if (intVar[j] == 1) {
	double value = elements[jj] ;
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
  If there's one positive coefficient, find it and reduce it, given a
  reasonable value for b<i>. If there's one negative coefficient, find it
  and increase it, given a reasonable value for blow<i>. Note that only one of
  these will execute, even if nPlus = nMinus = 1, because only one of b<i> or
  blow<i> will have a reasonable value. We just don't know which one.

  Be careful with rhs values here! Remember that shift up above.
*/
    double effectiveness = 0.0 ;
    if (nPlus == 1 && bi > 0.0 && bi < strengthenRhs) {
      for (CoinBigIndex jj = rstarti ; jj < rendi ; jj++) {
	double value = elements[jj] ;
	if (value > 0.0) {
	  elements[jj] -= bi ;
#         if CGL_DEBUG > 1
	  std::cout
	    << "      row " << i << ": decrease a<" << i << ","
	    << column[jj] << "> from " << value << " to "
	    << elements[jj] << "." << std::endl ;
#	  endif
	}
	effectiveness += fabs(value) ;
      }
      rowUpper[nKeep-1] = 0.0 ;
    }
    if (nMinus == 1 && blowi < 0.0 && blowi > -strengthenRhs) {
      for (CoinBigIndex jj = rstarti ; jj < rendi ; jj++) {
	double value = elements[jj] ;
	if (value < 0.0) {
	  elements[jj] -= blowi ;
#         if CGL_DEBUG > 1
	  std::cout
	    << "      row " << i << ": increase a<" << i << ","
	    << column[jj] << "> from " << value << " to "
	    << elements[jj] << "." << std::endl ;
#	  endif
	}
	effectiveness += fabs(value) ;
      }
      rowLower[nKeep-1] = 0.0 ;
    }
/*
  If we actually strengthened a coefficient, stash the strengthened row in
  strengthenRow for return to the caller.
*/
    if (effectiveness) {
      OsiRowCut *rc = new OsiRowCut() ;
      rc->setLb(rowLower[nKeep-1]) ;
      rc->setUb(rowUpper[nKeep-1]) ;
      rc->setRow(leni,column+rstarti,elements+rstarti,false) ;
      rc->setEffectiveness(effectiveness) ;
      assert (!info->strengthenRow[i]) ;
      info->strengthenRow[i] = rc ;
    }
  }
/*
  Deal with the objective, if we're using it. In particular, if we've deleted
  constraints, we need to shift the rhs entries for the objective.
*/
  if (useObj && nDelete)
  { rowLower[nKeep] = rowLower[nRealRows] ;
    rowUpper[nKeep] = rowUpper[nRealRows] ;
  }
/*
  Did we decide to delete anything (but not everything) in the processing
  loop?  If so, deal with it.  We need to physically delete the rows, of
  course. There may also be bookkeeping, if we're trying to strengthen rows
  in place.  In particular, strengthenRow will be passed back to the caller
  and is indexed using the original row index. We need a cross-reference,
  realRows, to keep track.

  realRows is a tradeoff -- it's bigger than necessary, by nDelete entries,
  but nDelete will typically be small and this approach to building the
  xref saves some work (we need a working array of size nRows, and we want
  to retain the block of indices at the front of which).

  But keep in mind that we can't strengthen the objective in place!
*/
  if (nDelete > 0 && nDelete < nRealRows) {
    if (info->strengthenRow) {
      realRows = new int [nRows] ;
      CoinZeroN(realRows,nRows) ;
      for (int k = 0 ; k < nDelete ; k++) {
        realRows[which[k]] = -1 ;
      }
      if (useObj) realRows[nRealRows] = -1 ;
      int inew = 0 ;
      for (int iold = 0 ; iold < nRows ; iold++) {
	if (!realRows[iold]) {
	  realRows[inew++] = iold ;
	}
      }
      assert(inew == nKeep) ;
    }
    rowCopy->deleteRows(nDelete,which) ;
    nRows = rowCopy->getNumRows() ;
    assert((useObj && (nRows == nKeep+1)) || (!useObj && (nRows == nKeep))) ;
  }
  delete[] which ;
/*
  Do we have any real constraints left to work with? (The objective doesn't
  count here.) If not, clean up and go home.
*/
  if (!nKeep) {
#   if CGL_DEBUG > 0
      std::cout
	<< "  All rows unsuitable for probing! Cleaning up." << std::endl ;
#   endif
    delete[] realRows ;
    return (false) ;
  }
/*
  So far, so good. We have something that looks like a reasonable collection
  of rows. Time to deal with fixed variables. We're going to walk the rows,
  compressing out the coefficients of the fixed variables, and we'll sort
  the remaining coefficients so that all negative coefficients precede
  all positive coefficients.  There will be no gaps when we're done.

  NOTE that columns are not physically removed.

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

  CoinBigIndex newSize = 0 ;
  int *columnPos = new int [nCols] ;
  double *elementsPos = new double [nCols] ;
  rowStartPos = new CoinBigIndex [nRows] ;
  for (int i = 0 ; i < nRows ; i++) {
    double offset = 0.0 ;
    const CoinBigIndex rstarti = rowStart[i] ;
    const CoinBigIndex rendi = rstarti+rowLength[i] ;
    rowStart[i] = newSize ;
    const CoinBigIndex saveNewStarti = newSize ;
    int nPos = 0 ;
    for (CoinBigIndex jj = rstarti ; jj < rendi ; jj++) { 
      int j = column[jj] ;
      double value = elements[jj] ;
      if (colUpper[j] > colLower[j]) {
	if (value < 0.0) {
	  elements[newSize] = value ;
	  column[newSize++] = j ;
	} else {
	  elementsPos[nPos] = value ;
	  columnPos[nPos++] = j ;
	}
      } else {
	offset += colUpper[j]*value ;
      }
    }
    rowStartPos[i] = newSize ;
    for (int k = 0 ; k < nPos ; k++) {
      elements[newSize] = elementsPos[k] ;
      column[newSize++] = columnPos[k] ;
    }
    rowLength[i] = newSize-saveNewStarti ;
    if (offset) {
      if (rowLower[i] > -weakRhs)
	rowLower[i] -= offset ;
      if (rowUpper[i] < weakRhs)
	rowUpper[i] -= offset ;
    }
  }
  delete [] columnPos ;
  delete [] elementsPos ;
  rowStart[nRows] = newSize ;
  rowCopy->setNumElements(newSize) ;

  return (true) ;
} 

}  // end file local namespace


//-------------------------------------------------------------------
// Generate column, disaggregation, and implication cuts
//------------------------------------------------------------------- 

/*
  The traditional const method. Row and column cuts will be returned in cs.
*/
void CglProbing::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
			      const CglTreeInfo info2) const
{

# if CGL_DEBUG > 0
  std::cout
    << "Entering CglProbing::generateCuts, matrix " << si.getNumRows()
    << " x " << si.getNumCols() << "." << std::endl ;
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
  working. Note that the row bound arrays are allocated with one extra
  position, in case we want to incorporate the objective as a constraint while
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
  Do the work. There's no guarantee about the state of the row and column
  bounds arrays if we come back infeasible.
*/
  int ninfeas =
	gutsOfGenerateCuts(si,cs,rowLower,rowUpper,colLower,colUpper,&info) ;
# if CGL_DEBUG > 0
  { int errs = 0 ;
    const CoinPackedMatrix *mtx = si.getMatrixByRow() ;
    errs = mtx->verifyMtx(2) ;
    if (errs > 0) {
      std::cout
        << "    generateCuts: verifyMtx failed, "
	<< errs << " errors." << std::endl ;
      assert (false) ;
    } else {
      std::cout << "    generateCuts: matrix verified." << std::endl ;
    }
  }
# endif
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

  TODO: Why are we fiddling with colLower_, colUpper_? Track this down.
        -- lh, 110212 --
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
  Note that this is not a const method (contrast with generateCuts).
  
  TODO: Note that the CglTreeInfo object is not const either. In generateCuts,
	it's passed in as a const and a nonconst clone is created. I need to
	look into that.  -- lh, 100924 --

	Possible reason is use of CglTreeInfo to pass back strengthened rows
	for row replacement during preprocessing.  -- lh, 110211 --

  Generally, check generateCuts for additional comment. I went through it first
  and didn't repeat comments on common code blocks here.
*/
int CglProbing::generateCutsAndModify (const OsiSolverInterface &si, 
				       OsiCuts &cs,
				       CglTreeInfo *info) 
{
# if CGL_DEBUG > 0
  std::cout
    << "Entering CglProbing::generateCutsAndModify, matrix " << si.getNumRows()
    << " x " << si.getNumCols() << "." << std::endl ;
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
  int saveMode = mode_ ;
/*
  Create working arrays for row and column bounds.

  We need nRows+1 because we may use the objective function as a constraint
  in the working system within gutsOfGenerateCuts.
*/
  int nRows = si.getNumRows() ; 
  double *rowLower = new double[nRows+1] ;
  double *rowUpper = new double[nRows+1] ;

  int nCols = si.getNumCols() ; 
  double *colLower = new double[nCols] ;
  double *colUpper = new double[nCols] ;
/*
  Do the work.
*/
# if CGL_DEBUG > 0
  { int errs = 0 ;
    const CoinPackedMatrix *mtx = si.getMatrixByRow() ;
    errs = mtx->verifyMtx(3) ;
    if (errs > 0) {
      std::cout
        << "    generateCutsAndModify: verifyMtx failed pre-cutgen, "
	<< errs << " errors." << std::endl ;
      assert (errs == 0) ;
    }
  }
# endif
/*
  Do the work. There's no guarantees about the state of the row and column
  bounds arrays if we come back infeasible.
*/
  int ninfeas = gutsOfGenerateCuts(si,cs,
				   rowLower,rowUpper,colLower,colUpper,info) ;
# if CGL_DEBUG > 0
  { int errs = 0 ;
    const CoinPackedMatrix *mtx = si.getMatrixByRow() ;
    errs = mtx->verifyMtx(3) ;
    if (errs > 0) {
      std::cout
        << "    generateCutsAndModify: verifyMtx failed post-cutgen, "
	<< errs << " errors." << std::endl ;
      assert (errs == 0) ;
    }
  }
# endif
/*
  Infeasibility is indicated by a stylized infeasible column cut.

  If we're debugging, check whether this is a plausible outcome.
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

  Don't replace bounds supplied by the client if we've gone infeasible.
  There's no guarantee of consistency in the values reported back from
  gutsOfGenerateCuts.
*/
  if (mode_ == 3 && !ninfeas) {
    delete [] rowLower_ ;
    delete [] rowUpper_ ;
    rowLower_ = rowLower ;
    rowUpper_ = rowUpper ;
  } else {
    delete [] rowLower ;
    delete [] rowUpper ;
  }
  if (!ninfeas) {
    delete [] colLower_ ;
    delete [] colUpper_ ;
    colLower_ = colLower ;
    colUpper_ = colUpper ;
  } else {
    delete [] colUpper ;
    delete [] colLower ;
  }

# if CGL_DEBUG > 0
  std::cout
    << "Leaving CglProbing::generateCutsAndModify, infeas "
    << ninfeas << std::endl ;
# endif

  return (ninfeas) ;
}


/*
  NOTE: It's assumed that the arrays passed in for rowLower and rowUpper
	have an extra slot that can be used for bounds on the objective,
	in the event that gutsOfGenerateCuts decides to put it in the matrix.

  rowLower, rowUpper, colLower, and colUpper should be allocated by the caller
  and will be filled in by gutsOfGenerateCuts. If the system is found to be
  infeasible, there's no particular guarantee of consistency for these arrays
  and the content should not be used.
*/
int CglProbing::gutsOfGenerateCuts (const OsiSolverInterface &si, 
				    OsiCuts &cs ,
				    double *rowLower, double *rowUpper,
                                    double *colLower, double *colUpper,
                                    CglTreeInfo *info) const
{

# if CGL_DEBUG > 0
  std::cout
    << "Entering CglProbing::gutsOfGenerateCuts, mode " << mode_ << ", " ;
  if (!info->inTree) std::cout << "not " ;
  std::cout << "in tree, pass " << info->pass << "." << std::endl ;
# endif

  int numberRowCutsBefore = cs.sizeRowCuts() ;

  int ninfeas = 0 ;
  bool feasible = true ;
  int mode = mode_ ;

  int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_ ;
  int maxRowLen = info->inTree ? maxElements_ : maxElementsRoot_ ;
/*
  Get a modifiable array of variable types. Force reevaluation of the
  integer classification (binary vs. general integer) based on current
  bounds.
*/
  int nCols = si.getNumCols(); 
  const char *intVarOriginal = si.getColType(true) ;
  char *intVar = CoinCopyOfArray(intVarOriginal,nCols) ;
/*
  Get modifiable bound arrays and the current solution.
*/
  CoinMemcpyN(si.getColLower(),nCols,colLower) ;
  CoinMemcpyN(si.getColUpper(),nCols,colUpper) ;
  const double *colsol = si.getColSolution() ;
/*
  Stage 1: Preliminary Processing (Common setup?)

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
/*
  If we're at the root of the search tree, scan the constraint system to see
  if there are naturally integer variables that are not declared as integer.
  Convert them to integer type.

  If we lose feasibility here, cut our losses and bail out while it's still
  uncomplicated.
*/
  if (!info->inTree && !info->pass) {
    feasible = analyze(&si,intVar,colLower,colUpper) ;
#   if CGL_DEBUG > 1
    std::cout
     << "  Analyze completed; " << ((feasible)?"feasible":"infeasible")
     << "." << std::endl ;
#   endif
    if (!feasible) {
      delete[] intVar ;
      return (1) ;
    }
  }
/*
  Get the objective offset and try to establish a nontrivial cutoff. The
  client can forbid use of the objective by setting usingObjective_ to -1.
  If the si claims a nontrivial cutoff, check that it's reasonable.

  TODO: Check through the code and make sure we're using cutoff and offset
	correctly. Should add the offset in one place and keep the value
	for later use.  -- lh, 110205 --
*/
  double offset = 0.0 ;
  si.getDblParam(OsiObjOffset,offset) ;
  double cutoff = 0.0 ;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff) ;
  if (!cutoff_available) {
    cutoff = si.getInfinity() ;
  } else {
    cutoff *= si.getObjSense() ;
    if (cutoff > 1.0e30) cutoff_available = false ;
  }
  bool useObj = (cutoff_available && (usingObjective_ >= 1)) ;
  bool useCutoff = (cutoff_available && (usingObjective_ >= 0)) ;
# if CGL_DEBUG > 1
  std::cout << "  cutoff " ;
  if (!cutoff_available)
    std::cout << "not available" ;
  else
    std::cout << cutoff ;
  std::cout << ", offset " << offset << "." << std::endl ;
  if (useCutoff)
    std::cout << "  using cutoff" ;
  else
    std::cout << "  not using cutoff" ;
  if (useObj)
    std::cout << ", using objective as constraint" ;
  else
    std::cout << ", not using objective as constraint" ;
  std::cout << "." << std::endl ;
# endif

/* EXTRA_PROBING_STUFF block (tighten using reduced cost) was here. */

/*
  Create the working row copy from the si's model or from a snapshot (only if
  the user specified mode 0 and a snapshot exists). rowLower and rowUpper are
  valid at completion, but if we initialise the local row copy from a snapshot
  the number of valid entries may not match the number of rows for the system
  held in the si.
*/
  CoinPackedMatrix *rowCopy = setupRowCopy(mode,useObj,cutoff,
  					   offset,rowCopy_,
					   numberRows_,numberColumns_,
					   rowLower_,rowUpper_,
					   si,rowLower,rowUpper) ;
  int nRows = rowCopy->getNumRows() ;
/*
  Groom the rowCopy: delete useless rows, substitute for fixed variables, sort
  coefficients so that negative precede positive, and try some specialised
  coefficient strengthening on the rows that remain. If groomModel comes
  back false, it's an indication that probing isn't worth pursuing. Clean
  up and go home while it's relatively easy.

  realRows is the cross-reference between rows in rowCopy and rows of the
  original model. rowStartPos is the first positive coefficient in a row.

  Note that the number of rows in the local copy is almost certainly not the
  same as the number of rows in the system held in the si. You can't just copy
  rowLower and rowUpper; they must be translated if realRows exists.
*/
  int *realRows = 0 ;
  int *rowStartPos = 0 ;
  feasible = groomModel(useObj,maxRowLen,si,intVar,rowCopy,rowLower,rowUpper,
  			colLower,colUpper,realRows,rowStartPos,info) ;
  nRows = rowCopy->getNumRows() ;
  if (!feasible) {
    delete[] intVar ;
    delete[] realRows ;
    delete[] rowStartPos ;
    delete rowCopy ;
    return (0) ;
  }
/*
  Make a column-major copy, after we've winnowed the rows to the set we want
  to work with. The break out the internal arrays for subsequent use.

  TODO: I'm asking myself, ``Why does the snapshot keep a colum-major
	copy? It doesn't seem to be used.  -- lh, 101021 --
*/
  CoinPackedMatrix * columnCopy = new CoinPackedMatrix(*rowCopy,0,0,true) ;
  const int *column = rowCopy->getIndices() ;
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *rowLength = rowCopy->getVectorLengths(); 
  const double *rowElements = rowCopy->getElements() ;
/*
  Allocate an array to hold indices of the columns we will probe.
*/
  if (!lookedAt_) {
    lookedAt_ = new int[nCols] ;
  }
  numberThisTime_ = 0 ;
/*
  Allocate the row processing record (markR) and the arrays that will hold
  row lhs bounds (minR, maxR).
*/
  int * markR = new int [nRows] ;
  double * minR = new double [nRows] ;
  double * maxR = new double [nRows] ;
/*
  End of setup. Start out by calling tighten to do bound propagation (two pass
  limit).
*/
  ninfeas = tighten(colLower,colUpper,column,rowElements,
		    rowStart,rowStartPos,rowLength,
		    rowLower,rowUpper,nRows,nCols,
		    intVar,2,primalTolerance_) ;
/*
  If we've lost feasibility, do nothing more. Otherwise, record the new bounds
  as column cuts and continue with setup and probing.
*/
  if (!ninfeas) {
    makeColCuts(nCols,cs,intVar,colsol,si.getColLower(),si.getColUpper(),
		si.getColSolution(),colLower,colUpper) ;
#   if CGL_DEBUG > 0
    std::cout
      << "    tighten found " << cs.sizeColCuts() << " column cuts."
      << std::endl ;
    if (cs.sizeColCuts()) {
      CglProbingDebug::checkBounds(si,cs.colCut((cs.sizeColCuts()-1))) ;
    }
#   endif
/*
  We've calculated changes in column bounds due to propagation. To go further,
  we need to be able to probe.
*/
    if (maxProbe > 0) {    // block extends to end of if (mode)
      numberThisTime_ = 0 ;
/*
  Calculate L<i> and U<i> for each row. If one or both are suitable to
  process further, markR[i] is set to -1. If both are too large to warrant
  further processing, markR[i] is set to -2.

  TODO: Why not do this as part of tighten()? Tighten must calculate row
	bounds in order to propagate column bounds. It's not like we're
	saving any work, other than transcribing the bounds into minR and
	maxR. On a purely mechanical level, tighten does not fill in minR,
	maxR, markR.  -- lh, 100923 --

*/
      calcRowBounds(colLower,colUpper,column,rowElements,
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
#       ifdef ZEROFAULT
	std::memset(array,0,sizeof(double_int_pair)*nCols) ;
#	endif
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
  lookedAt_ now contains a list of indices of variables to probe. Rows
  worth processing are tagged with -1 in markR. Do the probing.
*/
      ninfeas = probe(si,cs,colLower,colUpper,rowCopy,columnCopy,
		      rowStartPos,realRows,rowLower,rowUpper,
		      intVar,minR,maxR,markR,info,
		      useObj,useCutoff,cutoff) ;
    }
  }
/*
  Done! Time to clean up and go home.
*/
  delete [] markR ;
  delete [] minR ;
  delete [] maxR ;

  delete rowCopy ;
  delete columnCopy ;

  if (rowCopy_) {
    delete [] rowLower ;
    delete [] rowUpper ;
  }

  delete [] intVar ;
  delete [] rowStartPos ;
  delete [] realRows ;

/*
  jjf: put back unreasonable bounds on integer variables

  TODO: If we're infeasible, is this necessary? Better to just not use the
	bounds.   -- lh, 110217 --
*/
  const double *trueLower = si.getColLower() ;
  const double *trueUpper = si.getColUpper() ;
  if (!ninfeas) {
    for (int i = 0 ; i < nCols ; i++) {
      if (intVarOriginal[i] == 2) {
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
/*
  If we're requested to set the global cut flag, do it.
*/
  if (!info->inTree &&
      ((info->options&4) == 4 || ((info->options&8) && !info->pass))) {
    int numberRowCutsAfter = cs.sizeRowCuts() ;
    for (int i=numberRowCutsBefore;i<numberRowCutsAfter;i++)
      cs.rowCutPtr(i)->setGloballyValid() ;
  }

# if CGL_DEBUG > 0
  std::cout
    << "Leaving CglProbing::gutsOfGenerateCuts, "
    << cs.sizeRowCuts() << " row cuts, " << cs.sizeColCuts()
    << " col cuts, infeas " << ninfeas << std::endl ;
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  if (debugger) {
    debugger->validateCuts(cs,0,cs.sizeCuts()) ;
# if CGL_DEBUG > 2
    for (OsiCuts::iterator oneCut = cs.begin() ; oneCut != cs.end() ; oneCut++)
      (*oneCut)->print() ;
# endif
  } else if (si.getRowCutDebuggerAlways()) {
    std::cout << "  !! Not on optimal path." << std::endl ;
  }
# endif

  return ninfeas ;
}


/*
  Set mode

  Take care to change only the basic mode, without affecting higher bits.
*/
void CglProbing::setMode(int mode)
{
  if (mode >= 0 && mode < 3) {
    mode_ &= ~15 ;
    mode_ |= mode ;
  }
}

// Return only the basic mode
int CglProbing::getMode() const
{
  return mode_&15 ;
}


// Set maximum number of passes per node
void CglProbing::setMaxPass(int value)
{
  if (value>0)
    maxPass_ = value ;
}


// Get maximum number of passes per node
int CglProbing::getMaxPass() const
{
  return maxPass_ ;
}


// Set log level
void CglProbing::setLogLevel(int value)
{
  if (value >= 0)
    logLevel_ = value ;
}


// Get log level
int CglProbing::getLogLevel() const
{
  return logLevel_ ;
}


// Set maximum number of unsatisfied variables to look at
void CglProbing::setMaxProbe(int value)
{
  if (value >= 0)
    maxProbe_ = value ;
}


// Get maximum number of unsatisfied variables to look at
int CglProbing::getMaxProbe() const
{
  return maxProbe_ ;
}


// Set maximum number of variables to look at in one probe
void CglProbing::setMaxLook(int value)
{
  if (value >= 0)
    maxStack_ = value ;
}


// Get maximum number of variables to look at in one probe
int CglProbing::getMaxLook() const
{
  return maxStack_ ;
}


// Set maximum number of elements in row for scan
void CglProbing::setMaxElements(int value)
{
  if (value > 0)
    maxElements_ = value ;
}


// Get maximum number of elements in row for scan
int CglProbing::getMaxElements() const
{
  return maxElements_ ;
}


// Set maximum number of passes per node (root node)
void CglProbing::setMaxPassRoot(int value)
{
  if (value > 0)
    maxPassRoot_ = value ;
}


// Get maximum number of passes per node (root node)
int CglProbing::getMaxPassRoot() const
{
  return maxPassRoot_ ;
}


// Set maximum number of unsatisfied variables to look at (root node)
void CglProbing::setMaxProbeRoot(int value)
{
  if (value > 0)
    maxProbeRoot_ = value ;
}


// Get maximum number of unsatisfied variables to look at (root node)
int CglProbing::getMaxProbeRoot() const
{
  return maxProbeRoot_ ;
}


// Set maximum number of variables to look at in one probe (root node)
void CglProbing::setMaxLookRoot(int value)
{
  if (value > 0)
    maxStackRoot_ = value ;
}


// Get maximum number of variables to look at in one probe (root node)
int CglProbing::getMaxLookRoot() const
{
  return maxStackRoot_ ;
}


// Set maximum number of elements in row for scan (root node)
void CglProbing::setMaxElementsRoot(int value)
{
  if (value > 0)
    maxElementsRoot_ = value ;
}


// Get maximum number of elements in row for scan (root node)
int CglProbing::getMaxElementsRoot() const
{
  return maxElementsRoot_ ;
}


// Set whether to use objective
void CglProbing::setUsingObjective(int yesNo)
{
  usingObjective_ = yesNo ;
}


// Get whether objective is being used
int CglProbing::getUsingObjective() const
{
  return usingObjective_ ;
}


// Decide whether to do row cuts
void CglProbing::setRowCuts(int type)
{
  if (type > -5 && type < 5)
    rowCuts_ = type ;
}


// Returns row cuts generation type
int CglProbing::rowCuts() const
{
  return rowCuts_ ;
}


// Returns tight lower
const double * CglProbing::tightLower() const
{
  return colLower_ ;
}


// Returns tight upper
const double * CglProbing::tightUpper() const
{
  return colUpper_ ;
}


// Returns relaxed Row lower
const double * CglProbing::relaxedRowLower() const
{
  return rowLower_ ;
}


// Returns relaxed Row upper
const double * CglProbing::relaxedRowUpper() const
{
  return rowUpper_ ;
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

  numberRows_ = 0 ;
  numberColumns_ = 0 ;
  rowCopy_ = NULL ;
  columnCopy_ = NULL ;
  rowLower_ = NULL ;
  rowUpper_ = NULL ;
  colLower_ = NULL ;
  colUpper_ = NULL ;
  numberIntegers_ = 0 ;
  number01Integers_ = 0 ;
  numberThisTime_ = 0 ;
  totalTimesCalled_ = 0 ;
  lookedAt_ = NULL ;
  cutVector_ = NULL ;
  tightenBounds_ = NULL ;

# ifdef CLIQUE_ANALYSIS
  numberCliques_=0 ;
  cliqueType_=NULL ;
  cliqueStart_=NULL ;
  cliqueEntry_=NULL ;
  oneFixStart_=NULL ;
  zeroFixStart_=NULL ;
  endFixStart_=NULL ;
  whichClique_=NULL ;
  cliqueRow_=NULL ;
  cliqueRowStart_=NULL ;
# endif

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing (const CglProbing & rhs)
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
  numberRows_ = rhs.numberRows_ ;
  numberColumns_ = rhs.numberColumns_ ;
  if (rhs.rowCopy_) {
    rowCopy_= new CoinPackedMatrix(*(rhs.rowCopy_)) ;
    columnCopy_= new CoinPackedMatrix(*(rhs.columnCopy_)) ;
    rowLower_ = new double[numberRows_] ;
    CoinMemcpyN(rhs.rowLower_,numberRows_,rowLower_) ;
    rowUpper_ = new double[numberRows_] ;
    CoinMemcpyN(rhs.rowUpper_,numberRows_,rowUpper_) ;
    colLower_ = new double[numberColumns_] ;
    CoinMemcpyN(rhs.colLower_,numberColumns_,colLower_) ;
    colUpper_ = new double[numberColumns_] ;
    CoinMemcpyN(rhs.colUpper_,numberColumns_,colUpper_) ;
    int i ;
    numberIntegers_ = rhs.numberIntegers_ ;
    number01Integers_ = rhs.number01Integers_ ;
    cutVector_ = new disaggregation [number01Integers_] ;
    CoinMemcpyN(rhs.cutVector_,number01Integers_,cutVector_) ;
    for (i = 0 ; i < number01Integers_ ; i++) {
      if (cutVector_[i].index) {
	cutVector_[i].index =
	    CoinCopyOfArray(rhs.cutVector_[i].index,cutVector_[i].length) ;
      }
    }
  } else {
    rowCopy_ = NULL ;
    columnCopy_ = NULL ;
    rowLower_ = NULL ;
    rowUpper_ = NULL ;
    colLower_ = NULL ;
    colUpper_ = NULL ;
    numberIntegers_ = 0 ;
    number01Integers_ = 0 ;
    cutVector_ = NULL ;
  }
  numberThisTime_ = rhs.numberThisTime_ ;
  totalTimesCalled_ = rhs.totalTimesCalled_ ;
  if (numberColumns_)
    lookedAt_ = CoinCopyOfArray(rhs.lookedAt_,numberColumns_) ;
  else
    lookedAt_ = NULL ;
  if (rhs.tightenBounds_) {
    assert (numberColumns_) ;
    tightenBounds_ = CoinCopyOfArray(rhs.tightenBounds_,numberColumns_) ;
  } else {
    tightenBounds_ = NULL ;
  }
# ifdef CLIQUE_ANALYSIS
  numberCliques_=rhs.numberCliques_ ;
  if (numberCliques_) {
    cliqueType_ = new cliqueType [numberCliques_] ;
    CoinMemcpyN(rhs.cliqueType_,numberCliques_,cliqueType_) ;
    cliqueStart_ = new int [numberCliques_+1] ;
    CoinMemcpyN(rhs.cliqueStart_,(numberCliques_+1),cliqueStart_) ;
    int n = cliqueStart_[numberCliques_] ;
    cliqueEntry_ = new cliqueEntry [n] ;
    CoinMemcpyN(rhs.cliqueEntry_,n,cliqueEntry_) ;
    oneFixStart_ = new int [numberColumns_] ;
    CoinMemcpyN(rhs.oneFixStart_,numberColumns_,oneFixStart_) ;
    zeroFixStart_ = new int [numberColumns_] ;
    CoinMemcpyN(rhs.zeroFixStart_,numberColumns_,zeroFixStart_) ;
    endFixStart_ = new int [numberColumns_] ;
    CoinMemcpyN(rhs.endFixStart_,numberColumns_,endFixStart_) ;
    int n2=-1 ;
    for (int i=numberColumns_-1;i>=0;i--) {
      if (oneFixStart_[i]>=0) {
	n2=endFixStart_[i] ;
	break ;
      }
    }
    assert (n==n2) ;
    whichClique_ = new int [n] ;
    CoinMemcpyN(rhs.whichClique_,n,whichClique_) ;
    if (rhs.cliqueRowStart_) {
      cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_,numberRows_+1) ;
      n=cliqueRowStart_[numberRows_] ;
      cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_,n) ;
    } else {
      cliqueRow_=NULL ;
      cliqueRowStart_=NULL ;
    }
  } else {
    cliqueType_=NULL ;
    cliqueStart_=NULL ;
    cliqueEntry_=NULL ;
    oneFixStart_=NULL ;
    zeroFixStart_=NULL ;
    endFixStart_=NULL ;
    cliqueRow_=NULL ;
    cliqueRowStart_=NULL ;
    whichClique_=NULL ;
  }
# endif    // CLIQUE_ANALYSIS
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *CglProbing::clone() const
{
  return new CglProbing(*this) ;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglProbing::~CglProbing ()
{
  // free memory
  delete [] rowLower_ ;
  delete [] rowUpper_ ;
  delete [] colLower_ ;
  delete [] colUpper_ ;
  delete rowCopy_ ;
  delete columnCopy_ ;
  delete [] lookedAt_ ;
  if (cutVector_) {
    for (int i = 0 ; i < number01Integers_ ; i++) {
      delete [] cutVector_[i].index ;
    }
    delete [] cutVector_ ;
  }
  delete [] tightenBounds_ ;
# ifdef CLIQUE_ANALYSIS
  delete [] cliqueType_ ;
  delete [] cliqueStart_ ;
  delete [] cliqueEntry_ ;
  delete [] oneFixStart_ ;
  delete [] zeroFixStart_ ;
  delete [] endFixStart_ ;
  delete [] whichClique_ ;
  delete [] cliqueRow_ ;
  delete [] cliqueRowStart_ ;
# endif
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglProbing &
CglProbing::operator=(
                                         const CglProbing& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs) ;
    primalTolerance_ = rhs.primalTolerance_ ;
    numberRows_ = rhs.numberRows_ ;
    numberColumns_ = rhs.numberColumns_ ;
    delete [] rowLower_ ;
    delete [] rowUpper_ ;
    delete [] colLower_ ;
    delete [] colUpper_ ;
    delete rowCopy_ ;
    delete columnCopy_ ;
    delete [] lookedAt_ ;
    delete [] tightenBounds_ ;

#   ifdef CLIQUE_ANALYSIS
    delete [] cliqueType_ ;
    delete [] cliqueStart_ ;
    delete [] cliqueEntry_ ;
    delete [] oneFixStart_ ;
    delete [] zeroFixStart_ ;
    delete [] endFixStart_ ;
    delete [] whichClique_ ;
    delete [] cliqueRow_ ;
    delete [] cliqueRowStart_ ;
#   endif

    mode_ = rhs.mode_ ;
    rowCuts_ = rhs.rowCuts_ ;
    maxPass_ = rhs.maxPass_ ;
    logLevel_ = rhs.logLevel_ ;
    maxProbe_ = rhs.maxProbe_ ;
    maxStack_ = rhs.maxStack_ ;
    maxElements_ = rhs.maxElements_ ;
    maxPassRoot_ = rhs.maxPassRoot_ ;
    maxProbeRoot_ = rhs.maxProbeRoot_ ;
    maxStackRoot_ = rhs.maxStackRoot_ ;
    maxElementsRoot_ = rhs.maxElementsRoot_ ;
    usingObjective_ = rhs.usingObjective_ ;
    if (rhs.rowCopy_) {
      rowCopy_= new CoinPackedMatrix(*(rhs.rowCopy_)) ;
      columnCopy_= new CoinPackedMatrix(*(rhs.columnCopy_)) ;
      rowLower_ = new double[numberRows_] ;
      CoinMemcpyN(rhs.rowLower_,numberRows_,rowLower_) ;
      rowUpper_ = new double[numberRows_] ;
      CoinMemcpyN(rhs.rowUpper_,numberRows_,rowUpper_) ;
      colLower_ = new double[numberColumns_] ;
      CoinMemcpyN(rhs.colLower_,numberColumns_,colLower_) ;
      colUpper_ = new double[numberColumns_] ;
      CoinMemcpyN(rhs.colUpper_,numberColumns_,colUpper_) ;
      numberIntegers_ = rhs.numberIntegers_ ;
      number01Integers_ = rhs.number01Integers_ ;
      for (int i = 0 ; i < number01Integers_ ; i++) {
        delete [] cutVector_[i].index ;
      }
      delete [] cutVector_ ;
      cutVector_ = new disaggregation [number01Integers_] ;
      CoinMemcpyN(rhs.cutVector_,number01Integers_,cutVector_) ;
      for (int i = 0 ; i < number01Integers_ ; i++) {
        if (cutVector_[i].index) {
          cutVector_[i].index =
	      CoinCopyOfArray(rhs.cutVector_[i].index,cutVector_[i].length) ;
        }
      }
    } else {
      rowCopy_ = NULL ;
      columnCopy_ = NULL ;
      rowLower_ = NULL ;
      rowUpper_ = NULL ;
      colLower_ = NULL ;
      colUpper_ = NULL ;
      numberIntegers_ = 0 ;
      number01Integers_ = 0 ;
      cutVector_ = NULL ;
    }
    numberThisTime_ = rhs.numberThisTime_ ;
    totalTimesCalled_ = rhs.totalTimesCalled_ ;
    if (numberColumns_)
      lookedAt_ = CoinCopyOfArray(rhs.lookedAt_,numberColumns_) ;
    else
      lookedAt_ = NULL ;
    if (rhs.tightenBounds_) {
      assert (numberColumns_) ;
      tightenBounds_ = CoinCopyOfArray(rhs.tightenBounds_,numberColumns_) ;
    } else {
      tightenBounds_ = NULL ;
    }
#   ifdef CLIQUE_ANALYSIS
    numberCliques_ = rhs.numberCliques_ ;
    if (numberCliques_) {
      cliqueType_ = new cliqueType [numberCliques_] ;
      CoinMemcpyN(rhs.cliqueType_,numberCliques_,cliqueType_) ;
      cliqueStart_ = new int [numberCliques_+1] ;
      CoinMemcpyN(rhs.cliqueStart_,(numberCliques_+1),cliqueStart_) ;
      int n = cliqueStart_[numberCliques_] ;
      cliqueEntry_ = new cliqueEntry [n] ;
      CoinMemcpyN(rhs.cliqueEntry_,n,cliqueEntry_) ;
      oneFixStart_ = new int [numberColumns_] ;
      CoinMemcpyN(rhs.oneFixStart_,numberColumns_,oneFixStart_) ;
      zeroFixStart_ = new int [numberColumns_] ;
      CoinMemcpyN(rhs.zeroFixStart_,numberColumns_,zeroFixStart_) ;
      endFixStart_ = new int [numberColumns_] ;
      CoinMemcpyN(rhs.endFixStart_,numberColumns_,endFixStart_) ;
      int n2=-1 ;
      for (int i=numberColumns_-1;i>=0;i--) {
	if (oneFixStart_[i]>=0) {
	  n2=endFixStart_[i] ;
	  break ;
	}
      }
      assert (n==n2) ;
      whichClique_ = new int [n] ;
      CoinMemcpyN(rhs.whichClique_,n,whichClique_) ;
      if (rhs.cliqueRowStart_) {
        cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_,numberRows_+1) ;
        n=cliqueRowStart_[numberRows_] ;
        cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_,n) ;
      } else {
        cliqueRow_=NULL ;
        cliqueRowStart_=NULL ;
      }
    } else {
      cliqueType_=NULL ;
      cliqueStart_=NULL ;
      cliqueEntry_=NULL ;
      oneFixStart_=NULL ;
      zeroFixStart_=NULL ;
      endFixStart_=NULL ;
      whichClique_=NULL ;
      cliqueRow_=NULL ;
      cliqueRowStart_=NULL ;
    }
#   endif
  }
  return *this ;
}


/*
  Returns true if CglProbing thinks it's allowed to generate row cuts in the
  search tree (rather than just at the root node).  Provided so the client
  can know if the matrix will change in the tree.  Really meant so column cut
  generators can still be active without worrying code.  Default is true

  The above is enforced in generateCuts and generateCutsAndModify. Note that
  this method may well return an incorrect answer if called during cut
  generation.
*/
bool 
CglProbing::mayGenerateRowCutsInTree() const
{
  return rowCuts_ > 0 ;
}

// Mark variables to be tightened
void CglProbing::tightenThese (const OsiSolverInterface &solver, int number,
			       const int * which)
{
  delete [] tightenBounds_ ;
  int numberColumns = solver.getNumCols() ;
  if (numberColumns_)
    assert (numberColumns_==numberColumns) ;
  tightenBounds_ = new char [numberColumns] ;
  memset(tightenBounds_,0,numberColumns) ;
  for (int i=0;i<number;i++) {
    int k=which[i] ;
    if (k>=0&&k<numberColumns)
      tightenBounds_[k]=1 ;
  }
}


// Create C++ lines to get to current state
std::string
CglProbing::generateCpp( FILE * fp) 
{
  CglProbing other ;
  fprintf(fp,"0#include \"CglProbing.hpp\"\n") ;
  fprintf(fp,"3  CglProbing probing;\n") ;
  if (getMode()!=other.getMode())
    fprintf(fp,"3  probing.setMode(%d);\n",getMode()) ;
  else
    fprintf(fp,"4  probing.setMode(%d);\n",getMode()) ;
  if (getMaxPass()!=other.getMaxPass())
    fprintf(fp,"3  probing.setMaxPass(%d);\n",getMaxPass()) ;
  else
    fprintf(fp,"4  probing.setMaxPass(%d);\n",getMaxPass()) ;
  if (getLogLevel()!=other.getLogLevel())
    fprintf(fp,"3  probing.setLogLevel(%d);\n",getLogLevel()) ;
  else
    fprintf(fp,"4  probing.setLogLevel(%d);\n",getLogLevel()) ;
  if (getMaxProbe()!=other.getMaxProbe())
    fprintf(fp,"3  probing.setMaxProbe(%d);\n",getMaxProbe()) ;
  else
    fprintf(fp,"4  probing.setMaxProbe(%d);\n",getMaxProbe()) ;
  if (getMaxLook()!=other.getMaxLook())
    fprintf(fp,"3  probing.setMaxLook(%d);\n",getMaxLook()) ;
  else
    fprintf(fp,"4  probing.setMaxLook(%d);\n",getMaxLook()) ;
  if (getMaxElements()!=other.getMaxElements())
    fprintf(fp,"3  probing.setMaxElements(%d);\n",getMaxElements()) ;
  else
    fprintf(fp,"4  probing.setMaxElements(%d);\n",getMaxElements()) ;
  if (getMaxPassRoot()!=other.getMaxPassRoot())
    fprintf(fp,"3  probing.setMaxPassRoot(%d);\n",getMaxPassRoot()) ;
  else
    fprintf(fp,"4  probing.setMaxPassRoot(%d);\n",getMaxPassRoot()) ;
  if (getMaxProbeRoot()!=other.getMaxProbeRoot())
    fprintf(fp,"3  probing.setMaxProbeRoot(%d);\n",getMaxProbeRoot()) ;
  else
    fprintf(fp,"4  probing.setMaxProbeRoot(%d);\n",getMaxProbeRoot()) ;
  if (getMaxLookRoot()!=other.getMaxLookRoot())
    fprintf(fp,"3  probing.setMaxLookRoot(%d);\n",getMaxLookRoot()) ;
  else
    fprintf(fp,"4  probing.setMaxLookRoot(%d);\n",getMaxLookRoot()) ;
  if (getMaxElementsRoot()!=other.getMaxElementsRoot())
    fprintf(fp,"3  probing.setMaxElementsRoot(%d);\n",getMaxElementsRoot()) ;
  else
    fprintf(fp,"4  probing.setMaxElementsRoot(%d);\n",getMaxElementsRoot()) ;
  if (rowCuts()!=other.rowCuts())
    fprintf(fp,"3  probing.setRowCuts(%d);\n",rowCuts()) ;
  else
    fprintf(fp,"4  probing.setRowCuts(%d);\n",rowCuts()) ;
  if (getUsingObjective()!=other.getUsingObjective())
    fprintf(fp,"3  probing.setUsingObjective(%d);\n",getUsingObjective()) ;
  else
    fprintf(fp,"4  probing.setUsingObjective(%d);\n",getUsingObjective()) ;
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  probing.setAggressiveness(%d);\n",getAggressiveness()) ;
  else
    fprintf(fp,"4  probing.setAggressiveness(%d);\n",getAggressiveness()) ;
  return "probing" ;
}

