
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

namespace {

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
      if (up<e20Inf) {
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
      outside of CglProbing;
    * rowLower and rowUpper are allocated to the number of rows in the model
      plus 1 (in case we want the objective), so we need to be sure there's
      enough space;
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
// Generate disaggregation cuts
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
    errs = mtx->verifyMtx(2) ;
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
    errs = mtx->verifyMtx(2) ;
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
/*
  jjf: Setup information 

  Boil down detailed clique information into a form that's easily useable when
  calculating row bounds.
*/
  if (rowCliques && numberRows_ && numberColumns_)
    setupRowCliqueInformation(si) ;

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
/*
  If we don't have a snapshot, force mode 0 to mode 1.

  TODO: Ah ... maybe a `You're confused' message and an error return would be
	more appropriate? Seems like a serious algorithm error at a higher
	level.  -- lh, 100917 --

	Or maybe the client is within the documentation. After all, John's
	comment for mode 0 says `only available with snapshot; otherwise as
	mode 1'.  -- lh, 101021 --

	So, maybe I should try to figure out if requesting mode 0 results in
	creating a snapshot somewhere along the way.  -- lh, 110204 --
*/
  int mode = mode_ ;
  if (!rowCopy_ && mode == 0) {
    mode = 1 ;
#   if CGL_DEBUG > 0
    std::cout << "  forcing mode 0 to mode 1; no snapshot." << std::endl ;
#   endif
  }

  int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_ ;
  int maxRowLen = info->inTree ? maxElements_ : maxElementsRoot_ ;
  // Forcing probing to consider really dense constraints wasn't good,
  // even at the root.
  // if (!info->inTree && !info->pass) maxRowLen = nCols ;

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

  TODO: Down in probeClique, there's a clause that forces the cutoff to
  	infty for mode 0. Why not here?  -- lh, 110211 --
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
  Allocate a container to hold the cuts we'll generate locally. We'll only use
  this in mode 0.

  TODO: There's a similar calculation in probe() to size the local container
	used there. It allocates much more space. But note that the container
	passed to probe() is the OsiCuts container cs, given as a parameter to
	this method. rowCut is going to be used down in mode 0 for something.
*/
  int nRowsSafe = CoinMin(nRows,si.getNumRows()) ;
  int cutCapacity = info->inTree ? nRowsSafe/3 : nRowsSafe ;
  if (!info->inTree && !info->pass) cutCapacity *= 5 ;
  CglProbingRowCut rowCut(cutCapacity,!info->inTree) ;
/*
  Allocate the row processing record (markR) and the arrays that will hold
  row lhs bounds (minR, maxR).
*/
  int * markR = new int [nRows] ;
  double * minR = new double [nRows] ;
  double * maxR = new double [nRows] ;
/*
  End of Stage 1: Setup.

  Begin Stage 2: Mode-dependent Setup and Probing
*/
/*
  In modes 1 and 2, we're allowed to tighten bounds. Call tighten() to do bound
  propagation (two pass limit), then create `column cuts' to record the
  tightened bounds.

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
    if (!ninfeas) {   // block extends to end of if (mode)
/*
  Create column cuts where integer bounds have changed.
*/
      makeColCuts(nCols,cs,intVar,colsol,si.getColLower(),si.getColUpper(),
    		  si.getColSolution(),colLower,colUpper) ;
#     if CGL_DEBUG > 0
      CglProbingDebug::checkBounds(si,cs.getColCut(cs.sizeColCuts())) ;
#     endif
/*
  We've calculated changes in column bounds due to propagation. To go further,
  we need to be able to probe.
*/
      if (maxProbe > 0) {    // block extends to end of if (mode)
        numberThisTime_ = 0 ;
/*
  jjf: get min max etc for rows

  Calculate L<i> and U<i> for each row. If one or both are sufficiently small
  to process further, markR[i] is set to -1. If both are too large to warrant
  further processing, markR[i] is set to -2.

  TODO: Why do we need a second run for this? Why not do it as part of
	tighten()? My old notes indicate that tighten calculates the row
	bounds (and indeed it must to propagate column bounds). It's not
	like we're saving any work, other than transcribing the bounds into
	minR and maxR. On a purely mechanical level, tighten does not fill
	in minR, maxR, markR.  -- lh, 100923 --
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
  lookedAt_ now contains a list of indices of variables to probe. Rows
  worth processing are tagged with -1 in markR. Do the probing, with or
  without cliques.
*/
        if (!numberCliques_) {
          ninfeas = probe(si,cs,colLower,colUpper,rowCopy,columnCopy,
                          rowStartPos,realRows,rowLower,rowUpper,
                          intVar,minR,maxR,markR,info,
			  useObj,useCutoff,cutoff) ;
        } else {
          ninfeas = probeCliques(si,cs,colLower,colUpper,rowCopy,columnCopy,
                                 realRows,rowLower,rowUpper,
                                 intVar,minR,maxR,markR,info,
				 useObj,useCutoff,cutoff) ;
        }
      } // end maxProbe > 0
    } // end !nInfeas
  } // end mode != 0
/*
  End of block for mode = 1 or 2; begin mode = 0.
  
  Mode 0 works with the snapshot and does not do bound propagation. But note
  that much of the preparatory work that's done above for modes 1 and 2 in
  Stage 1 is done in snapshot() when the snapshot is created, including
  bound propagation. It really should be possible to unify a lot of this.

  TODO: Note that this code has been effectively disabled since 080428 (r629)
	when the only use (in CglPreProcess::modified) was guarded with 
	#if 0 / #endif. Expect oddities.  -- lh, 110212 --

  TODO: And now that I've worked through it, I have real questions. The key is
  	the table down at ACTIONS.
	
	We start out by calling probeClique and get back a private container
	of cuts. The code is clear that it's expecting only disaggregation
	cuts. We go through a lot of work to classify these cuts and process
	them into an alternate form which is stored in cutVector_. Then
	we work through cutVector_, looking for the situations identified
	in ACTIONS. Which are exactly the things that should be picked up
	by monotone actions, disaggregation cuts, and implication cuts. So
	if this code ever finds anything, it indicates an algorithm error
	someplace else. (Just possibly, we might pick up something across
	calls to this CglProbing object. cutVector_ is preserved for as
	long as the object is in existence.) Then we go through a lot
	more work to go back over cutVector_ one more time, checking for
	violations and reconstituting the cuts we processed away at the
	start. I'm still shaking my head.

	What's worth keeping? Well, the basic concept of the snapshot: We
	don't have to recreate the row and column copies with each call. And
	there's a simple placeholder where John's comment says `could augment
	cliques'. Which is true, but I wonder if that, too, isn't already done
	somewhere else, or could be done somewhere else more efficiently.

	So my conclusion, at this point, is that the only distinction worth
	keeping between mode 0 and modes 1 -- 3 is that mode 0 works from the
	snapshot. I'll see if I notice more along the way.

	-- lh, 110212 --

	Essentially confirmed in correspondence with JJF. Clique code
	imported from OSL but never fully functional.  -- lh, 110216 --
*/
  else
  if (maxProbe > 0) {
/*
  jjf: global cuts from previous calculations
       could check more thoroughly that integers are correct (sanity)
       make up list of new variables to look at 

  We're working from a snapshot.  In particular, we should have information
  about probing stashed in cutVector_ (for each x<p> = sequence, a vector of
  affected variables x<t> = index). So we're only going to probe variables
  where we have no information.

  When the snapshot is created, cutVector_ is initialised with all binary
  variables. So the snapshot documentation is a little deceptive. On first use
  after creation, we'll probe all binary variables, in order of distance from
  integrality, subject to whatever limits have been set.
*/
    assert(numberIntegers == numberIntegers_) ;
    numberThisTime_ = 0 ;
    const double *colsol = si.getColSolution() ;
    double_int_pair *array = new double_int_pair [number01Integers_] ;
#   ifdef ZEROFAULT
    std::memset(array,0,sizeof(double_int_pair)*number01Integers_) ;
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
    delete [] array ;
/*
  lookedAt_ now contains the indices of the binary variables for which we have
  no probing information. Calculate row lhs min and max, and call probeClique
  to do the probing (note that we will have clique information!)
*/
    calcRowBounds(colLower,colUpper,column,rowElements,rowStart,rowLength,
    		  rowLower,rowUpper,minR,maxR,markR,nRows) ;
/*
  jjf: don't do cuts at all if 0 (i.e. we are just checking bounds)

  TODO: Why don't we get a choice between probe and probeCliques? And why did
	we create a new cut collection (csNew) instead of using parameter cs?
	-- lh, 100924 --
*/
    OsiCuts csNew ;
    if (rowCuts_) {
      ninfeas = probeCliques(si,csNew,colLower,colUpper,
			     rowCopy,columnCopy,
			     realRows,rowLower,rowUpper,
			     intVar,minR,maxR,markR,info,
			     useObj,useCutoff,cutoff) ;
    }
/*
  TODO: I've made it to this point twice now (100924, 101021) and each time
  	bailed out to sort out other parts. Now it's time.

	The first thing that hits me is `If we have no snapshot to work from
	(rowCuts_ = 0), csNew is surely empty. Why are we even executing
	this next bit of code? For that matter, why are we here?! Mode 0
	should have been converted to mode 1.

	After reading ahead about 350 lines (look for END CLASSIFICATION) it's
	clear what we're doing in the next bit of code. We take the
	disaggregation cuts returned by probeCliques (and gods help us if
	there's anything else) and do a complicated case analysis to try and
	classify the type of cut we're looking at. We have about 150 lines of
	classification code whose sole purpose is to do a count, followed by
	the same 150 lines of classification code, this time filling in the
	entries and labelling them with the classification. This whole code
	block is highly suspect.

	-- lh, 110211 --

  AFTER PROBECLIQUE
*/
    if (!ninfeas) {
      int nCuts = csNew.sizeRowCuts() ;
      int iCut ;
/*
  backward will allow us to go from a probe index p to the corresponding in
  cutVector_.
  onList is an indicator that x<p> is in cutVector_

  Presumably onList will become more? It's not clear below why we need it.
*/
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
/*
  jjf: first do counts; we know initialized to zero

  We know this because we only probed variables that didn't have information.
  (Of course, if we guess wrong below while trying to determine the type of
  constraint, that goes out the window.)

  TODO: I've fixed up a few obvious errors (redundant tests, incorrect
	semantics) in the next two loops to classify cuts, but there are
	lots more. I'm going to wait 'til I have a grasp of the whole
	activity before trying to do more detail fixes. The possibility
	of classification error remains. Not to mention the `can do
	something?' comments.  -- lh, 110212 --

  TODO: There were a number of instances in the next two loops where there's a
  	test for backward[k] < 0 in one place and, close by, a test for
	onList[k]. Given initialisation immediately above, these tests should
	be equivalent. It'd be nice if the code made that apparent, and I've
	settled for backward[.] < 0.
	-- lh, 110212 --
*/
      for (iCut = 0 ; iCut < nCuts ; iCut++) {
	OsiRowCut rcut ;
	CoinPackedVector rpv ;
	rcut = csNew.rowCut(iCut) ;
	rpv = rcut.row() ;
/*
  Hmmm ... probeClique looks to be capable of generating coefficient
  strengthening cuts. The code is there, at least. But this assert is only
  guaranteed for disaggregation cuts.  -- lh, 110211 --

  Examination of the Cgl timeline says that mode 0 has been effectively
  disabled since 080428 r629, so this assert may well predate the introduction
  of coefficient strengthening.  -- lh, 110212 --
*/
	assert(rpv.getNumElements() == 2) ;
/*
  Get the index of the probing variable x<p> from the cut and use it to locate
  the entry in cutVector_. The other index will be x<t>.

  TODO: The probe variable x<p> could be either index, apparently. But
	it's entirely possible that x<p> in one cut could be x<t> in
	another. What happens then?  -- lh, 110211 --

  TODO: John allows that neither variable might be in cutVector_, because it's
  	possible for a variable to be converted to binary if its bounds are
	reduced to 0-1. So we have the freak accident of a disaggregation
	constraint where x<t> is general integer and x<p> was converted from
	general integer to binary.  -- lh, 110212 --
*/
	const int *indices = rpv.getIndices() ;
        int which = 0 ;
        int i = backward[indices[0]] ;
        if (i < 0) {
          which = 1 ;
          i = backward[indices[1]] ;
          if (i < 0) continue ;
        }
        int other = indices[1-which] ;
	double *elements = rpv.getElements() ;
	double lb = rcut.lb() ;
/*
  Sort out what kind of cut we're looking at. It's a bit hard to say
  whether this is correct for John's disaggregation cut generation code
  (I think not), and it's a cinch it's wrong for my disaggCuts method. But
  the details can be worked out.

  TODO: If this is going to be done easily and efficiently, it should be done
  	when the disaggregation cut is generated. The critical questions
	(identity of x<p>, probe direction, identity of x<t>, bound change)
	are all known at that point. It's far more difficult to recover that
	information here. Having said that, the current implementation always
	constructs disaggregation cuts in a canonical form, with x<p> in the
	second position.   -- lh, 110212 --

  In-code comments for this next bit are JJF.
*/
	if (lb == -DBL_MAX) {
          if (!rcut.ub()) {
            // UB
	    double elWhich = elements[which] ;
	    double elOther = elements[1-which] ;
            if (elWhich < 0.0) {
              // assert (elOther > 0.0) ;
              // delta to 0 => x to 0.0
              cutVector_[i].length++ ;
            } else {
              if (elOther < 0.0 &&
		  fabs(elWhich/elOther-colUpper[other]) < 1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++ ;
              } else {
                if (onList[other]) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1 == -1.0) {
                    // can do something ?
		    // (lh) is x<t> binary? if not, j < 0 here.
                    int j = backward[other] ;
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
              if (elements[0] == 1.0 &&
	          elements[1] == 1.0 && rcut.ub() == 1.0) {
                // can do something ?
		// (lh) is x<t> binary? if not, j < 0 here.
                int j = backward[other] ;
                cutVector_[i].length++ ;
                cutVector_[j].length++ ;
              } else {
                continue ;
              }
            }
          }
	} else {
/*
  End (stuff) <= b, begin b <= (stuff). Note that this case is not symmetric
  with the above case --- there's nothing for lb != 0.
*/
          assert(rcut.ub() == DBL_MAX) ;
          if (!lb) {
            // LB
	    double elWhich = elements[which] ;
	    double elOther = elements[1-which] ;
            if (elWhich > 0.0) {
              // assert (elOther < 0.0) ;
              // delta to 0 => x to 0.0
              // flip so same as UB
              cutVector_[i].length++; 
            } else {
              if (elOther < 0.0 &&
	      	  fabs(elWhich/elOther-colUpper[other]) < 1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++ ;
              } else {
                if (onList[other]) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1 == -1.0) {
                    // can do something ?
                    int j = backward[other] ;
                    cutVector_[i].length++ ;
                    cutVector_[j].length++ ;
                  } else {
                    continue ;
                  }
                }
              }
            }
          }        // end of b == 0
	}       // end of b <= (stuff)
      }      // end of counting loop
/*
  We have a set of counts attached to entries in cutVector_ that tell us how
  long the index arrays will need to be. Allocate the arrays.
*/
      for (int i = 0 ; i < number01Integers_ ; i++) {
	int j = cutVector_[i].sequence ;
	if (onList[j] && !cutVector_[i].index) {
	  disaggregation thisOne = cutVector_[i] ;
	  cutVector_[i].index = new disaggregationAction [thisOne.length] ;
          cutVector_[i].length = 0 ;
	}
      }
/*
  And load up the arrays.

  The classification structure here is cut-and-paste from above; all the
  comments about the case analysis apply again. The Agree/Disagree forms
  assume binary probe and target.   -- lh, 110212 --

  TODO: Note that the index stored in disaggregation.index[] can be either an
	index for cutVector_ (zeroOneInDisagg true) or the original
	column index (false). Presence in cutVector_ implies binary; it's
	still unclear to me if the reverse is always true. cutVector_
	may not be loaded with every binary variable.  -- lh, 110212 --

  TODO: What you can see, after a bit of thought, is that for a disaggregation
  	cut on two binary variables, there's a symmetry. For example,
	x<p> -> 1 ==> x<t> -> 0 gives x<p> <= (1-x<t>), which also arises from
	x<t> -> 1 ==> x<p> -> 0 (with the roles of x<p> and x<t> exchanged).
	It looks like the classification scheme is trying to exploit this, but
	not getting the decision tree right. We need to check for two binary
	variables right at the top.  -- lh, 110212 --

*/
      for (iCut = 0 ; iCut < nCuts ; iCut++) {
	OsiRowCut rcut ;
	CoinPackedVector rpv ;
	int iput ;
	rcut = csNew.rowCut(iCut) ;
	rpv = rcut.row() ;
	assert(rpv.getNumElements() == 2) ;
	const int *indices = rpv.getIndices() ;
	double *elements = rpv.getElements() ;
	double lb = rcut.lb() ;
	// (lh) assume x<p> is the first binary variable we find?
        int which = 0 ;
        int i = backward[indices[0]] ;
        if (i < 0) {
          which = 1 ;
          i = backward[indices[1]] ;
          if (i < 0) continue ;
        }
	disaggregation *cVi = &cutVector_[i] ;
        int other = indices[1-which] ;
        int j = backward[other] ;
	disaggregation *cVj = 0 ;
	if (j >= 0) cVj = &cutVector_[j] ;

	if (lb == -DBL_MAX) {	// <= constraint
          if (!rcut.ub()) {     // b = 0
	    double elWhich = elements[which] ;
	    double elOther = elements[1-which] ;
            if (elWhich < 0.0) {
	      // delta to 0 => x to 0.0
	      // Agree: -x<p>+x<t> <= 0
              iput = cVi->length ;
              if (j >= 0)
                setAffectedInDisagg(cVi->index[iput],j) ;
              else
                setAffectedInDisagg(cVi->index[iput],other) ;
              setWhenAtUBInDisagg(cVi->index[iput],false) ;
              setAffectedToUBInDisagg(cVi->index[iput],false) ;
	      setZeroOneInDisagg(cVi->index[iput],(j >= 0)) ;
              cVi->length++ ;
            } else { 
	      /*
	        Wrong classification check. Assume x<p> binary, x<t> general
		integer, l<p> = l<t> = 0. l'<p> = 1. Canonical form is
		(l<t>-l'<t>)x<p> + (l'<p>-l<p>)x<t> >= l<t>l'<p>-l'<t>l<p>
		reduces to
		l'<t>x<p> - x<t> <= 0
		l'<t> is tightened colLower, equivalent to colUpper only for
		binary variables or fixed general integer (shades of
		goingToTrueBound). Also, we've carefully tested to ensure
		that a<p> > 0 and a<t> < 0.  -- lh, 110212 --
	      */
              if (elOther < 0.0 &&
	          fabs(elWhich/elOther-colUpper[other]) < 1.0e-5) {
                // delta to 1 => x to upper bound
		// Agree: x<p>-x<t> <= 0
                iput = cVi->length ;
                if (j >= 0)
                  setAffectedInDisagg(cVi->index[iput],j) ;
                else
                  setAffectedInDisagg(cVi->index[iput],other) ;
                setWhenAtUBInDisagg(cVi->index[iput],true) ;
                setAffectedToUBInDisagg(cVi->index[iput],true) ;
                setZeroOneInDisagg(cVi->index[iput],(j >= 0)) ;
                cVi->length++ ;
              } else {
	        /*
		  a<p> > 0, a<t> < 0, b == 0, x<p>  binary. Now add x<t>
		  binary (j >= 0) and a<p>*a<t> = -1.0.  Effectively, we're
		  looking for x<p> - x<t> <= 0. Sure, we can do something. We
		  just did it, immediately above. Symmetry here says that
		  x<p> -> 1 ==> x<t> -> 1 and x<t> -> 0 ==> x<p> -> 0
		  -- lh, 110212 --
		*/
                if (j >= 0) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1 == -1.0) {
                    // can do something ?
                    // flip so value0 1.0
                    if (value1 == 1.0) {
                      j = i ;
                      i = backward[other] ;
                      value1 = value0 ;
                      value0 = 1.0 ;
                    }
                    assert (value0 == 1.0) ;
                    assert (value1 == -1.0) ;
                    iput = cVi->length ;
                    setAffectedInDisagg(cVi->index[iput],j) ;
                    setWhenAtUBInDisagg(cVi->index[iput],true) ;
                    setAffectedToUBInDisagg(cVi->index[iput],true) ;
                    setZeroOneInDisagg(cVi->index[iput],true) ;
                    cVi->length++ ;
                    iput = cVj->length ;
                    setAffectedInDisagg(cVj->index[iput],i) ;
                    setWhenAtUBInDisagg(cVj->index[iput],false) ;
                    setAffectedToUBInDisagg(cVj->index[iput],false) ;
                    setZeroOneInDisagg(cVj->index[iput],true) ;
                    cVj->length++ ;
                  }
                }
              }
            }
          } else {
	    /*
	      stuff <= b, b != 0. Additionally require x<t> binary,
	      a<p> = a<t> = 1.0, b = 1. So we're looking for x<p> + x<t> <= 1.
	      This is the binary disaggregation for x<p> to 1 ==> x<t> to 0,
	      and x<t> -> 1 ==> x<p> -> 0.
	    */
            if (onList[other]) {
              if (elements[0] == 1.0 &&
	          elements[1] == 1.0 && rcut.ub() == 1.0) {
                // can do something ?
                int j = backward[other] ;
                assert ( j >= 0) ;
                iput = cVi->length ;
                setAffectedInDisagg(cVi->index[iput],j) ;
                setWhenAtUBInDisagg(cVi->index[iput],true) ;
                setAffectedToUBInDisagg(cVi->index[iput],false) ;
                setZeroOneInDisagg(cVi->index[iput],true) ;
                cVi->length++ ;
                iput = cVj->length ;
                setAffectedInDisagg(cVj->index[iput],i) ;
                setWhenAtUBInDisagg(cVj->index[iput],true) ;
                setAffectedToUBInDisagg(cVj->index[iput],false) ;
                setZeroOneInDisagg(cVj->index[iput],true) ;
                cVj->length++ ;
              } else {
#ifdef COIN_DEVELOP
                abort() ;
#endif
		continue ;
              }
            }
          }
	} else {
/*
  More of the same, this time for constraints of the form b <= (stuff).
*/
          assert(rcut.ub() == DBL_MAX) ;
          if (!lb) {
            // LB
	    double elWhich = elements[which] ;
	    double elOther = elements[1-which] ;
            if (elWhich > 0.0) {
              iput = cVi->length ;
              if (j >= 0)
                setAffectedInDisagg(cVi->index[iput],j) ;
              else
                setAffectedInDisagg(cVi->index[iput],other) ;
              setWhenAtUBInDisagg(cVi->index[iput],false) ;
              setAffectedToUBInDisagg(cVi->index[iput],false) ;
              setZeroOneInDisagg(cVi->index[iput],onList[other]!=0) ;
              cVi->length++ ;
            } else { 
              if (elOther < 0.0 &&
	          fabs(elWhich/elOther-colUpper[other]) < 1.0e-5) {
                iput = cVi->length ;
                if (j >= 0)
                  setAffectedInDisagg(cVi->index[iput],j) ;
                else
                  setAffectedInDisagg(cVi->index[iput],other) ;
                setWhenAtUBInDisagg(cVi->index[iput],true) ;
                setAffectedToUBInDisagg(cVi->index[iput],true) ;
                setZeroOneInDisagg(cVi->index[iput],onList[other]!=0) ;
                cVi->length++ ;
              } else {
                if (onList[other]) {
                  double value0 = elements[0] ;
                  double value1 = elements[1] ;
                  if (value0*value1 == -1.0) {
                    // can do something ?
                    int j = backward[other] ;
                    assert (j>=0) ;
                    // flip so value0 -1.0
                    if (value1 == -1.0) {
                      j = i ;
                      i = backward[other] ;
                      value1 = value0 ;
                      value0 = -1.0 ;
                    }
                    assert (value0 == -1.0) ;
                    assert (value1 == 1.0) ;
                    iput = cVi->length ;
                    setAffectedInDisagg(cVi->index[iput],j) ;
                    setWhenAtUBInDisagg(cVi->index[iput],true) ;
                    setAffectedToUBInDisagg(cVi->index[iput],true) ;
                    setZeroOneInDisagg(cVi->index[iput],true) ;
                    cVi->length++ ;
                    iput = cVj->length ;
                    setAffectedInDisagg(cVj->index[iput],i) ;
                    setWhenAtUBInDisagg(cVj->index[iput],false) ;
                    setAffectedToUBInDisagg(cVj->index[iput],false) ;
                    setZeroOneInDisagg(cVj->index[iput],true) ;
                    cVj->length++ ;
                  }
                }
              }
            }
          }
	} 
      }
      delete [] backward ;
/*
  Note that onList is just a pointer into the block of space headed by
  backward.

  END CLASSIFICATION

  We know have a vector of indices attached to each entry in cutVector_,
  with bits set to indicate just what sort of disaggregation cut we're
  dealing with.

  jjf: Now sort and get rid of duplicates;  could also see if any are cliques

  To do the sort, compose a unique signature for each disaggregation entry, as
  (index)@(status bits), and sort the resulting vector. Then compare adjacent
  entries, squeezing out duplicates.
*/
      int longest = 0 ;
      for (int i = 0 ; i < number01Integers_ ; i++) 
        longest = CoinMax(longest,cutVector_[i].length) ;
      unsigned int *sortit = new unsigned int[longest] ;

      for (int i = 0 ; i < number01Integers_ ; i++) {

        disaggregation &thisOne = cutVector_[i] ;
        int k ;
        int number = thisOne.length ;
        for (k = 0 ; k < number ; k++) {
          int affected = affectedInDisagg(thisOne.index[k]) ;
          int zeroOne = zeroOneInDisagg(thisOne.index[k]) ? 1 : 0 ;
          int whenAtUB = whenAtUBInDisagg(thisOne.index[k]) ? 1 : 0 ;
          int affectedToUB = affectedToUBInDisagg(thisOne.index[k]) ? 1 : 0 ;
          sortit[k] = (affected<<3)|(zeroOne<<2)|(whenAtUB<<1)|affectedToUB ;
        }
        std::sort(sortit,sortit+number) ;

        int affectedLast = 0xffffffff ;
        int zeroOneLast = 0 ;
        int whenAtUBLast = 0 ;
        int affectedToUBLast = 0 ;
        int put = 0 ;
        for (k = 0 ; k < number ; k++) {
          int affected = sortit[k]>>3 ;
          int zeroOne = (sortit[k]&4)>>2 ;
          int whenAtUB = (sortit[k]&2)>>1 ;
          int affectedToUB = sortit[k]&1 ;
          disaggregationAction action ;
	  action.affected = 0 ;
          setAffectedInDisagg(action,affected) ;
          setZeroOneInDisagg(action,(zeroOne != 0)) ;
          setWhenAtUBInDisagg(action,(whenAtUB != 0)) ;
          setAffectedToUBInDisagg(action,(affectedToUB != 0)) ;
          if (affected != affectedLast || zeroOne != zeroOneLast) {
            // new variable
            thisOne.index[put++] = action ;
          } else if (whenAtUB != whenAtUBLast ||
	  	     affectedToUB != affectedToUBLast) {
            // new action - what can we discover
            thisOne.index[put++] = action ;
            int j = cutVector_[i].sequence ;
            int k = affected ;
            if (zeroOne) {
              k = cutVector_[k].sequence ;
              if (logLevel_ > 1)
	        std::cout
		  << "  for 0-1 pair " << j << " and "
		  << k << ":" << std::endl ;
            } else {
              if (logLevel_ > 1)
	        std::cout
		  << "  for pair " << j << " and " << k << ":" << std::endl ;
            }
            if (logLevel_ > 1)
	      std::cout
	        << "    previous: " << ((whenAtUBLast)?"up":"down") << " probe "
		<< "tightened " << ((affectedToUBLast)?"lower":"upper")
		<< " bound" << std::endl ;
	      std::cout
	        << "    current: " << ((whenAtUB)?"up":"down") << " probe "
		<< "tightened " << ((affectedToUB)?"lower":"upper")
		<< " bound" << std::endl ;
          }
          affectedLast = affected ;
          zeroOneLast = zeroOne ;
          whenAtUBLast = whenAtUB ;
          affectedToUBLast = affectedToUB ;
        }
        if (put < number) {
	  if (logLevel_ > 1)
	    std::cout
	      << "  actions for " << i << "reduced from " << number
	      << " to " << put << "." << std::endl ;
          thisOne.length = put ;
        }
      }
/*
  End sort. Now look over the lists of disaggregation cuts and see what we can
  learn from the ones where both x<p> and x<t> are binary.

  TODO: Unless I totally miss my guess here, we're about to derive implication
	cuts. And indeed, if I'd read ahead 20 lines, that's exactly what
	we're up to. Not to mention looking again for monotone actions and
	disaggregation cuts. So the real question is `Why aren't we doing all
	this by calling the appropriate methods in probeClique?

	What could be done here is to try and extend cliques. But there's only
	a placeholder case statement. Which is consistent with the comment up
	in CglPreProcess at the call to snapshot: `out for now - think about
	cliques'.
	
	-- lh, 110212 --
*/
      for (int i = 0 ; i < number01Integers_ ; i++) {
        disaggregation &thisOne = cutVector_[i] ;
        int number = thisOne.length ;
        for (int k = 0 ; k < number ; k++) {
          int affected = affectedInDisagg(thisOne.index[k]) ;
          bool zeroOne = zeroOneInDisagg(thisOne.index[k]) ;
	  /*
	    affected < i means we've already processed the equivalent cut
	    when we processed the list associated with affected, with the
	    roles of affected and i exchanged.
	  */
          if (zeroOne && affected > i) {
            bool whenAtUB = whenAtUBInDisagg(thisOne.index[k]) ;
            bool affectedToUB = affectedToUBInDisagg(thisOne.index[k]) ;
            disaggregation otherOne = cutVector_[affected] ;
            int numberOther = otherOne.length ;
            // Could do binary search if a lot
            int lastAction = -1 ;
            for (int j = 0 ; j < numberOther ; j++) {
              if (affectedInDisagg(otherOne.index[j])==i) {
                bool whenAtUBOther = whenAtUBInDisagg(otherOne.index[j]) ;
                bool affectedToUBOther = affectedToUBInDisagg(otherOne.index[j]) ;
                /* ACTIONS
		
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

		  (lh) 0 -- 3 are standard disaggregation forms. 10 -- 11 are
		  implication cuts. 20 -- 23 are caught as monotone actions or
		  implication cuts. We're just treading old ground. If this
		  code discovers anything, it indicates an algorithm error
		  elsewhere.  -- lh, 110212 --
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
    }		// end AFTER PROBECLIQUE
/*
  End of the block that starts with if (!ninfeas), right after the call to
  probeClique.

  Now have a look at the cuts in cutVector_ and see if any are violated. Only
  look at the entries for unfixed variables, and then only if the current
  solution value is not integral.

  TODO: The immediate question is ``Why are we doing this if we're
	infeasible?' I suppose you can argue that cutVector_ will be mostly
	empty (we haven't loaded in the disaggregation cuts), but it will be
	present. And didn't we check these cuts for violation when we created
	them? Apparently not -- a quick look ahead says we're finally going to
	add the cut to cs if we find it's violated.   -- lh, 110212 -- 

  TODO: Good grief! We're about to reconstitute the cuts that we processed
	down into the entries in cutVector_. And we're not even going to do
	a complete job of it. Surely there must be a better way.
	-- lh, 110212 --
*/
    if (cutVector_) {
      for (int i = 0 ; i < number01Integers_ ; i++) {
	int j = cutVector_[i].sequence ;
	double solInt = colsol[j] ;
	double  upper ;
	double solValue ;
	int icol = -1 ;
	int index[2] ;
	double element[2] ;
	if (colUpper[j]-colLower[j] > 1.0e-8) {
	  double away = fabs(0.5-(solInt-floor(solInt))) ;
	  if (away < 0.4999999) {
	    disaggregation thisOne = cutVector_[i] ;
	    int k ;
	    OsiRowCut rc ;
	    for (k = 0 ; k < thisOne.length ; k++) {
	      icol = affectedInDisagg(thisOne.index[k]) ;
              if (zeroOneInDisagg(thisOne.index[k]))
                icol = cutVector_[icol].sequence ;
	      solValue = colsol[icol] ;
	      upper = colUpper_[icol] ;
              double infeasibility = 0.0 ;
              if (!whenAtUBInDisagg(thisOne.index[k])) {
                if (!affectedToUBInDisagg(thisOne.index[k])) {
                  // delta -> 0 => x to lb (at present just 0)
                  infeasibility = solValue - upper * solInt ;
                  if (infeasibility > 1.0e-3) {
                    rc.setLb(-DBL_MAX) ;
                    rc.setUb(0.0) ;
                    index[0] = icol ;
                    element[0] = 1.0 ;
                    index[1] = j ;
                    element[1] = -upper ;
                  } else {
                    infeasibility = 0.0 ;
                  }
                } else {
                  // delta -> 0 => x to ub
                  abort() ;
                }
              } else {
                if (affectedToUBInDisagg(thisOne.index[k])) {
                  // delta -> 1 => x to ub (?)
                  icol = affectedInDisagg(thisOne.index[k]) ;
                  if (zeroOneInDisagg(thisOne.index[k]))
                    icol = cutVector_[icol].sequence ;
                  solValue = colsol[icol] ;
                  upper = colUpper_[icol] ;
                  if (!colLower[icol]) {
                    infeasibility = upper * solInt - solValue ;
                    if (infeasibility > 1.0e-3) {
                      rc.setLb(-DBL_MAX) ;
                      rc.setUb(0.0) ;
                      index[0] = icol ;
                      element[0] = -1.0 ;
                      index[1] = j ;
                      element[1] = upper ;
                    } else {
                      infeasibility = 0.0 ;
                    }
                  } else {
                    assert (upper == colLower[icol]) ;
                    infeasibility = 0.0 ;
                  }
                } else {
                  // delta + delta2 <= 1
                  assert (zeroOneInDisagg(thisOne.index[k])) ;
                  // delta -> 1 => delta2 -> 0
                  icol = affectedInDisagg(thisOne.index[k]) ;
                  icol = cutVector_[icol].sequence ;
                  // only do if icol > j
                  if (icol > j && colUpper[icol] ) {
                    solValue = colsol[icol] ;
                    if (!colLower[icol]) {
                      infeasibility = solInt + solValue - 1.0 ;
                      if (infeasibility > 1.0e-3) {
                        rc.setLb(-DBL_MAX) ;
                        rc.setUb(1.0) ;
                        index[0] = icol ;
                        element[0] = 1.0 ;
                        index[1] = j ;
                        element[1] = 1.0 ;
                      } else {
                        infeasibility = 0.0 ;
                      }
                    } else {
                      assert (upper == colLower[icol]) ;
                      infeasibility = 0.0 ;
                    }
                  }
                }
              }
              if (infeasibility) {
                rc.setEffectiveness(infeasibility) ;
                rc.setRow(2,index,element,false) ;
                if (logLevel_>1)
                  printf("%g <= %g * x%d + %g * x%d <= %g\n",
                         rc.lb(),element[0],index[0],element[1],
			 index[1],rc.ub()) ;
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
    << "Leaving CglProbing::gutsOfGenerateCuts, "
    << cs.sizeRowCuts() << " row cuts, " << cs.sizeColCuts()
    << " col cuts, infeas " << ninfeas << std::endl ;
  debugger = si.getRowCutDebugger() ;
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

/*
  Creates cliques for use by probing.
  Can also try and extend cliques as a result of probing (root node).
  Returns number of cliques found.

  minimumSize and maximumSize limit the size of the cliques reported.

  Mathematically, we look for constraints of the form
      L<i> <= SUM{P}x<j> - SUM{M}x<j> <= U<i>
  where all x<j> are binary, all x<j> j in P have a coefficient of +1, and
  all x<j>, j in M, have a coefficient of -1. The smallest value of the lhs
  is -|M|, the largest, |P|. We can make the following conclusions (numbered to
  match state codes below):

  3) -|M| > U<i> or |P| < L<i>  ==> infeasible

  2) -|M| = U<i> ==> all x<j>, j in M, at 1, all x<j>, j in P, at 0
      |P| = L<i> ==> all x<j>, j in P, at 1, all x<j>, j in M, at 0

  1) U<i> = -|M|+1 ==> any x<t>, t in P, at 1 ==> all x<j>, j in M, at 1
					      ==> all x<j>, j in P\t, at 0
		   ==> any x<t>, t in M, at 0 ==> all x<j>, j in P, at 0
					      ==> all x<j>, j in M\t, at 1
     L<i> =  |P|-1 ==> any x<t>, j in M, at 1 ==> all x<j>, j in P, at 1
					      ==> all x<j>, j in M\t, at 0
		   ==> any x<t>, j in P, at 0 ==> all x<j>, j in M, at 0
					      ==> all x<j>, j in P\t, at 1

  1) is the only case where we have an interesting clique to work with.

*/
int 
CglProbing::createCliques (OsiSolverInterface &si, 
			   int minimumSize, int maximumSize)
{
  // Remove existing clique information
  deleteCliques() ;
/*
  Generate cliques from the original matrix. Don't use a snapshot, even if
  we have one.
*/
  CoinPackedMatrix matrixByRow(*si.getMatrixByRow()) ;
  int numberRows = si.getNumRows() ;
  if (!rowCopy_)
    numberRows_ = numberRows ;
  numberColumns_ = si.getNumCols() ;
/*
  Initialise some bookkeeping. lookup will enable us to find the integer
  variable's sequence number given the column number.
*/
  numberCliques_ = 0 ;
  int numberEntries = 0 ;
  int numberIntegers = 0 ;
  int *lookup = new int[numberColumns_] ;
  int i ;
  for (i = 0 ; i < numberColumns_ ; i++) {
    if (si.isBinary(i))
      lookup[i] = numberIntegers++ ;
    else
      lookup[i] = -1 ;
  }
  int *which = new int[numberColumns_] ;
  int *whichRow = new int[numberRows] ;

  // Statistics
  int totalP1 = 0 ;
  int totalM1 = 0 ;
  int numberBig = 0 ;
  int totalBig = 0 ;
  int numberFixed = 0 ;

  // Expose the internal arrays of the matrix
  const double *elementByRow = matrixByRow.getElements() ;
  const int *column = matrixByRow.getIndices() ;
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts() ;
  const int *rowLength = matrixByRow.getVectorLengths() ;

  /*
    jjf: Column lengths for slacks
    
    See note below re. slacks. Vestigial  -- lh, 100924
  */
  const int *columnLength = si.getMatrixByCol()->getVectorLengths() ;

  const double *lower = si.getColLower() ;
  const double *upper = si.getColUpper() ;
  const double *rowLower = si.getRowLower() ;
  const double *rowUpper = si.getRowUpper() ;
/*
  Get down to the business of scanning rows.
*/
  int iRow ;
  for (iRow = 0 ; iRow < numberRows ; iRow++) {
    int numberP1 = 0 ;
    int numberM1 = 0 ;
    int j ;
    double upperValue = rowUpper[iRow] ;
    double lowerValue = rowLower[iRow] ;
    bool good = true ;
    int slack = -1 ;
/*
  Scan the columns of the current row. In order to generate cliques, the row
  must contain only binary variables with coefficients of +/- 1.0. Adjust the
  row bounds for fixed variables while we're here. If the row survives, the
  small end of which contains the indices of the positive coefficients, the
  large end the indices of the negative coefficients.

  TODO: We set slack when we run across a column of length 1, but don't use it
	for anything.   -- lh, 100924 --
*/
    for (j = rowStart[iRow] ; j < rowStart[iRow]+rowLength[iRow] ; j++) {
      int iColumn = column[j] ;
      int iInteger = lookup[iColumn] ;
      if (upper[iColumn]-lower[iColumn] < 1.0e-8) {
	upperValue -= lower[iColumn]*elementByRow[j] ;
	lowerValue -= lower[iColumn]*elementByRow[j] ;
	continue ;
      } else if (upper[iColumn] != 1.0 || lower[iColumn] != 0.0) {
	good = false ;
	break ;
      } else if (iInteger < 0) {
	good = false ;
	break ;
      } else {
	if (columnLength[iColumn] == 1)
	  slack = iInteger ;
      }
      if (fabs(elementByRow[j]) != 1.0) {
	good = false ;
	break ;
      } else if (elementByRow[j] > 0.0) {
	which[numberP1++] = iColumn ;
      } else {
	numberM1++ ;
	which[numberIntegers-numberM1] = iColumn ;
      }
    }
/*
  See if we have any strength. As explained at the head of the routine,
  3: infeasible; 2: all variables fixed ; 1: clique of some sort.
*/
    int iUpper = static_cast<int>(floor(upperValue+1.0e-5)) ;
    int iLower = static_cast<int>(ceil(lowerValue-1.0e-5)) ;
    int state = 0 ;
    if (upperValue < 1.0e6) {
      if (iUpper == 1-numberM1)
	state = 1 ;
      else if (iUpper == -numberM1)
	state = 2 ;
      else if (iUpper < -numberM1)
	state = 3 ;
    }
    if (!state && lowerValue > -1.0e6) {
      if (-iLower == 1-numberP1)
	state = -1 ;
      else if (-iLower == -numberP1)
	state = -2 ;
      else if (-iLower < -numberP1)
	state = -3 ;
    }
/*
  +/- 3 is infeasible; we're done.
*/
    if (good && state) {
      if (abs(state) == 3) {
	numberCliques_ = -99999 ;
	break ;
      }
/*
  +2 means all x<j>, j in M, must be 1, all x<j>, j in P, must be 0.
  -2 is the opposite.
*/
      else if (abs(state) == 2) {
	numberFixed += numberP1+numberM1 ;
	if (state > 0) {
	  for (i = 0 ; i < numberP1 ; i++)
	    si.setColUpper(which[i],0.0) ;
	  for (i = 0 ; i < numberM1 ; i++)
	    si.setColLower(which[numberIntegers-i-1],1.0) ;
	} else {
	  for (i = 0 ; i < numberP1 ; i++)
	    si.setColLower(which[i],1.0) ;
	  for (i = 0 ; i < numberM1 ; i++)
	    si.setColUpper(which[numberIntegers-i-1],0.0) ;
	}
      }
/*
  +1 means that any x<j>, j in P, at 1 implies all x<j>, j in M, at 1.
  -1 means that any x<j>, j in M, at 1 implies all x<j>, j in P, at 1.
*/
      else {
	int length = numberP1+numberM1 ;
        totalP1 += numberP1 ;
        totalM1 += numberM1 ;
	if (length >= minimumSize && length < maximumSize) {
	  whichRow[numberCliques_++] = iRow ;
	  numberEntries += length ;
	} else if (numberP1+numberM1 >= maximumSize) {
	  numberBig++ ;
	  totalBig += numberP1+numberM1 ;
	}
      }
    }
  }
/*
  Done with the row scan. Print a bit of information.
*/
  if (numberCliques_ < 0) {
    if (logLevel_)
      printf("*** Problem infeasible\n") ;
  } else {
    if (numberCliques_) {
      if (logLevel_)
        printf("%d cliques of average size %g found, %d P1, %d M1\n",
               numberCliques_,
               (static_cast<double>(totalP1+totalM1))/
	       (static_cast<double> (numberCliques_)),
               totalP1,totalM1) ;
    } else {
      if (logLevel_ > 1)
        printf("No cliques found\n") ;
    }
    if (numberBig) {
      if (logLevel_)
        printf("%d large cliques ( >= %d) found, total %d\n",
	     numberBig,maximumSize,totalBig) ;
    }
    if (numberFixed) {
      if (logLevel_)
        printf("%d variables fixed\n",numberFixed) ;
    }
  }
/*
  Deal with any cliques we found. The only information saved from the row scan
  is the row indices of the cliques. Start by setting up the bookkeeping arrays.

    cliqueType_  is a bit array; if set, the clique row is an equality
    cliqueStart_ is an index array; each entry points to the start of a clique
		 in cliqueEntry_
    cliqueEntry_ is an array with one entry per clique member; each entry
		 gives the column number of the variable and a single bit
		 indicating whether 1 or 0 is the strong direction.
    whichClique_ crossreferences a variable back to its clique; for each
		 variable there's a block of clique indices where the variable
		 is strong-1, followed by a block where the variable is
		 strong-0.
    oneFixStart_ points to the start of the strong-1 block for the variable
    zeroFixStart_ points to the start of the strong-0 block for the variable
    endFixStart_ is one past the end of the entire block for the variable
*/
  if (numberCliques_ > 0) {
    cliqueType_ = new cliqueType [numberCliques_] ;
    cliqueStart_ = new int [numberCliques_+1] ;
    cliqueEntry_ = new cliqueEntry [numberEntries] ;
    oneFixStart_ = new int [numberColumns_] ;
    zeroFixStart_ = new int [numberColumns_] ;
    endFixStart_ = new int [numberColumns_] ;
    whichClique_ = new int [numberEntries] ;
    numberEntries = 0 ;
    cliqueStart_[0] = 0 ;
    for (i = 0 ; i < numberColumns_ ; i++) {
      oneFixStart_[i] = -1 ;
      zeroFixStart_[i] = -1 ;
      endFixStart_[i] = -1 ;
    }
/*
  Process each clique into cliqueEntry storage format. Keep in mind that some
  cliques may now be vacant because their variables have been fixed, so the
  first thing we need to do is evaluate each row again to see if the clique
  conditions still hold, and to recover the index sets P and M. At least we
  know the row has the correct form.

  TODO: Is there any possibility that a type 1 row has converted to type 2? Is
	it worth checking here? -- lh, 101007 --
*/
    int iClique ;
    int numberCliques = numberCliques_ ;
    numberCliques_ = 0 ;
    for (iClique = 0 ; iClique < numberCliques ; iClique++) {
      int iRow = whichRow[iClique] ;
      whichRow[numberCliques_] = iRow ;
      int numberP1 = 0 ;
      int numberM1 = 0 ;
      int j ;
      double upperValue = rowUpper[iRow] ;
      double lowerValue = rowLower[iRow] ;
      for (j = rowStart[iRow] ; j < rowStart[iRow]+rowLength[iRow] ; j++) {
	int iColumn = column[j] ;
	if (upper[iColumn]-lower[iColumn] < 1.0e-8) {
	  upperValue -= lower[iColumn]*elementByRow[j] ;
	  lowerValue -= lower[iColumn]*elementByRow[j] ;
	  continue ;
	}
	if (elementByRow[j] > 0.0) {
	  which[numberP1++] = iColumn ;
	} else {
	  numberM1++ ;
	  which[numberIntegers-numberM1] = iColumn ;
	}
      }
      int iUpper = static_cast<int>(floor(upperValue+1.0e-5)) ;
      int iLower = static_cast<int>(ceil(lowerValue-1.0e-5)) ;
      int state = 0 ;
      if (upperValue < 1.0e6) {
	if (iUpper == 1-numberM1)
	  state = 1 ;
      }
      if (!state && lowerValue > -1.0e6) {
	state = -1 ;
      }
/*
  If the clique conditions no longer hold, move on to the next candidate.
*/
      if (abs(state) != 1)
	continue ;
/*
  Mark whether the clique row is an equality in cliqueType_. Record the
  column index of each clique member and the strong branching direction in
  cliqueEntry_.

  TODO: Unless there's some indication that other information will go into
	cliqueType, converting it to a simple boolean array, or even a struct
	with a simple boolean array, would be better than this single-bit
	setup.  -- lh, 101007 --
*/
      if (iLower == iUpper) {
	cliqueType_[numberCliques_].equality = 1 ;
      } else {
	cliqueType_[numberCliques_].equality = 0 ;
      }
/*
  U<i> = -|M|+1. 1 is the strong value for x<j> in P, 0 is the strong value for
  x<j> in M.
*/
      if (state > 0) {
	for (i = 0 ; i < numberP1 ; i++) {
	  int iColumn = which[i] ;
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn) ;
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],true) ;
	  numberEntries++ ;
	  oneFixStart_[iColumn] = 0 ;
	  zeroFixStart_[iColumn] = 0 ;
	}
	for (i = 0 ; i < numberM1 ; i++) {
	  int iColumn = which[numberIntegers-i-1] ;
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn) ;
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],false) ;
	  numberEntries++ ;
	  oneFixStart_[iColumn] = 0 ;
	  zeroFixStart_[iColumn] = 0 ;
	}
      }
/*
  L<i> = |P|-1. 0 is the strong value for x<j> in P, 1 is the strong value for
  x<j> in M.
*/
      else {
	for (i = 0 ; i < numberP1 ; i++) {
	  int iColumn = which[i] ;
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn) ;
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],false) ;
	  numberEntries++ ;
	  oneFixStart_[iColumn] = 0 ;
	  zeroFixStart_[iColumn] = 0 ;
	}
	for (i = 0 ; i < numberM1 ; i++) {
	  int iColumn = which[numberIntegers-i-1] ;
	  setSequenceInCliqueEntry(cliqueEntry_[numberEntries],iColumn) ;
	  setOneFixesInCliqueEntry(cliqueEntry_[numberEntries],true) ;
	  numberEntries++ ;
	  oneFixStart_[iColumn] = 0 ;
	  zeroFixStart_[iColumn] = 0 ;
	}
      }
      numberCliques_++ ;
      cliqueStart_[numberCliques_] = numberEntries ;
    }
/*
  cliqueEntry now contains a block for each clique, P variables followed by
  M variables. In a U<i> clique, the strong-1 entries precede the strong-0
  entries. In a L<i> clique, the strong-0 entries precede the strong-1
  entries. cliqueStart_ points to the start of the block for each clique.

  Now count the number of times each variable occurs as strong-1 or
  strong-0.  By virtue of previous initialisation, entries for variables not
  in any clique will be -1.

  TODO: I can't see any point in the nested loops here. Just walking
	cliqueEntry_ would achieve the same result.  -- lh, 101007 --
*/
    for (iClique = 0 ; iClique < numberCliques_ ; iClique++) {
      for (int j = cliqueStart_[iClique] ; j < cliqueStart_[iClique+1] ; j++) {
	int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]) ;
	if (oneFixesInCliqueEntry(cliqueEntry_[j]))
	  oneFixStart_[iColumn]++ ;
	else
	  zeroFixStart_[iColumn]++ ;
      }
    }
/*
  We want to set up whichClique_ with blocks of entries for each variable,
  pointing back to the cliques that use the variable.  The block of entries
  for a variable will be divided into two subblocks, first the strong-1
  cliques, then the strong-0 cliques.

  This first loop establishes the start of each strong-1 and strong-0
  subblock in oneFixStart_ and zeroFixStart_. The arrays which and
  endFixStart mirror oneFixStart_ and zeroFixStart_, respectively.
*/
    numberEntries = 0 ;
    for (int iColumn = 0 ; iColumn < numberColumns_ ; iColumn++) {
      if (oneFixStart_[iColumn] >= 0) {
	int n1 = oneFixStart_[iColumn] ;
	int n2 = zeroFixStart_[iColumn] ;
	oneFixStart_[iColumn] = numberEntries ;
	which[iColumn] = numberEntries ;
	numberEntries += n1 ;
	zeroFixStart_[iColumn] = numberEntries ;
	endFixStart_[iColumn] = numberEntries ;
	numberEntries += n2 ;
      }
    }
/*
  Now drop the clique indices into whichClique. which and endFixStart advance
  through each subblock as we process the cliques. When we're done, we can
  discard which (it's equal to zeroFixStart). endFixStart will be one past
  the end of the block for the variable. Note that endFixStart is not quite
  equal to oneFixStart shifted one index, because not all variables are in use
  in cliques --- all of oneFixStart, zeroFixStart, and endFixStart have entries
  of -1 interspersed with meaningful entries.
*/
    for (iClique=0 ;iClique<numberCliques_ ;iClique++) {
      for (int j=cliqueStart_[iClique] ;j<cliqueStart_[iClique+1] ;j++) {
	int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]) ;
	if (oneFixesInCliqueEntry(cliqueEntry_[j])) {
	  int put = which[iColumn] ;
	  which[iColumn]++ ;
	  whichClique_[put]=iClique ;
	} else {
	  int put = endFixStart_[iColumn] ;
	  endFixStart_[iColumn]++ ;
	  whichClique_[put]=iClique ;
	}
      }
    }
  }
  delete [] which ;
  delete [] whichRow ;
  delete [] lookup ;
  return numberCliques_ ;
}

// Delete all clique information
void 
CglProbing::deleteCliques()
{
  delete [] cliqueType_ ;
  delete [] cliqueStart_ ;
  delete [] cliqueEntry_ ;
  delete [] oneFixStart_ ;
  delete [] zeroFixStart_ ;
  delete [] endFixStart_ ;
  delete [] whichClique_ ;
  delete [] cliqueRow_ ;
  delete [] cliqueRowStart_ ;
  cliqueType_ = NULL ;
  cliqueStart_ = NULL ;
  cliqueEntry_ = NULL ;
  oneFixStart_ = NULL ;
  zeroFixStart_ = NULL ;
  endFixStart_ = NULL ;
  whichClique_ = NULL ;
  cliqueRow_ = NULL ;
  cliqueRowStart_ = NULL ;
  numberCliques_ = 0 ;
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


/*
  jjf: Sets up clique information for each row

  Boils down clique information into a form suitable for easy use when
  calculating row bounds U<i> and L<i>. To motivate the process, consider
  calculating U<i> for a row that is entangled with several cliques. If
  the cliques don't intersect, we can treat each of them as a super-variable
  and calculate

    U<i> = (upper bound for non-clique variables) +
		SUM{cliques}(upper bound clique contribution)

  A nontrivial amount of thinking will convince you that all we need to know
  to calculate the upper bound on the contribution of a clique is which
  variables are in a given clique, and the strong direction for the variable.
  See the written doc'n for details. But what if the cliques intersect? We
  don't want to double count. So this method partitions the variables. The
  largest clique goes in intact. Then the clique with the most unclaimed
  variables goes in. And so on, until the number of unclaimed variables is one
  or zero. We're not getting the full strength from cliques later in the
  sequence, but that's conservative.

  TODO: The current algorithm allows later cliques to poach variables that
	are in the intersection with earlier cliques. It's not clear to me
	that this is a good idea. It weakens the earlier (larger) cliques
	with no guarantee of major benefits in terms of strengthening the
	later (smaller) cliques.  -- lh, 101008 --
*/

void 
CglProbing::setupRowCliqueInformation (const OsiSolverInterface &si) 
{
  if (!numberCliques_)
    return ;
/*
  Use a snapshot if we have one, otherwise ask the solver for a copy of the
  constraint matrix. It can happen that CglProbing managed to eliminate rows
  when generating the snapshot, but the column count should be the same.
*/
  CoinPackedMatrix *rowCopy ;
  if (!rowCopy_) {
    numberRows_ = si.getNumRows() ; 
    numberColumns_ = si.getNumCols() ; 
    rowCopy = new CoinPackedMatrix(*si.getMatrixByRow()) ;
  } else {
    rowCopy = rowCopy_ ;
    assert(numberRows_ <= si.getNumRows()) ; 
    assert(numberColumns_ == si.getNumCols()) ; 
  }
  assert(numberRows_ && numberColumns_) ;
  const int *column = rowCopy->getIndices() ;
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *rowLength = rowCopy->getVectorLengths() ; 

  cliqueRowStart_ = new int [numberRows_+1] ;
  cliqueRowStart_[0] = 0 ;

  // Temporary array while building list
  cliqueEntry **array = new cliqueEntry* [numberRows_] ;

  // Clique ids
  int *which = new int [numberCliques_] ;
  // Number of times a given clique is encountered
  int *count = new int [numberCliques_] ;
  // Position of column in row (unfixed variables only; -1 otherwise).
  int *back = new int [numberColumns_] ;
  CoinZeroN(count,numberCliques_) ;
  CoinFillN(back,numberColumns_,-1) ;

  const double *lower = si.getColLower() ;
  const double *upper = si.getColUpper() ;

/*
  Open the main loop to scan rows. It won't be apparent for a while yet, but
  we're only interested in rows that a) were not used to form cliques, and
  b) contain variables that are members of cliques. We need to collect some
  data before we can answer the question.
*/
  int iRow ;
  for (iRow = 0 ; iRow < numberRows_ ; iRow++) {
    int j ;
/*
  Scan this row. Count the number of free variables and record the position
  of each variable in the row. Then update clique entanglement (the clique
  id, and the number of variables entangled with the clique).
*/
    int numberFree = 0 ;
    int numberUsed = 0 ;
    for (j = rowStart[iRow] ; j < rowStart[iRow]+rowLength[iRow] ; j++) {
      int iColumn = column[j] ;
      if (upper[iColumn] > lower[iColumn]) {
        back[iColumn] = j-rowStart[iRow] ;
        numberFree++ ;
        for (int k = oneFixStart_[iColumn] ; k < endFixStart_[iColumn] ; k++) {
          int iClique = whichClique_[k] ;
          if (!count[iClique]) {
            which[numberUsed++] = iClique ;
          }
          count[iClique]++ ;
        }
      }
    }
/*
  jjf: find largest cliques

  Recast the clique entanglement data into a form that's easily useable
  when calculating row bounds.
*/
    bool finished = false ;
    int numberInThis = 0 ;
    cliqueEntry *entries = NULL ;
    array[iRow] = entries ;
    while (!finished) {
/*
  Find the clique entangled with the largest number of variables from this
  row. As the loop progresses, a more accurate description will be `Find the
  clique with the largest number of entangled variables not yet claimed by
  cliques already processed.'
*/
      int largest = 1 ;
      int whichClique = -1 ;
      for (int i = 0 ; i < numberUsed ; i++) {
        int iClique = which[i] ;
        if (count[iClique] > largest) {
          largest = count[iClique] ;
          whichClique = iClique ;
        }
      }
/*
  jjf: Add in if >1 (but not if all as that means clique==row)

  If we formed a clique from this row then largest == numberFree by
  construction and we're done. We'll skip down to the else, set finished =
  true, and move on to the next row. So the rows we're processing here are
  rows that were not used to form cliques but contain 2 or more variables
  that are part of the same clique.

  TODO: This might be overly restrictive. You could picture a scenario where
	one constraint defines a clique, and another constraint, with the same
	variables but arbitrary coefficients, enforces some other restriction.
	(Admittedly, no example leaps to mind.)  -- lh, 101008 --

  TODO: It should be the case that if this is the row that formed the clique,
	count[which[0]] == numberFree. Because we have to suck in all the
	variables, and that means that the very first variable needs to be
	part of a clique, and that clique id will be in which[0]. We could make
	this check immediately after generating the entanglement data.
	-- lh, 101008 --
*/
      if (whichClique >= 0 && largest < numberFree) {
        if (!numberInThis) {
          int length = rowLength[iRow] ;
          entries = new cliqueEntry [length] ;
          array[iRow] = entries ;
          for (int i = 0 ; i < length ; i++) {
            setOneFixesInCliqueEntry(entries[i],false) ;
            setSequenceInCliqueEntry(entries[i],numberColumns_+1) ;
          }
        }
      
/*
  jjf: put in (and take out all counts)

  Scan the row again. For each unfixed variable, look to see if it's a member
  of the clique we're processing. If so, `remove' the variable from the count
  for other cliques where it's a member, and record the variable as a member
  of this clique, along with its strong direction.

  Presumably the rationale for scanning the block for this variable in
  whichClique_, as opposed to the block for this clique in cliqueEntry_,
  is that on average the number of variables in a clique is larger than the
  number of cliques using a variable, so the block in whichClique_ is smaller.
*/
        for (j = rowStart[iRow] ; j < rowStart[iRow]+rowLength[iRow] ; j++) {
          int iColumn = column[j] ;
          if (upper[iColumn] > lower[iColumn]) {
            bool found = false ;
            int k ;
            for (k = oneFixStart_[iColumn] ; k < endFixStart_[iColumn] ; k++) {
              int iClique = whichClique_[k] ;
              if (iClique == whichClique) {
                found = true ;
                break ;
              }
            }

            if (found) {
              for (k = oneFixStart_[iColumn] ;
		   k < endFixStart_[iColumn] ; k++) {
                int iClique = whichClique_[k] ;
                count[iClique]-- ;
              }
/*
  Note that the id given here is only locally significant (for this row). But
  that's ok, all we need to know to use this information is the variables that
  are part of a given clique.

  TODO: Now that I've figured out the motivation for this whole activity, the
	question that remains for me is this: Why do we allow later cliques
	to overwrite the entries for earlier cliques? It seems like we're
	weakening the earlier cliques.  -- lh, 101008 --
*/
              for (k = cliqueStart_[whichClique] ;
		   k < cliqueStart_[whichClique+1] ; k++) {
                if (sequenceInCliqueEntry(cliqueEntry_[k]) == iColumn) {
                  int iback = back[iColumn] ;
                  setSequenceInCliqueEntry(entries[iback],numberInThis) ;
                  setOneFixesInCliqueEntry(entries[iback],
				 oneFixesInCliqueEntry(cliqueEntry_[k])) ;
                  break ;
                }
              }
            }
          }
        }
        numberInThis++ ;
      } else {
        finished = true ;
      }
    }
/*
  Reserve space in the final output array.
*/
    if (numberInThis) 
      cliqueRowStart_[iRow+1] = cliqueRowStart_[iRow]+rowLength[iRow] ;
    else
      cliqueRowStart_[iRow+1] = cliqueRowStart_[iRow] ;
/*
  Clear the bookkeeping arrays in preparation for the next row.
*/
    for (int i = 0 ; i < numberUsed ; i++) {
      int iClique = which[i] ;
      count[iClique] = 0 ;
    }
    for (j = rowStart[iRow] ; j < rowStart[iRow]+rowLength[iRow] ; j++) {
      int iColumn = column[j] ;
      back[iColumn] = -1 ;
    }
  }
/*
  End main loop. Clean up, store the results, and go home. Rather than keep all
  the row-length fragments allocated above for individual entries arrays,
  copy them off into one big array, cliqueRow_, with the starting position for
  each row indexed through cliqueRowStart_.
*/
  delete [] which ;
  delete [] count ;
  delete [] back ;
  cliqueRow_ = new cliqueEntry [cliqueRowStart_[numberRows_]] ;
  for (iRow = 0 ; iRow < numberRows_ ; iRow++) {
    if (array[iRow]) {
      int start = cliqueRowStart_[iRow] ;
      CoinMemcpyN(array[iRow],rowLength[iRow],cliqueRow_+start) ;
      delete [] array[iRow] ;
    }
  }
  delete [] array ;
  if (rowCopy != rowCopy_)
    delete rowCopy ;
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

