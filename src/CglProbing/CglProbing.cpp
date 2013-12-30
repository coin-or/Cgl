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


#include "CoinPackedMatrix.hpp"
#include "CglProbing.hpp"
#include "CglProbingRowCut.hpp"
#include "CglProbingDebug.hpp"
#include "CglPhic.hpp"

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

}  // end file local namespace

/*
  Compare two sets of upper and lower bound arrays and generate column cuts.
*/

int CglProbing::makeColCuts (int nCols, OsiCuts &cs,
			     const char *const intVar,
			     const double *const origsol,
			     const double *const origlbs,
			     const double *const origubs,
			     const double *const newlbs,
			     const double *const newubs) const
{
# if CGL_DEBUG > 0
  if (verbosity_ >= 3)
    std::cout << "      " << "start makeColCuts." << std::endl ;
# endif
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
      if (newubs[i] < origubs[i]-1.0e-8) {
	if (newubs[i] < origsol[i]-1.0e-8) ifCut++ ;
#       if CGL_DEBUG > 0
	if (verbosity_ >= 6)
	  std::cout
	    << "          " << "u<" << i << "> " << origubs[i] << " --> "
	    << newubs[i] << ", diff " << origubs[i]-newubs[i]
	    << "." << std::endl ;
#       endif
	ubs.insert(i,newubs[i]) ;
	numberChanged++ ;
      }
      if (newlbs[i] > origlbs[i]+1.0e-8) {
	if (newlbs[i] > origsol[i]+1.0e-8) ifCut++ ;
#       if CGL_DEBUG > 0
	if (verbosity_ >= 6)
	  std::cout
	    << "          " << "l<" << i << "> " << origlbs[i] << " --> "
	    << newlbs[i] << ", diff " << newlbs[i]-origlbs[i]
	    << "." << std::endl ;
#       endif
	lbs.insert(i,newlbs[i]) ;
	numberChanged++ ;
      }
    }
  }
/*
  Stash the changes in a column cut. If the changes actually cut off some part
  of the current solution, boost the effectiveness.
*/
# if CGL_DEBUG > 0
  if (verbosity_ >= 2)
    std::cout
      << "      " << "improved " << numberChanged << " bounds; "
      << ifCut << " cuts." << std::endl ;
# endif
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
# if CGL_DEBUG > 0
  if (verbosity_ >= 3)
    std::cout << "      " << "end makeColCuts." << std::endl ;
# endif
  return (numberChanged) ;
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

bool CglProbing::analyze (const OsiSolverInterface &si,
			  char *const intVar,
			  double *const lbs, double *const ubs) const
{
# if CGL_DEBUG > 0
  if (verbosity_ >= 3)
    std::cout << "    " << "start analyze." << std::endl ;
# endif
  const double e20Inf = 1.0e20 ;
  const int changedToInt = 77 ;
/*
  Acquire various data structures for ease of work: 
*/
  const int n = si.getNumCols() ;
  const int m = si.getNumRows() ;
  const double *const objective = si.getObjCoefficients() ;
  const double objsense = si.getObjSense() ;

  // Row copy
  const CoinPackedMatrix *const matrixByRow = si.getMatrixByRow() ;
  const double *const elementByRow = matrixByRow->getElements() ;
  const int *const colIndices = matrixByRow->getIndices() ;
  const CoinBigIndex *const rowStarts = matrixByRow->getVectorStarts() ;
  const int *const rowLens = matrixByRow->getVectorLengths() ;

  // Column copy
  const CoinPackedMatrix *const matrixByCol = si.getMatrixByCol() ;
  const double *const elementByCol = matrixByCol->getElements() ;
  const int *const rowIndices = matrixByCol->getIndices() ;
  const CoinBigIndex *const colStarts = matrixByCol->getVectorStarts() ;
  const int *const colLens = matrixByCol->getVectorLengths() ;

  const double *const rowLower = si.getRowLower() ;
  const double *const rowUpper = si.getRowUpper() ;

  char *const ignore = new char [m] ;
  int *const which = new int[m] ;
  double *const changeRhs = new double[m] ;
  CoinZeroN(changeRhs,m) ;
  CoinZeroN(ignore,m) ;

  int numberChanged = 0 ;
  bool finished = false ;
/*
  Open main analysis loop. Keep going until nothing changes.
*/
  while (!finished) {
    int saveNumberChanged = numberChanged ;
/*
  Open loop to scan each constraint. With luck, we can conclude that some
  variable must be integer. With less luck, the constraint will not prevent
  the variable from being integer.
*/
    for (int i = 0 ; i < m ; i++) {
      int numberContinuous = 0 ;
      double value1 = 0.0,value2 = 0.0 ;
      bool allIntegerCoeff = true ;
      double sumFixed = 0.0 ;
      int jColumn1 = -1 ;
      int jColumn2 = -1 ;
/*
  Scan the coefficients of the constraint and accumulate some information:
    * Contribution of fixed variables.
    * Count of continuous variables, and column index and coefficients for
      the first and last continuous variable.
    * A boolean, allIntegerCoeff, true if all coefficients of unfixed
      integer variables are integer.
*/
      const CoinBigIndex startRowi = rowStarts[i] ;
      const CoinBigIndex endRowi = startRowi+rowLens[i] ;
      for (CoinBigIndex jj = startRowi ; jj < endRowi ; jj++) {
        int j = colIndices[jj] ;
        double aij = elementByRow[jj] ;
        if (ubs[j] > lbs[j]+1.0e-8) {
          if (!intVar[j]) {
            if (numberContinuous == 0) {
              jColumn1 = j ;
              value1 = aij ;
            } else {
              jColumn2 = j ;
              value2 = aij ;
            }
            numberContinuous++ ;
          } else {
            if (fabs(aij-floor(aij+0.5)) > 1.0e-12)
              allIntegerCoeff = false ;
          }
        } else {
          sumFixed += lbs[j]*aij ;
        }
      }
/*
  See if the row bounds are integer after adjusting for the contribution of
  fixed variables.
*/
      double blowi = rowLower[i] ;
      if (blowi > -e20Inf) {
        blowi -= sumFixed ;
        if (fabs(blowi-floor(blowi+0.5)) > 1.0e-12)
          allIntegerCoeff = false ;
      }
      double bi = rowUpper[i] ;
      if (bi < e20Inf) {
        bi -= sumFixed ;
        if (fabs(bi-floor(bi+0.5)) > 1.0e-12)
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
*/
      if (numberContinuous == 1) {
        if (blowi == bi) {
          if (fabs(value1) > 1.0e-3) {
            value1 = 1.0/value1 ;
            if (fabs(value1-floor(value1+0.5)) < 1.0e-12) {
              numberChanged++ ;
              intVar[jColumn1] = changedToInt ;
#	      if CGL_DEBUG > 0
	      if (verbosity_ > 2) {
	      	std::cout
		  << "    " << "x<" << jColumn1 << "> can be integer; "
		  << "Single continuous variable in integer equality."
		  << std::endl ;
	      }
#	      endif
            }
          }
        } else {
          if (fabs(value1) > 1.0e-3) {
            value1 = 1.0/value1 ;
            if (fabs(value1-floor(value1+0.5)) < 1.0e-12) {
              ignore[i] = 1 ;
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
        if (blowi == bi) {
          if (fabs(value1) == 1.0 && value1*value2 == -1.0 &&
	      !lbs[jColumn1] && !lbs[jColumn2] &&
	      colLens[jColumn1] == 1 && colLens[jColumn2] == 1) {
            int n = 0 ;
            int i ;
            double objChange =
		objsense*(objective[jColumn1]+objective[jColumn2]) ;
            double bound = CoinMin(ubs[jColumn1],ubs[jColumn2]) ;
            bound = CoinMin(bound,e20Inf) ;
	    // Since column lengths are 1, there is no second coeff.
            for (i = colStarts[jColumn1] ;
		 i < colStarts[jColumn1]+colLens[jColumn1] ; i++) {
              int jRow = rowIndices[i] ;
              double value = elementByCol[i] ;
              if (jRow != i) {
                which[n++] = jRow ;
                changeRhs[jRow] = value ;
              }
            }
            for (i = colStarts[jColumn2] ;
		 i < colStarts[jColumn2]+colLens[jColumn2] ; i++) {
              int jRow = rowIndices[i] ;
              double value = elementByCol[i] ;
              if (jRow != i) {
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
                  if (rowLens[jRow] == 1) {
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
                  } else if (rowLens[jRow] == 2) {
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
    for (int j = 0 ; j < n ; j++) {
      if (ubs[j] > lbs[j]+1.0e-8 && !intVar[j]) {
        double value ;
        value = ubs[j] ;
        if (value < e20Inf && fabs(value-floor(value+0.5)) > 1.0e-12) 
          continue ;
        value = lbs[j] ;
        if (value > -e20Inf && fabs(value-floor(value+0.5)) > 1.0e-12) 
          continue ;
        bool integer = true ;
        for (CoinBigIndex ii = colStarts[j] ;
	     ii < colStarts[j]+colLens[j] ; ii++) {
          if (!ignore[rowIndices[ii]]) {
            integer = false ;
            break ;
          }
        }
        if (integer) {
          numberChanged++ ;
          intVar[j] = changedToInt ;
#	  if CGL_DEBUG > 0
	  if (verbosity_ > 2) {
	    std::cout
	      << "    " << "x<" << j << "> can be integer; "
	      << "no constraint prevents it." << std::endl ;
	  }
#	  endif
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
  for (int j = 0 ; j < n ; j++) {
    if (intVar[j] == changedToInt) {
      if (ubs[j] > e20Inf) {
        ubs[j] = e20Inf ;
      } else {
        ubs[j] = floor(ubs[j]+1.0e-5) ;
      }
      if (lbs[j] < -e20Inf) {
        lbs[j] = -e20Inf ;
      } else {
        lbs[j] = ceil(lbs[j]-1.0e-5) ;
        if (lbs[j] > ubs[j])
          feasible = false ;
      }
      if (lbs[j] == 0.0 && ubs[j] == 1.0)
        intVar[j] = 1 ;
      else if (lbs[j] == ubs[j])
        intVar[j] = 0 ;
      else
        intVar[j] = 2 ;
    }
  }
/*
  Clean up and return.
*/
  delete [] which ;
  delete [] changeRhs ;
  delete [] ignore ;

# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout
      << "    " << "end analyze, " << ((feasible)?"feasible":"infeasible") ;
    if (numberChanged)
      std::cout
        << ", " << numberChanged << " variables could be made integer" ;
    std::cout << "." << std::endl ;
  }
# endif

  return (feasible) ;
}


/*
  Set up the working row copy.
  
  If we have a snapshot handy, most of the setup is done -- we just make
  a copy of the snapshot.  Otherwise, we need to create a row-major copy
  of the constraint matrix and the rhs and rhslow vectors.  If the client
  has specified mode 0 (work from snapshot), we silently correct it.
*/

CoinPackedMatrix *CglProbing::setupRowCopy (int mode, bool useObj,
					    double cutoff, double offset,
					    const OsiSolverInterface &si,
					    double *const rowLower,
					    double *const rowUpper) const
{
# if CGL_DEBUG > 0
  if (verbosity_ >= 3)
    std::cout
      << "    " << "start setupRowCopy, mode " << mode << "." << std::endl ;
# endif

  int nCols = si.getNumCols() ;
  int nRows = -1 ;
  CoinPackedMatrix *rowCopy = NULL ;
/*
  Make the working copy from the model in the si. If we're using the
  objective, add in the objective coefficients.
*/
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

# if CGL_DEBUG > 0
  if (paranoia_ > 0) {
    if (rowCopy) {
      int errs = 0 ;
      if (paranoia_ >= 3) errs = rowCopy->verifyMtx(1) ;
      if (errs > 0) {
	std::cout
	  << "    generateCutsAndModify: verifyMtx failed post-cutgen, "
	  << errs << " errors." << std::endl ;
	assert (errs == 0) ;
      }
    } else {
      std::cout << "ERROR! Failed to create working row copy!" << std::endl ;
      assert (false) ;
    }
  }
  if (verbosity_ >= 3)
    std::cout << "    " << "end setupRowCopy." << std::endl ;
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

  An alternate code block for deleting constraints on bounds would delete
  any constraint where an infinite (1.0e3 :-) bound on a variable would
  lead to an infinite row LB or UB. This is a recipe for deleting the
  entire model, when run prior to an initial round of bound propagation. It's
  been disabled since r675 (2008). I've chopped it entirely. Something useful
  could be done along this line, but only after an initial round of bound
  propagation. Arguably, an initial round of bound propagation would help all
  of CglProbing (not to mention CglPreProcess) and should be done in any event.

  Returns true if we're good to proceed, false if the groomed model is not
  worth pursuing.

  TODO: Remove the code that sorts coefficients in a row so that negative
	precedes positive. Used to be used in the propagation code over
	in probe(), but no longer. Keep it in place just to minimise the bugs
	as I do the conversion to CglPhic.  -- lh, 110405 --

*/

bool CglProbing::groomModel (bool useObj, int maxRowLen,
			     const OsiSolverInterface &si,
			     const char *const intVar,
			     CoinPackedMatrix *rowCopy,
			     double *const rowLower, double *const rowUpper,
			     const double *const colLower,
			     const double *const colUpper,
			     int *&realRows,
			     const CglTreeInfo *const info) const
{
# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout << "    " << "start groomModel, strengthening " ;
    if (info->strengthenRow && !info->pass) std::cout << "not " ;
    std::cout << "allowed, max coeffs " << maxRowLen << "." << std::endl ;
  }
# endif

  realRows = 0 ;

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
  int denseCoeffs = 0 ;
  int nWeak = 0 ;
  int nEmpty = 0 ;
  int nFixed = 0 ;
/*
  Preliminary assessment for dense rows, 'free' rows, empty rows, and rows
  with no unfixed variables.  We're looking for rows that are too long or
  have row bounds so weak as to be ineffective.  If the total coefficients
  in dense rows amounts to less than 1/10 of the matrix, just deal with it.

  At the end of the markup, which[i] <= 0 indicates the reason for deletion.
*/
  for (int i = 0 ; i < nRealRows ; i++) {
    const int &leni = rowLength[i] ;
    which[i] = leni ;
    if ((rowLower[i] < -weakRhs) && (rowUpper[i] > weakRhs)) {
      nWeak++ ;
      which[i] = -1 ;
    } else if (leni > maxRowLen) {
      nDense++ ;
      denseCoeffs += leni ;
      which[i] = -2 ;
    } else if (leni == 0) {
      nEmpty++ ;
    }
    // Not obviously bogus; scan for unfixed variables.
    const CoinBigIndex &rstarti = rowStart[i] ;
    const CoinBigIndex rendi = rstarti+leni ;
    bool allFixed = true ;
    for (CoinBigIndex jj = rstarti ; jj < rendi ; jj++) {
      int j = column[jj] ;
      if (colUpper[j]-colLower[j] > primalTolerance_) {
        allFixed = false ;
	break ;
      }
    }
    if (allFixed) {
      which[i] = -3 ;
      nFixed++ ;
    }
  }
  if (denseCoeffs*10 < nElements) maxRowLen = nCols ;
# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout << "      " << nDense << " dense " ;
    if (nDense) std::cout << " with " << denseCoeffs << " coeffs" ;
    std::cout
      << ", " << nWeak << " weak, " << nEmpty << " empty, "
      << nFixed << " fixed" ;
    if (maxRowLen == nCols) std::cout << "; keeping dense rows" ;
    std::cout << "." << std::endl ;
  }
# endif
/*
  Walk the rows again. Rows that should be deleted will come up positive
  when we subtract maxRowLen and can be summarily dismissed.  Note that we're
  converting which[] as we go, so that at the end it will contain the
  indices of the constraints to be deleted in a block at the beginning.
*/
  int nDelete = 0 ;
  int nKeep = 0 ;
  for (int i = 0 ; i < nRealRows ; i++) {
    if (which[i] <= 0) {
#     if CGL_DEBUG > 0
      if (verbosity_ >= 5) {
        char *why[] = {"empty","weak","dense","fixed"} ;
        std::cout
	  << "          " << "deleting r(" << i << "), " << why[-which[i]]
	  << std::endl ;
	const CoinBigIndex &jj = rowStart[i] ;
	CglProbingDebug::dump_row(i,rowLower[i],rowUpper[i],nan(""),nan(""),
				  0,true,true,10,rowLength[i],&column[jj],
				  &elements[jj],primalTolerance_,
				  colLower,colUpper) ;
      }
#     endif
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
    const CoinBigIndex &rstarti = rowStart[i] ;
    const int &leni = rowLength[i] ;
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
#         if CGL_DEBUG > 0
	  if (verbosity_ >= 2) {
	    std::cout
	      << "      r(" << i << "): decrease a<" << i << ","
	      << column[jj] << "> from " << value << " to "
	      << elements[jj] << "." << std::endl ;
	  }
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
#         if CGL_DEBUG > 0
	  if (verbosity_ >= 2) {
	    std::cout
	      << "      r(" << i << "): increase a<" << i << ","
	      << column[jj] << "> from " << value << " to "
	      << elements[jj] << "." << std::endl ;
	  }
#	  endif
	}
	effectiveness += fabs(value) ;
      }
      rowLower[nKeep-1] = 0.0 ;
    }
/*
  If we actually strengthened a coefficient, groom the row and stash the
  strengthened row in strengthenRow for return to the caller.
*/
    if (effectiveness) {
      OsiRowCut *rc = new OsiRowCut() ;
      rc->setLb(rowLower[nKeep-1]) ;
      rc->setUb(rowUpper[nKeep-1]) ;
      rc->setRow(leni,column+rstarti,elements+rstarti,false) ;
      CoinPackedVector &row = rc.mutableRow() ;
      double *elements = row.getElements() ;
      int k = 0 ;
      for ( ; k < leni ; k++)
        if (fabs(elements[k] < 1.0e-12) break ;
      if (k < leni) {
        int *columns = row.getIndices() ;
	for (int j = k+1 ; j < leni ; j++) {
	  if (fabs(elements[j] > 1.0e-12) {
	    elements[k] = elements[j] ;
	    columns[k++] = columns[j] ;
	  }
	}
	row.truncate(k) ;
      }
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
      if (verbosity_ >= 1)
	std::cout
	  << "  All rows unsuitable for probing! Cleaning up." << std::endl ;
#   endif
    delete[] realRows ;
    return (false) ;
  }
/*
  So far, so good. We have something that looks like a reasonable collection
  of rows. Time to deal with fixed variables. We're going to walk the rows,
  compressing out zeros and coefficients of fixed variables, and we'll sort
  the remaining coefficients so that all negative coefficients precede
  all positive coefficients.  There will be no gaps when we're done.

  FIXME: Figure out the appropriate tolerance to use to decide if a
	 coefficient is too small to retain in the matrix.

  NOTE that columns are not physically removed.

  We need to refresh the mutable pointers; deleteRows above may have rendered
  the previous pointers invalid.
*/
  elements = rowCopy->getMutableElements() ;
  column = rowCopy->getMutableIndices() ;
  rowStart = rowCopy->getMutableVectorStarts() ;
  rowLength = rowCopy->getMutableVectorLengths(); 

  CoinBigIndex newSize = 0 ;
  int *columnPos = new int [nCols] ;
  double *elementsPos = new double [nCols] ;
  CoinBigIndex *rowStartPos = new CoinBigIndex [nRows] ;
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
        // tolerance?
	if (value < 0.0) {
	  elements[newSize] = value ;
	  column[newSize++] = j ;
	} else if (value > 0.0) {
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
  delete [] rowStartPos ;
  rowStart[nRows] = newSize ;
  rowCopy->setNumElements(newSize) ;

# if CGL_DEBUG > 0
  if (verbosity_ >= 3)
    std::cout << "    " << "end groomModel." << std::endl ;
# endif

  return (true) ;
} 

//-------------------------------------------------------------------
// Generate column, disaggregation, and implication cuts
//------------------------------------------------------------------- 

/*
  The traditional const method. Row and column cuts will be returned in cs.

  But note Stefan went through and removed const/mutable in r1119.
*/
void CglProbing::generateCuts (const OsiSolverInterface &si, OsiCuts &cs,
			       const CglTreeInfo info2)
{

# if CGL_DEBUG > 0
  verbosity_ = 0 ;
  paranoia_ = 0 ;
  if (verbosity_ >= 3)
    std::cout
      << "CglProbing::generateCuts, matrix " << si.getNumRows()
      << " x " << si.getNumCols() << "." << std::endl ;
  int retval =
      CglProbingDebug::checkForRowCutDebugger(si,paranoia_,verbosity_) ;
  const bool optPathOnEntry = (retval == 2)?true:false ;
  if (verbosity_ >= 3) CglProbingDebug::dumpSoln(si,primalTolerance_) ;
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
*/
  int nRows = si.getNumRows() ; 
  double *rowLower = new double[nRows+1] ;
  double *rowUpper = new double[nRows+1] ;

  int nCols = si.getNumCols() ;
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
  bool feasible = gutsOfGenerateCuts(si,cs,rowLower,rowUpper,
  				     colLower,colUpper,&info) ;
# if CGL_DEBUG > 0
  if (paranoia_ > 0) {
    int errs = 0 ;
    const CoinPackedMatrix *mtx = si.getMatrixByRow() ;
    if (paranoia_ >= 4) {
      errs = mtx->verifyMtx(1) ;
      if (errs > 0) {
	std::cout
	  << "    generateCuts: verifyMtx failed, "
	  << errs << " errors." << std::endl ;
	assert (false) ;
      } else {
	std::cout << "    generateCuts: matrix verified." << std::endl ;
      }
    }
  }
# endif
/*
  Infeasibility is indicated by returning a stylized infeasible column cut.

  If we're debugging, do a sanity check that this cut really is infeasible.
*/
  if (!feasible) {
    OsiRowCut rc ;
    rc.setLb(DBL_MAX) ;
    rc.setUb(0.0) ;   
    cs.insert(rc) ;
#   if CGL_DEBUG > 0
    if (paranoia_ > 0) {
      const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
      if (debugger)
	assert(!debugger->invalidCut(rc)) ; 
    }
#   endif
  }
/*
  Delete the working arrays and restore the row cut mode. Setting colLower_
  and colUpper_ to NULL ensures that the client doesn't accidentally think it
  has access to the tightened bound arrays. Contrast with the end of
  generateCutsAndModify.
*/
  delete [] rowLower ;
  delete [] rowUpper ;
  delete [] colLower ;
  delete [] colUpper ;
  colLower_ = NULL ;
  colUpper_ = NULL ;

  rowCuts_ = saveRowCuts ;

# if CGL_DEBUG > 0
  if (verbosity_ >= 1) {
    std::cout
      << "  " << "end probing: " << ((feasible)?"feasible":"infeasible") ;
    if (cs.sizeCuts() > 0)
      std::cout
        << ", " << cs.sizeRowCuts() << " row cuts"
	<< ", " << cs.sizeColCuts() << " col cuts" ;
    std::cout << "." << std::endl ;
  }
  if (paranoia_ > 0 && optPathOnEntry && !feasible) {
    std::cout
      << "ERROR: entered on optimal path but result is infeasible.!"
      << std::endl ;
    assert(false) ;
  }
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
  verbosity_ = 0 ;
  paranoia_ = 0 ;
  if (verbosity_ >= 3)
    std::cout
      << "CglProbing::generateCutsAndModify, matrix " << si.getNumRows()
      << " x " << si.getNumCols() << "." << std::endl ;
  int retval =
      CglProbingDebug::checkForRowCutDebugger(si,paranoia_,verbosity_) ;
  const bool optPathOnEntry = (retval == 2)?true:false ;
  if (verbosity_ >= 3) CglProbingDebug::dumpSoln(si,primalTolerance_) ;
# endif

/*
  Enforce the rule that a negative value prohibits row cut generation in the
  search tree. (Code 0x4 is `column cuts only'.)
*/
  int saveRowCuts = rowCuts_ ;
  if (rowCuts_ < 0) {
    if (info->inTree)
     rowCuts_ = 0x4 ;
    else
      rowCuts_ = -rowCuts_ ;
  }
  int saveMode = mode_ ;
/*
  Create working arrays for row and column bounds.

  We need nRows+1 in rowLower, rowUpper because we may use the objective
  function as a constraint in the working system within gutsOfGenerateCuts.

  Clear out old variable bounds arrays now, so we don't have to worry about it
  later.
*/
  int nRows = si.getNumRows() ; 
  double *rowLower = new double[nRows+1] ;
  double *rowUpper = new double[nRows+1] ;

  int nCols = si.getNumCols() ;
  if (colLower_) {
    delete[] colLower_ ;
    colLower_ = NULL ;
  }
  if (colUpper_) {
    delete[] colUpper_ ;
    colUpper_ = NULL ;
  }
  double *colLower = new double[nCols] ;
  double *colUpper = new double[nCols] ;

# if CGL_DEBUG > 0
  if (paranoia_ > 0) {
    int errs = 0 ;
    const CoinPackedMatrix *mtx = si.getMatrixByRow() ;
    if (paranoia_ >= 4) errs = mtx->verifyMtx(1) ;
    if (errs > 0) {
      std::cout
        << "    generateCutsAndModify: verifyMtx failed pre-cutgen, "
	<< errs << " errors." << std::endl ;
      assert (errs == 0) ;
    }
  }
# endif

/*
  Do the work. There are no guarantees about the state of the row and column
  bounds arrays if we come back infeasible.
*/
  bool feasible = gutsOfGenerateCuts(si,cs,rowLower,rowUpper,
  				     colLower,colUpper,info) ;
# if CGL_DEBUG > 0
  if (paranoia_ > 0) {
    int errs = 0 ;
    const CoinPackedMatrix *mtx = si.getMatrixByRow() ;
    if (paranoia_ >= 4) errs = mtx->verifyMtx(1) ;
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
  if (!feasible) {
    OsiRowCut rc ;
    rc.setLb(DBL_MAX) ;
    rc.setUb(0.0) ;   
    cs.insert(rc) ;
#   if CGL_DEBUG > 0
    if (paranoia_ > 0) {
      const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
      if (debugger)
	assert(!debugger->invalidCut(rc)) ;
    }
#   endif
  }
/*
  Restore modes.
*/
  rowCuts_ = saveRowCuts ;
  mode_ = saveMode ;
/*
  Don't replace bounds supplied by the client if we've gone infeasible.
  There's no guarantee of consistency in the values reported back from
  gutsOfGenerateCuts.
*/
  delete[] rowLower ;
  delete[] rowUpper ;

  if (feasible) {
    colLower_ = colLower ;
    colUpper_ = colUpper ;
  } else {
    delete[] colUpper ;
    delete[] colLower ;
  }

# if CGL_DEBUG > 0
  if (verbosity_ >= 1) {
    std::cout
      << "  " << "end probing: " << ((feasible)?"feasible":"infeasible") ;
    if (cs.sizeCuts() > 0)
      std::cout
        << ", " << cs.sizeRowCuts() << " row cuts"
	<< ", " << cs.sizeColCuts() << " col cuts" ;
    std::cout << "." << std::endl ;
  }
  if (paranoia_ > 0 && optPathOnEntry && !feasible) {
    std::cout
      << "ERROR: entered on optimal path but result is infeasible.!"
      << std::endl ;
    assert(false) ;
  }
# endif

  return (static_cast<int>(feasible)) ;
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
bool CglProbing::gutsOfGenerateCuts (const OsiSolverInterface &si, 
				     OsiCuts &cs ,
				     double *rowLower, double *rowUpper,
                                     double *colLower, double *colUpper,
                                     CglTreeInfo *info)
{

# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout
      << "  start gutsOfGenerateCuts, mode " << mode_ << ", " ;
    if (!info->inTree) std::cout << "not " ;
    std::cout << "in tree, pass " << info->pass << "." << std::endl ;
  }
# endif

  int numberRowCutsBefore = cs.sizeRowCuts() ;

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
  const double *const colsol = si.getColSolution() ;
/*
  Stage 1: Preliminary Processing

  During Stage 1 we'll get rid of constraints that are judged unsuitable for
  further probing (too dense, or not amenable to propagation or implication),
  and do one bit of really simple coefficient tightening. At the end of Stage
  1, we'll have a row- and column-major copy of the reduced and strengthened
  constraint system.

  If we're at the root of the search tree, scan the constraint system to see
  if there are naturally integer variables that are not declared as integer.
  Convert them to integer type.

  If we lose feasibility here, cut our losses and bail out while it's still
  uncomplicated.
*/
  if (!info->inTree && !info->pass) {
    feasible = analyze(si,intVar,colLower,colUpper) ;
    if (!feasible) {
      delete[] intVar ;
      delete[] colLower ;
      colLower_ = 0 ;
      delete[] colUpper ;
      return (feasible) ;
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
  if (verbosity_ > 3) {
    std::cout << "      " << "cutoff " ;
    if (!cutoff_available)
      std::cout << "not available" ;
    else
      std::cout << cutoff ;
    std::cout << ", offset " << offset << "." << std::endl ;
    std::cout << "      " ;
    if (useCutoff)
      std::cout << "using cutoff" ;
    else
      std::cout << "  not using cutoff" ;
    if (useObj)
      std::cout << ", using objective as constraint" ;
    else
      std::cout << ", not using objective as constraint" ;
    std::cout << "." << std::endl ;
  }
# endif

/*
  Create the working row copy from the si's model or from a snapshot (only if
  the user specified mode 0 and a snapshot exists). rowLower and rowUpper are
  valid at completion, but if we initialise the local row copy from a snapshot
  the number of valid entries may not match the number of rows for the system
  held in the si.
*/
  CoinPackedMatrix *rowCopy = setupRowCopy(mode,useObj,cutoff,
  					   offset,si,rowLower,rowUpper) ;
  int nRows = rowCopy->getNumRows() ;
/*
  Groom the rowCopy: delete useless rows, substitute for fixed variables, sort
  coefficients so that negative precede positive, and try some specialised
  coefficient strengthening on the rows that remain. If groomModel comes
  back false, it's an indication that probing isn't worth pursuing. Clean
  up and go home while it's relatively easy.

  realRows is the cross-reference between rows in rowCopy and rows of the
  original model.

  Note that the number of rows in the local copy is almost certainly not the
  same as the number of rows in the system held in the si. You can't just copy
  rowLower and rowUpper; they must be translated if realRows exists.
*/
  int *realRows = 0 ;
  bool worthProbing = groomModel(useObj,maxRowLen,si,intVar,
				 rowCopy,rowLower,rowUpper,
				 colLower,colUpper,realRows,info) ;
  nRows = rowCopy->getNumRows() ;
  if (!worthProbing) {
    delete[] intVar ;
    delete[] realRows ;
    delete rowCopy ;
    return (feasible) ;
  }
/*
  Make a column-major copy, after we've winnowed the rows to the set we want
  to work with.
*/
  CoinPackedMatrix *columnCopy = new CoinPackedMatrix(*rowCopy,0,0,true) ;
/*
  We now have row- and column-major copies of the working constraint system.
  Set up and invoke bound propagation to tighten up the row and column bounds.
*/
  CglPhic phic(rowCopy,columnCopy,rowLower,rowUpper) ;
  phic.setVerbosity(verbosity_) ;
  // phic.setVerbosity(4) ;
  phic.setParanoia(0) ;
  phic.setRowPropTol(1.0e-3) ;
  phic.setColPropTol(1.0e-4) ;
  phic.setDisturbTol(1.0e-1) ;
  phic.setPropType(CglPhic::PropGenInt|CglPhic::PropBinary) ;
  phic.loanColBnds(colLower,colUpper) ;
  phic.loanColType(intVar) ;
  phic.initLhsBnds() ;
  phic.initPropagation() ;
  phic.tightenAbInitio(feasible) ;
/*
  Cut back on propagation now that we've tightened the system.
*/
  phic.setRowPropTol(1.0) ;
  phic.setPropType(CglPhic::PropBinary) ;
  // phic.setPropType(CglPhic::PropGenInt|CglPhic::PropBinary) ;
  // phic.setPropType(CglPhic::PropCon|CglPhic::PropGenInt|CglPhic::PropBinary) ;
  // phic.setVerbosity(6) ;
/*
  If we've lost feasibility, do nothing more. Otherwise, record the new bounds
  as column cuts and continue with setup and probing.
*/
  if (feasible) {
    phic.clearPropagation() ;
    int numTightened = makeColCuts(nCols,cs,intVar,colsol,
    				   si.getColLower(),si.getColUpper(),
				   colLower,colUpper) ;
#   if CGL_DEBUG > 0
    if (verbosity_ > 3 && numTightened) {
      std::cout
	<< "      " << "initial processing improved " << numTightened
	<< " variable bounds." << std::endl ;
    }
    if (cs.sizeColCuts() && paranoia_ > 0) {
      assert(CglProbingDebug::checkBounds(si,
			cs.colCut((cs.sizeColCuts()-1)),verbosity_) == 0) ;
    }
#   endif
/*
  We've tightened up the initial system as best we can. To go further,
  we need to be able to probe.
*/
    if (maxProbe > 0) {
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
      numberThisTime_ = 0 ;
      if (!lookedAt_) {
	lookedAt_ = new int[nCols] ;
      }
      if (mode == 1) {
	double_int_pair *array = new double_int_pair [nCols] ;
#       ifdef ZEROFAULT
	CoinZeroN(array,nCols) ;
#	endif
	const double multiplier =
	    ((info->inTree || (info->pass&1) != 0))?1.0:-1.0 ;
/*
  Examine unfixed integer variables and record the ones that are unsatisfied
  (i.e., not at an integral value). The multiplier serves to reverse the sort
  order on alternate calls.

  TODO: the constant .49999 is actually a hardwired integer feasibility
	tolerance of 1e-5. This should be a parameter; better, an attribute
	of the CglProbing object. Presumably the game with .49999-(*) and
	multiplier is a convenience for sorting (most infeasible are close
	to zero, nearly integral are large and negative). But the comparison
	method (double_int_pair_compare) isn't used elsewhere, so maybe this
	can be simplified.  -- lh, 100924 --
*/
	for (int j = 0 ; j < nCols ; j++) {
	  if (intVar[j] && (colUpper[j]-colLower[j]) > 1.0e-8) {
	    double away = fabs(0.5-(colsol[j]-floor(colsol[j]))) ;
	    if (away < 0.49999 || !info->inTree) {
	      //array[numberThisTime_].infeasibility = away ;
	      array[numberThisTime_].infeasibility = away*multiplier ;
	      //array[numberThisTime_].infeasibility = -columnLength[i] ;
	      array[numberThisTime_++].sequence = j ;
	    }
	  }
	}
	std::sort(array,array+numberThisTime_,double_int_pair_compare()) ;
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
#     if CGL_DEBUG > 0
      if (verbosity_ >= 3) {
        std::cout
	  << "      " << "probing " << numberThisTime_ ;
	if (mode == 1 && info->inTree) std::cout << " unsatisfied" ;
	std::cout << " variables." << std::endl ;
      }
#     endif
/*
  lookedAt_ now contains a list of indices of variables to probe. Do the
  probing.
*/
      feasible = probe(si,cs,phic,realRows,info,useObj,useCutoff,cutoff) ;
    }
  }
/*
  Done! Time to clean up and go home.
*/
  delete rowCopy ;
  delete columnCopy ;

  delete[] intVar ;
  delete[] realRows ;
/*
  If we're requested to set the global cut flag, do it.
*/
  if (!info->inTree &&
      ((info->options&4) == 4 || ((info->options&8) && !info->pass))) {
    int numberRowCutsAfter = cs.sizeRowCuts() ;
    for (int i = numberRowCutsBefore ; i < numberRowCutsAfter ; i++)
      cs.rowCutPtr(i)->setGloballyValid() ;
  }

# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout
      << "      " << "end gutsOfGenerateCuts, "
      << ((feasible)?"feasible":"infeasible")
      << "; " << cs.sizeRowCuts() << " row cuts, " << cs.sizeColCuts()
      << " col cuts." << std::endl ;
  }
  if (paranoia_ > 0) {
    const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
    if (debugger) {
      debugger->validateCuts(cs,0,cs.sizeCuts()) ;
    }
  }
  if (verbosity_ >= 2 && (cs.sizeCuts() > 0)) {
    for (OsiCuts::iterator cut = cs.begin() ; cut != cs.end() ; cut++)
      (*cut)->print() ;
  }
# endif

  return (feasible) ;
}


/*
  Set mode
*/
void CglProbing::setMode(int mode)
{
  assert(1 <= mode && mode <= 2) ;
    mode_ = mode ;
}

// Return the mode
int CglProbing::getMode() const
{
  return mode_ ;
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



//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing ()
  : CglCutGenerator(),
    numberRows_(0),
    numberColumns_(0),
    primalTolerance_(1.0e-07),
    mode_(1),
    rowCuts_(1),
    maxPass_(3),
    logLevel_(0),
    verbosity_(0),
    paranoia_(0),
    maxProbe_(100),
    maxStack_(50),
    maxElements_(1000),
    maxPassRoot_(3),
    maxProbeRoot_(100),
    maxStackRoot_(50),
    maxElementsRoot_(10000),
    usingObjective_(0)
{

  colLower_ = NULL ;
  colUpper_ = NULL ;
  numberIntegers_ = 0 ;
  number01Integers_ = 0 ;
  numberThisTime_ = 0 ;
  totalTimesCalled_ = 0 ;
  lookedAt_ = NULL ;

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing (const CglProbing & rhs)
  : CglCutGenerator(rhs),
    numberRows_(rhs.numberRows_),
    numberColumns_(rhs.numberColumns_),
    primalTolerance_(rhs.primalTolerance_),
    mode_(rhs.mode_),
    rowCuts_(rhs.rowCuts_),
    maxPass_(rhs.maxPass_),
    logLevel_(rhs.logLevel_),
    verbosity_(rhs.verbosity_),
    paranoia_(rhs.paranoia_),
    maxProbe_(rhs.maxProbe_),
    maxStack_(rhs.maxStack_),
    maxElements_(rhs.maxElements_),
    maxPassRoot_(rhs.maxPassRoot_),
    maxProbeRoot_(rhs.maxProbeRoot_),
    maxStackRoot_(rhs.maxStackRoot_),
    maxElementsRoot_(rhs.maxElementsRoot_),
    usingObjective_(rhs.usingObjective_)
{  
  numberThisTime_ = rhs.numberThisTime_ ;
  totalTimesCalled_ = rhs.totalTimesCalled_ ;
  if (numberColumns_) {
    lookedAt_ = CoinCopyOfArray(rhs.lookedAt_,numberColumns_) ;
    if (rhs.colLower_)
      colLower_ = CoinCopyOfArray(rhs.colLower_,numberColumns_) ;
    else
      colLower_ = NULL ;
    if (rhs.colUpper_)
      colUpper_ = CoinCopyOfArray(rhs.colUpper_,numberColumns_) ;
    else
      colUpper_ = NULL ;
  } else {
    lookedAt_ = NULL ;
    colLower_ = NULL ;
    colUpper_ = NULL ;
  }
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
  delete [] colLower_ ;
  delete [] colUpper_ ;
  delete [] lookedAt_ ;
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
    delete [] colLower_ ;
    delete [] colUpper_ ;
    delete [] lookedAt_ ;

    mode_ = rhs.mode_ ;
    rowCuts_ = rhs.rowCuts_ ;
    maxPass_ = rhs.maxPass_ ;
    logLevel_ = rhs.logLevel_ ;
    verbosity_ = rhs.verbosity_ ;
    paranoia_ = rhs.paranoia_ ;
    maxProbe_ = rhs.maxProbe_ ;
    maxStack_ = rhs.maxStack_ ;
    maxElements_ = rhs.maxElements_ ;
    maxPassRoot_ = rhs.maxPassRoot_ ;
    maxProbeRoot_ = rhs.maxProbeRoot_ ;
    maxStackRoot_ = rhs.maxStackRoot_ ;
    maxElementsRoot_ = rhs.maxElementsRoot_ ;
    usingObjective_ = rhs.usingObjective_ ;
    numberThisTime_ = rhs.numberThisTime_ ;
    totalTimesCalled_ = rhs.totalTimesCalled_ ;
    if (numberColumns_)
      lookedAt_ = CoinCopyOfArray(rhs.lookedAt_,numberColumns_) ;
    else
      lookedAt_ = NULL ;
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

