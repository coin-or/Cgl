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
#include "CglProbingDebug.hpp"


const double CglProbing::CGL_REASONABLE_INTEGER_BOUND = 1.23456789e10 ;
const double CglProbing::CGL_BOGUS_1E60_BOUND = 1.0e60 ;

/* 
  jjf: This tightens column bounds (and can declare infeasibility)
  It may also declare rows to be redundant

  Returns the number of infeasible variables.

  Very brute-force. All constraints are reexamined with every pass. L(i) and
  U(i) are recomputed each time a constraint is processed.
  
  Updates colLower and colUpper. Can detect infeasibility due to constraint or
  variable bounds.

  'The paper' referred to in the original comments seems likely to be
  Andersen, Erling D. and Knud D. Andersen, "Presolving in Linear
  Programming," Mathematical Programming, Vol. 71, pp. 221–245, 1995.
*/
int CglProbing::tighten(double *const colLower, double *const colUpper,
                    const int *const column, const double *const rowElements, 
                    const CoinBigIndex *const rowStart, 
		    const CoinBigIndex *const rowStartPos,
		    const int *const rowLength,
                    double *const rowLower, double *const rowUpper, 
                    int nRows, int nCols, const char *const intVar,
		    int maxpass, double tolerance) const
{
  int dolrows ;
  int nchange = 1 ;
  int jpass = 0 ;
  int ninfeas = 0 ;

  assert (rowStartPos) ;
/*
  Main loop: propagate until we converge or until the pass limit (maxpass) is
  reached.
*/
  while (nchange && jpass < maxpass) {
    nchange = 0 ; 
    jpass++ ;
    // Appears to force an alternation between working off bi and blowi.
    dolrows = (jpass & 1) == 1 ;
/*
  Open a loop to process each row. We need at least one finite bound to work
  with.
*/
    for (int i = 0 ; i < nRows ; ++i) {
      if (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20) continue ;

      int iflagu = 0 ;
      int iflagl = 0 ;
      double dmaxup = 0.0 ;
      double dmaxdown = 0.0 ;
      double dbound = 0.0 ;
      const CoinBigIndex starti = rowStart[i] ;
      const CoinBigIndex startPosi = rowStartPos[i] ;
      const CoinBigIndex endi = rowStart[i]+rowLength[i] ;
/*
  Compute L(i) and U(i) for the row. Remember that the constraint coefficients
  are sorted, negative first, then positive.
*/
      for (CoinBigIndex jj = starti ; jj < startPosi ; ++jj) {
	const double &value = rowElements[jj] ;
	const int &j = column[jj] ;
	if (colUpper[j] < 1.0e12) 
	  dmaxdown += colUpper[j]*value ;
	else
	  ++iflagl ;
	if (colLower[j] > -1.0e12) 
	  dmaxup += colLower[j]*value ;
	else
	  ++iflagu ;
      }
      for (CoinBigIndex jj = startPosi ; jj < endi ; ++jj) {
	const double &value = rowElements[jj] ;
	const int &j = column[jj] ;
	if (colUpper[j] < 1.0e12) 
	  dmaxup += colUpper[j] * value ;
	else
	  ++iflagu ;
	if (colLower[j] > -1.0e12) 
	  dmaxdown += colLower[j] * value ;
	else
	  ++iflagl ;
      }
      if (iflagu) dmaxup = 1.0e31 ;
      if (iflagl) dmaxdown = -1.0e31 ;
/*
  If bilow < L(i) <= U(i) < bi, the constraint is redundant. Make sure we
  don't ever process it by forcing rowLower and rowUpper to infinity. 
  If U(i) < bilow or L(i) > bi, we're infeasible.
*/
      if (dmaxup < rowLower[i]-tolerance ||
	  dmaxdown > rowUpper[i]+tolerance) {
	ninfeas++ ;
	break ;
      }
      if (dmaxup <= rowUpper[i]+tolerance &&
	  dmaxdown >= rowLower[i]-tolerance) {
	++nchange ;
	rowLower[i] = -DBL_MAX ;
	rowUpper[i] = DBL_MAX ;
	continue ;
      }
/*
  See what kind of progress we can make. If the constraint has finite L(i) and
  U(i), chose one or the other.

  The reasonableness condition on rowUpper makes sense, but I don't see any
  reason why rowLower must be positive. Particularly in view of the test
  immediately below, where all we require is > -1e15.   -- lh, 110330 --
*/
      if (iflagu == 0 && rowLower[i] > 0.0 &&
          iflagl == 0 && rowUpper[i] < 1e15) {
	if (dolrows) {
	  iflagu = 1 ;
	} else {
	  iflagl = 1 ;
	}
      }
/*
  Work against blow<i>-U(i\t). Doesn't recognise an opportunity to move from
  an infinite bound to a finite bound. Check for infeasibility based on
  variable bounds.
*/
      if (iflagu == 0 && rowLower[i] > -1e15) {
	for (CoinBigIndex jj = starti ; jj < endi ; ++jj) {
	  const double &value = rowElements[jj] ;
	  const int &j = column[jj] ;
	  if (value > 0.0) {
	    if (colUpper[j] < 1.0e12) {
	      // tighten lower bound
	      dbound = colUpper[j]+(rowLower[i]-dmaxup)/value ;
	      if (dbound > colLower[j]+1.0e-8) {
		colLower[j] = dbound ;
		++nchange ;
		if (colUpper[j]-colLower[j] <= tolerance) {
		  if (colUpper[j]-colLower[j] < -100.0*tolerance) {
		    ninfeas++ ;
		  }
		}
	      }
	    }
	  } else {
	    if (colLower[j] > -1.0e12) {
	      // tighten upper bound
	      dbound = colLower[j]+(rowLower[i]-dmaxup)/value ;
	      if (dbound < colUpper[j] - 1.0e-8) {
		colUpper[j] = dbound ;
		++nchange ;
		if (colUpper[j]-colLower[j] <= tolerance) {
		  if (colUpper[j]-colLower[j] < -100.0*tolerance) {
		    ninfeas++ ;
		  }
		}
	      }
	    }
	  }
	}
      }
/*
  Work against b<i>-L(i\t). Doesn't recognise an opportunity to move from
  an infinite bound to a finite bound. Check for infeasibility based on
  variable bounds.
*/
      if (iflagl == 0 && rowUpper[i] < 1e15) {
	for (CoinBigIndex jj = starti; jj < endi; ++jj) {
	  const double &value = rowElements[jj] ;
	  const int &j = column[jj] ;
	  if (value < 0.0) {
	    if (colUpper[j] < 1.0e12) {
	      dbound = colUpper[j]+(rowUpper[i]-dmaxdown)/value ;
	      if (dbound > colLower[j]+1.0e-8) {
		colLower[j] = dbound ;
		++nchange ;
		// Interesting variant from previous stanza
		if (!(colUpper[j]-colLower[j] > tolerance)) {
		  if (colUpper[j]-colLower[j] < -100.0*tolerance) {
		    ninfeas++ ;
		  }
		}
	      }
	    } 
	  } else {
	    if (colLower[j] > -1.0e12) {
	      dbound = colLower[j]+(rowUpper[i]-dmaxdown)/value ;
	      if (dbound < colUpper[j]-1.0e-8) {
		colUpper[j] = dbound ;
		++nchange ;
		if (!(colUpper[j]-colLower[j] > tolerance)) {
		  if (colUpper[j]-colLower[j] < -100.0*tolerance) {
		    ninfeas++ ;
		  }
		}
	      }
	    }
	  }
	}
      }
    }    // end of loop to scan all rows
/*
  Scan the columns. For integer variables, force the bounds to integral
  values. Only count this as a change if we're just slightly below the integral
  value. We could well end up changing 1.9 to 1.0, which arguably is the wrong
  direction.

  Seems like we could skip this if we're already infeasible. The utility of
  the test for infeasibility here is questionable, but I can see where it
  could happen.
*/
    for (int j = 0 ; j < nCols ; ++j) {
      if (intVar[j]) {
	if (colUpper[j]-colLower[j] > 1.0e-8) {
	  if (floor(colUpper[j]+1.0e-4) < colUpper[j]) nchange++ ;
	  colUpper[j] = floor(colUpper[j]+1.0e-4) ;
	  if (ceil(colLower[j]-1.0e-4) > colLower[j]) nchange++ ;
	  colLower[j] = ceil(colLower[j]-1.0e-4) ;
	  if (colUpper[j] < colLower[j]) {
	    /*printf("infeasible\n");*/
	    ninfeas++ ;
	  }
	}
      }
    }
    if (ninfeas) break ;
  }    // end of 'propagate til done'
  return (ninfeas) ;
}


/*
  Calculate row lower and upper bounds L<i> and U<i>, respectively, and load
  them into minR and maxR. markR is loaded with codes indicating whether it's
  worth processing the row when propagating bounds.

  TODO: There doesn't seem to be any good reason to pass the guts of the row
	copy as parameters. It should be easy enough to unpack them.
	-- lh, 100924 --

  TODO: minR and maxR will have artificial values of +/- 1e60 installed as
	yet another arbitrary infinity.  -- lh, 100924 --

  TODO: Confirm the meaning of the codes in markR. Tentatively,
	-1: one or both of L<i> and U<i> small enough to process further
	-2: both L<i> and U<i> too large to process further (here, both
	    larger than 1e10).
	-- lh, 100924 --

  TODO: At the end of this method is the root of the 1e60 bug. Instead of
	marking constraints with infinite L<i> and U<i> as unsuitable, an
	arbitrary 1e60 is set and the constraint is flagged as ok for
	further processing. This is the root of the `improved' 1e60 bounds
	on variables. Look for CGL_BOGUS_1E60_BOUND.  -- lh, 100924 --

  TODO: Dusty! Uses int instead of CoinBigIndex.  -- lh, 110330 --
*/

void CglProbing::calcRowBounds (const double *const colLower,
				const double *const colUpper,
			        const int *const column,
				const double *const rowElements, 
			        const CoinBigIndex *const rowStart, 
			        const int *const rowLength,
			        const double *const rowLower,
				const double *const rowUpper, 
			        double *const minR, double *const maxR,
				int *const markR, int nRows) const
{
  int i, j, k, kre ;
  int krs ;
  int iflagu, iflagl ;
  double dmaxup, dmaxdown ;

# if CGL_DEBUG > 1
  std::cout << "Entering CglProbing::calcRowBounds." << std::endl ;
# endif
/*
  Open a loop to step through the rows, calculating L<i> and U<i> for each row
  based on variable (column) bounds l<j> and u<j>.
*/
  for (i = 0 ; i < nRows ; ++i) {
    if (rowLower[i] > -1.0e20 || rowUpper[i] < 1.0e20) {
      iflagu = 0 ;
      iflagl = 0 ;
      dmaxup = 0.0 ;
      dmaxdown = 0.0 ;
      krs = rowStart[i] ;
      kre = rowStart[i]+rowLength[i] ;
/*
   Walk the row and compute L<i> and U<i>. Arbitrarily, declare 1e12 as
   effective infinity and too large to process further.
   
   TODO: 1e12 should be a parameter or defined constant. Same for 1e60 below.
	 -- lh, 100924 --
*/
      for (k = krs ; k < kre ; ++k) {
	double value = rowElements[k] ;
	j = column[k] ;
	if (value > 0.0) {
	  if (colUpper[j] < 1.0e12) 
	    dmaxup += colUpper[j] * value ;
	  else
	    ++iflagu ;
	  if (colLower[j] > -1.0e12) 
	    dmaxdown += colLower[j] * value ;
	  else
	    ++iflagl ;
	} else if (value < 0.0) {
	  if (colUpper[j] < 1.0e12) 
	    dmaxdown += colUpper[j] * value ;
	  else
	    ++iflagl ;
	  if (colLower[j] > -1.0e12) 
	    dmaxup += colLower[j] * value ;
	  else
	    ++iflagu ;
	}
      }
      if (iflagu)
	maxR[i] = CglProbing::CGL_BOGUS_1E60_BOUND ;
      else
	maxR[i] = dmaxup ;
      if (iflagl) 
	minR[i] = -CglProbing::CGL_BOGUS_1E60_BOUND ;
      else
	minR[i] = dmaxdown ;
/*
  Mark the constraint as worth further processing (-1) or not (-2).

  TODO: This bit of code is one of the roots of the 1e60 bug. By arbitrarily
	installing 1e60 and processing all constraints, CglProbing ends up
	exporting bogus `improved' variable bounds of 1e60 which must then be
	filtered out by all clients. Either this code goes, or a filter has
	to be installed at the point where generateCuts is about to return.
	-- lh, 100924 --
*/

#     if 1
      markR[i] = -1 ;
#     else
      if (minR[i] < -1.0e10 && maxR[i] > 1.0e10) {
	markR[i] = -2 ;
      } else {
	markR[i] = -1 ;
      }
#     endif
    } else {
      minR[i] = -CglProbing::CGL_BOGUS_1E60_BOUND ;
      maxR[i] = CglProbing::CGL_BOGUS_1E60_BOUND ;
#     if 1
      markR[i] = -1 ;
#     else
      markR[i] = -2 ;
#     endif
    }
#   if CGL_DEBUG > 2
    // Just in case we executed the else and didn't set these ...
    krs = rowStart[i] ;
    kre = rowStart[i]+rowLength[i] ;
    CglProbingDebug::dump_row(i,rowLower[i],rowUpper[i],minR[i],maxR[i],
    			      0,false,false,0,
			      rowLength[i],&column[krs],&rowElements[krs],
			      primalTolerance_,colLower,colUpper) ;
#   endif
  }

# if CGL_DEBUG > 1
  std::cout << "Leaving CglProbing::calcRowBounds." << std::endl ;
# endif

}


