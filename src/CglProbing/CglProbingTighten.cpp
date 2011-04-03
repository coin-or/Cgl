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


