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
*/
int 
CglProbing::tighten(double *colLower, double * colUpper,
                    const int *column, const double *rowElements, 
                    const CoinBigIndex *rowStart, 
		    const CoinBigIndex * rowStartPos,const int * rowLength,
                    double *rowLower, double *rowUpper, 
                    int nRows,int nCols,char * intVar,int maxpass,
                    double tolerance) const
{
  int i,j,k ;
  int dolrows;
  int nchange = 1 ;
  int jpass = 0 ;
  int ninfeas = 0 ;

  // do without cliques and using sorted version
  assert (rowStartPos);
  while(nchange) {
    nchange = 0; 
    if (jpass==maxpass) break;
    jpass++;
    dolrows = (jpass & 1) == 1;
    
    for (i = 0; i < nRows; ++i) {
      if (rowLower[i]>-1.0e20||rowUpper[i]<1.0e20) {
	int iflagu = 0;
	int iflagl = 0;
	double dmaxup = 0.0;
	double dmaxdown = 0.0;
	double dbound = 0.0 ;
	int krs = rowStart[i];
	int krs2 = rowStartPos[i];
	int kre = rowStart[i]+rowLength[i];
	
	/* ------------------------------------------------------------*/
	/* Compute L(i) and U(i) */
	/* ------------------------------------------------------------*/
	for (k = krs; k < krs2; ++k) {
	  double value=rowElements[k];
	  int j = column[k];
	  if (colUpper[j] < 1.0e12) 
	    dmaxdown += colUpper[j] * value;
	  else
	    ++iflagl;
	  if (colLower[j] > -1.0e12) 
	    dmaxup += colLower[j] * value;
	  else
	    ++iflagu;
	}
	for (k = krs2; k < kre; ++k) {
	  double value=rowElements[k];
	  int j = column[k];
	  if (colUpper[j] < 1.0e12) 
	    dmaxup += colUpper[j] * value;
	  else
	    ++iflagu;
	  if (colLower[j] > -1.0e12) 
	    dmaxdown += colLower[j] * value;
	  else
	    ++iflagl;
	}
	if (iflagu)
	  dmaxup=1.0e31;
	if (iflagl)
	  dmaxdown=-1.0e31;
	if (dmaxup <= rowUpper[i] + tolerance && dmaxdown >= rowLower[i] - tolerance) {
	  /*
	   * The sum of the column maxs is at most the row ub, and
	   * the sum of the column mins is at least the row lb;
	   * this row says nothing at all.
	   * I suspect that this corresponds to
	   * an implied column singleton in the paper (viii, on p. 325),
	   * where the singleton in question is the row slack.
	   */
	  ++nchange;
	  rowLower[i]=-DBL_MAX;
	  rowUpper[i]=DBL_MAX;
	} else {
	  if (dmaxup < rowLower[i] -tolerance || dmaxdown > rowUpper[i]+tolerance) {
	    ninfeas++;
	    break;
	  }
	  /*        Finite U(i) */
	  /* -------------------------------------------------------------*/
	  /* below is deliberate mistake (previously was by chance) */
	  /*        never do both */
	  if (iflagu == 0 && rowLower[i] > 0.0 && iflagl == 0 && rowUpper[i] < 1e15) {
	    if (dolrows) {
	      iflagu = 1;
	    } else {
	      iflagl = 1;
	    }
	  }
	  if (iflagu == 0 && rowLower[i] > -1e15) {
	    for (k = krs; k < kre; ++k) {
	      double value=rowElements[k];
	      j = column[k];
	      if (value > 0.0) {
		if (colUpper[j] < 1.0e12) {
		  dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
		  if (dbound > colLower[j] + 1.0e-8) {
		    /* we can tighten the lower bound */
		    /* the paper mentions this as a possibility on p. 227 */
		    colLower[j] = dbound;
		    ++nchange;
		    
		    /* this may have fixed the variable */
		    /* I believe that this roughly corresponds to a
		     * forcing constraint in the paper (p. 226).
		     * If there is a forcing constraint (with respect
		     * to the original, untightened bounds), then in this 
		     * loop we will go through all the columns and fix
		     * each of them to their implied bound, rather than
		     * determining that the row as a whole is forced
		     * and just fixing them without doing computation for
		     * each column (as in the paper).
		     * By doing it this way, we can tighten bounds and
		     * get forcing constraints for free.
		     */
		    if (colUpper[j] - colLower[j] <= tolerance) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      } else {
		if (colLower[j] > -1.0e12) {
		  dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
		  if (dbound < colUpper[j] - 1.0e-8) {
		    colUpper[j] = dbound;
		    ++nchange;
		    if (colUpper[j] - colLower[j] <= tolerance) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  /* ----------------------------------------------------------------*/
	  /*        Finite L(i) */
	  /* ----------------------------------------------------------------*/
	  if (iflagl == 0 && rowUpper[i] < 1e15) {
	    for (k = krs; k < kre; ++k) {
	      double value=rowElements[k];
	      j = column[k];
	      if (value < 0.0) {
		if (colUpper[j] < 1.0e12) {
		  dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
		  if (dbound > colLower[j] + 1.0e-8) {
		    colLower[j] = dbound;
		    ++nchange;
		    if (! (colUpper[j] - colLower[j] > tolerance)) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		} 
	      } else {
		if (colLower[j] > -1.0e12) {
		  dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
		  if (dbound < colUpper[j] - 1.0e-8) {
		    colUpper[j] = dbound;
		    ++nchange;
		    if (! (colUpper[j] - colLower[j] > tolerance)) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    for (j = 0; j < nCols; ++j) {
      if (intVar[j]) {
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  if (floor(colUpper[j]+1.0e-4)<colUpper[j]) 
	    nchange++;
	  // clean up anyway
	  colUpper[j]=floor(colUpper[j]+1.0e-4);
	  if (ceil(colLower[j]-1.0e-4)>colLower[j]) 
	    nchange++;
	  // clean up anyway
	  colLower[j]=ceil(colLower[j]-1.0e-4);
	  if (colUpper[j]<colLower[j]) {
	    /*printf("infeasible\n");*/
	    ninfeas++;
	  }
	}
      }
    }
    if (ninfeas) break;
  }
  return (ninfeas);
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
*/

void CglProbing::calcRowBounds (double *colLower, double * colUpper,
			        const int *column, const double *rowElements, 
			        const CoinBigIndex *rowStart, 
			        const int * rowLength,
			        double *rowLower, double *rowUpper, 
			        double * minR, double * maxR, int * markR,
			        int nRows) const
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


