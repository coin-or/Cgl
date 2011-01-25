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
  Separated tighten() and tightenClique() as the two really didn't share any
  code other than the variable declarations.  -- lh, 101124 --
*/
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
	and column copies as parameters. It should be easy enough to unpack
	them.  -- lh, 100924 --

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

void 
CglProbing::calcRowBounds(double *colLower, double * colUpper,
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

# if CGL_DEBUG > 0
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
#   if CGL_DEBUG > 1
    dump_row(i,rowLower[i],rowUpper[i],minR[i],maxR[i],0,false,false,0,
	     rowLength[i],&column[krs],&rowElements[krs],
	     primalTolerance_,colLower,colUpper) ;
#   endif
  }

# if CGL_DEBUG > 0
  std::cout << "Leaving CglProbing::calcRowBounds." << std::endl ;
# endif

}


/*
  I've chopped the portion of the method that doesn't use clique information.
  Easy enough to do, it really should have been separated from the start, as it
  shares nothing with the part that uses clique information.  -- lh, 101007 --

  And subsequently discovered that tighten(Clique) is called from exactly one
  point, so moved the conditional test for valid clique information there and
  renamed this version to tightenClique. I haven't yet gone through to prune
  useless parameters and/or variables. -- lh, 101124 --
*/

/*
  jjf: This tightens column bounds (and can declare infeasibility)
  It may also declare rows to be redundant
*/
int 
CglProbing::tightenClique(double *colLower, double * colUpper,
                    const int *column, const double *rowElements, 
                    const CoinBigIndex *rowStart, 
		    const CoinBigIndex * rowStartPos,const int * rowLength,
                    double *rowLower, double *rowUpper, 
                    int nRows,int nCols,char * intVar,int maxpass,
                    double tolerance) const
{
  int i, j, k, kre ;
  int krs ;
  int dolrows ;
  int iflagu, iflagl ;
  int ntotal = 0 ;
  int nchange = 1 ;
  int jpass = 0 ;
  double dmaxup, dmaxdown, dbound ;
  int ninfeas = 0 ;
/*
  We have non-trivial clique information from previous processing, or we
  wouldn't be here.
*/
  assert(cliqueRowStart_ && numberRows_ && cliqueRowStart_[numberRows_]) ;

  // For clique stuff
  double *cliqueMin = new double[nCols] ;
  double *cliqueMax = new double[nCols] ;
  // And second best ones
  double *cliqueMin2 = new double[nCols] ;
  double *cliqueMax2 = new double[nCols] ;
/*
  Lopped out about 200 lines of non-clique code here which is completely
  independent of the code that follows. This is the only part of CglProbing
  that uses the cliqueRow_ information constructed by
  setupRowCliqueInformation.
*/
  while (nchange) {
    int ilbred = 0 ; /* lower bounds improved */
    int iubred = 0 ; /* upper bounds improved */
    int nrwdrp = 0 ; /* redundant rows */
    if (jpass == maxpass) break ;
    jpass++ ;
    dolrows = (jpass & 1) == 1 ;
    
    for (i = 0 ; i < nRows ; ++i) {
      bool cliqueChanges = false ;
      if (rowLower[i] > -1.0e20 || rowUpper[i] < 1.0e20) {

	iflagu = 0 ;		// infinite U<i>
	iflagl = 0 ;		// infinite L<i>
	dmaxup = 0.0 ;		// finite U<i>
	dmaxdown = 0.0 ;	// finite U<i>

	krs = rowStart[i] ;
	kre = rowStart[i]+rowLength[i] ;
/*
  Recompute L<i> and U<i>. It could be that we have no clique entanglement in
  this row. In that case, calculating L<i> and U<i> is a straightforward matter
  of accumulating the row bounds based on upper and lower bounds of each
  variable.

  TODO: Why are we checking cliqueMin? If it's null, we've concluded above that
	there's no clique information, period. For that matter, how can it be
	that i >= numberRows_? This implies nRows > numberRows_. On second
	thought, this would be a case where the constraint system has grown
	since we took the snapshot.
	-- lh, 101008 --
*/
        if (!cliqueMin || i >= numberRows_ ||
	    cliqueRowStart_[i] == cliqueRowStart_[i+1]) {
          for (k = krs ; k < kre ; ++k) {
            double value = rowElements[k] ;
            j = column[k] ;
            if (value > 0.0) {
              if (colUpper[j] < 1.0e12) 
                dmaxup += colUpper[j]*value ;
	      else
                ++iflagu ;
              if (colLower[j] > -1.0e12) 
                dmaxdown += colLower[j]*value ;
	      else
                ++iflagl ;
            } else if (value < 0.0) {
              if (colUpper[j] < 1.0e12) 
                dmaxdown += colUpper[j]*value ;
	      else
                ++iflagl ;
              if (colLower[j] > -1.0e12) 
                dmaxup += colLower[j]*value ;
	      else
                ++iflagu ;
            }
          }
/*
  TODO: Yet another magic number for infinity.  -- lh, 101008 --
*/
	  if (iflagu)
	    dmaxup = 1.0e31 ;
	  if (iflagl)
	    dmaxdown = -1.0e31 ;
        }
/*
  We have clique entanglement. Try to use the clique information to get a
  stronger bound. Recall that the cliqueRow_ block for this row is indexed
  by pos'n in row. So we need to look at entry cliqueRowStart_[i]+(k-krs);
  hence the bias.
*/
	else {
          int nClique = 0 ;
          int bias = cliqueRowStart_[i]-krs ;
          double dmaxup2 = 0.0 ;
          double dmaxdown2 = 0.0 ;
          double sumZeroFixes = 0.0 ;
          for (k = krs ; k < kre ; ++k) {
            double value = rowElements[k] ;
            j = column[k] ;
	    // iClique is the locally significant clique id
            int iClique = sequenceInCliqueEntry(cliqueRow_[k+bias]) ;
            bool oneFixes = oneFixesInCliqueEntry(cliqueRow_[k+bias]) ;
/*
  Variable is not a member of a clique, or is fixed (hence not in clique by
  definition). Accumulate contribution for this variable in the normal way.

  TODO: Notice that even though we're doing the identical calculation to the
	previous code block, it's slightly different, just to confuse the
	innocent. Why do we bother to set dmax to infinity?  -- lh, 101008 --
*/
            if (iClique >= numberColumns_ || colUpper[j] == colLower[j]) {
              if (value > 0.0) {
                if (colUpper[j] >= 1.0e12) {
                  dmaxup = 1e31 ;
                  ++iflagu ;
                } else {
                  dmaxup += colUpper[j]*value ;
                }
                if (colLower[j] <= -1.0e12) {
                  dmaxdown = -1e31 ;
                  ++iflagl ;
                } else {
                  dmaxdown += colLower[j]*value ;
                }
              } else if (value < 0.0) {
                if (colUpper[j] >= 1.0e12) {
                  dmaxdown = -1e31 ;
                  ++iflagl ;
                } else {
                  dmaxdown += colUpper[j] * value ;
                }
                if (colLower[j] <= -1.0e12) {
                  dmaxup = 1e31 ;
                  ++iflagu ;
                } else {
                  dmaxup += colLower[j] * value ;
                }
              }
            }
/*
  Ok, the interesting bit --- try to use clique information. Have we
  encountered this clique id yet? If not, zero the contribution entries for
  all ids up to and including this one (lazy initialisation; do it as far as
  the max id encountered to date).

  TODO: Now that I've worked out what the calculation for a clique should be,
	I'm certain this is not it. We are correctly accumulating the sum of
	the coefficients for variables in strong-0 (sumZeroFixes), but we're
	not correctly determining max and min a<ij> for variables in strong-1
	and max and min a<ij> for variables in strong-0. The latter four are
	needed to get the correct max and min clique values. Seems like this
	was on the author's mind, but got corrupted along the way. We have the
	necessary four arrays.  -- lh, 101008 --
*/
	    else {
              if (iClique >= nClique) {
                for (int j = nClique ; j <= iClique ;j++) {
                  cliqueMin[j] = 0.0 ;
                  cliqueMax[j] = 0.0 ;
                  cliqueMin2[j] = 0.0 ;
                  cliqueMax2[j] = 0.0 ;
                }
                nClique = iClique+1 ;
              }
              //  Update best and second best
              if (oneFixes) {
                if (value > 0.0) {
                  dmaxup2 += value ;
                  cliqueMax2[iClique] = cliqueMax[iClique] ;
                  cliqueMax[iClique] = CoinMax(cliqueMax[iClique],value) ;
                } else if (value < 0.0) {
                  dmaxdown2 +=  value ;
                  cliqueMin2[iClique] = cliqueMin[iClique] ;
                  cliqueMin[iClique] = CoinMin(cliqueMin[iClique],value) ;
                }
              } else {
                sumZeroFixes += value ;
                if (value > 0.0) {
                  dmaxup2 += value ;
                  cliqueMin2[iClique] = cliqueMin[iClique] ;
                  cliqueMin[iClique] = CoinMin(cliqueMin[iClique],-value) ;
                } else if (value < 0.0) {
                  dmaxdown2 +=  value ;
                  cliqueMax2[iClique] = cliqueMax[iClique] ;
                  cliqueMax[iClique] = CoinMax(cliqueMax[iClique],-value) ;
                }
              }
            }
          }
/*
  Clearly this is a work-in-progress. clique[Min,Max]2 is never referenced.
  But this is the bones of the correct calculation. We can lump together the
  sums of strong-0 coefficients over all cliques into one variable,
  sumZeroFixes. But we can't lump the max & min calculations in an individual
  clique over strong-1 and strong-0 partitions.
*/
          double dmaxup3 = dmaxup + sumZeroFixes ;
          double dmaxdown3 = dmaxdown + sumZeroFixes ;
          for (int iClique = 0 ; iClique < nClique ; iClique++) {
            dmaxup3 += cliqueMax[iClique] ;
            dmaxdown3 += cliqueMin[iClique] ;
          }
/*
  Add in the contribution of the clique variables calculated as if they were
  completely independent.
*/
          dmaxup += dmaxup2 ;
          dmaxdown += dmaxdown2 ;
/*
  Sanity. The bound augmented with clique information should never be worse
  than the bound calculated without it.
*/
          assert (dmaxup3 <= dmaxup+1.0e-8) ;
          assert (dmaxdown3 >= dmaxdown-1.0e-8) ;
/*
  Did taking the cliques into account help us?
*/
          if (dmaxup3 < dmaxup-1.0e-8 || dmaxdown3 > dmaxdown+1.0e-8) {
            cliqueChanges = true ;
            //printf("normal min/max %g , %g clique %g , %g\n",
            //     dmaxdown,dmaxup,dmaxdown3,dmaxup3) ;
            dmaxdown = dmaxdown3 ;
            dmaxup = dmaxup3 ;
          }
        }
/*
  What do U<i> and L<i> tell us? If they're both within the bounds set in the
  constraint system, then this constraint is redundant. Mark it as such by
  setting infinite bounds.

  TODO: With yet another finite infinity. That makes around eight different
	values I've encountered so far.  -- lh, 101008 --

  Previous comment (author?; not jjf style):  The sum of the column maxs is at
  most the row ub, and the sum of the column mins is at least the row lb ;
  this row says nothing at all.  I suspect that this corresponds to an
  implied column singleton in the paper (viii, on p. 325), where the
  singleton in question is the row slack.

  TODO: WHAT PAPER? Sheesh. But I have the vague recollection of asking this
	question in the past. Maybe I should search my archives. `Implied
	column singleton' is continuous preprocessing terminology.
	-- lh, 101008 --
*/
	if (dmaxup <= rowUpper[i]+tolerance &&
	    dmaxdown >= rowLower[i]-tolerance) {
	  ++nrwdrp ;
	  rowLower[i] = -DBL_MAX ;
	  rowUpper[i] = DBL_MAX ;
	}
/*
  Check for U<i> < rowLower or L<i> > rowUpper. This is infeasibility.

  TODO: Clearly we don't yet trust the clique calculation (for good reason!)
	Which strongly implies that cliques must not be exercised. I'll need
	to check that.  -- lh, 101008 --
*/
	else {
	  if (dmaxup < rowLower[i]-tolerance ||
	      dmaxdown > rowUpper[i]+tolerance) {
	    ninfeas++ ;
            assert (!cliqueChanges) ;
	    break ;
	  }
/*
   Original comment: Finite U(i)
       below is deliberate mistake (previously was by chance)
       never do both

   Time to tighten variable bounds based on the lhs bounds.  If we have the
   ability to tighten against both bounds, work against the lower bound on
   odd-numbered passes, upper bound on even numbered passes (by pretending
   that the opposite bound has an infinity).
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
  Cliques are not present or didn't help. Tighten variable bounds without
  complications.

  TODO: We don't seem to be entertaining the possibility of putting a finite
	bound on a variable with a currently infinite bound, which is odd
	because it works perfectly well. Maybe that's just hidden by some
	earlier assignment of an arbitrary value for infinity.
	-- lh, 101008 --
*/
          if (!cliqueChanges) {
            if (iflagu == 0 && rowLower[i] > -1e15) {
              for (k = krs ; k < kre ; ++k) {
                double value = rowElements[k] ;
                j = column[k] ;
                if (value > 0.0) {
                  if (colUpper[j] < 1.0e12) {
                    dbound = colUpper[j] + (rowLower[i] - dmaxup) / value ;
                    if (dbound > colLower[j] + 1.0e-8) {
                      /* we can tighten the lower bound */
                      /* the paper mentions this as a possibility on p. 227 */
                      colLower[j] = dbound ;
                      ++ilbred ;
                      
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
                          ninfeas++ ;
                        }
                      }
                    }
                  }
                } else {
                  if (colLower[j] > -1.0e12) {
                    dbound = colLower[j] + (rowLower[i] - dmaxup) / value ;
                    if (dbound < colUpper[j] - 1.0e-8) {
                      colUpper[j] = dbound ;
                      ++iubred ;
                      if (colUpper[j] - colLower[j] <= tolerance) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                          ninfeas++ ;
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
              for (k = krs ; k < kre ; ++k) {
                double value=rowElements[k] ;
                j = column[k] ;
                if (value < 0.0) {
                  if (colUpper[j] < 1.0e12) {
                    dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value ;
                    if (dbound > colLower[j] + 1.0e-8) {
                      colLower[j] = dbound ;
                      ++ilbred ;
                      if (! (colUpper[j] - colLower[j] > tolerance)) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                          ninfeas++ ;
                        }
                      }
                    }
                  } 
                } else {
                  if (colLower[j] > -1.0e12) {
                    dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value ;
                    if (dbound < colUpper[j] - 1.0e-8) {
                      colUpper[j] = dbound ;
                      ++iubred ;
                      if (! (colUpper[j] - colLower[j] > tolerance)) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                          ninfeas++ ;
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            // with cliques
            int bias = cliqueRowStart_[i]-krs ;
            if (iflagu == 0 && rowLower[i] > -1e15) {
              for (k = krs ; k < kre ; ++k) {
                double value=rowElements[k] ;
                j = column[k] ;
                int iClique = sequenceInCliqueEntry(cliqueRow_[k+bias]) ;
                //bool oneFixes = (cliqueRow_[k+bias].oneFixes!=0) ;
                if (iClique>=numberColumns_) {
                  if (value > 0.0) {
                    if (colUpper[j] < 1.0e12) {
                      dbound = colUpper[j] + (rowLower[i] - dmaxup) / value ;
                      if (dbound > colLower[j] + 1.0e-8) {
                        /* we can tighten the lower bound */
                        /* the paper mentions this as a possibility on p. 227 */
                        colLower[j] = dbound ;
                        ++ilbred ;
                        
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
                            ninfeas++ ;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowUpper[i]>1.0e20) {
                        // can we modify coefficient
                        if (dmaxdown+value>rowLower[i]+1.0e-8) {
                          assert (dmaxdown<rowLower[i]+1.0e-8) ;
                          double change = dmaxdown+value - rowLower[i] ;
                          double newValue = value - change ;
                          if (newValue<1.0e-12)
                            newValue=0.0 ;
                          printf("Could change value from %g to %g\n",
                                 value,newValue) ;
                          // dmaxup -= change ; 
                        }
#endif
                      }
                    }
                  } else {
                    if (colLower[j] > -1.0e12) {
                      dbound = colLower[j] + (rowLower[i] - dmaxup) / value ;
                      if (dbound < colUpper[j] - 1.0e-8) {
                        colUpper[j] = dbound ;
                        ++iubred ;
                        if (colUpper[j] - colLower[j] <= tolerance) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                            ninfeas++ ;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowUpper[i]>1.0e20) {
                        // can we modify coefficient
                        if (dmaxdown-value>rowLower[i]+1.0e-8) {
                          assert (dmaxdown<rowLower[i]+1.0e-8) ;
                          double change = dmaxdown-value-rowLower[i] ;
                          double newValue = value+change ;
                          double newLower = rowLower[i]+change ;
                          if (newValue<1.0e-12)
                            newValue=0.0 ;
                          printf("Could change value from %g to %g and lorow from %g to %g\n",
                                 value,newValue,rowLower[i],newLower) ;
                         // dmaxdown += change                          
                        }
#endif
                      }
                    }
                  }
                } else if (colUpper[j]>colLower[j]) {
                  // in clique
                  // adjustment
                  double dmaxup2=dmaxup ;
                  assert (cliqueMax[iClique]>=0) ;
                  assert (cliqueMax2[iClique]>=0) ;
                  /* get max up if at other bound
                     May not go down at all but will not go up */
                  if (fabs(value)==fabs(cliqueMax[iClique]))
                    dmaxup2 -= cliqueMax[iClique]-cliqueMax2[iClique] ;
                  if (dmaxup2<rowLower[i]-1.0e-8) {
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if ( dmaxup<rowLower[i]-1.0e-8) {
                      ninfeas++ ;
                    } else {
                      if (value > 0.0) {
                        colLower[j] = 1.0 ;
                        ++ilbred ;
                      } else {
                        colUpper[j] = 0.0 ;
                        ++iubred ;
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
              for (k = krs ; k < kre ; ++k) {
                double value=rowElements[k] ;
                j = column[k] ;
                int iClique = sequenceInCliqueEntry(cliqueRow_[k+bias]) ;
                //bool oneFixes = (cliqueRow_[k+bias].oneFixes!=0) ;
                if (iClique>=numberColumns_) {
                  if (value < 0.0) {
                    if (colUpper[j] < 1.0e12) {
                      dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value ;
                      if (dbound > colLower[j] + 1.0e-8) {
                        colLower[j] = dbound ;
                        ++ilbred ;
                        if (! (colUpper[j] - colLower[j] > tolerance)) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                            ninfeas++ ;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowLower[i]<-1.0e20) {
                        // can we modify coefficient
                        if (dmaxup+value<rowUpper[i]-1.0e-8) {
                          assert (dmaxup>rowUpper[i]-1.0e-8) ;
                          double change = dmaxup+value - rowUpper[i] ;
                          double newValue = value - change ;
                          if (newValue<1.0e-12)
                            newValue=0.0 ;
                          printf("Could change value from %g to %g b\n",
                                 value,newValue) ;
                          // dmaxdown -= change ; 
                        }
#endif
                      }
                    } 
                  } else {
                    if (colLower[j] > -1.0e12) {
                      dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value ;
                      if (dbound < colUpper[j] - 1.0e-8) {
                        colUpper[j] = dbound ;
                        ++iubred ;
                        if (! (colUpper[j] - colLower[j] > tolerance)) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                            ninfeas++ ;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowLower[i]<-1.0e20) {
                        // can we modify coefficient
                        if (dmaxup-value<rowUpper[i]-1.0e-8) {
                          assert (dmaxup>rowUpper[i]-1.0e-8) ;
                          double change = dmaxup-value-rowUpper[i] ;
                          double newValue = value+change ;
                          double newUpper = rowUpper[i]+change ;
                          if (newValue<1.0e-12)
                            newValue=0.0 ;
                          printf("Could change value from %g to %g and uprow from %g to %g b\n",
                                 value,newValue,rowLower[i],newUpper) ;
                         // dmaxup += change                          
                        }
#endif
                      }
                    }
                  }
                } else if (colUpper[j]>colLower[j]) {
                  // in clique
                  // adjustment
                  double dmaxdown2=dmaxdown ;
                  assert (cliqueMin[iClique]<=0) ;
                  assert (cliqueMin2[iClique]<=0) ;
                  /* get max down if this is at other bound
                     May not go up at all but will not go down */
                  if (fabs(value)==fabs(cliqueMin[iClique]))
                    dmaxdown2 -= cliqueMin[iClique]-cliqueMin2[iClique] ;
                  if (dmaxdown2>rowUpper[i]+1.0e-8) {
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if ( dmaxdown>rowUpper[i]+1.0e-8) {
                      ninfeas++ ;
                    } else {
                      if (value < 0.0) {
                        colLower[j] = 1.0 ;
                        ++ilbred ;
                      } else {
                        colUpper[j] = 0.0 ;
                        ++iubred ;
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
    for (j = 0 ; j < nCols ; ++j) {
      if (intVar[j]) {
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  if (floor(colUpper[j]+1.0e-4)<colUpper[j]) 
	    nchange++ ;
	  // clean up anyway
	  colUpper[j]=floor(colUpper[j]+1.0e-4) ;
	  if (ceil(colLower[j]-1.0e-4)>colLower[j]) 
	    nchange++ ;
	  // clean up anyway
	  colLower[j]=ceil(colLower[j]-1.0e-4) ;
	  if (colUpper[j]<colLower[j]) {
	    /*printf("infeasible\n") ;*/
	    ninfeas++ ;
	  }
	}
      }
    }
    nchange=ilbred+iubred+nrwdrp ;
    ntotal += nchange ;
    if (ninfeas) break ;
  }
  delete [] cliqueMin ;
  delete [] cliqueMax ;
  delete [] cliqueMin2 ;
  delete [] cliqueMax2 ;
  return (ninfeas) ;
}

