
/*
  Copyright (C) 2010, International Business Machines Corporation and others.
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

#include "CoinPackedMatrix.hpp"
#include "CglProbing.hpp"
#include "CglProbingRowCut.hpp"

// To enable debugging, set CGL_DEBUG in CglProbingDebug.hpp
#include "CglProbingDebug.hpp"

namespace {

/*
  Allocating one big block for data arrays is efficient but makes it
  difficult for debugging tools to detect accesses outside the boundary of an
  array. The assert checks that a double is larger than an int and that the
  ratio is a power of two.
*/
#if CGL_DEBUG == 0

const unsigned int DIratio = sizeof(double)/sizeof(int) ;

#define ONE_ARRAY

#endif

/*
  Bit constants used to specify bound information.

   0x0: no bounds tightened
   0x1: u<j> tightened
   0x2: l<j> tightened
   0x4: l<j> < -1e10 (-infty)
   0x8: u<j> >  1e10 (+infty)
*/
const unsigned int tightenUpper = 0x1 ;
const unsigned int tightenLower = 0x2 ;
const unsigned int infLower = 0x4 ;
const unsigned int infUpper = 0x8 ;

/*
  Integer and bit constants used to specify probing.

  The second down pass is used to restore the down probe state when down
  is feasible but up is not.

  The first set (downProbe, upProbe, downProbe2, oneProbeTooMany) *must*
  match the iteration number for the corresponding pass in the PROBE
  loop. The symbolic names make the code a bit easier to read.
*/
const unsigned int downProbe = 0 ;
const unsigned int upProbe = 1 ;
const unsigned int downProbe2 = 2 ;
const unsigned int oneProbeTooMany = 3 ;

const unsigned int probeDown = tightenUpper ;
const unsigned int probeUp = tightenLower ;

const unsigned int downProbeFeas = 0x1 ;
const unsigned int upProbeFeas = 0x2 ;
const unsigned int downProbe2Feas = 0x4 ;

}  // end file-local namespace

// ===============================================


/*
  Run through stackC and create disaggregation cuts. See the typeset
  documentation for the full derivation. Suppose the probe variable is x<p>
  and we're trying to generate a disaggregation cut for target x<t>.

  For a down probe, we have orig_u<p> reduced to u<p> = down = colUpper[p].
  It may be that this causes orig_u<t> to be reduced to u<t>. The upper
  disaggregation cut is
    -(orig_u<t>-u<t>)x<p> + x<t> <= u<t> - u<p>(orig_u<t>-u<t>)

  It may be that this causes orig_l<t> to be increased to l<t>. The lower
  disaggregation cut is
    -(orig_l<t>-l<t>)x<p> + x<t> >= l<t> - u<p>(orig_l<t>-l<t>)

  For an up probe, we have orig_l<p> increased to l<p> = up = colLower[p].
  It may be that this causes orig_u<t> to be reduced to u<t>. The upper
  disaggregation cut is
     (orig_u<t>-u<t>)x<p> + x<t> <= u<t> + l<p>(orig_u<t>-u<t>)

  It may be that this causes orig_l<t> to be increased to l<t>. The lower
  disaggregation cut is
     (orig_l<t>-l<t>)x<p> + x<t> >= l<t> + l<p>(orig_l<t>-l<t>)

  Note that stackC[0] = pndx, the index of the probed variable, so we
  don't want to process stackC[0] in the main loop.

  It must be true that goingToTrueBound == 2 and justReplace == false when
  this method is called.

  TODO: If you look at the math, the critical thing to get the cut form
	used here is that the old and new bounds on x<p> be separated by
	one. E.g., for a down probe, (orig_u<p>-u<p>) = 1.  It looks to
	me like goingToTrueBound might have been an attempt to guarantee this
	condition, but the test is wrong, hence it's patched up with tweak
	that essentially cuts this cut back to binary probing variables.
	-- lh, 101214 --

  probeDir should be coded using probeDown, probeUp.
*/
void disaggCuts (int nstackC, int probeDir,
		 double primalTolerance_, double disaggEffectiveness,
		 const OsiSolverInterface &si,
		 CglProbingRowCut &rowCut, const int *const stackC,
		 const double *const colsol,
		 const double *const colUpper, const double *const colLower,
		 const double *const saveU, const double *const saveL,
		 int *const index, double *const element)
{ 
  int pndx = stackC[0] ;

# if CGL_DEBUG > 0
  std::cout
    << "Entering disaggCuts, probed variable "
    << si.getColName(pndx) << "(" << pndx << ")." << std::endl ;
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  // down or up are the only possibilities
  assert((probeDir&(~(probeDown|probeUp))) == 0) ;
# endif

  int plusMinus = 0 ;
  double lu_p = 0.0 ;
  double origsol_p = colsol[pndx] ;
  double origdelta_p = 0.0 ;
/*
  Set up the quantities that vary depending on whether this is an up or down
  probe:
  * lu_p is the new upper (down probe) or lower (up probe) bound on x<p>
  * plusMinus will be -1 (down probe) or +1 (up probe).
  * origdelta_p will be x*<p>-u<p> (down probe) or l<p>-x*<p> (up probe)
*/
  if (probeDir == probeDown) {
    lu_p = colUpper[pndx] ;
    plusMinus = -1 ;
    origdelta_p = origsol_p-lu_p ;
  } else {
    lu_p = colLower[pndx] ;
    plusMinus = 1 ;
    origdelta_p = lu_p-origsol_p ;
  }

  for (int istackC = nstackC-1 ; istackC > 0 ; istackC--) {
    int icol = stackC[istackC] ;
    double origu_t = saveU[istackC] ;
    double origl_t = saveL[istackC] ;
    double u_t = colUpper[icol] ;
    double sol_t = colsol[icol] ;
    double l_t = colLower[icol] ;
    double boundChange = origu_t-u_t ;
    double coeff = plusMinus*boundChange ;
/*
  Generate the upper disaggregation cut, if it's violated and the coefficients
  will be reasonable. The effectiveness is the improvement (change) in
  x<pndx>.
*/
    if (boundChange > 0.0 && origu_t < 1.0e10 &&
	(sol_t > u_t+boundChange*origdelta_p+primalTolerance_)) {
      OsiRowCut rc ;
      rc.setLb(-DBL_MAX) ;
      rc.setUb(u_t+coeff*lu_p) ;
      index[0] = icol ;
      element[0] = 1.0 ;
      index[1] = pndx ;
      element[1] = coeff ;
/*
  Effectiveness is how far x<p> moves. We can save effort by transforming
  sol_p = lu_p + ((u_t-sol_t)/coeff) to delta_p = ((u_t-sol_t)/coeff).
*/
      double delta_p = ((u_t-sol_t)/coeff) ;
      assert(delta_p > origdelta_p) ;
      rc.setEffectiveness(delta_p-origdelta_p) ;
      if (rc.effectiveness() > disaggEffectiveness) {
	rc.setRow(2,index,element,false) ;
#	if CGL_DEBUG > 0
	if (debugger) assert(!debugger->invalidCut(rc)); 
#	endif
	rowCut.addCutIfNotDuplicate(rc) ;
      }
    }
/*
  Generate the upper disaggregation cut.
*/
    boundChange = origl_t-l_t ;
    coeff = plusMinus*boundChange ;
    if (boundChange < 0.0 && origl_t > -1.0e10 &&
	(sol_t < l_t+boundChange*origdelta_p-primalTolerance_)) {
      OsiRowCut rc ;
      rc.setLb(l_t+coeff*lu_p) ;
      rc.setUb(DBL_MAX) ;
      index[0] = icol ;
      element[0] = 1.0 ;
      index[1] = pndx ;
      element[1] = coeff ;
      double delta_p = ((sol_t-l_t)/coeff) ;
      assert(delta_p > origdelta_p) ;
      rc.setEffectiveness(delta_p-origdelta_p) ;
      if (rc.effectiveness() > disaggEffectiveness) {
	rc.setRow(2,index,element,false) ;
#	if CGL_DEBUG > 0
	if (debugger) assert(!debugger->invalidCut(rc)); 
#	endif
	rowCut.addCutIfNotDuplicate(rc) ;
#if 0
	printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,origl_t,origu_t,colLower[icol],colsol[icol], pndx, colsol[pndx]) ;
	printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
	       icol,boundChange,pndx,rc.ub()) ;
#endif
      }
    }
  }

# if CGL_DEBUG > 0
  std::cout
    << "Leaving disaggCuts, "
    << rowCut.numberCuts() << " cuts." << std::endl ;
# endif
  return ;
}




/*
  Walk the column for this variable and propagate a bound change on x<j> to
  the row bounds for rows referencing x<j>. Recall that the row bounds are

    L<i> = SUM{P}(l<j>) + SUM{M}(u<j>)
    U<i> = SUM{P}(u<j>) + SUM{M}(l<j>)

  for negative coefficients M and positive coefficients P.
  The update cases are:

    column bound    a<ij>    row bound
	l<j>         >0        L<i>
	l<j>         <0        U<i>
	u<j>         >0        U<i>
	u<j>         <0        L<i>
  
  Movement is used to determine the column bound that's being tightened. It
  should be > 0 to tighten l<j>, movement < 0 to tighten u<j>.  Once we've
  selected the row bound to update, the update is always += value*movement.
*/

bool updateRowBounds (int j, double movement,
		      const int *const columnStart,
		      const int *const columnLength,
		      const int *const row,
		      const double *const columnElements, 
		      const double *const rowUpper,
		      const double *const rowLower,
		      int &nstackR, int *const stackR, int *const markR,
		      double *const minR, double *const maxR,
		      double *const saveMin, double *const saveMax)
{
  bool feasible = true ;

  for (int k = columnStart[j] ; k < columnStart[j]+columnLength[j] ; k++) {
    int irow = row[k] ;
    double value = columnElements[k] ;
    double delta = value*movement ;
/*
  markR is initialised by tighten2. A code of -1 indicates this row should be
  processed; -2 says ignore. Other codes should not happen.

  TODO: This assert, and the disabled continue, go back to the change in
	tighten2 where rows with infinite bounds are handed arbitrary
	bounds of 1e60 and labelled suitable for further processing. The
	assert should be removed and the continue clause reinstated.
	-- lh, 101127 --

  TODO: Replace the original assert (!= -2) with a more effective one.
	But don't remove the rest of the code structure. I still think I'll
	end up reintroducing the `ignore' code.  -- lh, 101203 --
*/
#   if CGL_DEBUG > 0
    if (markR[irow] < -1 || markR[irow]%10000 >= nstackR) {
      std::cout
	<< "Row " << irow << " marked " << markR[irow]
	<< "; expected between -1 and " << nstackR << "." << std::endl ;
      
      dump_row(irow,rowLower[irow],rowUpper[irow],minR[irow],maxR[irow],
	       0,false,false,0,0,0,0,0,0,0) ;
      std::cout << std::endl ;
    }
    assert (markR[irow] >= -1 && markR[irow]%10000 < nstackR) ;
#   endif
    if (markR[irow] == -1) {
      stackR[nstackR] = irow ;
      markR[irow] = nstackR ;
      saveMin[nstackR] = minR[irow] ;
      saveMax[nstackR] = maxR[irow] ;
      nstackR++ ;
#   if 0
    } else if (markR[irow] == -2) {
      continue ;
#   endif
    }
/*
  LOU_DEBUG: clear count as row bound will change.
*/
    markR[irow] = markR[irow]%10000 ;
/*
  Case analysis to tighten the appropriate bound.
*/
    if (movement > 0.0) {
      if (value > 0.0) {
	if (minR[irow] > -1.0e10)
	  minR[irow] += delta ;
	if (minR[irow] > rowUpper[irow]+1.0e-5) {
	  feasible = false ;
	  break ;
	}
      } else {
	if (maxR[irow] < 1.0e10)
	  maxR[irow] += delta ;
	if (maxR[irow] < rowLower[irow]-1.0e-5) {
	  feasible = false ;
	  break ;
	}
      }
    } else {
      if (value > 0.0) {
	if (maxR[irow] < 1.0e10)
	  maxR[irow] += delta ;
	if (maxR[irow] < rowLower[irow]-1.0e-5) {
	  feasible = false ;
	  break ;
	}
      } else {
	if (minR[irow] > -1.0e10)
	  minR[irow] += delta ;
	if (minR[irow] > rowUpper[irow]+1.0e-5) {
	  feasible = false ;
	  break ;
	}
      }
    }
  }
  return (feasible) ;
}



/*
  jjf: clean up djs and solution

  Calculate domains (u<j>-l<j>) for variables, and potential swings in rows. Do
  some cleanup of reduced costs.

  If CGL_DEBUG > 0, we'll also do some matrix consistency checks.
*/
void groomSoln (double direction, double primalTolerance, double *const djs,
		const double *const colLower, double *const colsol,
		const double *const colUpper, double *const columnGap,
		const CoinPackedMatrix *const rowCopy,
		const int *const rowStartPos,
		double *const largestNegativeInRow,
		double *const largestPositiveInRow)
{
  int nRows = rowCopy->getNumRows() ;
  int nCols = rowCopy->getNumCols() ;
  const int *rowStart = rowCopy->getVectorStarts() ;
  const int *column = rowCopy->getIndices() ;
  const double *rowElements = rowCopy->getElements() ;

# ifndef NDEBUG
  // Note that rowLength is only used in the assert
  const int *rowLength = rowCopy->getVectorLengths() ;
# endif

/*
  Find the largest potential postive and negative swing in a row, where
  swing is defined as a<ij>(u<j>-l<j>). Check correctness of rowStartPos if
  we're debugging.
*/
  for (int i = 0 ; i < nRows ; i++) {

    assert (rowStart[i]+rowLength[i] == rowStart[i+1]) ;

    int kk ;
    double value ;

#if CGL_DEBUG > 0
    // Check correctness of rowStartPos
    for (kk = rowStart[i] ; kk < rowStart[i+1] ; kk++) {
      value = rowElements[kk] ;
      if (value > 0.0) 
	break ;
    }
    assert (rowStartPos[i] == kk) ;
#endif

    value = 0.0 ;
    for (kk = rowStart[i] ; kk < rowStartPos[i] ; kk++) {
      int iColumn = column[kk] ;
      double gap = CoinMin(1.0e100,colUpper[iColumn]-colLower[iColumn]) ;
      value = CoinMin(value,gap*rowElements[kk]) ;
    }
    largestNegativeInRow[i] = value ;
    value = 0.0 ;
    for ( ; kk < rowStart[i+1] ; kk++) {
      int iColumn = column[kk] ;
      double gap = CoinMin(1.0e100,colUpper[iColumn]-colLower[iColumn]) ;
      value = CoinMax(value,gap*rowElements[kk]) ;
    }
    largestPositiveInRow[i] = value ;
  }
/*
  Convert reduced costs to minimisation and clean up any that are on the wrong
  side of zero.

  TODO: Seems like we could move this ahead of the previous loops and then
	use columnGap[i] to calculate largest[Negative,Positive]InRow.
	-- lh, 101125 --
*/
  for (int i = 0 ; i < nCols ; ++i) {
    double djValue = djs[i]*direction ;
    double gap = colUpper[i]-colLower[i] ;
    if (gap > 1.0e-8) {
      if (colsol[i] < colLower[i]+primalTolerance) {
        colsol[i] = colLower[i] ;
        djs[i] = CoinMax(0.0,djValue) ;
      } else if (colsol[i] > colUpper[i]-primalTolerance) {
        colsol[i] = colUpper[i] ;
        djs[i] = CoinMin(0.0,djValue) ;
      } else {
        djs[i] = 0.0 ;
      }
    }
    columnGap[i] = gap-primalTolerance ;
  }

  return ;
}

/*
  Given that the variable being probed has been discovered to be monotone,
  this method will save the resulting bound changes for integer variables as
  column cuts. 

  One could do the same for continuous variables, but at present we don't.
  Presumably we need to test against origColLower, origColUpper (the
  bounds held in the si) to decide if we have real cuts on the original
  problem. One could argue that we could test against saveL and saveU and
  save the trouble of repeatedly recording the same column cut.

  Index and elements are scratch arrays.

  Returns true if any of the bound improvements are actually cuts, false
  otherwise.
*/
bool monotoneActions (double primalTolerance_,
		      const OsiSolverInterface &si,
		      OsiCuts &cs,
		      int nstackC, const int *const stackC,
		      const char *const intVar,
		      const double *const colLower,
		      const double *const colsol,
		      const double *const colUpper,
		      int *const index, double *const element)
{
  OsiColCut cc ;
  bool anyColumnCuts = false ;
  int nTot = 0 ;
  int nFix = 0 ;
  const double *origColLower = si.getColLower() ;
  const double *origColUpper = si.getColUpper() ;
# if CGL_DEBUG > 0
  int probedVar = stackC[0] ;
  std::cout
    << "    " << si.getColName(probedVar) << " (" << probedVar
    << ") monotone [" << colLower[probedVar] << ", "
    << colUpper[probedVar] << "] from " << colsol[probedVar]
    << "." << std::endl ;
# endif
/*
  Stash the improved column bounds, u<j> first, then l<j>. If any of them
  actually cut off the current solution, indicate that we have cuts.
*/
  bool ifCut = false ;
  for (int istackC = 0 ; istackC < nstackC ; istackC++) {
    int icol = stackC[istackC] ;
    if (intVar[icol]) {
      if (colUpper[icol] < origColUpper[icol]-1.0e-4) {
	element[nFix] = colUpper[icol] ;
	index[nFix++] = icol ;
	if (colsol[icol] > colUpper[icol]+primalTolerance_) {
	  ifCut = true ;
	  anyColumnCuts = true ;
	}
      }
    }
  }
  if (nFix) {
    nTot = nFix ;
    cc.setUbs(nFix,index,element) ;
    nFix = 0 ;
  }
  for (int istackC = 0 ; istackC < nstackC ; istackC++) {
    int icol = stackC[istackC] ;
    if (intVar[icol]) {
      if (colLower[icol] > origColLower[icol]+1.0e-4) {
	element[nFix] = colLower[icol] ;
	index[nFix++] = icol ;
	if (colsol[icol] < colLower[icol]-primalTolerance_) {
	  ifCut = true ;
	  anyColumnCuts = true ;
	}
      }
    }
  }
  if (nFix) {
    nTot += nFix ;
    cc.setLbs(nFix,index,element) ;
  }
  if (nTot) {
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

  return (anyColumnCuts) ;
}

void clearStacks (double primalTolerance_,
		  int nstackC, int *const stackC, int *const markC,
		  const double *const colLower, const double *const colUpper,
		  int nstackR, int *const stackR, int *const markR)
{
/*
  Clear off the bound change markers, except for fixed variables. We'll be
  moving on to another variable.
*/
  for (int istackC = 0 ; istackC < nstackC ; istackC++) {
    int icol = stackC[istackC] ;
    if (colUpper[icol]-colLower[icol] > primalTolerance_) {
      markC[icol] &= ~(tightenLower|tightenUpper) ;
    } else {
      markC[icol] = tightenLower|tightenUpper ;
    }
  }
  for (int istackR = 0 ; istackR < nstackR ; istackR++) {
    int irow = stackR[istackR] ;
    markR[irow] = -1 ;
  }
}


// =========================================================



/*
  jjf: Does probing and adding cuts

  Note that this method is heavily commented and instrumented in CbcAnn.

  It would appear that probe() has received more attention that probeClique or
  probeSlack. Neither of them has the ONE_ARRAY option.
*/

/*
  We're going to probe integer variables, and we're interested in propagating
  the effect of each probe (force a variable up or down).

  The bulk of the method consists of three nested loops:
    * PASS LOOP: makes multiple passes, probing all variables in lookedAt_
    * LOOKEDAT LOOP: iterates over the variables in lookedAt_
    * PROBE LOOP: probes down/up/down for each variable.

  The body of the probe loop runs a bit over 2600 lines and propagates
  the effect of forcing the probed variable down or up.

    * The start of the body updates the row bounds of affected rows, then
      initialises stackC with the probed variable.

    * The middle is a 1,000 line loop that does the bulk of propagation.
      At the start there's some dead code (PROBING100) associated with
      CglTreeProbingInfo. The purpose remains to be determined. This is
      followed by six specialised replications of the same functionality:
      walk the column of a variable from stackC, and for each affected
      row, try to tighten the bounds of other variables in the row. If
      we successfully tighten bounds on a variable, walk the column of
      that variable, tighten the bounds of affected rows, and place the
      variable on stackC.

    * The end of the body deals with the result of the probe. At the end
      of each iteration, there's a decision: have we proven infeasibility
      (INFEASIBILITY) or monotonicity (MONOTONE), or do we need another
      iteration (ITERATE). Monotonicity and iteration are large blocks
      in themselves.

  There was a fair bit of preamble and postamble around the PASS loop.
  The preamble code is protected by PROBING_EXTRA_STUFF and has been disabled
  since 2006. It appears to find disaggregation cuts.  The postamble code
  purports to find bigM constraints. It's clearly a work in progress. I've
  moved both out to tmp files.  -- lh, 101203 --

  Comments on the data structures:

  markC	 0: not fixed
	 1: tightened upper bound
	 2: tightened lower bound
	 3: fixed

    or perhaps the above should be interpreted as bit codes, along with the
    bits used to indicate infinite bounds:

       0x0: no bounds tightened
       0x1: u<j> tightened
       0x2: l<j> tightened
       0x4: l<j> < -1e10 (-infty)
       0x8: u<j> >  1e10 (+infty)

    Note that we'll only tighten a bound once --- this is enforced during
    bounds propagation.

  stackC      Records columns to be processed. This record is started anew
	      each time we probe, pushing a variable in a given direction.
	      The variable being probed is always entry 0 on stackC.
	      saveL and saveU are correlated with stackC.

	      This isn't a stack --- it's a record of every variable
	      processed. nstackC is the tail of the queue, and istackC is
	      the head. At the end of a probing pass it's scanned to recover
	      the actual implications that were generated. It's not clear
	      to me why it's allocated at 2*ncols. Perhaps, at one time,
	      there were separate entries for upper and lower bound changes?

	      No, just as I'm thinking I'm almost done, I notice there's some
	      fancy footwork in the portion of the code that handles the case
	      where both the up and down probes were feasible. Seems likely
	      that we're allowing for each direction to process all variables
	      at most once.

	      And as I'm working through the code that handles the case
	      where down and up probes were both feasible, there in the
	      move singleton code is an assert that nstackC < nCols.

	      And by the time I finally reach the end, I'm pretty well
	      convinced that it's an error if we exceed nCols on stackC.
	      The only thing to worry about is whether singleton processing
	      could add a variable that's already been added due to
	      propagation, and a check of the code says the answer is 'no'.

       TODO:  When I rewrite this, allocate stackC as nCols and put in some
	      asserts to catch any violations. There are already asserts
	      in place at the point where a variable is placed on stackC
	      during propagation. (Look for `Is there any change to
	      propagate'.)   -- lh, 101128 --


  stackR      Stack rows to be processed. This stack is started anew each
	      time we probe pushing a variable in a given direction.

  minR, maxR, markR  Initialised externally by tighten2. For row i, minR[i] =
		     LB(i), maxR[i] = UB(i) (row lhs lower and upper bounds,
		     respectively).  markR is set to -1 if there's at least
		     one finite lhs bound, -2 if we shouldn't process this
		     row (no finite bounds, or too many coefficients).

  Then, as we work, markR[i] will be set to point to the entry in stackR for
  row i. saveMin and saveMax are correlated with stackR, and hold the
  original values for LB(i) and UB(i) (from minR and maxR) while we're
  working. As it turns out, however, the the back pointer from markR is never
  used.

  largestPositiveInRow	for row i, max{j in P} a<ij>(u<i>-l<j>), where P is
			the set of positive coefficients a<ij>
  largestNegativeInRow	for row i, min{j in M} a<ij>(u<i>-l<j>), where M is
			the set of negative coefficients a<ij>

  element and index are scratch arrays used to construct cuts.

  realRows is the translation array created if we modified the constraint
  system during preprocessing in gutsOfGenerateCuts.

*/
/*
  The #ifdef around the formal parameter apparently satisfies the method
  signature while clearly indicating the parameter is unused in the body of
  the method.
*/

int CglProbing::probe( const OsiSolverInterface & si, 
		       OsiCuts & cs, 
		       double * colLower, double * colUpper, 
		       CoinPackedMatrix *rowCopy,
		       CoinPackedMatrix *columnCopy,
		       const CoinBigIndex * rowStartPos,const int * realRows, 
		       const double * rowLower, const double * rowUpper,
		       const char * intVar, double * minR, double * maxR, 
		       int * markR, 
                       CglTreeInfo * info) const

{
# if CGL_DEBUG > 0
  std::cout << "Entering CglProbing::probe." << std::endl ;

  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
# endif
/*
  PREPARATION

  All the code down through `PASS LOOP: HEAD' is preparation. Do all of the
  setup that will not change over the nested loops that do the work of probing.
  Note that rowCopy may have considerably fewer rows than the original system.
*/

  int nRows = rowCopy->getNumRows() ;
  int nCols = rowCopy->getNumCols() ;
  int maxStack = info->inTree ? maxStack_ : maxStackRoot_ ;

/*
  Fiddle maxStack. Odd sort of function:
     1 < maxStack <= 25 => 2*maxStack
    26 < maxStack <= 50 => 50
    51 < maxStack       => maxStack

  TODO: Grepping through the code says there's no way that totalTimesCalled_
	can become negative. Perhaps in some derived class? More likely a
	magic number for something, or else a quick way to disable this
	code.  -- lh, 101021 --
*/
  if ((totalTimesCalled_%10) == -1) {
    int newMax = CoinMin(2*maxStack,50) ;
    maxStack = CoinMax(newMax,maxStack) ;
  }
  totalTimesCalled_++ ;

# ifdef ONE_ARRAY
  assert((DIratio != 0) && ((DIratio&(DIratio-1)) == 0)) ;
  int nSpace = 8*nCols+4*nRows+2*maxStack ;
  nSpace += (4*nCols+nRows+maxStack+DIratio-1)>>(DIratio-1) ;
  double * colsol = new double[nSpace] ;
  double * djs = colsol + nCols ;
  double * columnGap = djs + nCols ;
  double * saveL = columnGap + nCols ;
  double * saveU = saveL + 2*nCols ;
  double * saveMin = saveU + 2*nCols ;
  double * saveMax = saveMin + nRows ;
  double * largestPositiveInRow = saveMax + nRows ;
  double * largestNegativeInRow = largestPositiveInRow + nRows ;
  double * element = largestNegativeInRow + nRows ;
  double * lo0 = element + nCols ;
  double * up0 = lo0 + maxStack ;
  int * markC = reinterpret_cast<int *> (up0+maxStack) ;
  int * stackC = markC + nCols ;
  int * stackR = stackC + 2*nCols ;
  int * index = stackR + nRows ;
  int * stackC0 = index + nCols ;
# else 
  double * colsol = new double[nCols] ;
  double * djs = new double[nCols] ;
  double * columnGap = new double [nCols] ;
  double * saveL = new double [2*nCols] ;
  double * saveU = new double [2*nCols] ;
  double * saveMin = new double [nRows] ;
  double * saveMax = new double [nRows] ;
  double * largestPositiveInRow = new double [nRows] ;
  double * largestNegativeInRow = new double [nRows] ;
  double * element = new double[nCols] ;
  double * lo0 = new double[maxStack] ;
  double * up0 = new double[maxStack] ;
  int * markC = new int [nCols] ;
  int * stackC = new int [2*nCols] ;
  int * stackR = new int [nRows] ;
  int * index = new int[nCols] ;
  int * stackC0 = new int[maxStack] ;
# endif

/*
  Create a local container to hold the cuts we generate. At the end we'll
  transfer these to the container passed as a parameter.

  Typically, rowCopy (nRows) will be smaller than the original system in the
  si, but it could be one larger (all original constraints plus the objective).
  At the root, allow lots of room; even more if this is the first call at the
  root. In the tree, much less.

  If all we're doing is strengthening rows in place (option 0x40 indicates
  we're in preprocessing, a good place to do this), we can only work with
  what we have in rowCopy (and we need a translation array).

  TODO: There's a problem with realRows here --- if we didn't delete any
	rows during preprocessing in gutsOfGenerateCuts, realRows is not
	allocated. So we can't strengthen in place in that case. Not likely
	the intent.  -- lh, 101209 --
*/
  bool justReplace = ((info->options&0x40) != 0) && (realRows != NULL) ;
  int nRowsFake = 0 ;
  if (justReplace) {
    nRowsFake = nRows ;
  } else {
    int nRowsSafe = CoinMin(nRows,si.getNumRows()) ;
    nRowsFake = info->inTree ? nRowsSafe/3 : nRowsSafe*10 ;
    if (!info->inTree && info->pass == 0) nRowsFake *= 10 ;
  }
  CglProbingRowCut rowCut(nRowsFake,!info->inTree) ;


  // Unpack matrices
  const int * column = rowCopy->getIndices() ;
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts() ;
  const double * rowElements = rowCopy->getElements() ;

  const int * row = columnCopy->getIndices() ;
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts() ;
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements() ;

/*
  This won't be used until we're sorting out the results of propagating a probe
  (look for the next occurrence of columnLength2). We need to make sure the
  column we're looking at is a singleton in the original system.
*/
#define MOVE_SINGLETONS
#ifdef MOVE_SINGLETONS
  const double * objective = si.getObjCoefficients() ;
  const int * columnLength2 = si.getMatrixByCol()->getVectorLengths(); 
#endif
/*
  Grab the column solution and reduced costs and groom them.
*/
  CoinMemcpyN(si.getReducedCost(),nCols,djs) ;
  CoinMemcpyN(si.getColSolution(),nCols,colsol) ;
  double direction = si.getObjSense() ;
  groomSoln(direction,primalTolerance_,djs,colLower,colsol,colUpper,columnGap,
	    rowCopy,rowStartPos,largestNegativeInRow,largestPositiveInRow) ;
/*
  Determine objective cutoff (minimisation convention).

  TODO: Seems to me that this bit of code is a guaranteed fail for a
	maximisation problem with no cutoff. It also repeats in a number
	of places.  -- lh, 101125 --
*/
  double cutoff ;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff) ;
  if (!cutoff_available || usingObjective_ < 0) {
    cutoff = si.getInfinity() ;
  }
  cutoff *= direction ;
  if (fabs(cutoff) > 1.0e30)
    assert (cutoff > 1.0e30) ;
  double current = si.getObjValue() ;
  current *= direction ;
/*
  Scan the variables, noting the ones that are fixed or have a bound too large
  to be useful.
*/
  //int nFix = 0 ;
  for (int i = 0 ; i < nCols ; i++) {
    if (colUpper[i]-colLower[i] < 1.0e-8) {
      markC[i] = tightenLower|tightenUpper ;
      //nFix++ ;
    } else {
      markC[i] = 0 ;
      if (colUpper[i] > 1.0e10)
	markC[i] |= infUpper ;
      if (colLower[i] < -1.0e10)
	markC[i] |= infLower ;
    }
  }
  //printf("PROBE %d fixed out of %d\n",nFix,nCols) ;

/*
  jjf: If we are going to replace coefficient then we don't need to be
       effective

  Seems like this comment does not agree with CglProbingRowCut.addCuts(),
  which checks effectiveness before adding a cut to the OsiCuts collection or
  entering it into the strenghtenRow array.

  But it does agree with the way coefficient strengthening cuts are handled
  down in the end-of-iteration processing for the down/up/down probe loop.

  It looks like the idea of strengthenRow is that a cut that's simply
  strengthening an existing row i is entered in slot i of strengthenRow. It's
  left to the client to deal with it on return. I need to get a better idea of
  how justReplace is handled, here and in the client.  -- lh, 101125 --
*/
  //double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3 ;
  double needEffectiveness = info->strengthenRow ? 1.0e-8 : 1.0e-3 ;
  if (justReplace && ((info->pass&1) != 0))
    needEffectiveness = -1.0e10 ;

  double tolerance = 1.0e1*primalTolerance_ ;

  /* for both way coding */
  int nstackC0 = -1 ;
  int nstackR,nstackC ;
/*
  So, let's assume that the real business starts here, and the previous block
  is irrelevant.
*/
/*
  All the initialisation that occurs here is specific to CglTreeProbingInfo.
  Sort of begs the question `Why is this handled as a derived class of
  CglTreeInfo?' And why (what?) is CglTreeInfo, in the grand scheme of things?

  Some grepping about in Cbc says that we might see a CglTreeProbingInfo
  object during processing of the root node. My tentative hypothesis is that
  it's a vehicle to pass implications back to cbc, where they will be stashed
  in a CglImplication object for further use. It's worth mentioning that the
  existence of a CglTreeProbingInfo object depends on an independent symbol
  (CLIQUE_ANALYSIS) being defined in CbcModel::branchAndBound. It's also worth
  mentioning that the code here will core dump if the info object is not a
  CglTreeProbingInfo object.

  And finally, it's worth mentioning that this code is disabled --- PROBING100
  is defined to 0 at the head of CglProbing --- since 080722.

  The trivial implementation CglTreeInfo::initializeFixing() always returns 0,
  so in effect this block ensures saveFixingInfo is false.

  -- lh, 101126 --
*/
  bool saveFixingInfo = false ;
#if PROBING100
  CglTreeProbingInfo *probingInfo = dynamic_cast<CglTreeProbingInfo *> (info) ;
  const int *backward = NULL ;
  const int *integerVariable = NULL ;
  const int *toZero = NULL ;
  const int *toOne = NULL ;
  const fixEntry *fixEntries = NULL ;
#endif
  if (info->inTree) {
#if PROBING100
    backward = probingInfo->backward() ;
    integerVariable = probingInfo->integerVariable() ;
    toZero = probingInfo->toZero() ;
    toOne = probingInfo->toOne() ;
    fixEntries = probingInfo->fixEntries() ;
#endif
  } else {
    saveFixingInfo = (info->initializeFixing(&si) > 0) ;
  }
/*
  PASS LOOP: HEAD

  Major block #2: Main pass loop.
*/
  bool anyColumnCuts = false ;
  int ninfeas = 0 ;
  int rowCuts = rowCuts_ ;
  double disaggEffectiveness = 1.0e-3 ;
  int ipass = 0 ;
  int nfixed = -1 ;
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_ ;

  while (ipass < maxPass && nfixed) {
    int iLook ;
    ipass++ ;
    //printf("pass %d\n",ipass) ;
    nfixed = 0 ;
/*
  We went through a fair bit of trouble in gutsOfGenerateCut to determine
  the set of variables to be probed and loaded the indices into lookedAt_;
  numberThisTime_ reflects that.

  The net effect here is that the root gets special treatment
  (maxProbeRoot_) and the first pass at the root gets extra special treatment
  (numberThisTime_).

  See CbcCutGenerator::generateCuts for the origin of magic number 123. 
  Comments there indicate it's intended to take effect deep in the tree, but
  the code here will also work at the root.

  Special cases aside, the net effect of the loops is to walk lookedAt_,
  promoting every nth variable (n = cutDown) to the front of lookedAt_. When
  the loop finishes, the rest of the variables (which were copied off to
  stackC) are appended to lookedAt_. If XXXXXX is not defined, we have an
  expensive noop.  -- lh, 101126 --
*/
    int justFix = (!info->inTree && !info->pass) ? -1 : 0 ;
    int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_ ;
    if (justFix < 0)
      maxProbe = numberThisTime_ ;
    if (maxProbe == 123) {
      maxProbe = 0 ;
      if (!info->inTree) {
	if (!info->pass || numberThisTime_ < -100) {
	  maxProbe = numberThisTime_ ;
	} else {
	  int cutDown = 4 ;
	  int offset = info->pass%cutDown ;
	  int i ;
	  int k = 0 ;
	  int kk = offset ;
	  for (i = 0 ; i < numberThisTime_ ; i++) {
	    if (!kk) {
#define XXXXXX
#ifdef XXXXXX
	      lookedAt_[maxProbe] = lookedAt_[i] ;
#endif
	      maxProbe++ ;
	      kk = cutDown-1 ;
	    } else {
	      stackC[k++] = lookedAt_[i] ;
	      kk-- ;
	    }
	  }
#ifdef XXXXXX
	  memcpy(lookedAt_+maxProbe,stackC,k*sizeof(int)) ;
#endif
	}
      } else {
	if (numberThisTime_ < 200) {
	  maxProbe = numberThisTime_ ;
	} else {
	  int cutDown = CoinMax(numberThisTime_/100,4) ;
	  int offset = info->pass%cutDown ;
	  int i ;
	  int k = 0 ;
	  int kk = offset ;
	  for (i = 0 ; i < numberThisTime_ ; i++) {
	    if (!kk) {
#ifdef XXXXXX
	      lookedAt_[maxProbe] = lookedAt_[i] ;
#endif
	      maxProbe++ ;
	      kk = cutDown-1 ;
	    } else {
	      stackC[k++] = lookedAt_[i] ;
	      kk-- ;
	    }
	  }
#ifdef XXXXXX
	  memcpy(lookedAt_+maxProbe,stackC,k*sizeof(int)) ;
#endif
	}
      }
    }  // end maxProbe == 123

/*
  This looks to be an overall limit on probing. It's decremented every time a
  variable is popped off stackC for processing.

  TODO: PROBING5 would disable this limit; currently inactive. -- lh, 101127 --
*/
    int leftTotalStack = maxStack*CoinMax(200,maxProbe) ;
#ifdef PROBING5
    if (!info->inTree&&!info->pass)
      leftTotalStack = 1234567890 ;
#endif
    //printf("maxStack %d maxPass %d numberThisTime %d info pass %d\n",
    //   maxStack,maxPass,numberThisTime_,info->pass) ;
/*
  LOOKEDAT LOOP: HEAD

  Main loop to probe each variable in lookedAt_

  We're finally ready to get down to the business of probing. Open a loop to
  probe each variable in lookedAt_.
*/
    for (iLook = 0 ; iLook < numberThisTime_ ; iLook++) {
      double solval ;
      double down ;
      double up ;
/*
  Too successful? Consider bailing out.

  If we're generating row cuts, but haven't fixed any variables or we're in
  the tree, break from probing.

  JustFix < 0 says this is the first pass at the root, or (??). If it is the
  first pass at the root, turn off row cuts, keep on fixing variables, and
  stop at the end of this pass. If (??) and we haven't fixed any variables,
  break.

  TODO: Otherwise keep going with no changes? That doesn't seem right. The
  	logic here does not cover all situations.  -- lh, 101126 --
*/
      if (rowCut.outOfSpace() || leftTotalStack <= 0) {
	if (!justFix && (!nfixed || info->inTree)) {
#ifdef COIN_DEVELOP
	  if (!info->inTree)
	    printf("Exiting a on pass %d, maxProbe %d\n",
		   ipass,maxProbe) ;
#endif	  
	  break ;
	} else if (justFix <= 0) {
	  if (!info->inTree) {
	    rowCuts = 0 ;
	    justFix = 1 ;
	    disaggEffectiveness = COIN_DBL_MAX ;
	    needEffectiveness = COIN_DBL_MAX ;
	    //maxStack=10 ;
	    maxPass = 1 ;
	  } else if (!nfixed) {
#ifdef COIN_DEVELOP
	    printf("Exiting b on pass %d, maxProbe %d\n",
		   ipass,maxProbe) ;
#endif	  
	    break ;
	  }
	}
      }
/*
  awayFromBound isn't quite accurate --- we'll also set awayFromBound = 0
  for a general integer variable at an integer value strictly within bounds.
*/
      int awayFromBound = 1 ;
/*
  Have a look at the current variable. We're not interested in processing it
  if either or both bounds have been improved (0x02|0x01). If the variable has
  become fixed, claim both bounds have improved.

  TODO: if x<j> is not an integer variable we have big problems. This should
	be an assert. -- lh, 101126 --
*/
      int j = lookedAt_[iLook] ;
      solval = colsol[j] ;
      down = floor(solval+tolerance) ;
      up = ceil(solval-tolerance) ;
      if (columnGap[j] < 0.0) markC[j] = tightenUpper|tightenLower ;
      if ((markC[j]&(tightenUpper|tightenLower)) != 0 || !intVar[j]) continue ;
      double saveSolval = solval ;
/*
  `Normalize' variables that are near (or past) their bounds.

  In normal circumstances, u<j> - l<j> is at least 1, while tol will be
  down around 10e-5 (given a primal tolerance of 10e-6, about as large as
  you'd normally see). Then

  l<j> < l<j>+tol < l<j>+2tol << u<j>-2tol < u<j>-tol < u<j>

  For x<j> > u<j>-tol, (x<j>,u<j>,u<j>) => (u<j>-.1, l<j>, u<j>)
  For x<j> < l<j>+tol, (x<j>,l<j>,l<j>) => (l<j>+.1, l<j>, u<j>)
  For up == down,      (x<j>,v,v)       => (v+.1,v,v+1)

  So we're forcing a spread of 1 between up and down, and moving x<j>
  to something far enough (.1) from integer to avoid accuracy problems when
  calculating movement. A side effect is that primal infeasibility is
  corrected.

  TODO: Why this structure? Another approach would be to test using up and
  	down calculated just above. It works out pretty much the same if
	you stick with variables of type double. Moving to int would likely
	make it faster.  -- lh, 101127 --

  TODO: Things will get seriously wierd if (u<j> - l<j>) <= 3tol. For
        u<j>-tol < x<j> < l<j>+2tol, values will change. But the symmetric
	case u<j>-2tol < x<j> < l<j>+tol will never change values because
	of the order of the if statements. Clearly the code is not designed
	to cope with this level of pathology.  -- lh, 101127 --
  
  TODO: It's worth noting that awayFromBound is never used.  -- lh, 101127 --

  TODO: It's worth noting that the values set for solval are never used except
	in a debug print.  -- lh, 101213 --

  TODO: The final case below corresponds to l<j>+1 <= down == up < u<j>-1,
	e.g., an interior integer value. Seems like we'd want to probe +1 and
	-1. Presumably the code isn't set up to do this, so we arbitrarily
	pick the up probe. Could the code be augmented to do the stronger
	probe?  -- lh, 101127 --
*/
      if (solval >= colUpper[j]-tolerance ||
          solval <= colLower[j]+tolerance || up == down) {
	awayFromBound = 0 ;
	if (solval <= colLower[j]+2.0*tolerance) {
	  solval = colLower[j]+1.0e-1 ;
	  down = colLower[j] ;
	  up = down+1.0 ;
	} else if (solval >= colUpper[j]-2.0*tolerance) {
	  solval = colUpper[j]-1.0e-1 ;
	  up = colUpper[j] ;
	  down = up-1 ;
	} else {
          up = down+1.0 ;
          solval = down+1.0e-1 ;
        }
      }
      assert (up <= colUpper[j]) ;
      assert (down >= colLower[j]) ;
      assert (up > down) ;
#     if CGL_DEBUG > 0
      const double lj = colLower[j] ;
      const double uj = colUpper[j] ;
      const bool downIsLower = (fabs(down-colLower[j]) < 1.0e-7) ;
      const bool upIsUpper = (fabs(up-colUpper[j]) < 1.0e-7) ;
#     endif
/*
  Set up to probe each variable down (1), up (2), down (1).

  The notion is that we'll set a bit (1, 2, 4) in the result record if we're
  feasible for a particular way. Way is defined as:
    1: we're looking at movement between floor(x*<j>) and x*<j>, u<j>
    2: we're looking at movement between ceil(x*<j>) and x*<j>, l<j>
  As defined here, movement will be negative for way = 1.

  Why do we have a second pass with way = 1? Glad you asked. About 1200
  lines farther down, it finally becomes apparent that we execute the
  final pass only when the initial down probe was feasible, but the up probe
  was not. In a nutshell, we're restoring the state produced by the down
  probe, which we'll then nail down in the section of code labelled `keep',
  a mere 600 lines farther on.
*/
      unsigned int iway ;
      unsigned int way[] = { probeDown, probeUp, probeDown } ;
      unsigned int feasValue[] =
          { downProbeFeas, upProbeFeas, downProbe2Feas } ;
      unsigned int feasRecord = 0 ;
      bool notFeasible = false ;
      int istackC = 0 ;
      int istackR = 0 ;
/*
  PROBE LOOP: HEAD

  Open a loop to probe up and down for the current variable. Get things
  started by placing the variable on the propagation queue (stackC).

  As with the previous loops (variables in lookedAt_, passes) this loop
  extends to the end of major block #2).
*/
      for (iway = downProbe ; iway < oneProbeTooMany ; iway++) {
        stackC[0] = j ;
        markC[j] = way[iway] ;
/*
  Calculate movement given the current direction. Given the setup above,
  movement should be at least +/-1.

  TODO: Notice we reset up or down yet again. Really, this shouldn't be
  	necessary. For that matter, why did we bother with awayFromBound
	immediately above? -- lh, 101127 --

  TODO: The more I think about this, the less happy I am. Seems to me that
  	neither colUpper or colLower can change during iterations of this
	loop. If we discover monotone, we save bounds and break. If we
	don't discover monotone, we restore. Let's test that.
	-- lh, 110105 --
*/
        double solMovement ;
        double movement ;
        int goingToTrueBound = 0 ;

#	if CGL_DEBUG > 0
	assert(lj == colLower[j]) ;
	assert(uj == colUpper[j]) ;
	assert(downIsLower == (fabs(down-colLower[j]) < 1.0e-7)) ;
	assert(upIsUpper == (fabs(up-colUpper[j]) < 1.0e-7)) ;
#	endif

        if (way[iway] == probeDown) {
          movement = down-colUpper[j] ;
          solMovement = down-colsol[j] ;
          assert(movement < -0.99999) ;
          if (fabs(down-colLower[j]) < 1.0e-7) {
            goingToTrueBound = 2 ;
            down = colLower[j] ;
          }
        } else {
          movement = up-colLower[j] ;
          solMovement = up-colsol[j] ;
          assert(movement > 0.99999) ;
          if (fabs(up-colUpper[j]) < 1.0e-7) {
            goingToTrueBound = 2 ;
            up = colUpper[j] ;
          }
        }
/*
  Coding for goingToTrueBound is:
    0: We're somewhere in the interior of a general integer.
    1: We're heading for one of the bounds on a general integer.
    2: We're processing a binary variable (u<j>-l<j> = 1 and l<j> = 0).
  You can view this next test as correcting for the fact that the previous
  code assumes binary variables are all there is.

  TODO: Now that I've figured out the math for the constraints generated
	by disaggCuts (see typeset documentation), I see what's wrong
	(?) here. The disagg cut that's generated relies on the new probe
	bound being distance one from the original bound. But (for example)
	that's (colUpper-down) = 1 for a down probe.  Not quite what's
	being checked for here. But this next test ensures binary variables,
	and then the tests are equivalent. Before I willy-nilly change this, I
	need to sort out a couple of other places where goingToTrueBound
	controls behaviour.  -- lh, 101214 --

  TODO: Apropos the above, for example, the test for coefficient strengthening
  	cuts includes goingToTrueBound != 0  -- lh, 110105 --
*/
        if (goingToTrueBound && (colUpper[j]-colLower[j] > 1.5 || colLower[j]))
          goingToTrueBound = 1 ;
/*
  If the user doesn't want disaggregation cuts, pretend we're in the interior
  of a general integer variable.

  TODO: Why is row strengthening limited to binary variables? If we can
  	figure out what to do, the limitation seems artificial.
	-- lh, 101127 --
	The test here also says that we can't have coefficient strengthening
	without disaggregation. I don't see any reason for this dominance.
	-- lh, 110105 --
*/
        if ((rowCuts&1) == 0)
          goingToTrueBound = 0 ;
	bool canReplace = info->strengthenRow && (goingToTrueBound == 2) ;

#ifdef PRINT_DEBUG
	// Alert for general integer with nontrivial range.
	// Also last use of solval
        if (fabs(movement)>1.01) {
          printf("big %d %g %g %g\n",j,colLower[j],solval,colUpper[j]) ;
        }
#endif
/*
  Recall that we adjusted the reduced costs (djs) and current objective
  (current) to the minimisation sign convention. objVal will accumulate
  objective change due to forced bound changes.
*/
        double objVal = current ;
        if (solMovement*djs[j] > 0.0)
          objVal += solMovement*djs[j] ;
        nstackC = 1 ;
        nstackR = 0 ;
        saveL[0] = colLower[j] ;
        saveU[0] = colUpper[j] ;
        assert (saveU[0] > saveL[0]) ;
        notFeasible = false ;
/*
  Recalculate the upper (down probe) or lower (up probe) bound. You can see
  from the debug print that arbitrary integers are an afterthought in this
  code.

  TODO: And why not just set the bounds to down (colUpper) or up (colLower)?
        We set movement = down-colUpper just above, and we've taken a lot of
	care to groom up and down. This is pointless effort.  -- lh, 101214 --
*/
        if (movement < 0.0) {
          colUpper[j] += movement ;
          colUpper[j] = floor(colUpper[j]+0.5) ;
	  columnGap[j] = colUpper[j]-colLower[j]-primalTolerance_ ;
#ifdef PRINT_DEBUG
          printf("** Trying %d down to 0\n",j) ;
#endif
        } else {
          colLower[j] += movement ;
          colLower[j] = floor(colLower[j]+0.5) ;
	  columnGap[j] = colUpper[j]-colLower[j]-primalTolerance_ ;
#ifdef PRINT_DEBUG
          printf("** Trying %d up to 1\n",j) ;
#endif
        }
/*
  Is this variable now fixed?
*/
        if (fabs(colUpper[j]-colLower[j]) < 1.0e-6)
          markC[j] = tightenUpper|tightenLower ;
/*
  Clear the infinite bound bits (infUpper|infLower) and check again. Bounds
  may have improved.
*/
	markC[j] &= ~(infUpper|infLower) ;
	if (colUpper[j] > 1.0e10)
	  markC[j] |= infUpper ;
	if (colLower[j] < -1.0e10)
	  markC[j] |= infLower ;
        istackC = 0 ;
/*
  Update row bounds to reflect the change in variable bound.

  TODO: Replacing code to update row bounds with updateRowBounds(). We'll
	see how it looks.  -- lh, 101203 --
*/

  if (!updateRowBounds(j,movement,columnStart,columnLength,row,columnElements,
		       rowUpper,rowLower,nstackR,stackR,markR,
		       minR,maxR,saveMin,saveMax)) {
    notFeasible = true ;
    istackC = 1 ;
  }
/*
  PROBE LOOP: BEGIN PROPAGATION

  Row bounds are now adjusted for all rows with a<ij> != 0. Time to consider
  the effects. nstackC is incremented each time we add a variable, istackC is
  incremented each time we finish processing it.

  If we lost feasibility above, istackC = nstackC = 1 and this loop will not
  execute.

  TODO: Is stackC really a stack? Or a queue? And what is the role of
	maxStack? stackC is allocated at 2*ncols, but maxStack defaults to
	50. -- lh, 101127 --

  TODO: stackC is clearly a queue. Allocation strategy is still unclear.
        -- lh, 101128 --

  jjf:  could be istackC<maxStack?
*/
        while (istackC < nstackC && nstackC < maxStack) {
	  leftTotalStack-- ;
          int jway ;
          int jcol = stackC[istackC] ;
          jway = markC[jcol] ;
/*
  jjf: If not first and fixed then skip

  It's clear that x<j>, the probing variable, can be marked as fixed at this
  point (happens just above, when the upper or lower bound is adjusted). And
  there's an earlier check that avoids probing a variable that's fixed when
  plucked from lookedAt_. But ... why leave the if with an empty body?

  TODO: Disabled since r118 050225, but seems like the right thing to do. Is it
	somehow guaranteed that fixed variables are not stacked?
	-- lh, 101127 --

  TODO: Based on what's happening down below, it looks like we can encounter
	variables here which have had both bounds tightened and are then
	pushed onto istackC to propagate the effect of that. This bit of code
	should simply go away.  -- lh, 101127 --
*/
          if (((jway&(tightenLower|tightenUpper)) ==
	       (tightenLower|tightenUpper)) &&
	      istackC) {
            //istackC++ ;
            //continue ;
            //printf("fixed %d on stack\n",jcol) ;
          }
#if PROBING100
/*
  This block of code has been disabled since r661 080722. I have no notes on
  it from my old annotated CglProbing. Given the use of backward, it's tied to
  CglTreeProbingInfo.  -- lh, 101127 --
*/
	  if (backward) {
	    int jColumn = backward[jcol] ;
	    if (jColumn>=0) {
	      int nFix=0 ;
	      // 0-1 see what else could be fixed
	      if (jway==tightenUpper) {
		// fixed to 0
		int j ;
		for (j = toZero_[jColumn];j<toOne_[jColumn];j++) {
		  int kColumn=fixEntry_[j].sequence ;
		  kColumn = integerVariable_[kColumn] ;
		  bool fixToOne = fixEntry_[j].oneFixed ;
		  if (fixToOne) {
		    if (colLower[kColumn]==0.0) {
		      if (colUpper[kColumn]==1.0) {
			// See if on list
			if (!(markC[kColumn]&(tightenLower|tightenUpper))) {
			  if(nStackC<nCols) {
			    stackC[nstackC]=kColumn ;
			    saveL[nstackC]=colLower[kColumn] ;
			    saveU[nstackC]=colUpper[kColumn] ;
			    assert (saveU[nstackC]>saveL[nstackC]) ;
			    assert (nstackC<nCols) ;
			    nstackC++ ;
			    markC[kColumn]|=tightenLower ;
			    nFix++ ;
			  }
			} else if ((markC[kColumn]&(tightenLower|tightenUpper))==tightenUpper) {
			  notFeasible = true ;
			}
		      } else {
			// infeasible!
			notFeasible = true ;
		      }
		    }
		  } else {
		    if (colUpper[kColumn]==1.0) {
		      if (colLower[kColumn]==0.0) {
			// See if on list
			if (!(markC[kColumn]&(tightenLower|tightenUpper))) {
			  if(nStackC<nCols) {
			    stackC[nstackC]=kColumn ;
			    saveL[nstackC]=colLower[kColumn] ;
			    saveU[nstackC]=colUpper[kColumn] ;
			    assert (saveU[nstackC]>saveL[nstackC]) ;
			    assert (nstackC<nCols) ;
			    nstackC++ ;
			    markC[kColumn]|=tightenUpper ;
			    nFix++ ;
			  }
			} else if ((markC[kColumn]&(tightenLower|tightenUpper))==tightenLower) {
			  notFeasible = true ;
			}
		      } else {
			// infeasible!
			notFeasible = true ;
		      }
		    }
		  }
		}
	      } else if (jway==tightenLower) {
		int j ;
		for ( j=toOne_[jColumn];j<toZero_[jColumn+1];j++) {
		  int kColumn=fixEntry_[j].sequence ;
		  kColumn = integerVariable_[kColumn] ;
		  bool fixToOne = fixEntry_[j].oneFixed ;
		  if (fixToOne) {
		    if (colLower[kColumn]==0.0) {
		      if (colUpper[kColumn]==1.0) {
			// See if on list
			if (!(markC[kColumn]&(tightenLower|tightenUpper))) {
			  if(nStackC<nCols) {
			    stackC[nstackC]=kColumn ;
			    saveL[nstackC]=colLower[kColumn] ;
			    saveU[nstackC]=colUpper[kColumn] ;
			    assert (saveU[nstackC]>saveL[nstackC]) ;
			    assert (nstackC<nCols) ;
			    nstackC++ ;
			    markC[kColumn]|=tightenLower ;
			    nFix++ ;
			  }
			} else if ((markC[kColumn]&(tightenLower|tightenUpper))==tightenUpper) {
			  notFeasible = true ;
			}
		      } else {
			// infeasible!
			notFeasible = true ;
		      }
		    }
		  } else {
		    if (colUpper[kColumn]==1.0) {
		      if (colLower[kColumn]==0.0) {
			// See if on list
			if (!(markC[kColumn]&(tightenLower|tightenUpper))) {
			  if(nStackC<nCols) {
			    stackC[nstackC]=kColumn ;
			    saveL[nstackC]=colLower[kColumn] ;
			    saveU[nstackC]=colUpper[kColumn] ;
			    assert (saveU[nstackC]>saveL[nstackC]) ;
			    assert (nstackC<nCols) ;
			    nstackC++ ;
			    markC[kColumn]|=tightenUpper ;
			    nFix++ ;
			  }
			} else if ((markC[kColumn]&(tightenLower|tightenUpper))==tightenLower) {
			  notFeasible = true ;
			}
		      } else {
			// infeasible!
			notFeasible = true ;
		      }
		    }
		  }
		}
	      }
	    }
	  }
#endif  // PROBING100 (disabled)
/*
  PROBE LOOP: WALK COLUMN

  Loop to walk the column of a variable popped off stackC.

  We've pulled x<jcol> off stackC. We're going to walk the column and process
  each row where a<i,jcol> != 0. Note that we've already updated the row
  bounds to reflect the change in bound for x<jcol>. In effect, we've stacked
  this variable as a surrogate for stacking the rows whose bounds were
  changed by the change to bounds on this variable. Now we're going to examine
  those rows and see if any look promising for further propagation.

*/
          for (int k = columnStart[jcol] ;
	       k < columnStart[jcol]+columnLength[jcol] ; k++) {
            if (notFeasible)
              break ;
            int irow = row[k] ;
	    // Row already processed with these bounds.
	    if (markR[irow]/10000 > 0) continue ;

	    int rStart = rowStart[irow] ;
	    int rEnd = rowStartPos[irow] ;
	    double rowUp = rowUpper[irow] ;
	    double rowUp2 = 0.0 ;
	    bool doRowUpN ;
	    bool doRowUpP ;
/*
  So, is this row promising?

  If our current L<i> (minR) is larger than the original b<i> (rowUp), we're
  infeasible.

  Otherwise, a bit of linear algebra brings the conclusion that if the
  largest negative gap (min{j in N} a<ij>(u<j>-l<j>)) in the row is less than
  -(b<i>-L<i>), we aren't going to be able to make any progress shrinking the
  domains of variables with negative coefficients.  Similarly, if the largest
  positive gap (max{j in P} a<ij>(u<j>-l<j>) is greater than b<i>-L<i>, we
  can't shrink the domain of any variable with a positive coefficient.

  There are analogous conditions using U<i> and blow<i>.

  Note that we don't update the largest negative and positive gaps, so as
  propagation grinds on, this simple test becomes less accurate. On the other
  hand, the gaps (b<i>-L<i>) and (U<i>-blow<i>) are also shrinking. Still, it's
  another justification for the strict limits on propagation.

  Summary:

  minR = SUM{P}(a<ij>l<j>) + SUM{N}(a<ij>u<j>)
  maxR = SUM{P}(a<ij>u<j>) + SUM{N}(a<ij>l<j>)

  doRowUpN: a<ij> < 0, minR can violate rowUp => increase l<j> (more negative)
  doRowLoN: a<ij> < 0, maxR can violate rowLo => decrease u<j> (more positive)
  doRowUpP: a<ij> > 0, minR can violate rowUp => decrease u<j> (less positive)
  doRowLoP: a<ij> > 0, maxR can violate rowLo => increase l<j> (less negative)

  Note that the bound (l<j>, u<j>) that can be tightened in any given
  situation is *not* the bound involved in minR or maxR.

  For binary variables, of course, `shrink the domain' corresponds to fixing
  the variable.
*/

	    if (rowUp < 1.0e10) {
	      doRowUpN = true ;
	      doRowUpP = true ;
	      rowUp2 = rowUp-minR[irow] ;
	      if (rowUp2 < -primalTolerance_) {
		notFeasible = true ;
		break ;
	      } else {
		if (rowUp2+largestNegativeInRow[irow] > 0)
		  doRowUpN = false ;
		if (rowUp2-largestPositiveInRow[irow] > 0)
		  doRowUpP = false ;
	      }
	    } else {
	      doRowUpN = false ;
	      doRowUpP = false ;
	      rowUp2 = COIN_DBL_MAX ;
	    }
	    double rowLo = rowLower[irow] ;
	    double rowLo2 = 0.0 ;
	    bool doRowLoN ;
	    bool doRowLoP ;
	    if (rowLo > -1.0e10) {
	      doRowLoN = true ;
	      doRowLoP = true ;
	      rowLo2 = rowLo-maxR[irow] ;
	      if (rowLo2 > primalTolerance_) {
		notFeasible = true ;
		break ;
	      } else {
		if (rowLo2-largestNegativeInRow[irow] < 0)
		  doRowLoN = false ;
		if (rowLo2+largestPositiveInRow[irow] < 0)
		  doRowLoP = false ;
	      }
	    } else {
	      doRowLoN = false ;
	      doRowLoP = false ;
	      rowLo2 = -COIN_DBL_MAX ;
	    }
	    markR[irow] += 10000 ;
/*
  LOU_DEBUG: see if we're processing again without a bound change.
	    if (markR[irow]/10000 > 1)
	      std::cout
		<< "REDUNDANT: processing " << irow
		<< " for the " << markR[irow]/10000 << " time." << std::endl ;
*/
/*
  PROBE LOOP: WALK ROW

  Given the above analysis, work on the columns with negative coefficients.

  doRowUpN: a<ij> < 0, minR can violate rowUp => increase l<j> (more negative)
  doRowLoN: a<ij> < 0, maxR can violate rowLo => decrease u<j> (more positive)

  TODO: The asserts here should be compiled out unless we're debugging. 
	-- lh, 101210 --
*/
	    if (doRowUpN && doRowLoN) {
	      //doRowUpN=doRowLoN=false ;
	      // Start neg values loop
	      for (int kk = rStart ; kk < rEnd ; kk++) {
		int kcol = column[kk] ;
		int markIt = markC[kcol] ;
		// Skip columns with fixed variables.
		if ((markIt&(tightenLower|tightenUpper)) != (tightenLower|tightenUpper)) {
		  double value2 = rowElements[kk] ;
		  if (colUpper[kcol] <= 1e10)
		    assert ((markIt&infUpper) == 0) ;
		  else
		    assert ((markIt&infUpper) != 0) ;
		  if (colLower[kcol] >= -1e10)
		    assert ((markIt&infLower) == 0) ;
		  else
		    assert ((markIt&infLower) != 0) ;
		  assert (value2 < 0.0) ;
/*
  Not every column will be productive. Can we do anything with this one?
*/
		  double gap = columnGap[kcol]*value2 ;
		  bool doUp = (rowUp2+gap < 0.0) ;
		  bool doDown = (rowLo2-gap > 0.0) ;
		  if (doUp || doDown) {
		    double moveUp = 0.0 ;
		    double moveDown = 0.0 ;
		    double newUpper = -1.0 ;
		    double newLower = 1.0 ;
/*
  doUp attempts to increase the lower bound. The math looks like this:
    a<irow,kcol>x<kcol> + (minR[irow]-a<irow,kcol>u<kcol>) <= b<irow>
    value2*x<kcol> - value2*u<kcol> <= (b<irow>-minR[irow]) = rowUp2
    x<kcol> - u<kcol> >= rowUp2/value2
    l<kcol> = u<kcol> + rowUp2/value2
  Hence we need u<kcol> to be finite (0x8 not set). We also refuse to tighten
  l<kcol> a second time (0x2 not set).

  TODO: So, why don't we allow tightening a second time?  Clearly we might
	process a row on some other iteration that imposes a tighter bound.
	Is this an optimisation for binary variables, or is there some
	deeper issue at work in terms of correctness?  -- lh, 101127 --

  TODO: Also, it's clear that this bit of code would like to handle
	continuous variables, but it's also clear that it does it incorrectly.
	When the new lower bound straddles the existing upper bound, that
	looks to me to be sufficient justification to fix the variable. But
	the code here does just the opposite: merely tightening the bound is
	flagged as fixing the variable, with no consideration that we may have
	gone infeasible. But if we prove we straddle the upper bound, all we
	claim is we've tightened the upper bound!  -- lh, 101127 --
*/
		    if (doUp && ((markIt&(tightenLower|infUpper)) == 0)) {
		      double dbound = colUpper[kcol]+rowUp2/value2 ;
		      if (intVar[kcol]) {
			markIt |= tightenLower ;
			newLower = ceil(dbound-primalTolerance_) ;
		      } else {
			newLower = dbound ;
			if (newLower+primalTolerance_ > colUpper[kcol] &&
			    newLower-primalTolerance_ <= colUpper[kcol]) {
			  newLower = colUpper[kcol] ;
			  markIt |= tightenLower ;
			  //markIt = (tightenUpper|tightenLower) ;
			} else {
			  // avoid problems - fix later ?
			  markIt = (tightenUpper|tightenLower) ;
			}
		      }
		      moveUp = newLower-colLower[kcol] ;
		    }
/*
  The same, but attempt to decrease the upper bound by comparison with the
  blow for the row.
*/
		    if (doDown && ((markIt&(tightenUpper|infLower)) == 0)) {
		      double dbound = colLower[kcol] + rowLo2/value2 ;
		      if (intVar[kcol]) {
			markIt |= tightenUpper ;
			newUpper = floor(dbound+primalTolerance_) ;
		      } else {
			newUpper = dbound ;
			if (newUpper-primalTolerance_ < colLower[kcol] &&
			    newUpper+primalTolerance_ >= colLower[kcol]) {
			  newUpper = colLower[kcol] ;
			  markIt |= tightenUpper ;
			  //markIt = (tightenUpper|tightenLower) ;
			} else {
			  // avoid problems - fix later ?
			  markIt = (tightenUpper|tightenLower) ;
			}
		      }
		      moveDown = newUpper-colUpper[kcol] ;
		    }
/*
  Is there any change to propagate? Assuming we haven't exceeded the limit
  for propagating, queue the variable. (If the variable's already on the list,
  that's not a problem.)

  Note that markIt is initially loaded from markC and whatever bounds are
  changed are then or'd in. So changes accumulate.
*/
		    if (!moveUp && !moveDown)
		      continue ;
		    bool onList = ((markC[kcol]&(tightenLower|tightenUpper)) != 0) ;
		    if (nstackC < 2*maxStack) {
		      markC[kcol] = markIt ;
		    }
		    if (moveUp && (nstackC < 2*maxStack)) {
		      if (!onList) {
			stackC[nstackC] = kcol ;
			saveL[nstackC] = colLower[kcol] ;
			saveU[nstackC] = colUpper[kcol] ;
			assert (saveU[nstackC] > saveL[nstackC]) ;
			assert (nstackC < nCols) ;
			nstackC++ ;
			onList = true ;
		      }
/*
  The first assert here says `If we're minimising, and d<j> < 0, then
  x<j> must be NBUB, so if the new l<j> > x*<j> = u<j>, we're infeasible.'

  TODO: Why haven't we detected infeasibility? Put another way, why isn't
	there an assert(notFeasible)?  -- lh, 101127 --

  TODO: Seems like the change in obj should be ((new l<j>) - x*<j>)*d<j>, but
  	moveUp is (newLower-colLower). Was there some guarantee that x*<j>
	equals l<j>?  -- lh, 101127 --
  TODO: I think this is another symptom of incomplete conversion for general
	integers. For binary variables, it's not an issue.  -- lh, 101203 --
*/
		      if (newLower > colsol[kcol]) {
			if (djs[kcol] < 0.0) {
			  assert (newLower > colUpper[kcol]+primalTolerance_) ;
			} else {
			  objVal += moveUp*djs[kcol] ;
			}
		      }
/*
  TODO: Up above we've already set newLower = ceil(dbound-primalTolerance_).
        It's highly unlikely that ceil(newLower-1.0e4) will be different.
	-- lh, 101127 --
  TODO: My first inclination is to claim that newLower should never be worse
        than the existing bound. But I'd better extend the benefit of the
	doubt for a bit. Some rows will produce stronger bounds than others.
	-- lh, 101127 --
*/
		      if (intVar[kcol]) 
			newLower = CoinMax(colLower[kcol],
					   ceil(newLower-1.0e-4)) ;
		      colLower[kcol] = newLower ;
		      columnGap[kcol] = colUpper[kcol]-newLower-primalTolerance_ ;
/*
  Propagate the column bound change.

  TODO: Notice that we're not setting a variable here when we discover
        infeasibility; rather, we're setting the column bound to a ridiculous
	value which will be discovered farther down.  -- lh, 101127 --
  TODO: Replacing code to update row bounds with updateRowBounds(). We'll
	see how it looks.  -- lh, 101203 --
*/
		      if (fabs(colUpper[kcol]-colLower[kcol]) < 1.0e-6) {
			markC[kcol] = (tightenLower|tightenUpper) ;
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol] > 1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol] < -1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveUp,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colLower[kcol] = 1.0e50 ;
		      }
		    }
/*
  Repeat the whole business, in the down direction. Stack the variable, if
  it's not already there, do the infeasibility check / objective update, and
  update minR / maxR for affected constraints.
*/
		    if (moveDown&&nstackC<2*maxStack) {
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
			assert (saveU[nstackC]>saveL[nstackC]) ;
			assert (nstackC<nCols) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_) ;
			} else {
			  objVal += moveDown*djs[kcol] ;
			}
		      }
		      if (intVar[kcol])
			newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4)) ;
		      colUpper[kcol]=newUpper ;
		      columnGap[kcol] = newUpper-colLower[kcol]-primalTolerance_ ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol] = (tightenLower|tightenUpper); // say fixed
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol]>1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol]<-1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveDown,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colUpper[kcol] = -1.0e50 ;
		      }
		    }
/*
  We've propagated the bound changes for this column. Check to see if we
  discovered infeasibility. Abort by claiming we've cleared the propagation
  stack.
*/
		    if (colLower[kcol] > colUpper[kcol]+primalTolerance_) {
		      notFeasible = true ;
		      k = columnStart[jcol]+columnLength[jcol] ;
		      istackC = nstackC+1 ;
		      break ;
		    }
		  }   // end if (doUp || doDown) (productive column)
		}  // end column not fixed
	      } // end loop on negative coefficients of this row
	    } else if (doRowUpN) {
/*
  The above block propagated change for a variable with a negative coefficient,
  where we could work against both the row upper and lower bounds. Now we're
  going to do it again, but specialised for the case where we can only work
  against the upper bound.
*/
	      // Start neg values loop
	      for (int kk = rStart ; kk < rEnd ; kk++) {
		int kcol =column[kk] ;
		int markIt=markC[kcol] ;
		if ((markIt&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[kk] ;
		  double gap = columnGap[kcol]*value2 ;
		  if (!(rowUp2 + gap < 0.0))
		    continue ;
		  double moveUp=0.0 ;
		  double newLower=1.0 ;
		  if ((markIt&(tightenLower|infUpper))==0) {
		    double dbound = colUpper[kcol]+rowUp2/value2 ;
		    if (intVar[kcol]) {
		      markIt |= tightenLower ;
		      newLower = ceil(dbound-primalTolerance_) ;
		    } else {
		      newLower=dbound ;
		      if (newLower+primalTolerance_>colUpper[kcol]&&
			  newLower-primalTolerance_<=colUpper[kcol]) {
			newLower=colUpper[kcol] ;
			markIt |= tightenLower ;
			//markIt= (tightenLower|tightenUpper) ;
		      } else {
			// avoid problems - fix later ?
			markIt= (tightenLower|tightenUpper) ;
		      }
		    }
		    moveUp = newLower-colLower[kcol] ;
		    if (!moveUp)
		      continue ;
		    bool onList = ((markC[kcol]&(tightenLower|tightenUpper))!=0) ;
		    if (nstackC<2*maxStack) {
		      markC[kcol] = markIt ;
		    }
		    if (moveUp&&nstackC<2*maxStack) {
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
			assert (saveU[nstackC]>saveL[nstackC]) ;
			assert (nstackC<nCols) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newLower>colsol[kcol]) {
			if (djs[kcol]<0.0) {
			  /* should be infeasible */
			  assert (newLower>colUpper[kcol]+primalTolerance_) ;
			} else {
			  objVal += moveUp*djs[kcol] ;
			}
		      }
		      if (intVar[kcol]) 
			newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4)) ;
		      colLower[kcol]=newLower ;
		      columnGap[kcol] = colUpper[kcol]-newLower-primalTolerance_ ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]= (tightenLower|tightenUpper); // say fixed
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol]>1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol]<-1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveUp,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colLower[kcol] = 1.0e50 ;
		      }
		    }
		    if (colLower[kcol] > colUpper[kcol]+primalTolerance_) {
		      notFeasible = true ;
		      k = columnStart[jcol]+columnLength[jcol] ;
		      istackC = nstackC+1 ;
		      break ;
		    }
		  }
		}
	      } // end big loop rStart->rPos
	    } else if (doRowLoN) {
/*
  And yet again, for the case where we can only work against the lower bound.
*/
	      // Start neg values loop
	      for (int kk = rStart ; kk < rEnd ; kk++) {
		int kcol =column[kk] ;
		if ((markC[kcol]&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double moveDown=0.0 ;
		  double newUpper=-1.0 ;
		  double value2=rowElements[kk] ;
		  int markIt=markC[kcol] ;
		  assert (value2<0.0) ;
		  double gap = columnGap[kcol]*value2 ;
		  bool doDown = (rowLo2 -gap > 0.0) ;
		  if (doDown && ((markIt&(tightenUpper|infLower)) == 0)) {
		    double dbound = colLower[kcol] + rowLo2/value2 ;
		    if (intVar[kcol]) {
		      markIt |= tightenUpper ;
		      newUpper = floor(dbound+primalTolerance_) ;
		    } else {
		      newUpper=dbound ;
		      if (newUpper-primalTolerance_<colLower[kcol]&&
			  newUpper+primalTolerance_>=colLower[kcol]) {
			newUpper=colLower[kcol] ;
			markIt |= tightenUpper ;
			//markIt= (tightenLower|tightenUpper) ;
		      } else {
			// avoid problems - fix later ?
			markIt= (tightenLower|tightenUpper) ;
		      }
		    }
		    moveDown = newUpper-colUpper[kcol] ;
		    if (!moveDown)
		      continue ;
		    bool onList = ((markC[kcol]&(tightenLower|tightenUpper))!=0) ;
		    if (nstackC<2*maxStack) {
		      markC[kcol] = markIt ;
		    }
		    if (moveDown&&nstackC<2*maxStack) {
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
			assert (saveU[nstackC]>saveL[nstackC]) ;
			assert (nstackC<nCols) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_) ;
			} else {
			  objVal += moveDown*djs[kcol] ;
			}
		      }
		      if (intVar[kcol])
			newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4)) ;
		      colUpper[kcol]=newUpper ;
		      columnGap[kcol] = newUpper-colLower[kcol]-primalTolerance_ ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]= (tightenLower|tightenUpper); // say fixed
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol]>1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol]<-1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveDown,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colUpper[kcol] = -1.0e50 ;
		      }
		    }
		    if (colLower[kcol] > colUpper[kcol]+primalTolerance_) {
		      notFeasible = true ;
		      k = columnStart[jcol]+columnLength[jcol] ;
		      istackC = nstackC+1 ;
		      break ;
		    }
		  }
		}
	      }  // end big loop rStart->rPos
	    }
/*
  We've finished working on the negative coefficients of the row. Advance
  rStart and rEnd to cover the positive coefficients and repeat the previous
  500 lines.
*/
	    rStart = rEnd ;
	    rEnd = rowStart[irow+1] ;
	    if (doRowUpP&&doRowLoP) {
	      //doRowUpP=doRowLoP=false ;
	      // Start pos values loop
	      for (int kk =rStart;kk<rEnd;kk++) {
		int kcol=column[kk] ;
		int markIt=markC[kcol] ;
		if ((markIt&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[kk] ;
		  assert (value2 > 0.0) ;
		  /* positive element */
		  double gap = columnGap[kcol]*value2 ;
		  bool doDown = (rowLo2 + gap > 0.0) ;
		  bool doUp = (rowUp2 - gap < 0.0) ;
		  if (doDown||doUp) {
		    double moveUp=0.0 ;
		    double moveDown=0.0 ;
		    double newUpper=-1.0 ;
		    double newLower=1.0 ;
		    if (doDown && ((markIt&(tightenLower|infUpper)) == 0)) {
		      double dbound = colUpper[kcol] + rowLo2/value2 ;
		      if (intVar[kcol]) {
			markIt |= tightenLower ;
			newLower = ceil(dbound-primalTolerance_) ;
		      } else {
			newLower=dbound ;
			if (newLower+primalTolerance_>colUpper[kcol]&&
			    newLower-primalTolerance_<=colUpper[kcol]) {
			  newLower=colUpper[kcol] ;
			  markIt |= tightenLower ;
			  //markIt= (tightenLower|tightenUpper) ;
			} else {
			  // avoid problems - fix later ?
			  markIt= (tightenLower|tightenUpper) ;
			}
		      }
		      moveUp = newLower-colLower[kcol] ;
		    }
		    if (doUp && ((markIt&(tightenUpper|infLower)) == 0)) {
		      double dbound = colLower[kcol] + rowUp2/value2 ;
		      if (intVar[kcol]) {
			markIt |= tightenUpper ;
			newUpper = floor(dbound+primalTolerance_) ;
		      } else {
			newUpper=dbound ;
			if (newUpper-primalTolerance_<colLower[kcol]&&
			    newUpper+primalTolerance_>=colLower[kcol]) {
			  newUpper=colLower[kcol] ;
			  markIt |= tightenUpper ;
			  //markIt= (tightenLower|tightenUpper) ;
			} else {
			  // avoid problems - fix later ?
			  markIt= (tightenLower|tightenUpper) ;
			}
		      }
		      moveDown = newUpper-colUpper[kcol] ;
		    }
		    if (!moveUp&&!moveDown)
		      continue ;
		    bool onList = ((markC[kcol]&(tightenLower|tightenUpper))!=0) ;
		    if (nstackC<2*maxStack) {
		      markC[kcol] = markIt ;
		    }
		    if (moveUp&&nstackC<2*maxStack) {
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
			assert (saveU[nstackC]>saveL[nstackC]) ;
			assert (nstackC<nCols) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newLower>colsol[kcol]) {
			if (djs[kcol]<0.0) {
			  /* should be infeasible */
			  assert (newLower>colUpper[kcol]+primalTolerance_) ;
			} else {
			  objVal += moveUp*djs[kcol] ;
			}
		      }
		      if (intVar[kcol]) 
			newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4)) ;
		      colLower[kcol]=newLower ;
		      columnGap[kcol] = colUpper[kcol]-newLower-primalTolerance_ ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]= (tightenLower|tightenUpper); // say fixed
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol]>1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol]<-1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveUp,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colLower[kcol] = 1.0e50 ;
		      }
		    }
		    if (moveDown&&nstackC<2*maxStack) {
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
			assert (saveU[nstackC]>saveL[nstackC]) ;
			assert (nstackC<nCols) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_) ;
			} else {
			  objVal += moveDown*djs[kcol] ;
			}
		      }
		      if (intVar[kcol])
			newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4)) ;
		      colUpper[kcol]=newUpper ;
		      columnGap[kcol] = newUpper-colLower[kcol]-primalTolerance_ ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]= (tightenLower|tightenUpper); // say fixed
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol]>1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol]<-1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveDown,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colUpper[kcol] = -1.0e50 ;
		      }
		    }
		    if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
		      notFeasible = true; ;
		      k=columnStart[jcol]+columnLength[jcol] ;
		      istackC=nstackC+1 ;
		      break ;
		    }
		  }
		}
	      } // end big loop rPos->rEnd
	    } else if (doRowUpP) {
	      // Start pos values loop
	      for (int kk =rStart;kk<rEnd;kk++) {
		int kcol =column[kk] ;
		int markIt=markC[kcol] ;
		if ((markIt&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[kk] ;
		  assert (value2 > 0.0) ;
		  /* positive element */
		  double gap = columnGap[kcol]*value2 ;
		  bool doUp = (rowUp2 - gap < 0.0) ;
		  if (doUp && ((markIt&(tightenUpper|infLower)) == 0)) {
		    double newUpper=-1.0 ;
		    double dbound = colLower[kcol] + rowUp2/value2 ;
		    if (intVar[kcol]) {
		      markIt |= tightenUpper ;
		      newUpper = floor(dbound+primalTolerance_) ;
		    } else {
		      newUpper=dbound ;
		      if (newUpper-primalTolerance_<colLower[kcol]&&
			  newUpper+primalTolerance_>=colLower[kcol]) {
			newUpper=colLower[kcol] ;
			markIt |= tightenUpper ;
			//markIt= (tightenLower|tightenUpper) ;
		      } else {
			// avoid problems - fix later ?
			markIt= (tightenLower|tightenUpper) ;
		      }
		    }
		    double moveDown = newUpper-colUpper[kcol] ;
		    if (!moveDown)
		      continue ;
		    bool onList = ((markC[kcol]&(tightenLower|tightenUpper))!=0) ;
		    if (nstackC<2*maxStack) {
		      markC[kcol] = markIt ;
		    }
		    if (nstackC<2*maxStack) {
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
			assert (saveU[nstackC]>saveL[nstackC]) ;
			assert (nstackC<nCols) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_) ;
			} else {
			  objVal += moveDown*djs[kcol] ;
			}
		      }
		      if (intVar[kcol])
			newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4)) ;
		      colUpper[kcol]=newUpper ;
		      columnGap[kcol] = newUpper-colLower[kcol]-primalTolerance_ ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]= (tightenLower|tightenUpper); // say fixed
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol]>1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol]<-1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveDown,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colUpper[kcol] = -1.0e50 ;
		      }
		    }
		    if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
		      notFeasible = true ;
		      k=columnStart[jcol]+columnLength[jcol] ;
		      istackC=nstackC+1 ;
		      break ;
		    }
		  }
		}
	      } // end big loop rPos->rEnd
	    } else if (doRowLoP) {
	      // Start pos values loop
	      for (int kk =rStart;kk<rEnd;kk++) {
		int kcol =column[kk] ;
		if ((markC[kcol]&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[kk] ;
		  int markIt=markC[kcol] ;
		  assert (value2 > 0.0) ;
		  /* positive element */
		  double gap = columnGap[kcol]*value2 ;
		  bool doDown = (rowLo2 +gap > 0.0) ;
		  if (doDown&&(markIt&(tightenLower|infUpper))==0) {
		    double newLower=1.0 ;
		    double dbound = colUpper[kcol] + rowLo2/value2 ;
		    if (intVar[kcol]) {
		      markIt |= tightenLower ;
		      newLower = ceil(dbound-primalTolerance_) ;
		    } else {
		      newLower=dbound ;
		      if (newLower+primalTolerance_>colUpper[kcol]&&
			  newLower-primalTolerance_<=colUpper[kcol]) {
			newLower=colUpper[kcol] ;
			markIt |= tightenLower ;
			//markIt= (tightenLower|tightenUpper) ;
		      } else {
			// avoid problems - fix later ?
			markIt= (tightenLower|tightenUpper) ;
			}
		    }
		    double moveUp = newLower-colLower[kcol] ;
		    if (!moveUp)
		      continue ;
		    bool onList = ((markC[kcol]&(tightenLower|tightenUpper))!=0) ;
		    if (nstackC<2*maxStack) {
		      markC[kcol] = markIt ;
		    }
		    if (nstackC<2*maxStack) {
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
			assert (saveU[nstackC]>saveL[nstackC]) ;
			assert (nstackC<nCols) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newLower>colsol[kcol]) {
			if (djs[kcol]<0.0) {
			  /* should be infeasible */
			  assert (newLower>colUpper[kcol]+primalTolerance_) ;
			} else {
			  objVal += moveUp*djs[kcol] ;
			}
		      }
		      if (intVar[kcol]) 
			newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4)) ;
		      colLower[kcol]=newLower ;
		      columnGap[kcol] = colUpper[kcol]-newLower-primalTolerance_ ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]= (tightenLower|tightenUpper); // say fixed
		      }
		      markC[kcol] &= ~(infLower|infUpper) ;
		      if (colUpper[kcol]>1.0e10)
			markC[kcol] |= infUpper ;
		      if (colLower[kcol]<-1.0e10)
			markC[kcol] |= infLower ;

		      if (!updateRowBounds(kcol,moveUp,
			       columnStart,columnLength,row,columnElements,
			       rowUpper,rowLower,nstackR,stackR,markR,
			       minR,maxR,saveMin,saveMax)) {
			colLower[kcol] = 1.0e50 ;
		      }
		    }
		    if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
		      notFeasible = true ;
		      k=columnStart[jcol]+columnLength[jcol] ;
		      istackC=nstackC+1 ;
		      break ;
		    }
		  }
		}
	      }      // end big loop rPos->rEnd
	    }    // end processing of positive coefficients of row.
          }   // end loop to walk the column of a variable popped off stackC
          istackC++ ;
        }  // end stackC processing loop
/*
  PROBE LOOP: END PROPAGATION

  End propagation of probe in one direction for a single variable. Hard to
  believe 1,000 lines of code can achieve so little.
*/
/*
  PROBE LOOP: INFEASIBILITY

  Feasibility check. Primal infeasibility is noted in notFeasible; we still
  need to test against the objective bound. If we have shown up and down
  infeasibility, we're infeasible, period. Break from the probe loop and
  terminate the lookedAt and pass loops by forcing their iteration counts past
  the maximum.

  If we're not feasible, then we surely can't drive the variable to bound,
  so reset goingToTrueBound.

  TODO: I'm coming to think forcing goingToTrueBound = 0 is really an overload
	to suppress various post-probe activities (cut generation, singleton
	motion, etc). It might be that this is not the best approach.
	-- lh, 101216 --
*/
	if (notFeasible || (objVal > cutoff)) {
#	  if CGL_DEBUG > 1
	  std::cout << "  Column " << j << " infeasible: " ;
	  if (notFeasible && (objVal > cutoff))
	    std::cout << "primal bounds and objective cutoff" ;
	  else if (notFeasible)
	    std::cout << "primal bounds" ;
	  else
	    std::cout << "objective cutoff" ;
	  std::cout << "." << std::endl ;
#	  endif
	  notFeasible = true ;
	  if (iway == upProbe && feasRecord == 0) {
	    ninfeas = 1 ;
	    j = nCols-1 ;
	    iLook = numberThisTime_ ;
	    ipass = maxPass ;
	    break ;
	  }
	  goingToTrueBound = 0 ;
	} else {
	  feasRecord |= feasValue[iway] ; 
	}
/*
  Save the implications? Currently restricted to binary variables. If info
  were a CglTreeProbing object, fixes() would save the implication in the
  arrays kept there, returning false only when space is exceeded.
  
  TODO: This would be one of the places that will fail if fixes() and
  	initializeFixing() are removed from CglTreeInfo.  -- lh, 101127 --
*/
	if (!notFeasible && saveFixingInfo) {
	  // save fixing info
	  assert (j == stackC[0]) ;
	  int toValue = (way[iway] == probeDown) ? -1 : +1 ;
	  for (istackC = 1 ; istackC < nstackC ; istackC++) {
	    int icol = stackC[istackC] ;
	    // for now back to just 0-1
	    if (colUpper[icol]-colLower[icol]<1.0e-12 &&
		!saveL[istackC] && saveU[istackC]==1.0) {
	      assert(saveL[istackC] == colLower[icol] ||
		     saveU[istackC] == colUpper[icol]) ;
	      saveFixingInfo =
		  info->fixes(j,toValue,icol,
			      (colLower[icol] == saveL[istackC])) ;
	    }
	  }
	}
/*
  PROBE LOOP: MONOTONE (BREAK AND KEEP)

  We can't prove infeasibility. What's left? If we're at the end of the up
  probe, we may know enough to prove monotonicity:

  * If this is the up probe and we were feasible on this probe but not on
    the down pass, we're monotone.
  * If this is the second down probe, we were feasible on down, infeasible on
    up (monotone) and repeated the down probe to recreate the bounds.

  Either way, execute monotoneActions to capture the state and move on to the
  next variable. Terminate the probe loop by forcing iway = oneProbeTooMany.

  This if-then-else extends to the end of the down/up/down probe loop.

  TODO: There's something not right about this `second down pass' business.
        Why do we copy stackC to stackC0 on a feasible first pass? Surely that
	could be used to restore the state of the first down pass.
	-- lh, 101128 --
*/
        if (iway == downProbe2 ||
	    (iway == upProbe && feasRecord == upProbeFeas)) {
          nfixed++ ;
	  bool retVal = monotoneActions(primalTolerance_,si,cs,
	  				nstackC,stackC,intVar,
					colLower,colsol,colUpper,
					index,element) ;
	  if (retVal) anyColumnCuts = true ;
	  clearStacks(primalTolerance_,nstackC,stackC,markC,colLower,colUpper,
	  	      nstackR,stackR,markR) ;
	  break ;
        }
/*
  PROBE LOOP: ITERATE

  If we reach here, we will iterate the probe loop.

  Cases remaining include

  a) End of the first down probe; in this case we will iterate regardless of
     feasibility.
  b) End of the up probe, up and down feasible; in this case we will not
     iterate (we're done, cannot prove monotonicity).
  c) End of the up probe, down feasible, up infeasible; in this case we will
     iterate to restore the down feasible state.

  The next block deals with preparing to iterate at the end of the first down
  probe.  Whatever the result, we'll need to go to the up probe. If the down
  probe was feasible, save stackC and the new bounds.

  TODO: Unless I miss my guess, large chunks of this block of code will be
	replicated in the code that handles case b), as we generate cuts
	related to forcing the variable up.  Further inspection confirms
	this. -- lh, 101128 --

  TODO: I'm becoming more convinced that the second down probe iteration is
	unnecessary. Given that we have stackC0, lo0, and up0, seems like
	calling monotoneActions with lo0 = colLower, up0 = colUpper, and
	stackC0 = stackC should do the trick.  -- lh, 101216 --

  jjf: is it worth seeing if can increase coefficients or maybe better see
       if it is a cut

*/
	if (iway == downProbe) {
	  nstackC0 = CoinMin(nstackC,maxStack) ;
	  if (notFeasible) {
	    nstackC0 = 0 ;
	  } else {
	    for (istackC = 0 ; istackC < nstackC0 ; istackC++) {
	      int icol = stackC[istackC] ;
	      stackC0[istackC] = icol ;
	      lo0[istackC] = colLower[icol] ;
	      up0[istackC] = colUpper[icol] ;
	    }
	  }
/*
  Run through stackC and make disaggregation cuts, if possible. At the end of
  each iteration, restore the bounds for the variable and reset markC.

  goingToTrueBound was set to 0 way back in initial setup unless the user
  requested disagg cuts (rowCuts&0x01).
*/
	  assert (iway == downProbe) ;
	  if (goingToTrueBound == 2 && !justReplace) {
	    disaggCuts(nstackC,probeDown,primalTolerance_,
		       disaggEffectiveness,si,
		       rowCut,stackC,colsol,colUpper,colLower,saveU,saveL,
		       index,element) ;
	  }
	  for (istackC = nstackC-1 ; istackC >= 0 ; istackC--) {
	    int icol = stackC[istackC] ;
	    double oldU = saveU[istackC] ;
	    double oldL = saveL[istackC] ;
	    colUpper[icol] = oldU ;
	    colLower[icol] = oldL ;
	    columnGap[icol] = oldU-oldL-primalTolerance_ ;
	    markC[icol] = 0 ;
	    if (oldU > 1.0e10)
	      markC[icol] |= infUpper ;
	    if (oldL < -1.0e10)
	      markC[icol] |= infLower ;
	  }
/*
  And now for the rows. Remember, we're still finishing off the first down
  probe and readying for the next iteration (up probe). The final three lines
  (of about 300) actually restore the row bounds. The rest is cut generation.
*/
	  for (istackR = 0 ; istackR < nstackR ; istackR++) {
	    int irow = stackR[istackR] ;
/*
  Generate coefficient strengthening cuts if the user has asked for them,
  and move singletons if possible.

  TODO: It looks like `move singletons' amounts to this: If all variables
	in the row except the probing variable are column singletons, we
	should be able to simply transfer them into the objective and let
	them bang against their implicit bounds. And now that I write that,
	I'm suspicious. More investigation is required.  -- lh, 101127 --

  TODO: Figure out why goingToTrueBound is relevant here, other than it's
	set to zero if we're infeasible.  -- lh, 101127 --
	I think that might be all there is to it.   -- lh, 101128 --
	Well, maybe not. Check the math for the contained cuts to see if
	it's valid only for binary variables; when within distance 1 of
	the upper or lower bound of a general integer; or more generally.
	GoingToTrueBound indicates the first of these two. -- lh, 110105 --

  TODO: Figure out what the third and fourth comment lines below (`Taken out
	...') mean. Because MOVE_SINGLETONS is active. If I were to hazard
	a guess, the problem would be that the constraint system being probed
	may have been preprocessed from the original and a column that looks
	like a singleton here may not be a singleton on the original system.
	Which would explain columnLength2, and why this code is again active.

	But in the end, I'm inclined towards `should be found elsewhere', as
	in `do this as a preprocessing step.'
	-- lh, 101127 --

  generate coefficient strengthening cuts (1)
*/
	    if ((rowCuts&2) != 0 && goingToTrueBound) {
	      bool ifCut = anyColumnCuts ;
/*
  Start out by working against the row upper bound U(i). Note that we're
  calculating row activity (sum) and will use it to decide on the viability of
  generating a cut. So if we've generated column cuts (colsol now violates a
  bound for some x<j>) then the calculation of row activity will be incorrect.

  TODO: Sort out why anyColumnCuts is an obstacle. Theoretical or practical?
	Could we use the new bound instead? The guard predates the singleton
	business.   -- lh, 101128 --
*/
	      double gap = rowUpper[irow]-maxR[irow] ;
	      double sum = 0.0 ;
	      if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
		// see if the strengthened row is a cut
		// also see if singletons can go to good objective
		// Taken out as should be found elsewhere
		// and has to be original column length
#ifdef MOVE_SINGLETONS
		bool moveSingletons = true ;
#endif
		for (int kk = rowStart[irow] ; kk < rowStart[irow+1] ; kk++) {
		  int iColumn = column[kk] ;
		  double value = rowElements[kk] ;
		  sum += value*colsol[iColumn] ;
#ifdef MOVE_SINGLETONS
		  if (moveSingletons && j != iColumn) {
		    if (colUpper[iColumn] > colLower[iColumn]) {
		      if (columnLength2[iColumn] != 1) {
			moveSingletons = false ;
		      }
		    }
		  }
#endif
		}
#ifdef MOVE_SINGLETONS
/*
  jjf: can fix any with good costs.

  And room on stackC0. Comparing with similar code in the segment where down
  and up probes both show feasible, we're pretending that this was already
  on stackC and hence was copied to stackC0 above.

  TODO: Seems to me that the checks here for column length and bounds are
	unnecessary. If they passed above, they'll pass here. -- lh, 101127 --

  TODO: Work out the math for this. In my previous passes over this bit of
	code, I keep coming to the same conclusion: the math is wrong. I
	need to take a serious run at it. See notes in CglProbingAnn.
	-- lh, 101128 --
*/
		if (moveSingletons) {
		  for (int kk = rowStart[irow] ;
		       kk < rowStart[irow+1] ; kk++) {
		    int iColumn = column[kk] ;
		    if (j != iColumn) {
		      if (colUpper[iColumn] > colLower[iColumn]) {
			if (columnLength2[iColumn] == 1) {
			  double value = rowElements[kk] ;
			  if (direction*objective[iColumn]*value < 0.0 &&
			      !(markC[iColumn]&(tightenLower|tightenUpper))) {
			    if (nstackC0+1 < maxStack) {
			      stackC0[nstackC0] = iColumn ;
			      if (value > 0.0) {
				lo0[nstackC0] = colUpper[iColumn] ;
				up0[nstackC0] = colUpper[iColumn] ;
			      } else {
				lo0[nstackC0] = colLower[iColumn] ;
				up0[nstackC0] = colLower[iColumn] ;
			      }
			      nstackC0++ ;
			    }
			  }
			}
		      }
		    }
		  }
		}
#endif
/*
  jjf: can be a cut
       subtract gap from upper and integer coefficient
       saveU and saveL spare

  Given b<i> - U<i> > 0 as a result of probing x<j>, we can reduce the
  coefficient of x<j> and reduce b<i>.

  Borrow saveU and saveL as working arrays; their usefulness ended when stackC
  processing finished, above.

  TODO: Just what sort of cut are we making here? Look this up when next up
	at SFU. Must be in the presentation notes.  -- lh, 101128 --

  TODO: Why are we borrowing saveU and saveL? There are method-global element
	and index arrays used everywhere else.  -- lh, 101128 --
*/
		if (sum-gap*colsol[j] > maxR[irow]+primalTolerance_ ||
		    (info->strengthenRow && rowLower[irow] < -1.0e20)) {
		  int *index = reinterpret_cast<int *>(saveL) ;
		  double *element = saveU ;
		  int n = 0 ;
		  bool coefficientExists = false ;
		  double sum2 = 0.0 ;
/*
  Generate the coefficients.
*/
		  for (int kk = rowStart[irow] ;
		       kk < rowStart[irow+1] ; kk++) {
		    int kColumn = column[kk] ;
		    double el = rowElements[kk] ;
		    if (kColumn != j) {
		      index[n] = kColumn ;
		      element[n++] = el ;
		    } else {
		      el = el-gap ;
		      if (fabs(el) > 1.0e-12) {
			index[n] = kColumn ;
			element[n++] = el ;
		      }
		      coefficientExists = true ;
		    }
		    sum2 += colsol[kColumn]*el ;
		  }
		  if (!coefficientExists) {
		    index[n] = j ;
		    element[n++] = -gap ;
		    sum2 -= colsol[j]*gap ;
		  }
/*
  Check effectiveness. Add the cut only if it's sufficiently effective.

  If we're simply replacing an existing constraint, this is a little bit
  misleading. Effectiveness is set extremely low.
*/
		  OsiRowCut rc ;
		  rc.setLb(-DBL_MAX) ;
		  double ub = rowUpper[irow]-gap*(colLower[j]+1.0) ;
		  rc.setUb(ub) ;
		  double effectiveness = sum2-ub ;
		  effectiveness =
		      CoinMax(effectiveness,
			      (sum-gap*colsol[j]-maxR[irow])/gap) ;
		  if (!coefficientExists)
		    effectiveness = CoinMax(1.0e-7,effectiveness) ;
		  rc.setEffectiveness(effectiveness) ;
		  //rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap) ;
		  if (rc.effectiveness() > needEffectiveness) {
		    rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
		    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		    // If strengthenRow point to row
		    //if(info->strengthenRow)
		    //printf("a point to row %d\n",irow) ;
#define STRENGTHEN_PRINT
#ifdef STRENGTHEN_PRINT
/*
Replace this code with the debugging print method from CglProbingAnn.
*/
		    if (canReplace && rowLower[irow] < -1.0e20) {
		      printf("1Cut %g <= ",rc.lb()) ;
		      int k ;
		      //printf("original row %d - %g <= <= %g - j = %d\n",iow,rowLower[irow],rowUpper[irow],j) ;
		      //for (int kk=rowStart[irow];kk<rowStart[irow+1];kk++) 
		      //printf("(%d,%g) ",column[kk],rowElements[kk]) ;
		      //printf("\n") ;
		      for ( k=0;k<n;k++) {
			int iColumn = index[k] ;
			printf("%g*",element[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rc.ub()) ;
		      printf("Row %g <= ",rowLower[irow]) ;
		      for (k=rowStart[irow];k<rowStart[irow+1];k++) {
			int iColumn = column[k] ;
			printf("%g*",rowElements[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rowUpper[irow]) ;
		    }
#endif
/*
Can we simply replace the existing constraint?

realRows comes in as a parameter. This is the translation array created if
we modified the constraint system during preprocessing in gutsOfGenerateCuts.

TODO: Seems a bit disingenuous to fail at replacement now, given that
effectiveness is artificially low. Notice again the inconsistent use of
finite infinities on the row lb.  -- lh, 101128 --
*/
		    int realRow =
		      (canReplace && rowLower[irow] < -1.0e20)?(irow):(-1) ;
		    if (realRows && realRow >= 0)
		      realRow = realRows[realRow] ;
		    if (!justReplace) {
		      rowCut.addCutIfNotDuplicate(rc,realRow) ;
		    } else if (realRow >= 0) {
		      double effectiveness = 0.0 ;
		      for (int i = 0 ; i < n ; i++)
			effectiveness += fabs(element[i]) ;
		      if (!info->strengthenRow[realRow] ||
			  info->strengthenRow[realRow]->effectiveness() > effectiveness) {
			delete info->strengthenRow[realRow] ;
			rc.setEffectiveness(effectiveness) ;
			info->strengthenRow[realRow]=rc.clone() ;
		      }
		    }
		  }
		}
	      }
/*
And repeat working against the lower bound L(i). As in so many other places
in this code, it's the identical functionality, with some subtle edits that
distinguish the L(i) math from the U(i) math.
*/
	      gap = minR[irow]-rowLower[irow] ;
	      if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		// see if the strengthened row is a cut
		sum =0.0 ;
		// also see if singletons can go to good objective
#ifdef MOVE_SINGLETONS
		bool moveSingletons=true ;
#endif
		for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		     kk++) {
		  int iColumn = column[kk] ;
		  double value = rowElements[kk] ;
		  sum += value*colsol[iColumn] ;
#ifdef MOVE_SINGLETONS
		  if (moveSingletons&&j!=iColumn) {
		    if (colUpper[iColumn]>colLower[iColumn]) {
		      if (columnLength2[iColumn]!=1) {
			moveSingletons=false ;
		      }
		    }
		  }
#endif
		}
#ifdef MOVE_SINGLETONS
		if (moveSingletons) {
		  // can fix any with good costs
		  for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		       kk++) {
		    int iColumn = column[kk] ;
		    if (j!=iColumn) {
		      if (colUpper[iColumn]>colLower[iColumn]) {
			if (columnLength2[iColumn]==1) {
			  double value = rowElements[kk] ;
			  if
			  (direction*objective[iColumn]*value>0.0&&!(markC[iColumn]&(tightenLower|tightenUpper))) {
			    // Fix
			    if (nstackC0+1<maxStack) {
			      stackC0[nstackC0]=iColumn ;
			      if (value<0.0) {
				lo0[nstackC0]=colUpper[iColumn] ;
				up0[nstackC0]=colUpper[iColumn] ;
			      } else {
				lo0[nstackC0]=colLower[iColumn] ;
				up0[nstackC0]=colLower[iColumn] ;
			      }
			      nstackC0++ ;
			    }
			  }
			}
		      }
		    }
		  }
		}
#endif
		if (sum+gap*colsol[j]<minR[irow]-primalTolerance_||(info->strengthenRow&&rowUpper[irow]>1.0e20)) {
		  // can be a cut
		  // add gap to lower and integer coefficient
		  // saveU and saveL spare
		  int * index = reinterpret_cast<int *>(saveL) ;
		  double * element = saveU ;
		  int n=0 ;
		  bool coefficientExists=false ;
		  double sum2=0.0 ;
		  for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		       kk++) {
		    int kColumn = column[kk] ;
		    double el = rowElements[kk] ;
		    if (kColumn!=j) {
		      index[n]=kColumn ;
		      element[n++]=el ;
		    } else {
		      el=el+gap ;
		      if (fabs(el)>1.0e-12) {
			index[n]=kColumn ;
			element[n++]=el ;
		      }
		      coefficientExists=true ;
		    }
		    sum2 += colsol[kColumn]*el ;
		  }
		  if (!coefficientExists) {
		    index[n]=j ;
		    element[n++]=gap ;
		    sum2 += colsol[j]*gap ;
		  }
		  OsiRowCut rc ;
		  double lb = rowLower[irow]+gap*(colLower[j]+1.0) ;
		  rc.setLb(lb) ;
		  rc.setUb(DBL_MAX) ;
		  // effectiveness
		  double effectiveness=lb-sum2 ;
		  effectiveness = CoinMax(effectiveness,
					  (minR[irow]-
					   sum-gap*colsol[j])/gap) ;
		  if (!coefficientExists)
		    effectiveness=CoinMax(1.0e-7,
					  effectiveness) ;
		  rc.setEffectiveness(effectiveness) ;
		  if (rc.effectiveness()>needEffectiveness) {
		    rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
		    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		    //if(info->strengthenRow)
		    //printf("b point to row %d\n",irow) ;
#ifdef STRENGTHEN_PRINT
		    if (canReplace&&rowUpper[irow]>1.0e20) {
		      printf("2Cut %g <= ",rc.lb()) ;
		      int k ;
		      for ( k=0;k<n;k++) {
			int iColumn = index[k] ;
			printf("%g*",element[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rc.ub()) ;
		      printf("Row %g <= ",rowLower[irow]) ;
		      for (k=rowStart[irow];k<rowStart[irow+1];k++) {
			int iColumn = column[k] ;
			printf("%g*",rowElements[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rowUpper[irow]) ;
		    }
#endif
		    int realRow = (canReplace&&rowUpper[irow]>1.0e20) ? irow : -1 ;
		    if (realRows&&realRow>=0)
		      realRow=realRows[realRow] ;
		    if (!justReplace) {
		      rowCut.addCutIfNotDuplicate(rc,realRow) ;
		    } else if (realRow>=0) {
		      double effectiveness=0.0 ;
		      for (int i=0;i<n;i++)
			effectiveness+=fabs(element[i]) ;
		      if (!info->strengthenRow[realRow]||info->strengthenRow[realRow]->effectiveness()>effectiveness) {
			delete info->strengthenRow[realRow] ;
			rc.setEffectiveness(effectiveness) ;
			info->strengthenRow[realRow]=rc.clone() ;
		      }
		    }
		  }
		}
	      }
	    }    // end generate coefficient strengthening cuts (1)
/*
Restore row bounds.
*/
	    minR[irow] = saveMin[istackR] ;
	    maxR[irow] = saveMax[istackR] ;
	    markR[irow] = -1 ;
	  }
	}    // end of code to iterate after first down pass
/*
PROBE LOOP: DOWN AND UP FEASIBLE

Again, not quite the full story. Cases remaining:
b) End of up probe, up and down feasible; in this case we will not
   iterate.
c) End of up probe, down feasible, up infeasible; in this case we will
   iterate to restore the down feasible state.

The next block deals with case b). We don't want to iterate the
down/up/down loop, so force iway to oneProbeTooMany. As usual, the move
singleton code consists of two nearly identical blocks, one working off
U(i), the next off L(i).

TODO: The conditions guarding the move singleton code here are different than
      those for pass 0: here there's an added guard using strengthenRow.
      Just a failure to change all instances? Nor do I see how John's
      comment (`just at root') is enforced here, unless it's from outside
      via rowCuts or strengthenRow.

      The check for any column cuts is not present, so presumably we can
      do singleton motion in the presence of column cuts. The check for gap
      does not have the high end (gap < 1.0e8) test.

      Note also that here the test for cut generation is moved outside the
      loop that iterates over rows. Presumably that's because in the
      previous case, the loop also restored markR, etc.

      Note that the code that handles the bound change is considerably
      different than the previous instance in pass 0. We're adding the
      changes to stackC rather than stackC0.
      -- lh, 101128 --
*/
	  else {
	  if (iway == upProbe &&
	      feasRecord == (downProbeFeas|upProbeFeas)) {
	    iway = oneProbeTooMany ;
#ifdef MOVE_SINGLETONS
	    // look for singletons that can move (just at root)
	    if ((rowCuts&2) != 0 &&
		 goingToTrueBound && info->strengthenRow) {
	      for (istackR = 0 ; istackR < nstackR ; istackR++) {
		int irow = stackR[istackR] ;
		double gap = rowUpper[irow]-maxR[irow] ;
		if (gap > primalTolerance_) {
		  // also see if singletons can go to good objective
		  bool moveSingletons = true ;
		  for (int kk = rowStart[irow] ;
		       kk < rowStart[irow+1] ; kk++) {
		    int iColumn = column[kk] ;
		    if (moveSingletons && j != iColumn) {
		      if (colUpper[iColumn] > colLower[iColumn]) {
			if (columnLength2[iColumn] != 1) {
			  moveSingletons = false ;
			}
		      }
		    }
		  }
		  if (moveSingletons) {
		    // can fix any with good costs
		    for (int kk = rowStart[irow] ;
			 kk < rowStart[irow+1] ; kk++) {
		      int iColumn = column[kk] ;
		      if (j != iColumn) {
			if (colUpper[iColumn] > colLower[iColumn]) {
			  if (columnLength2[iColumn] == 1) {
			    double value = rowElements[kk] ;
			    if (direction*objective[iColumn]*value < 0.0 &&
				!(markC[iColumn]&(tightenLower|tightenUpper))) {
			      stackC[nstackC] = iColumn ;
			      saveL[nstackC] = colLower[iColumn] ;
			      saveU[nstackC] = colUpper[iColumn] ;
			      assert(saveU[nstackC] > saveL[nstackC]) ;
			      if (value > 0.0) {
				colLower[iColumn] = colUpper[iColumn] ;
			      } else {
				colUpper[iColumn] = colLower[iColumn] ;
			      }
			      columnGap[iColumn] = -primalTolerance_ ;
			      assert (nstackC < nCols) ;
			      nstackC++ ;
			    }
			  }
			}
		      }
		    }
		  }
		}
/*
Repeat move singleton checks, against L(i).
*/
		gap = minR[irow]-rowLower[irow] ;
		if (gap > primalTolerance_) {
		  // also see if singletons can go to good objective
		  bool moveSingletons=true ;
		  for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		       kk++) {
		    int iColumn = column[kk] ;
		    if (moveSingletons&&j!=iColumn) {
		      if (colUpper[iColumn]>colLower[iColumn]) {
			if (columnLength2[iColumn]!=1) {
			  moveSingletons=false ;
			}
		      }
		    }
		  }
		  if (moveSingletons) {
		    // can fix any with good costs
		    for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		       kk++) {
		      int iColumn = column[kk] ;
		      if (j!=iColumn) {
			if (colUpper[iColumn]>colLower[iColumn]) {
			  if (columnLength2[iColumn]==1) {
			    double value = rowElements[kk] ;
			    if
			    (direction*objective[iColumn]*value>0.0&&!(markC[iColumn]&(tightenLower|tightenUpper))) {
			      // Fix
			      stackC[nstackC]=iColumn ;
			      saveL[nstackC]=colLower[iColumn] ;
			      saveU[nstackC]=colUpper[iColumn] ;
			      assert (saveU[nstackC]>saveL[nstackC]) ;
			      if (value<0.0) {
				colLower[iColumn]=colUpper[iColumn] ;
			      } else {
				colUpper[iColumn]=colLower[iColumn] ;
			      }
			      columnGap[iColumn] = -primalTolerance_ ;
			      assert (nstackC<nCols) ;
			      nstackC++ ;
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
#endif
/*
jjf: point back to stack

We're done playing with singletons. Get to the real work here. We don't need
markC any more; coopt it for a cross-reference, var index -> stackC index.
The +1000 offset allows us to distinguish invalid entries (not all variables
are on stackC).
*/
	    for (istackC = nstackC-1 ; istackC >= 0 ; istackC--) {
	      int icol = stackC[istackC] ;
	      markC[icol] = istackC+1000 ;
	    }
	    OsiColCut cc ;
	    int nTot = 0, nFix = 0, nInt = 0 ;
	    bool ifCut = false ;
/*
See if we have bounds improvement. Given a lower bound ld<j> from the
down probe, a bound lu<j> from the up probe, and the original bound lo<j>,
the bound we want is max(min(ld<j>,lu<j>),lo<j>). Use a conservative
tolerance.
*/
	    for (istackC = 1 ; istackC < nstackC0 ; istackC++) {
	      int icol = stackC0[istackC] ;
	      int istackC1 = markC[icol]-1000 ;
	      if (istackC1 >= 0) {
		if (CoinMin(lo0[istackC],colLower[icol]) >
					      saveL[istackC1]+1.0e-4) {
		  saveL[istackC1] = CoinMin(lo0[istackC],colLower[icol]) ;
		  if (intVar[icol] /*||!info->inTree*/ ) {
		    element[nFix] = saveL[istackC1] ;
		    index[nFix++] = icol ;
		    nInt++ ;
		    if (colsol[icol] < saveL[istackC1]-primalTolerance_)
		      ifCut = true ;
		  }
		  nfixed++ ;
		} 
	      }
	    }
	    if (nFix) {
	      nTot = nFix ;
	      cc.setLbs(nFix,index,element) ;
	      nFix = 0 ;
	    }
/*
Repeat for upper bounds. There's an odd bit of asymmetry here; this loop
tries to generate a cut if there's no bound improvement.

TODO: What, pray tell, is a `two cut'?  Binary variables only, it seems.
      Judging from the conditions, we're looking for a variable that's
      fixed at 0 on a down probe and fixed at 1 on an up probe, or vice
      versa.
      -- lh, 101128 --
*/
	    for (istackC = 1 ; istackC < nstackC0 ; istackC++) {
	      int icol = stackC0[istackC] ;
	      int istackC1 = markC[icol]-1000 ;
	      if (istackC1 >= 0) {
		if (CoinMax(up0[istackC],colUpper[icol]) <
					      saveU[istackC1]-1.0e-4) {
		  saveU[istackC1] = CoinMax(up0[istackC],colUpper[icol]) ;
		  if (intVar[icol] /*||!info->inTree*/ ) {
		    element[nFix] = saveU[istackC1] ;
		    index[nFix++] = icol ;
		    nInt++ ;
		    if (colsol[icol] > saveU[istackC1]+primalTolerance_)
		      ifCut = true ;
		  }
		  nfixed++ ;
		} else if (!info->inTree &&
			   saveL[0] == 0.0 && saveU[0] == 1.0) {
		  // See if can do two cut
		  double upperWhenDown = up0[istackC] ;
		  double lowerWhenDown = lo0[istackC] ;
		  double upperWhenUp = colUpper[icol] ;
		  double lowerWhenUp = colLower[icol] ;
		  double upperOriginal = saveU[istackC1] ;
		  double lowerOriginal = saveL[istackC1] ;
		  if (upperWhenDown < lowerOriginal+1.0e-12 &&
		      lowerWhenUp > upperOriginal-1.0e-12) {
		    OsiRowCut rc ;
		    rc.setLb(lowerOriginal) ;
		    rc.setUb(lowerOriginal) ;
		    rc.setEffectiveness(1.0e-5) ;
		    int index[2] ;
		    double element[2] ;
		    index[0] = j ;
		    index[1] = icol ;
		    element[0] = -(upperOriginal-lowerOriginal) ;
/*
jjf: If zero then - must have been fixed without noticing!

TODO: Uh, fixed without noticing?! That's an algorithm problem.
      -- lh, 101128 --
*/
		    if (fabs(element[0]) > 1.0e-8) {
		      element[1] = 1.0 ;
		      rc.setRow(2,index,element,false) ;
		      cs.insert(rc) ;
		    }
		  } else if (upperWhenUp < lowerOriginal+1.0e-12 &&
			     lowerWhenDown > upperOriginal-1.0e-12) {
		    OsiRowCut rc ;
		    rc.setLb(upperOriginal) ;
		    rc.setUb(upperOriginal) ;
		    rc.setEffectiveness(1.0e-5) ;
		    int index[2] ;
		    double element[2] ;
		    index[0] = j ;
		    index[1] = icol ;
		    element[0] = upperOriginal-lowerOriginal ;
		    element[1] = 1.0 ;
		    rc.setRow(2,index,element,false) ;
		    cs.insert(rc) ;
		  } 
		}
	      }
	    }    // end loop to update upper bounds (w/ two cuts)
	    if (nFix) {
	      nTot+=nFix ;
	      cc.setUbs(nFix,index,element) ;
	    }
	    // could tighten continuous as well
	    if (nInt) {
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
	  }    // end of code to handle up and down feasible.
/*
  PROBE LOOP: DOWN FEASIBLE, UP INFEASIBLE

  The only remaining case is down feasible, up infeasible, in which case we
  need to reset to the initial state and run the down probe again. Surely
  there must be a better way?
*/
	    else {
	    goingToTrueBound = 0 ;
	  }
/*
  PROBE LOOP: RESTORE

  And now the code to reset to the initial state. (But it seems to be doing a
  bit more than that. Looks like some of the actions are related to
  end-of-pass processing for the up probe, similar to end-of-pass processing
  for the down probe.
*/
	  if (goingToTrueBound == 2 && !justReplace) {
	    disaggCuts(nstackC,probeUp,primalTolerance_,
		       disaggEffectiveness,si,
		       rowCut,stackC,colsol,colUpper,colLower,saveU,saveL,
		       index,element) ;
	  }
	  /* restore all */
	  for (istackC = nstackC-1 ; istackC >= 0 ; istackC--) {
	    int icol = stackC[istackC] ;
	    double oldU = saveU[istackC] ;
	    double oldL = saveL[istackC] ;
/*
TODO: This next bit differs in subtle ways from the restore at the end
      of pass 0. Here we force markC to show fixed of the bounds are
      within 1.0e-4 of one another, a fairly broad tolerance; code for
      the pass 0 restore does not restore fixed status.  Also,
      we're not distinguishing here between actions for down feasible,
      up feasible, vs. down feasible, up infeasible.  -- lh, 101128 --

*/
	    colUpper[icol] = oldU ;
	    colLower[icol] = oldL ;
	    columnGap[icol] = oldU-oldL-primalTolerance_ ;
	    if (oldU > oldL+1.0e-4) {
	      markC[icol] = 0 ;
	      if (oldU > 1.0e10)
		markC[icol] |= infUpper ;
	      if (oldL < -1.0e10)
		markC[icol] |= infLower ;
	    } else {
	      markC[icol] = (tightenUpper|tightenLower) ;
	    }
	  }
/*
End of column restore. On to rows. Looks much like the end of pass 0 except
that the singleton stuff was done up in the down feasible, up feasible
block. Presumably we don't want to do it in the down feasible, up infeasible
case. And we won't generate cuts for up infeasible because goingToTrueBound
has been forced to zero.

TODO: Note that here we run a separate loop to calculate sum, because we
      can't piggyback on the singleton loop.  -- lh, 101128 --

TODO: Again, subtle but important differences between this code and the code
      that follows pass 0.   -- lh, 101128 --

TODO: Recalculation of canReplace seems completely pointless here.
      -- lh, 101128 --

begin generate coefficient strengthening cuts (2)
*/
	  for (istackR = 0 ; istackR < nstackR ; istackR++) {
	    int irow = stackR[istackR] ;
	    // switch off strengthening if not wanted
	    if ((rowCuts&2) != 0 && goingToTrueBound) {
	      bool canReplace = info->strengthenRow && (goingToTrueBound==2) ;
	      bool ifCut = anyColumnCuts ;
	      double gap = rowUpper[irow]-maxR[irow] ;
	      double sum = 0.0 ;
	      if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
		// see if the strengthened row is a cut
		for (int kk = rowStart[irow] ;
		     kk < rowStart[irow+1] ; kk++) {
		  sum += rowElements[kk]*colsol[column[kk]] ;
		}
/*
Same commentary as previous code in pass 0 finishing block.

TODO: This condition differs from the pass 0 condition. canReplace replaces
      a simple info->strengthenRow.  -- lh, 101128 --
*/
		if (sum+gap*colsol[j] > rowUpper[irow]+primalTolerance_ ||
		    (canReplace && rowLower[irow]<-1.e20)) {
		  int * index = reinterpret_cast<int *>(saveL) ;
		  double * element = saveU ;
		  int n=0 ;
		  bool coefficientExists=false ;
		  double sum2=0.0 ;
		  for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		       kk++) {
		    int kColumn = column[kk] ;
		    double el = rowElements[kk] ;
		    if (kColumn!=j) {
		      index[n]=kColumn ;
		      element[n++]=el ;
		    } else {
		      el=el+gap ;
		      if (fabs(el)>1.0e-12) {
			index[n]=kColumn ;
			element[n++]=el ;
		      }
		      coefficientExists=true ;
		    }
		    sum2 += colsol[kColumn]*el ;
		  }
		  if (!coefficientExists) {
		    index[n]=j ;
		    element[n++]=gap ;
		    sum2 += colsol[j]*gap ;
		  }
		  OsiRowCut rc ;
		  rc.setLb(-DBL_MAX) ;
		  double ub = rowUpper[irow]+gap*(colUpper[j]-1.0) ;
		  rc.setUb(ub) ;
		  // effectiveness
		  double effectiveness=sum2-ub ;
		  effectiveness = CoinMax(effectiveness,
					  (sum+gap*colsol[j]-
					   rowUpper[irow])/gap) ;
		  if (!coefficientExists)
		    effectiveness=CoinMax(1.0e-7,
					  effectiveness) ;
		  rc.setEffectiveness(effectiveness) ;
		  if (rc.effectiveness()>needEffectiveness) {
		    rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
		    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		    //if(canReplace)
		    //printf("c point to row %d\n",irow) ;
#ifdef STRENGTHEN_PRINT
		    if (canReplace&&rowLower[irow]<-1.0e20) {
		      printf("3Cut %g <= ",rc.lb()) ;
		      int k ;
		      for ( k=0;k<n;k++) {
			int iColumn = index[k] ;
			printf("%g*",element[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rc.ub()) ;
		      dump_row(irow,rc.lb(),rc.ub(),
		      	       nan(""),nan(""),&si,true,false,4,
			       n,index,element,
			       1.0e-10,colLower,colUpper) ;
		      printf("Row %g <= ",rowLower[irow]) ;
		      for (k=rowStart[irow];k<rowStart[irow+1];k++) {
			int iColumn = column[k] ;
			printf("%g*",rowElements[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rowUpper[irow]) ;
		      k = rowStart[irow] ;
		      dump_row(irow,rowLower[irow],rowUpper[irow],
		      	       minR[irow],maxR[irow],&si,true,false,4,
			       rowStart[irow+1]-k,&column[k],&rowElements[k],
			       1.0e-10,colLower,colUpper) ;
		    }
#endif
		    int realRow = (canReplace&&rowLower[irow]<-1.0e20) ? irow : -1 ;
		    if (realRows&&realRow>=0)
		      realRow=realRows[realRow] ;
		    if (!justReplace) {
		      rowCut.addCutIfNotDuplicate(rc,realRow) ;
		    } else if (realRow>=0) {
		      double effectiveness=0.0 ;
		      for (int i=0;i<n;i++)
			effectiveness+=fabs(element[i]) ;
		      if (!info->strengthenRow[realRow]||info->strengthenRow[realRow]->effectiveness()>effectiveness) {
			delete info->strengthenRow[realRow] ;
			rc.setEffectiveness(effectiveness) ;
			info->strengthenRow[realRow]=rc.clone() ;
		      }
		    }
		  }
		}
	      }
/*
As with pass 0 code, repeat working against row lower bound.
*/
	      gap = minR[irow]-rowLower[irow] ;
	      if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		// see if the strengthened row is a cut
		if (!sum) {
		  for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		       kk++) {
		    sum += rowElements[kk]*colsol[column[kk]] ;
		  }
		}
		if (sum-gap*colsol[j]<rowLower[irow]-primalTolerance_||(canReplace&&rowUpper[irow]>1.0e20)) {
		  // can be a cut
		  // subtract gap from integer coefficient
		  // saveU and saveL spare
		  int * index = reinterpret_cast<int *>(saveL) ;
		  double * element = saveU ;
		  int n=0 ;
		  bool coefficientExists=false ;
		  double sum2=0.0 ;
		  for (int kk =rowStart[irow];kk<rowStart[irow+1] ;
		       kk++) {
		    int kColumn = column[kk] ;
		    double el = rowElements[kk] ;
		    if (kColumn!=j) {
		      index[n]=kColumn ;
		      element[n++]=el ;
		    } else {
		      el=el-gap ;
		      if (fabs(el)>1.0e-12) {
			index[n]=kColumn ;
			element[n++]=el ;
		      }
		      coefficientExists=true ;
		    }
		    sum2 += colsol[kColumn]*el ;
		  }
		  if (!coefficientExists) {
		    index[n]=j ;
		    element[n++]=-gap ;
		    sum2 -= colsol[j]*gap ;
		  }
		  OsiRowCut rc ;
		  double lb = rowLower[irow]-gap*(colUpper[j]-1) ;
		  rc.setLb(lb) ;
		  rc.setUb(DBL_MAX) ;
		  double effectiveness=lb-sum2 ;
		  effectiveness = CoinMax(effectiveness,
					  (rowLower[irow]-
					   sum+gap*colsol[j])/gap) ;
		  if (!coefficientExists)
		    effectiveness=CoinMax(1.0e-7,
					  effectiveness) ;
		  rc.setEffectiveness(effectiveness) ;
		  if (rc.effectiveness()>needEffectiveness) {
		    rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
		    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		    //if(canReplace)
		    //printf("d point to row %d\n",irow) ;
#ifdef STRENGTHEN_PRINT
		    if (canReplace&&rowUpper[irow]>1.0e20) {
		      printf("4Cut %g <= ",rc.lb()) ;
		      int k ;
		      for ( k=0;k<n;k++) {
			int iColumn = index[k] ;
			printf("%g*",element[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rc.ub()) ;
		      printf("Row %g <= ",rowLower[irow]) ;
		      for (k=rowStart[irow];k<rowStart[irow+1];k++) {
			int iColumn = column[k] ;
			printf("%g*",rowElements[k]) ;
			if (si.isInteger(iColumn))
			  printf("i%d ",iColumn) ;
			else
			  printf("x%d ",iColumn) ;
		      }
		      printf("<= %g\n",rowUpper[irow]) ;
		    }
#endif
		    int realRow = (canReplace&&rowUpper[irow]>1.0e20) ? irow : -1 ;
		    if (realRows&&realRow>=0)
		      realRow=realRows[realRow] ;
		    if (!justReplace) {
		      rowCut.addCutIfNotDuplicate(rc,realRow) ;
		    } else if (realRow>=0) {
		      double effectiveness=0.0 ;
		      for (int i=0;i<n;i++)
			effectiveness+=fabs(element[i]) ;
		      if (!info->strengthenRow[realRow]||info->strengthenRow[realRow]->effectiveness()>effectiveness) {
			delete info->strengthenRow[realRow] ;
			rc.setEffectiveness(effectiveness) ;
			info->strengthenRow[realRow]=rc.clone() ;
		      }
		    }
		  }
		}
	      }
	    }    // end generate coefficient strengthening cuts (2)
	    minR[irow]=saveMin[istackR] ;
	    maxR[irow]=saveMax[istackR] ;
	    markR[irow]=-1 ;
	  }
	}       // end of PROBE: DOWN AND UP FEASIBLE
      }     // PROBE LOOP: END
    }    // LOOKEDAT LOOP: END
  }   // PASS LOOP: END
/*
  Free up the space we've been using.
*/
#ifndef ONE_ARRAY
  delete [] stackC0 ;
  delete [] lo0 ;
  delete [] up0 ;
  delete [] columnGap ;
  delete [] markC ;
  delete [] stackC ;
  delete [] stackR ;
  delete [] saveL ;
  delete [] saveU ;
  delete [] saveMin ;
  delete [] saveMax ;
  delete [] index ;
  delete [] element ;
  delete [] djs ;
  delete [] largestPositiveInRow ;
  delete [] largestNegativeInRow ;
#endif
  delete [] colsol ;
/*
  jjf:  Add in row cuts

  Transfer row cuts from the local container into something that can escape
  this method. If we're not in `just replace' mode, simply copy them over to
  the container (cs) passed as a parameter using addCuts. If we're in `just
  replace' mode, ? need to explore further ?  -- lh, 101125 --

  TODO: If there's any possibility of realRows[i] < 0, this code will break
  	trying to read strengthenRow[]!
        -- lh, 101128 --
*/
  if (!ninfeas) {
    if (!justReplace) {
      rowCut.addCuts(cs,info->strengthenRow,info->pass) ;
    } else {
      for (int i = 0 ; i < nRows ; i++) {
	int realRow = realRows[i] ;
	OsiRowCut *cut = info->strengthenRow[realRow] ;
	if (realRow >= 0 && cut) {
#ifdef CLP_INVESTIGATE
	  printf("Row %d, real row %d effectiveness %g\n",
	  	 i,realRow,cut->effectiveness()) ;
#endif
	  cs.insert(cut) ;
	}
      }
    }
  }
#if 0
/*
  Debug print.
*/
  {
    int numberRowCutsAfter = cs.sizeRowCuts() ;
    int k ;
    for (k = 0;k<numberRowCutsAfter;k++) {
      OsiRowCut thisCut = cs.rowCut(k) ;
      printf("Cut %d is %g <=",k,thisCut.lb()) ;
      int n=thisCut.row().getNumElements() ;
      const int * column = thisCut.row().getIndices() ;
      const double * element = thisCut.row().getElements() ;
      assert (n) ;
      for (int i=0;i<n;i++) {
	printf(" %g*x%d",element[i],column[i]) ;
      }
      printf(" <= %g\n",thisCut.ub()) ;
    }
  }
#endif

# if CGL_DEBUG > 0
  std::cout
    << "Leaving CglProbing::probe, ninfeas " << ninfeas
    << "." << std::endl ;
# endif

  return (ninfeas) ;
}
