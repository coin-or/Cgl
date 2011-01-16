
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
void disaggCuts (int nstackC, unsigned int probeDir,
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



/*
  This method examines the rows on stackR looking for redundant rows that
  consist solely of singleton variables (i.e., variables that appear in just
  this constraint). It then looks to see if any of the variables in the row
  can be fixed at bound.

  Consider the case where U<i> < b<i> and all unfixed variables in the row
  are singletons (occur in no other constraint). Given x<j> is a singleton,
  the reduced cost is c<j> - ya<j> = c<j> - ye<i> = c<j> - y<i>.  But if
  U<i> < b<i>, the constraint is loose by definition, hence y<i> = 0 and
  the reduced cost is c<j>. If a<ij> > 0, x<j> should be u<j> in U<i>. If
  c<j> < 0, minimising will tend to drive x<j> to u<j>. This cannot violate
  constraint i (U<i> < b<i>) and (since x<j> is a singleton) will have no
  effect elsewhere. Hence we can fix x<j> at u<j>.

  Do the case analysis, and what you find is that against U<i> we want
    a<ij> > 0  ==>  c<j> < 0  ==>  push x<j> to u<j>
    a<ij> < 0  ==>  c<j> > 0  ==>  push x<j> to l<j>
  and against L<i> we want
    a<ij> > 0  ==>  c<j> > 0  ==>  push x<j> to l<j>
    a<ij> < 0  ==>  c<j> < 0  ==>  push x<j> to u<j>

  Extend it one more time by taking the objective direction (max/min) into
  account (dir = 1.0 for min, -1.0 for max) and you have
    against U<i> ==> a<ij>*c<j>*dir < 0
    against L<i> ==> a<ij>*c<j>*dir > 0

  John's original comment for this code was
    // also see if singletons can go to good objective
    // Taken out as should be found elsewhere
    // and has to be original column length
  but he reinstated it. Apparently there were cases where fixing the probing
  variable was required to satisfy the condition that all unfixed variables be
  singletons. Enough of them to justify reinstating this code.
*/

void singletonRows (int jProbe, double primalTolerance_,
		    const OsiSolverInterface &si,
		    const CoinPackedMatrix *rowCopy,
		    int *const markC,
		    int &nstackC, int *const stackC,
		    double *const saveL, double *const saveU,
		    double *const colUpper, double *const colLower,
		    double *const colGap,
		    const int nstackR, const int *const stackR,
		    const double *const rowUpper,
		    const double *const rowLower,
		    const double *const maxR, const double *const minR)
{
/*
  Unpack a few vectors from the row-major matrix.
*/
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *column = rowCopy->getIndices() ;
  const double *rowElements = rowCopy->getElements() ;
  const int nCols = rowCopy->getNumCols() ;
/*
  `Singleton' must be based on the column length in the original system.
*/
  const double *objective = si.getObjCoefficients() ;
  const int *columnLengths = si.getMatrixByCol()->getVectorLengths() ;
  const double objSense = si.getObjSense() ;
/*
  Open a loop to work through the rows on stackR.
*/
  for (int istackR = 0 ; istackR < nstackR ; istackR++) {
    int i = stackR[istackR] ;
/*
  Check the gaps. If the constraint is potentially tight in both directions,
  there's nothing more to do here.
*/
    const double uGap = rowUpper[i]-maxR[i] ;
    const double lGap = minR[i]-rowLower[i] ;
    if (uGap < primalTolerance_ && lGap < primalTolerance_) continue ;
/*
  Scan the row to see if it meets the `all singletons' condition. Again,
  if this fails, there's nothing more to be done.

  Note that the original code didn't check the probing variable x<p>,
  because if you're probing a binary variable it's fixed. But for general
  integer variables a probe does not in general fix the variable.  So we
  check all variables.

  We should not be executing this method if we have prima facie infeasibility.
*/
    bool canFix = true ;
    for (int kk = rowStart[i] ; kk < rowStart[i+1] ; kk++) {
      int j = column[kk] ;
      assert(colUpper[j]-colLower[j] > -primalTolerance_) ;
      if (colUpper[j] > colLower[j] && columnLengths[j] != 1) {
	canFix = false ;
	break ;
      }
    }
    if (!canFix) continue ;
/*
  If we've passed the tests, we can look for variables suitable to drive to
  bound. Work against U<i> first. We're looking for variables with a<ij> > 0
  that will be naturally driven to u<j>, or variables with a<ij> < 0 that will
  be naturally driven to l<j>.
  
  Don't overwrite the saved bounds if we've
  tightened this variable already!
*/
    if (uGap > primalTolerance_) {
      for (int kk = rowStart[i] ; kk < rowStart[i+1] ; kk++) {
	int j = column[kk] ;
	const double lj = colLower[j] ;
	const double uj = colUpper[j] ;
	if (uj > lj) {
	  double value = rowElements[kk] ;
	  if (objSense*objective[j]*value < 0.0)
	  { assert(jProbe != j) ;
	    if (!(markC[j]&(tightenLower|tightenUpper))) {
	      stackC[nstackC] = j ;
	      saveL[nstackC] = lj ;
	      saveU[nstackC] = uj ;
	    }
	    assert(saveU[nstackC] > saveL[nstackC]) ;
	    if (value > 0.0) {
	      colLower[j] = uj ;
	    } else {
	      colUpper[j] = lj ;
	    }
	    markC[j] |= tightenLower|tightenUpper ;
	    colGap[j] = -primalTolerance_ ;
	    assert(nstackC < nCols) ;
	    nstackC++ ;
	  }
	}
      }
    }
/*
  And now the identical code, except that we're working against L<i>, hence
  the sense of the elibigility test is reversed and we want variables with
  a<ij> > 0 that will be naturally driven to l<j>, or variables with
  a<ij> < 0 that will be naturally driven to u<j>.
*/
    if (lGap > primalTolerance_) {
      for (int kk = rowStart[i] ; kk < rowStart[i+1] ; kk++) {
	int j = column[kk] ;
	const double lj = colLower[j] ;
	const double uj = colUpper[j] ;
	if (uj > lj) {
	  double value = rowElements[kk] ;
	  if (objSense*objective[j]*value > 0.0)
	  { assert(jProbe != j) ;
	    if (!(markC[j]&(tightenLower|tightenUpper))) {
	      stackC[nstackC] = j ;
	      saveL[nstackC] = lj ;
	      saveU[nstackC] = uj ;
	    }
	    assert(saveU[nstackC] > saveL[nstackC]) ;
	    if (value < 0.0) {
	      colLower[j] = uj ;
	    } else {
	      colUpper[j] = lj ;
	    }
	    markC[j] |= tightenLower|tightenUpper ;
	    colGap[j] = -primalTolerance_ ;
	    assert(nstackC < nCols) ;
	    nstackC++ ;
	  }
	}
      }
    }
  }
  return ;
}


/*
  Generate coefficient strengthening cuts.
  
  Assume we have a binary probe variable x<p> with a<ip> > 0. Assume a
  down probe.  Assume that U<i> > b<i> before the probe, but now U'<i> < b<i>
  (i.e., the constraint is redundant against b<i> after the down probe forces
  x<p> to 0 and reduces U<i> to U'<i>). We would like to have the following in
  a strengthened constraint a'<ip>x<p> + (stuff) <= b'<i>:
    * When x<p> = 0, b'<i> = U'<i>  (the lhs can't be larger than U'<i>)
    * When x<p> = 1, b'<i> = b<i>   (the original value)
  Define delta = b<i> - U'<i>, b'<i> = b<i>-delta, and  a'<ip> = a<ip>-delta.
  When x<p> = 1, the delta term on each side cancels and we're left with the
  original constraint. When x<p> = 0, the rhs tightens to
  b'<i> = b<i>+delta<i> = U'<i>.

  For an honest derivation that works for binary and general integer
  variables, see the typeset documentation. As you'd expect, there are four
  cases (a<ip> < 0, a<ip> > 0) X (up probe, down probe). Assume that a down
  probe bounds x<p> as [l<p>, u'<p>] and an up probe bounds x<p> as [l'<p>,
  u<p>]. Then the summary (omitting subscripts to fit the line) is:

  a<ip>  dir    delta             blow'       a'<ip>      b'

    >0    d    b - U'                         a-delta   b-(u'+1)delta
    >0    u    blow - L'   blow+(l'-1)delta   a+delta
    <0    d    b - L'      blow-(u'+1)delta   a-delta
    <0    u    blow - U'                      a+delta   b+(l'-1)delta

  In the code, negating delta for down probes unifies the up and down math.

  TODO: In the original code, this method will not execute if column cuts
	are present (i.e., bounds were tightened so that some x*<j> may not
	be feasible). The current code retains those guards.  Clearly it's
	a problem if x*<p> isn't feasible --- we can't generate a useful
	cut over floor(x*<p>), ceil(x*<p>) if the interval is outside the
	polytope. It's not so clear this is an obstacle for other variables,
	except that the lhs activity may be inaccurate.
*/

#define STRENGTHEN_PRINT
void strengthenCoeff (
		      int jProbe, unsigned int probeDir,
		      double primalTolerance_,
		      bool justReplace, bool canReplace,
		      double needEffectiveness,
		      const OsiSolverInterface &si,
		      CglProbingRowCut &rowCut,
		      const CoinPackedMatrix *rowCopy,
		      double *const colUpper, double *const colLower,
		      const double *const colsol,
		      const int nstackR, const int *const stackR,
		      const double *const rowUpper,
		      const double *const rowLower,
		      const double *const maxR,
		      const double *const minR,
		      const int *const realRows,
		      double *const element,
		      int *const index,
		      CglTreeInfo *const info
		     )
{
# if CGL_DEBUG > 0
  std::cout
    << "Entering strengthenCoeff, probing "
    << si.getColName(jProbe) << "(" << jProbe << ") "
    << ((probeDir == probeDown)?"down":"up") << "." << std::endl ;
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
# endif
/*
  Define magic numbers up front.
    * Limits on `reasonable' coefficients. We don't want coefficients
      outside this range.
    * `Infinity' for the rhs. This'll get fixed later.
    * Effectiveness in the case that a<ip> = 0 and a'<ip> != 0.
  All copied from the original code.
*/
  const double bigCoeff = 1.0e8 ;
  const double tinyCoeff = 1.0e-12 ;
  const double rhsInf = 1.0e20 ;
  const double aipZeroEffectiveness = 1.0e-7 ;
/*
  Unpack a few vectors from the row-major matrix.
*/
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *column = rowCopy->getIndices() ;
  const double *rowElements = rowCopy->getElements() ;
/*
  Open up an outer loop to walk stackR and look for interesting constraints.
*/
  for (int istackR = 0 ; istackR < nstackR ; istackR++) {
    int irow = stackR[istackR] ;
    double bi = rowUpper[irow] ;
    double blowi = rowLower[irow] ;
/*
  We can't get anywhere unless probing has made one end or the other of the
  constraint redundant. If there's no gap, we can't do anything.
*/
    double uGap = bi-maxR[irow] ;
    double lGap = blowi-minR[irow] ;
    if (uGap < primalTolerance_ && -lGap < primalTolerance_) continue ;
    bool isRangeCon = ((blowi  > -rhsInf) && (bi < rhsInf)) ;
/*
  We'll need the lhs activity, excluding the probed variable, to evaluate the
  effectiveness of the cut.  Extract the coefficient of the probe variable
  while we're here; we need it to determine which case we're working.

  TODO: As noted at the head of the routine, we could correct for column cuts
        by checking bounds when calculating the sum below, if it was worth the
	effort.   -- lh, 110113 --
*/
    double sum = 0.0 ;
    double aip = 0.0 ;
    bool aipNonZero = false ;
    for (CoinBigIndex kk = rowStart[irow] ; kk < rowStart[irow+1] ; kk++) {
      int k = column[kk] ;
      double aik = rowElements[kk] ;
      if (k == jProbe) {
	aipNonZero = true ;
        aip = aik ;
      } else {
	sum += aik*colsol[k] ;
      }
    }
/*
  Now figure out which case we're in and do some setup. At the end, see if the
  coefficient will have reasonable magnitude. If not, move on to the next
  iteration.
*/
    double delta ;
    double deltaMul ;
    bool revisebi = true ;
    if (probeDir == probeDown) {
      deltaMul = colUpper[jProbe]+1.0 ; 
      if (aip >= 0) {
        delta = -uGap ;
      } else {
        delta = -lGap ;
	revisebi = false ;
      }
    } else {
      deltaMul = colLower[jProbe]-1.0 ; 
      if (aip >= 0) {
        delta = lGap ;
	revisebi = false ;
      } else {
        delta = uGap ;
      }
    }
    if (CoinAbs(delta) < primalTolerance_ || CoinAbs(delta) > bigCoeff)
      continue ;
/*
  Now decide if we have something that cuts off the current solution. Augment
  the lhs activity with the contribution for a'<ip>x*<p>, calculate the new
  rhs value, and compare.
  As an alternate measure of effectiveness, consider the gap between the
  current activity and the revised lhs bound, normalised by the gap between
  the original rhs and the revised lhs bound. 

  We'll also generate the cut if we can strengthen it in place and we're not
  dealing with a range constraint. (The reason for excluding range constraints
  is that we might try to alter a<ip> against both b<i> and blow<i>. That way
  lies madness.)

  TODO: It's not clear to me how the first two effectiveness calculations
	relate to one another. Think more on this -- lh, 110115 --
*/
    double aipPrime = aip+delta ;
    bool aipPrimeNonZero = true ;
    if (CoinAbs(aipPrime) < tinyCoeff) {
      aipPrime = 0.0 ;
      aipPrimeNonZero = false ;
    }
    double xp = colsol[jProbe] ;
    sum += aipPrime*xp ;
    double biPrime = DBL_MAX ;
    double blowiPrime = -DBL_MAX ;
    double effectiveness = 0.0 ;
    bool genCut = false ;
    if (revisebi) {
      biPrime = rowUpper[irow]+deltaMul*delta ;
      effectiveness = sum-biPrime ;
      effectiveness = CoinMax(effectiveness,(sum-maxR[irow])/uGap) ;
    } else {
      blowiPrime = rowLower[irow]+deltaMul*delta ;
      effectiveness = blowiPrime-sum ;
      effectiveness = CoinMax(effectiveness,(minR[irow]-sum)/lGap) ;
    }
    if (!aipNonZero && aipPrimeNonZero)
      effectiveness = CoinMax(effectiveness,aipZeroEffectiveness) ;
    if (effectiveness > needEffectiveness) genCut = true ;
    if (info->strengthenRow && !isRangeCon) genCut = true ;
/*
  Are we going to generate the cut? If not, on to the next iteration.
*/
    if (!genCut) continue ;
/*
  Generate the coefficients. Copy whatever's present and plug in aipPrime at
  the end.
*/
    int n = 0 ;
    int aip_n = -1 ;
    for (CoinBigIndex kk = rowStart[irow] ; kk < rowStart[irow+1] ; kk++) {
      int k = column[kk] ;
      double aik = rowElements[kk] ;
      if (k == jProbe)
	aip_n = k ;
      index[n] = k ;
      element[n++] = aik ;
    }
    if (aip_n < 0) {
      aip_n = n ;
      n++ ;
    }
    index[aip_n] = jProbe ;
    element[aip_n] = aipPrime ;
/*
  Fill in a cut structure with the cut.
*/
    OsiRowCut rc ;
    rc.setLb(blowiPrime) ;
    rc.setUb(biPrime) ;
    rc.setEffectiveness(effectiveness) ;
    rc.setRow(n,index,element,false) ;
#   if CGL_DEBUG > 0
    if (debugger) assert(!debugger->invalidCut(rc)); 
#   endif
#   ifdef STRENGTHEN_PRINT
    if (canReplace && !isRangeCon) {
      printf("Strengthen Cut (1):\n") ;
      dump_row(irow,rc.lb(),rc.ub(),
	       nan(""),nan(""),&si,true,false,4,
	       n,index,element,
	       1.0e-10,colLower,colUpper) ;
      printf("Original Row:\n") ;
      int k = rowStart[irow] ;
      dump_row(irow,rowLower[irow],rowUpper[irow],
	       minR[irow],maxR[irow],&si,true,false,4,
	       rowStart[irow+1]-k,&column[k],&rowElements[k],
	       1.0e-10,colLower,colUpper) ;
    }
#   endif
/*
  If we're in preprocessing, we might try to simply replace the existing
  constraint (justReplace = true; preprocessing is the only place this will
  happen). Otherwise, drop the cut into the cut set.

  realRows comes in as a parameter. This is the translation array created if
  we modified the constraint system during preprocessing in gutsOfGenerateCuts.

  Effectiveness, if we're strengthening in place, seems to be absolute size of
  coefficients; smaller is better. (Why?)

  TODO: It seems to me that we can get here with
	justReplace = true and canReplace = false, and this condition is
	constant over an iteration of the way loop. In other words, we've done
	all the work to generate this cut and we knew before we started that
	we would discard it here.  -- lh, 110107 --
	Put in an assert to see if we ever do all the work, only to reject
	the cut.  -- lh, 110115 --
*/
      int realRow = (canReplace && !isRangeCon)?(irow):(-1) ;
      if (realRows && realRow >= 0)
	realRow = realRows[realRow] ;
      if (!justReplace) {
	rowCut.addCutIfNotDuplicate(rc,realRow) ;
      } else if (realRow >= 0) {
	double effectiveness = 0.0 ;
	for (int i = 0 ; i < n ; i++)
	  effectiveness += CoinAbs(element[i]) ;
	if (!info->strengthenRow[realRow] ||
	    info->strengthenRow[realRow]->effectiveness() > effectiveness) {
	  delete info->strengthenRow[realRow] ;
	  rc.setEffectiveness(effectiveness) ;
	  info->strengthenRow[realRow] = rc.clone() ;
	} else {
#       if CGL_DEBUG > 0
	  std::cout
	    << "INPLACE: rejected on cut coeff magnitude." << std::endl ;
#       endif
	}
      } else {
#       if CGL_DEBUG > 0
	  std::cout
	    << "INPLACE: rejected because no real row." << std::endl ;
#       endif
      }
    }

# if CGL_DEBUG > 0
  std::cout
    << "Leaving strengthenCoeff, "
    << rowCut.numberCuts() << " cuts." << std::endl ;
# endif

  return ;
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
int CglProbing::probe( const OsiSolverInterface &si, 
		       OsiCuts &cs, 
		       double *colLower, double *colUpper, 
		       CoinPackedMatrix *rowCopy,
		       CoinPackedMatrix *columnCopy,
		       const CoinBigIndex *rowStartPos, const int *realRows, 
		       const double *rowLower, const double *rowUpper,
		       const char *intVar, double *minR, double *maxR, 
		       int *markR, 
                       CglTreeInfo *info) const

{
# if CGL_DEBUG > 0
  std::cout << "Entering CglProbing::probe." << std::endl ;
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
  if (CoinAbs(cutoff) > 1.0e30)
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

  From CglProbingAnn: anyColumnCuts is set only in the case that we've
  fixed a variable by probing (i.e., one of the up or down probe resulted
  in infeasibility) and that probe entailed column cuts. Once set, it is
  never rescinded. In the reworked code, it's set as the return value of
  monotoneActions().
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

  JustFix < 0 says this is the first pass at the root; for any other
  iteration of the pass loop, it'll be initialised to 0.  If it is the first
  pass at the root, turn off row cuts, keep on fixing variables, and stop
  at the end of this pass. Otherwise, if we haven't fixed any variables, break.

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
      const bool downIsLower = (CoinAbs(down-colLower[j]) < 1.0e-7) ;
      const bool upIsUpper = (CoinAbs(up-colUpper[j]) < 1.0e-7) ;
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
	assert(downIsLower == (CoinAbs(down-colLower[j]) < 1.0e-7)) ;
	assert(upIsUpper == (CoinAbs(up-colUpper[j]) < 1.0e-7)) ;
#	endif

        if (way[iway] == probeDown) {
          movement = down-colUpper[j] ;
          solMovement = down-colsol[j] ;
          assert(movement < -0.99999) ;
          if (CoinAbs(down-colLower[j]) < 1.0e-7) {
            goingToTrueBound = 2 ;
            down = colLower[j] ;
          }
        } else {
          movement = up-colLower[j] ;
          solMovement = up-colsol[j] ;
          assert(movement > 0.99999) ;
          if (CoinAbs(up-colUpper[j]) < 1.0e-7) {
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
        if (CoinAbs(movement)>1.01) {
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
        if (CoinAbs(colUpper[j]-colLower[j]) < 1.0e-6)
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol]) < 1.0e-6) {
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
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
		      if (CoinAbs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
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
	motion, etc) that should not be performed when the probe shows
	infeasible. It might be that this is not the best approach.
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

  If we reach here, we may iterate the probe loop.

  Cases remaining include

  a) End of the first down probe; in this case we will iterate regardless of
     feasibility.
  b) End of the up probe, up and down feasible; in this case we will not
     iterate (we're done, cannot prove monotonicity).
  c) End of the up probe, down feasible, up infeasible; in this case we will
     iterate to restore the down feasible state.

  TODO: Unless I miss my guess, large chunks of this block of code will be
	replicated in the code that handles case b), as we generate cuts
	related to forcing the variable up.  Further inspection confirms
	this. -- lh, 101128 --

  TODO: I'm becoming more convinced that the second down probe iteration is
	unnecessary. Given that we have stackC0, lo0, and up0, seems like
	calling monotoneActions with lo0 = colLower, up0 = colUpper, and
	stackC0 = stackC should do the trick.  -- lh, 101216 --

  TODO: As I get the control flow sorted out a bit better, I should try to
  	separate b) (which does not iterate) from a) and c).

  The next block deals with the end of the first down probe. If the probe
  is feasible, attempt to tighten column bounds with singletonRows, process
  stackC to generate disaggregation cuts, and copy the tightened bounds
  to stackC0, lo0, and up0 for use after the up probe.  Feasible or not,
  we'll need to go on to the up probe, so restore the original bounds.

  Note that GoingToTrueBound was set to 0 way back in initial setup unless
  the user requested disaggregation cuts (rowCuts&0x01).

  I've rearranged the control flow here, moving the call to singletonRows
  ahead of the code that captures column bound changes to stackC0, et al.
  This makes it possible to use the same code at the end of case b) (down and
  up probe both feasible). Also, singletonRows is just another way to tighten
  column bounds. In its original position, it was subject to the same
  restrictions as coefficient strengthening (coefficient strengthening enabled
  and goingToTrueBound for the probing variable). I don't see any reason to
  keep them (but it might be good to add a separate control bit precisely for
  singletonRows).  -- lh, 110113 --
*/
	if (iway == downProbe) {
          if (!notFeasible) {
/*
  Attempt to tighten singleton bounds and generate disaggregation cuts.
*/
	    singletonRows(j,primalTolerance_,si,rowCopy,markC,
			  nstackC,stackC,saveL,saveU,
			  colUpper,colLower,columnGap,
			  nstackR,stackR,rowUpper,rowLower,maxR,minR) ;
	    if (goingToTrueBound == 2 && !justReplace) {
	      disaggCuts(nstackC,probeDown,primalTolerance_,
			 disaggEffectiveness,si,rowCut,stackC,colsol,
			 colUpper,colLower,saveU,saveL,index,element) ;
	    }
/*
  Copy off the bound changes from the down probe -- we'll need them after the
  up probe.
*/
	    nstackC0 = CoinMin(nstackC,maxStack) ;
	    for (istackC = 0 ; istackC < nstackC0 ; istackC++) {
	      int icol = stackC[istackC] ;
	      stackC0[istackC] = icol ;
	      lo0[istackC] = colLower[icol] ;
	      up0[istackC] = colUpper[icol] ;
	    }
	  } else {
	    nstackC0 = 0 ;
	  }
/*
  Now reset column bounds to their original values in preparation for the
  up probe.
*/
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
	  if ((rowCuts&0x02) != 0 && goingToTrueBound && !anyColumnCuts)
	    strengthenCoeff(j,iway,primalTolerance_,justReplace,canReplace,
	    		    needEffectiveness,si,rowCut,
	    		    rowCopy,colUpper,colLower,colsol,nstackR,stackR,
			    rowUpper,rowLower,maxR,minR,realRows,
			    element,index,info) ;
/*
  Restore row bounds.
*/
	  for (istackR = 0 ; istackR < nstackR ; istackR++) {
	    int irow = stackR[istackR] ;
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

  Moved the call to disaggCuts up from the common code for b) and c), because
  c) (up infeasible) sets goingToTrueBound to 0 and we can't execute
  disaggCuts. Similarly for strengthenCoeff.
*/
	  else {
	  if (iway == upProbe &&
	      feasRecord == (downProbeFeas|upProbeFeas)) {
	    iway = oneProbeTooMany ;

	    singletonRows(j,primalTolerance_,si,rowCopy,markC,
			  nstackC,stackC,saveL,saveU,
			  colUpper,colLower,columnGap,
			  nstackR,stackR,rowUpper,rowLower,maxR,minR) ;
	    if (goingToTrueBound == 2 && !justReplace) {
	      disaggCuts(nstackC,probeUp,primalTolerance_,
			 disaggEffectiveness,si,
			 rowCut,stackC,colsol,colUpper,colLower,saveU,saveL,
			 index,element) ;
	    if ((rowCuts&0x02) != 0 && goingToTrueBound && !anyColumnCuts)
	      strengthenCoeff(j,iway,primalTolerance_,justReplace,canReplace,
			      needEffectiveness,si,rowCut,
			      rowCopy,colUpper,colLower,colsol,nstackR,stackR,
			      rowUpper,rowLower,maxR,minR,realRows,
			      element,index,info) ;
	    }
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
		    if (CoinAbs(element[0]) > 1.0e-8) {
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

  And now the code to reset to the initial state.
*/
	  /* restore all */
	  for (istackC = nstackC-1 ; istackC >= 0 ; istackC--) {
	    int icol = stackC[istackC] ;
	    double oldU = saveU[istackC] ;
	    double oldL = saveL[istackC] ;
/*
  TODO: This next bit differs in subtle ways from the restore at the end
	of the down probe. Here we force markC to show fixed if the bounds are
	within 1.0e-4 of one another, a fairly broad tolerance; code for
	the pass 0 restore does not restore fixed status.  Also,
	we're not distinguishing here between actions for down feasible,
	up feasible, vs. down feasible, up infeasible.  -- lh, 101128 --

	I still don't see any reason for the restore here to differ from
	the restore after the down probe.  -- lh, 110113 --
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
  End of column restore. On to rows.
*/
	  for (istackR = 0 ; istackR < nstackR ; istackR++) {
	    int irow = stackR[istackR] ;
	    minR[irow] = saveMin[istackR] ;
	    maxR[irow] = saveMax[istackR] ;
	    markR[irow] = -1 ;
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
