
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
  array.
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

  The first set (downIter, upIter, oneIterTooMany) *must*
  match the iteration number for the corresponding pass in the PROBE
  loop. The symbolic names make the code a bit easier to read.
*/
const unsigned int downIter = 0 ;
const unsigned int upIter = 1 ;
const unsigned int oneIterTooMany = 2 ;

const unsigned int probeDown = tightenUpper ;
const unsigned int probeUp = tightenLower ;

const unsigned int downIterFeas = 0x1 ;
const unsigned int upIterFeas = 0x2 ;

}  // end file-local namespace

// ===============================================

/*
  This method restores column and row bounds after a probe. We need to do this
  in several circumstances:
    * After a feasible down probe. We need to restore the original state prior
      to the up probe.
    * After feasible down and up probes. We need to restore the original state
      so we can probe the next variable.
    * After a feasible down probe followed by an infeasible up probe. We're
      monotone down, but before we install the down state we need to erase the
      up state.
*/

void restoreRowColBounds (double tol, int nColStack,
			  const int *const colStack, int *const colMark,
			  const double *const savedColL,
			  const double *const savedColU,
			  double *const colL, double *const colU,
			  double *const colGap,
			  int nRowStack,
			  const int *const rowStack, int *const rowMark,
			  const double *const savedRowL,
			  const double *const savedRowU,
			  double *const rowL, double *const rowU
			 )

{
/*
  Reset column bounds to their original values and reset the mark flags.
*/
  for (int iColStack = nColStack-1 ; iColStack >= 0 ; iColStack--) {
    int icol = colStack[iColStack] ;
    double oldU = savedColU[iColStack] ;
    double oldL = savedColL[iColStack] ;
    double oldGap = oldU-oldL-tol ;
    colU[icol] = oldU ;
    colL[icol] = oldL ;
    colGap[icol] = oldGap ;
    if (oldGap < 1.0e-4) {
      colMark[icol] = (tightenUpper|tightenLower) ;
    } else {
      colMark[icol] = 0 ;
      if (oldU > 1.0e10) colMark[icol] |= infUpper ;
      if (oldL < -1.0e10) colMark[icol] |= infLower ;
    }
  }
/*
  And now the row bounds.
*/
  for (int iRowStack = 0 ; iRowStack < nRowStack ; iRowStack++) {
    int irow = rowStack[iRowStack] ;
    rowL[irow] = savedRowL[iRowStack] ;
    rowU[irow] = savedRowU[iRowStack] ;
    rowMark[irow] = -1 ;
  }

  return ;
}




/*
  This method looks for implication cuts: implications that can be generated
  given a binary probing variable x<p>. For an arbitrary binary variable x<j>,
  we're looking for things like the following:
    1) x<p> -> 0 ==> x<j> -> 0 and x<p> -> 1 ==> x<j> -> 0
       In this case, we can fix x<j> to 0. Similarly if were driven to 1.
    2) x<p> -> 0 ==> x<j> -> 0 and x<p> -> 1 ==> x<j> -> 1
       In this case, we can write x<j> = x<p>. If the sense is reversed, we
       get x<j> = (1-x<p>).
  Form 2) are actual cuts; form 1) just amounts to a bound change.

  You can extend this to arbitrary bounds on x<j>. Let ld<j> be the lower
  bound for the down probe and ud<j> be the upper bound for the down probe.
  Similarly, lu<j> and uu<j> are lower and upper bounds from the up probe, and
  l<j> and u<j> are the original bounds.
    1) The equivalent is to look for max(min(ld<j>,lu<j>),l<j>) and
       min(max(ud<j>,uu<j>),u<j>).
    2) The equivalent constraints are x<j> = (u<j>-l<j>)x<p> + l<j> and
       x<j> = -(u<j>-l<j>)x<p> + u<j>
  As with binary variables, form 2) are actual cuts; form 1) just amounts to
  a bound change.

  1) is applied to continuous variables as well as integer variables,
  and the result is captured in vec_l and vec_u. Normally, these
  parameters will be saveL and saveU, hence the result will propagate
  out from CglProbing::genCutsAndModify when the `original' bounds are
  restored. Only integer variables are recorded as column cuts.

  The method uses jProbe to decide if it should generate form 2) cuts. If
  jProbe >= 0, it's assumed to be a binary variable and form 2) cuts are
  attempted. Recall that stackC[0] and stackC0[0] are the probing variable,
  that's why we don't process them.

  Bounds held in vec_lu and vec_uu are assumed to be indexed by column
  number.  Bounds held in vec_l, and vec_u are assumed to be correlated
  with entries on stackC. Bounds held in vec_ld and vec_ud are assumed to
  be correlated with entries on stackC0. This makes sense when you consider
  the context: we're calling implicationCuts after discovering that the up and
  down probes were both feasible. The bounds from the down probe (ld, ud) are
  in lo0, up0; the bounds from the up probe (lu, uu) are in colLower,
  colUpper, and the original bounds (l, u) are in saveL, saveU.

  MarkC, element, and index are used as work arrays. It's relatively rare that
  any given probe will change bounds on half the variables, so we're going to
  bet on that. We'll store lower bounds from [0] of index and element and
  upper bounds from [nCols-1].

  TODO: Check that the above comment is correct --- do improvements to
	continuous bounds really propagate back to the caller?
	-- lh, 110201 --
*/

void implicationCuts (int jProbe, double primalTol, int nCols,
		      OsiCuts &cs, const OsiSolverInterface &si,
		      const char *const intVar,
		      const double *const colsol,
		      int nC, const int *const stackC, int *const markC,
		      const double *const vec_uu, const double *const vec_lu,
		      double *const vec_u, double *const vec_l,
		      int nC0, const int *const stackC0,
		      const double *const vec_ud, const double *const vec_ld,
		      double *const element, int *const index
		     )
{
# if CGL_DEBUG > 1
  std::cout
    << "Entering implicationCuts, probe " << jProbe << "." << std::endl ;
  bool hitTheWall = false ;
  int origRowCutCnt = cs.sizeRowCuts() ;
  int origColCutCnt = cs.sizeColCuts() ;
/*
  If we have a binary probe variable, it should not be fixed.
*/
  if (jProbe >= 0) assert(vec_l[0] == 0.0 && vec_u[0] == 1.0) ;
# endif
  const int stkOffset = CoinMax(1000,(2*nC)) ;
/*
  Construct a cross-reference so that we can correlate entries on stackC0 and
  stackC. markC is no longer of use, so load entry j with the position of x<j>
  on stackC. Keep in mind that stackC0 and stackC will not necessarily hold
  the same variables; the offset guarantees we can recognise the entries
  made here.
*/
  for (int iC = nC-1 ; iC >= 0 ; iC--) {
    int j = stackC[iC] ;
    markC[j] = iC+stkOffset ;
  }
/*
  See if we have bounds improvement. Given a lower bound ld<j> from the
  down probe, a bound lu<j> from the up probe, and the original bound lo<j>,
  we can surely impose the lesser of ld<j> and lu<j> if this is more than the
  original lo<j>.  In math, max(min(ld<j>,lu<j>),lo<j>). Use a conservative
  tolerance. Stash the improved bound in saveL so that it'll be restored
  when we restore the `original' (pre-probing) bounds. For integer variables,
  add a column cut.
*/
  int nTot = 0 ;
  int nInt = 0 ;
  int ilb = 0 ;
  int iub = nCols-1 ;
  bool ifCut = false ;
  const double minChange = 1.0e-4 ;
  int implIndices[2] ;
  double implCoeffs[2] ;
/*
  Scan stackC0 looking for variables that had bound improvement on both the up
  and down probe. If the variable was changed by only one probe, move on.
*/
  for (int iC0 = 1 ; iC0 < nC0 ; iC0++) {
    int j = stackC0[iC0] ;
    int iC = markC[j]-stkOffset ;
    if (iC < 0) continue ;
    double udj = vec_ud[iC0] ;
    double ldj = vec_ld[iC0] ;
    double uuj = vec_uu[j] ;
    double luj = vec_lu[j] ;
    double &uj = vec_u[iC] ;
    double &lj = vec_l[iC] ;
/*
  Obtain the consensus lower bound.
*/
    double bestProbelj = CoinMin(ldj,luj) ;
    if (bestProbelj > lj+minChange) {
#     if CGL_DEBUG > 1
      std::cout << "      consensus lb for " ;
      if (intVar[j]) std::cout << "(int) " ;
      std::cout
        << " x<" << j << "> improved from "
	<< lj << " to " << bestProbelj << "." << std::endl ;
#     endif
      lj = bestProbelj ;
      nTot++ ;
      if (intVar[j] && nInt < nCols) {
	element[ilb] = bestProbelj ;
	index[ilb++] = j ;
	nInt++ ;
	if (colsol[j] < bestProbelj-primalTol)
	  ifCut = true ;
      }
#     if CGL_DEBUG > 1
      else {
        if (intVar[j] && nInt >= nCols && hitTheWall == false) {
	  std::cout
	    << "CglProbing::implicationCuts: (lb) exceeded column cut space!."
	    << std::endl ;
	  hitTheWall = true ;
	}
      }
#     endif
    } 
/*
  Obtain the consensus upper bound.
*/
    double bestProbeuj = CoinMax(udj,uuj) ;
    if (bestProbeuj < uj-minChange) {
#     if CGL_DEBUG > 1
      std::cout << "      consensus ub for " ;
      if (intVar[j]) std::cout << "(int) " ;
      std::cout
        << " x<" << j << "> improved from "
	<< uj << " to " << bestProbeuj << "." << std::endl ;
#     endif
      uj = bestProbeuj ;
      nTot++ ;
      if (intVar[j] && nInt < nCols) {
	element[iub] = bestProbeuj ;
	index[iub--] = j ;
	nInt++ ;
	if (colsol[j] > bestProbeuj+primalTol)
	  ifCut = true ;
      }
#     if CGL_DEBUG > 1
      else {
        if (intVar[j] && nInt >= nCols && hitTheWall == false) {
	  std::cout
	    << "CglProbing::implicationCuts: (ub) exceeded column cut space!."
	    << std::endl ;
	  hitTheWall = true ;
	}
      }
#     endif
    } 
/*
  Now look for honest implication cuts (form 2 in the introductory comments).
  Recall that jProbe >= 0 is taken as an indication that the probe variable
  is binary and these cuts should be attempted.  It's possible (miracles
  happen) that the above code has managed to fix x<j> by pulling together
  the bounds from the down and up probe (but we should not lose feasibility
  here). In the case that x<j> is fixed, move on.
*/
    if (jProbe >= 0) {
      assert(uj >= lj) ;
      if (CoinAbs(uj-lj) < primalTol) continue ;
/*
  Check for x<j> = (u<j>-l<j>)x<p> + l<j>. In words, check that the down probe
  forced x<j> to l<j> and the up probe forced x<j> to u<j>.
*/
      if (udj < lj+primalTol && luj > uj-primalTol) {
	OsiRowCut rc ;
	rc.setLb(lj) ;
	rc.setUb(lj) ;
	rc.setEffectiveness(1.0e-5) ;
	implIndices[0] = jProbe ;
	implIndices[1] = j ;
	implCoeffs[0] = -(uj-lj) ;
	implCoeffs[1] = 1.0 ;
	rc.setRow(2,implIndices,implCoeffs,false) ;
#       if CGL_DEBUG > 1
	std::cout << "      implication:" ;
	rc.print() ;
#       endif
	cs.insert(rc) ;
      } else if (uuj < lj+primalTol && ldj > uj-primalTol) {
/*
  Repeat for x<j> = -(u<j>-l<j>)x<p> + u<j>
*/
	OsiRowCut rc ;
	rc.setLb(uj) ;
	rc.setUb(uj) ;
	rc.setEffectiveness(1.0e-5) ;
	implIndices[0] = jProbe ;
	implIndices[1] = j ;
	implCoeffs[0] = uj-lj ;
	implCoeffs[1] = 1.0 ;
	rc.setRow(2,implIndices,implCoeffs,false) ;
#       if CGL_DEBUG > 1
	std::cout << "      inverse implication:" ;
	rc.print() ;
#       endif
	cs.insert(rc) ;
      } 
    }
  }
/*
  If we have column cuts, load them into the cut collection passed in as a
  parameter.
*/
  if (nInt > 0) {
    OsiColCut cc ;
    if (ilb > 0) {
      cc.setLbs(ilb,index,element) ;
    }
    if (iub < nCols-1) {
      cc.setUbs((nCols-1)-iub,&index[iub+1],&element[iub+1]) ;
    }
    if (ifCut) {
      cc.setEffectiveness(100.0) ;
    } else {
      cc.setEffectiveness(1.0e-5) ;
    }
#   if CGL_DEBUG > 1
    std::cout
      << "    tightened " << ilb << " lb, " << (nCols-1)-iub << " ub" ;
    if (ifCut)
      std::cout << "; cut." << std::endl ;
    else
      std::cout << "; no cut." << std::endl ;
    CglProbingDebug::checkBounds(si,cc) ;
#   endif
    cs.insert(cc) ;
  }
# if CGL_DEBUG > 1
  std::cout
    << "Leaving implicationCuts, " << cs.sizeRowCuts()-origRowCutCnt
    << " row cuts, " << cs.sizeColCuts()-origColCutCnt << " col cuts."
    << std::endl ;
# endif

  return ;
}


/*
  Run through stackC and create disaggregation cuts. See the typeset
  documentation for the full derivation. Note that the constraints generated
  here are a specialised form and require that the change in bound
  (u<p> to u'<p>, or l<p> to l'<p>) is distance 1.
  
  Suppose the probe variable is x<p> and we're trying to generate a
  disaggregation cut for target x<t>.

  For a down probe, we have u<p> reduced to u'<p> = colUpper[p].
  It may be that this causes u<t> to be reduced to u'<t>. The upper
  disaggregation cut is
    -(u<t>-u'<t>)x<p> + x<t> <= u'<t> - u'<p>(u<t>-u'<t>)

  It may be that this causes l<t> to be increased to l'<t>. The lower
  disaggregation cut is
    -(l<t>-l'<t>)x<p> + x<t> >= l'<t> - u'<p>(l<t>-l'<t>)

  For an up probe, we have l<p> increased to l'<p> = colLower[p].
  It may be that this causes u<t> to be reduced to u'<t>. The upper
  disaggregation cut is
     (u<t>-u'<t>)x<p> + x<t> <= u'<t> + l'<p>(u<t>-u'<t>)

  It may be that this causes l<t> to be increased to l'<t>. The lower
  disaggregation cut is
     (l<t>-l'<t>)x<p> + x<t> >= l'<t> + l'<p>(l<t>-l'<t>)

  Note that stackC[0] = pndx, the index of the probed variable, so we
  don't want to process stackC[0] in the main loop.

  These cuts do not cast the probe value in concrete, so we can apply them
  whether or not the probe cuts off portions of the current solution.
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

# if CGL_DEBUG > 1
  assert((probeDir == probeDown) || (probeDir == probeUp)) ;
  std::cout
    << "Entering disaggCuts, "
    << ((probeDir == probeDown)?"down":"up") << " probe on "
    << "x<" << pndx << ">." << std::endl ;
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  int cutsAtStart = rowCut.numberCuts() ;
# endif

  int plusMinus = 0 ;
  double luPrime_p = 0.0 ;
  double x_p = colsol[pndx] ;
  double deltaProbe_p = 0.0 ;
/*
  Set up the quantities that vary depending on whether this is an up or down
  probe:
  * lu'<p> is the new upper (down probe) or lower (up probe) bound on x<p>.
  * plusMinus will be -1 (down probe) or +1 (up probe).

  deltaProbe_p is then x*<p>-lu'<p>; we'll use it below to check that the
  cut is violated and effective.

  (d,u)  -(u<t>-u'<t>)x<p> + x<t> <= u'<t> - u'<p>(u<t>-u'<t>)
  (d,l)  -(l<t>-l'<t>)x<p> + x<t> >= l'<t> - u'<p>(l<t>-l'<t>)
  (u,u)   (u<t>-u'<t>)x<p> + x<t> <= u'<t> + l'<p>(u<t>-u'<t>)
  (u,l)   (l<t>-l'<t>)x<p> + x<t> >= l'<t> + l'<p>(l<t>-l'<t>)
*/
  if (probeDir == probeDown) {
    luPrime_p = colUpper[pndx] ;
    plusMinus = -1 ;
  } else {
    luPrime_p = colLower[pndx] ;
    plusMinus = 1 ;
  }
  deltaProbe_p = x_p-luPrime_p ;

# if CGL_DEBUG > 2
  if (probeDir == probeDown) {
    std::cout
      << "  u<" << pndx << "> = " << saveU[0]
      << " reduced to " << luPrime_p << "." << std::endl ;
  } else {
    std::cout
      << "  l<" << pndx << "> = " << saveL[0]
      << " increased to " << luPrime_p << "." << std::endl ;
  }
  std::cout
    << "  x<" << pndx << "> = " << x_p << "; deltaProbe = "
    << deltaProbe_p << std::endl ;
# endif

  for (int istackC = nstackC-1 ; istackC > 0 ; istackC--) {
    int icol = stackC[istackC] ;
    double u_t = saveU[istackC] ;
    double l_t = saveL[istackC] ;
    double x_t = colsol[icol] ;
    double uPrime_t = colUpper[icol] ;
    double lPrime_t = colLower[icol] ;
    double deltaProbe_t = u_t-uPrime_t ;
    double a_p = plusMinus*deltaProbe_t ;
/*
  Generate the upper disaggregation cut, if it's violated and the coefficients
  will be reasonable.

    -(u<t>-u'<t>)x<p> + x<t> <= u'<t> - u'<p>(u<t>-u'<t>)  d,u
     (u<t>-u'<t>)x<p> + x<t> <= u'<t> + l'<p>(u<t>-u'<t>)  u,u

  The effectiveness is the improvement in x<p>.  For a down probe, x<p>
  should decrease (move closer to u'<p>); for an up probe, x<p> should
  increase (move closer to l'<p>). Given deltaProbe_p = (x*<p>-lu'<p>),
  check against deltaCut_p = ((u'<t>-x*<t>)/a<p>).

*/
    if (deltaProbe_t > 0.0 && u_t < 1.0e10 &&
	(x_t > uPrime_t-deltaProbe_p*a_p+primalTolerance_)) {
      double deltaCut_p = (uPrime_t-x_t)/a_p ;
      double effectiveness = plusMinus*(deltaProbe_p-deltaCut_p) ;
#     if CGL_DEBUG > 2
      std::cout
	<< "    Generating upper disaggregation cut for x<"
	<< icol << "> (" << istackC << ")." << std::endl ;
      std::cout
	<< "    u<" << icol << "> = " << u_t
	<< " reduced to " << uPrime_t << "." << std::endl ;
      std::cout
	<< "    x<" << icol << "> = " << x_t
	<< " reduced to " << (uPrime_t-deltaProbe_p*a_p) << "." << std::endl ;
      std::cout
	<< "    x<" << pndx << "> = " << x_p
	<< " changed to " << deltaCut_p+luPrime_p
	<< "; effectiveness = " << effectiveness
	<< "; required " << disaggEffectiveness << "." << std::endl ;
#     endif
      assert(effectiveness+primalTolerance_ > 0) ;
      if (effectiveness > disaggEffectiveness) {
	OsiRowCut rc ;
	rc.setEffectiveness(effectiveness) ;
	rc.setLb(-DBL_MAX) ;
	rc.setUb(uPrime_t+a_p*luPrime_p) ;
	index[0] = icol ;
	element[0] = 1.0 ;
	index[1] = pndx ;
	element[1] = a_p ;
	rc.setRow(2,index,element,false) ;
#	if CGL_DEBUG > 1
	if (debugger) assert(!debugger->invalidCut(rc)); 
	std::cout << "    Adding upper disaggregation cut." << std::endl ;
#	endif
	rowCut.addCutIfNotDuplicate(rc) ;
      }
    }
/*
  Generate the lower disaggregation cut.

    -(l<t>-l'<t>)x<p> + x<t> >= l'<t> - u'<p>(l<t>-l'<t>)  d,l
     (l<t>-l'<t>)x<p> + x<t> >= l'<t> + l'<p>(l<t>-l'<t>)  u,l
*/
    deltaProbe_t = l_t-lPrime_t ;
    a_p = plusMinus*deltaProbe_t ;
    if (deltaProbe_t < 0.0 && l_t > -1.0e10 &&
	(x_t < lPrime_t-deltaProbe_p*a_p-primalTolerance_)) {
      
      double deltaCut_p = (lPrime_t-x_t)/a_p ;
      double effectiveness = plusMinus*(deltaProbe_p-deltaCut_p) ;
#     if CGL_DEBUG > 2
      std::cout
	<< "    Generating lower disaggregation cut for x<"
	<< icol << "> (" << istackC << ")." << std::endl ;
      std::cout
	<< "    l<" << icol << "> = " << l_t
	<< " increased to " << lPrime_t << "." << std::endl ;
      std::cout
	<< "    x<" << icol << "> = " << x_t
	<< " increased to " << (lPrime_t-deltaProbe_p*a_p) << "." << std::endl ;
      std::cout
	<< "    x<" << pndx << "> = " << x_p
	<< " changed to " << deltaCut_p+luPrime_p
	<< "; effectiveness = " << effectiveness
	<< "; required " << disaggEffectiveness << "." << std::endl ;
#     endif
      assert(effectiveness+primalTolerance_ > 0) ;
      if (effectiveness > disaggEffectiveness) {
	OsiRowCut rc ;
	rc.setEffectiveness(effectiveness) ;
	rc.setLb(lPrime_t+a_p*luPrime_p) ;
	rc.setUb(DBL_MAX) ;
	index[0] = icol ;
	element[0] = 1.0 ;
	index[1] = pndx ;
	element[1] = a_p ;
	rc.setRow(2,index,element,false) ;
#	if CGL_DEBUG > 1
	if (debugger) assert(!debugger->invalidCut(rc)); 
	std::cout << "    Adding lower disaggregation cut." << std::endl ;
#	endif
	rowCut.addCutIfNotDuplicate(rc) ;
      }
    }
  }

# if CGL_DEBUG > 1
  std::cout << "Leaving disaggCuts" ;
  if (rowCut.numberCuts() > cutsAtStart)
    std::cout << ", " << rowCut.numberCuts()-cutsAtStart << " cuts" ;
  std::cout << "." << std::endl ;
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
		      const CoinBigIndex *const columnStart,
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

  const CoinBigIndex kkstart = columnStart[j] ;
  const CoinBigIndex kkend = kkstart+columnLength[j] ;

  for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
    int irow = row[kk] ;
    double value = columnElements[kk] ;
    double delta = value*movement ;
/*
  markR is initialised by calcRowBounds. A code of -1 indicates this row
  should be processed; -2 says ignore. Other codes should not happen.

  TODO: This assert, and the disabled continue, go back to the change in
	calcRowBounds where rows with infinite bounds are handed arbitrary
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
      
      CglProbingDebug::dump_row(irow,rowLower[irow],rowUpper[irow],
	  minR[irow],maxR[irow],0,false,false,0,0,0,0,0,0,0) ;
      std::cout << std::endl ;
    }
    assert (markR[irow] >= -1 && markR[irow]%10000 < nstackR) ;
#   endif
#   if CGL_DEBUG > 2
    double minRstart = minR[irow] ;
    double maxRstart = maxR[irow] ;
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
#   if CGL_DEBUG > 2
    if (minRstart != minR[irow]) {
      std::cout
	<< "      row " << irow << " min "
	<< minRstart << " -> " << minR[irow] << std::endl ;
    }
    if (maxRstart != maxR[irow]) {
      std::cout
	<< "      row " << irow << " max "
	<< maxRstart << " -> " << maxR[irow] << std::endl ;
    }
#   endif
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
# if CGL_DEBUG > 1
  std::cout << "Entering groomSoln." << std::endl ;
  std::cout << "  verifying matrix." ;
  rowCopy->verifyMtx(3) ;
# endif
  int nRows = rowCopy->getNumRows() ;
  int nCols = rowCopy->getNumCols() ;
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *rowLength = rowCopy->getVectorLengths() ;
  const int *column = rowCopy->getIndices() ;
  const double *rowElements = rowCopy->getElements() ;

/*
  Find the largest potential postive and negative swing in a row, where
  swing is defined as a<ij>(u<j>-l<j>). Check correctness of rowStartPos if
  we're debugging.
*/
  for (int i = 0 ; i < nRows ; i++) {
    const CoinBigIndex kkstart = rowStart[i] ;
    const CoinBigIndex kkendneg = rowStartPos[i] ;
    const CoinBigIndex kkend = kkstart+rowLength[i] ;
    double value = 0.0 ;
    CoinBigIndex kk = 0 ;

#if CGL_DEBUG > 0
    for (kk = kkstart ; kk < kkend ; kk++) {
      value = rowElements[kk] ;
      if (value > 0.0) 
	break ;
    }
    if (rowStartPos[i] != kk) {
      std::cout
        << "  Row " << i << ": " << kkstart << ".." << kkend << ", len "
	<< rowLength[i] << ", pos at " << rowStartPos[i] << "." << std::endl ;
    }
    assert (rowStartPos[i] == kk) ;
    value = 0.0 ;
#endif

    for (kk = kkstart ; kk < kkendneg ; kk++) {
      int iColumn = column[kk] ;
      double gap = CoinMin(1.0e100,colUpper[iColumn]-colLower[iColumn]) ;
      value = CoinMin(value,gap*rowElements[kk]) ;
    }
    largestNegativeInRow[i] = value ;
    value = 0.0 ;
    for ( ; kk < kkend ; kk++) {
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
# if CGL_DEBUG > 1
  std::cout << "Leaving groomSoln." << std::endl ;
# endif
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
  int anyColumnCuts = 0 ;
  int nTot = 0 ;
  int nFix = 0 ;
  const double *origColLower = si.getColLower() ;
  const double *origColUpper = si.getColUpper() ;
# if CGL_DEBUG > 1
  int probedVar = stackC[0] ;
  std::cout << "Entering monotoneActions." ;
# if CGL_DEBUG > 1
  std::cout
    << std::endl
    << "    " << " x(" << probedVar
    << ") monotone [" << colLower[probedVar] << ", "
    << colUpper[probedVar] << "] from " << colsol[probedVar]
    << "." << std::endl ;
# endif
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
	  anyColumnCuts++ ;
	}
#       if CGL_DEBUG > 1
        std::cout
	  << "    " << " x(" << icol
	  << ") ub " << origColUpper[icol] << " ==> " << colUpper[icol] ;
	if (ifCut) std::cout << " (cut)" ;
	std::cout << "." << std::endl ;
#       endif
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
	  anyColumnCuts++ ;
	}
#       if CGL_DEBUG > 1
        std::cout
	  << "    " << " x(" << icol
	  << ") lb " << origColLower[icol] << " ==> " << colLower[icol] ;
	if (ifCut) std::cout << " (cut)" ;
	std::cout << "." << std::endl ;
#       endif
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
#   if CGL_DEBUG > 0
    CglProbingDebug::checkBounds(si,cc) ;
#   endif
    cs.insert(cc) ;
  }

# if CGL_DEBUG > 1
  std::cout << "  Leaving monotoneActions" ;
  if (nTot > 0)
    std::cout << ", " << nTot << " monotone" ;
  if (anyColumnCuts > 0)
    std::cout << ", " << anyColumnCuts << " column cuts" ;
  std::cout << "." << std::endl ;
# endif
  return (anyColumnCuts > 0) ;
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

  Note that we cannot process the objective here --- that'd be nonsense.
  It's a minor pain to exclude it.

  John's original comment for this code was
    // also see if singletons can go to good objective
    // Taken out as should be found elsewhere
    // and has to be original column length
  but he reinstated it. Apparently there were cases where fixing the probing
  variable was required to satisfy the condition that all unfixed variables be
  singletons. Enough of them to justify reinstating this code.
*/

void singletonRows (int jProbe, double primalTolerance_, bool useObj,
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
# if CGL_DEBUG > 1
  std::cout << "Entering singletonRows." << std::endl ;
  int nstackC_orig = nstackC ;
# endif
/*
  Unpack a few vectors from the row-major matrix.
*/
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *rowLength = rowCopy->getVectorLengths() ;
  const int *column = rowCopy->getIndices() ;
  const double *rowElements = rowCopy->getElements() ;
  const int nCols = rowCopy->getNumCols() ;
  const int nRows = rowCopy->getNumRows() ;
/*
  `Singleton' must be based on the column length in the original system.
*/
  const double *objective = si.getObjCoefficients() ;
  const int *columnLengths = si.getMatrixByCol()->getVectorLengths() ;
  const double objSense = si.getObjSense() ;
/*
  Open a loop to work through the rows on stackR. Don't process the objective
  (if present, it's the last row).
*/
  for (int istackR = 0 ; istackR < nstackR ; istackR++) {
    int i = stackR[istackR] ;
    if (useObj && i == (nRows-1)) continue ;
/*
  Check the gaps. If the constraint is potentially tight in both directions,
  there's nothing more to do here.
*/
    const double uGap = rowUpper[i]-maxR[i] ;
    const double lGap = minR[i]-rowLower[i] ;
    if (uGap < primalTolerance_ && lGap < primalTolerance_) continue ;
/*
  Scan the row to see if it meets the `all singletons' condition and
  contains at least one unfixed singleton. Again, if this fails, there's
  nothing more to be done.

  Note that the original code didn't check the probing variable x<p>,
  because if you're probing a binary variable it's fixed. But for general
  integer variables a probe does not in general fix the variable.  So we
  check all variables.

  We should not be executing this method if we have prima facie infeasibility.
*/
    bool canFix = true ;
    bool target = false ;
    const CoinBigIndex kkstart = rowStart[i] ;
    const CoinBigIndex kkend = kkstart+rowLength[i] ;
    for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
      int j = column[kk] ;
      assert(colUpper[j]-colLower[j] > -primalTolerance_) ;
      if (colUpper[j] > colLower[j]) {
        if (columnLengths[j] != 1) {
	  canFix = false ;
	  break ;
	} else {
	  target = true ;
	}
      }
    }
    if (!canFix || !target) continue ;
#   if CGL_DEBUG > 2
    CglProbingDebug::dump_row(i,rowLower[i],rowUpper[i],minR[i],maxR[i],
	&si,true,true,4,rowLength[i],&column[rowStart[i]],
        &rowElements[rowStart[i]],primalTolerance_,colLower,colUpper) ;
#   endif
/*
  If we've passed the tests, we can look for variables suitable to drive to
  bound. Work against U<i> first. We're looking for variables with a<ij> > 0
  that will be naturally driven to u<j>, or variables with a<ij> < 0 that will
  be naturally driven to l<j>.
  
  Don't overwrite the saved bounds if we've tightened this variable already!
*/
    if (uGap > primalTolerance_) {
      for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
	int j = column[kk] ;
	const double lj = colLower[j] ;
	const double uj = colUpper[j] ;
	if (uj > lj) {
	  double value = rowElements[kk] ;
	  if (objSense*objective[j]*value < 0.0)
	  { if (!(markC[j]&(tightenLower|tightenUpper))) {
	      assert(nstackC < nCols) ;
	      stackC[nstackC] = j ;
	      saveL[nstackC] = lj ;
	      saveU[nstackC] = uj ;
	      nstackC++ ;
	    }
	    if (value > 0.0) {
	      colLower[j] = uj ;
	    } else {
	      colUpper[j] = lj ;
	    }
	    markC[j] |= tightenLower|tightenUpper ;
	    colGap[j] = -primalTolerance_ ;
#	    if CGL_DEBUG > 1
	    std::cout
	      << "    row " << i << " uGap " << rowUpper[i] << "-"
	      << maxR[i] << " = " << uGap
	      << " " << " x(" << j
	      << ") c = " << objective[j] << ", a = " << value ;
	    if (colLower[j] == uj)
	      std::cout << " l " << lj << " -> " << uj ;
	    else
	      std::cout << " u " << uj << " -> " << lj ;
	    std::cout << std::endl ;
#	    endif
	  }
	}
      }
    }
/*
  And now the identical code, except that we're working against L<i>, hence
  the sense of the eligibility test is reversed and we want variables with
  a<ij> > 0 that will be naturally driven to l<j>, or variables with
  a<ij> < 0 that will be naturally driven to u<j>.
*/
    if (lGap > primalTolerance_) {
      for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
	int j = column[kk] ;
	const double lj = colLower[j] ;
	const double uj = colUpper[j] ;
	if (uj > lj) {
	  double value = rowElements[kk] ;
	  if (objSense*objective[j]*value > 0.0)
	  { if (!(markC[j]&(tightenLower|tightenUpper))) {
	      assert(nstackC < nCols) ;
	      stackC[nstackC] = j ;
	      saveL[nstackC] = lj ;
	      saveU[nstackC] = uj ;
	      nstackC++ ;
	    }
	    if (value < 0.0) {
	      colLower[j] = uj ;
	    } else {
	      colUpper[j] = lj ;
	    }
	    markC[j] |= tightenLower|tightenUpper ;
	    colGap[j] = -primalTolerance_ ;
#	    if CGL_DEBUG > 1
	    std::cout
	      << "    row " << i << " lGap " << minR[i] << "-"
	      << rowLower[i] << " = " << lGap
	      << " " << " x(" << j
	      << ") c = " << objective[j] << ", a = " << value ;
	    if (colLower[j] == uj)
	      std::cout << " l " << lj << " -> " << uj ;
	    else
	      std::cout << " u " << uj << " -> " << lj ;
	    std::cout << std::endl ;
#	    endif
	  }
	}
      }
    }
  }
# if CGL_DEBUG > 1
  std::cout << "Leaving singletonRows" ;
  if (nstackC > nstackC_orig)
    std::cout << ", stacked " << nstackC-nstackC_orig << " vars" ;
  std::cout << "." << std::endl ;
# endif
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

  a<ip>  dir    delta          blow'          a'<ip>      b'

    >0    d    b - U'                         a-delta   b-(u'+1)delta
    >0    u    blow - L'   blow+(l'-1)delta   a+delta
    <0    d    blow - L'   blow-(u'+1)delta   a-delta
    <0    u    b - U'                         a+delta   b+(l'-1)delta

  In the code, negating delta for down probes unifies the up and down math.

  Note that coefficient strengthening isn't conditional; it embeds the probe
  value in the constraint system. We can use the results under two conditions:
    1) The probe is monotone, hence we know we'll be keeping this value.
    2) The bound changes associated with the probe didn't cut off any of the
       current optimal solution.
  Condition 2) is conservative, because most often the solution at hand is a
  solution to a continuous relaxation. You can argue that 1) is a special case
  of 2).

  There's another, more subtle consideration when we're in `just replace'
  mode. We're tilting a linear constraint. Making it tighter to one side of a
  pivot makes it looser to the other side. If we're going to outright replace
  a constraint, there cannot be any polytope past the active side of the
  dichotomy.
*/

void strengthenCoeff (
		      int jProbe, unsigned int probeDir,
		      double primalTolerance_,
		      bool justReplace,
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
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
# endif
# if CGL_DEBUG > 1
  std::cout
    << "Entering strengthenCoeff, probing "
    << "x<" << jProbe << "> "
    << ((probeDir == probeDown)?"down":"up") << "." << std::endl ;
  int newCuts = 0 ;
  int cutsAtStart = rowCut.numberCuts() ;
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
  const int *rowLen = rowCopy->getVectorLengths() ;
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
  constraint redundant (but not too redundant --- the gap value will become
  the new coefficient a<ip>', so it must be reasonable).

  Range constraints pose the danger that we'll want to assign two different
  values to a<ip>. If we're in `justReplace' mode, this is a nonstarter.
*/
    double uGap = bi-maxR[irow] ;
    double lGap = blowi-minR[irow] ;
    bool useableUGap = ((uGap > primalTolerance_) && (uGap < bigCoeff)) ;
    bool useableLGap = ((-lGap > primalTolerance_) && (-lGap < bigCoeff)) ;
    if (!(useableUGap || useableLGap)) continue ;
    bool isRangeCon = ((blowi  > -rhsInf) && (bi < rhsInf)) ;
    if (isRangeCon && justReplace) continue ;
/*
  We'll need the lhs activity, excluding the probed variable, to evaluate the
  effectiveness of the cut.  Extract the coefficient of the probe variable
  while we're here; we need it to determine which case we're working.
*/
    double sum = 0.0 ;
    double aip = 0.0 ;
    bool aipNonZero = false ;
    const CoinBigIndex kkstart = rowStart[irow] ;
    const CoinBigIndex kkend = kkstart+rowLen[irow] ;
    for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
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
  Now that we know a<ip> and the probe direction, we can finally decide
  which gap is important. Set up deltaMul while we're at it, it's little
  enough additional work.
*/
    double delta ;
    double deltaMul ;
    bool revisebi = true ;
    if (probeDir == probeDown) {
      if (aip >= 0) {
        if (!useableUGap) continue ;
        delta = -uGap ;
      } else {
        if (!useableLGap) continue ;
        delta = -lGap ;
	revisebi = false ;
      }
      deltaMul = colUpper[jProbe]+1.0 ; 
    } else {
      if (aip >= 0) {
        if (!useableLGap) continue ;
        delta = lGap ;
	revisebi = false ;
      } else {
        if (!useableUGap) continue ;
        delta = uGap ;
      }
      deltaMul = colLower[jProbe]-1.0 ; 
    }
/*
  Now decide if we have something that cuts off the current solution. Augment
  the lhs activity with the contribution for a'<ip>x*<p>, calculate the new
  rhs value, and compare.
  As an alternate measure of effectiveness, consider the gap between the
  current activity and the revised lhs bound, normalised by the gap between
  the original rhs and the revised lhs bound. 

  We'll also generate the cut if we can strengthen it in place.

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
    double violation = 0.0 ;
    bool genCut = false ;
    if (revisebi) {
      biPrime = rowUpper[irow]+deltaMul*delta ;
      violation = sum-biPrime ;
    } else {
      blowiPrime = rowLower[irow]+deltaMul*delta ;
      violation = blowiPrime-sum ;
    }
    double effectiveness = violation ;
    if (!aipNonZero && aipPrimeNonZero)
      effectiveness = CoinMax(effectiveness,aipZeroEffectiveness) ;
    if (effectiveness > needEffectiveness) genCut = true ;
    if (info->strengthenRow) genCut = true ;
/*
  Are we going to generate the cut? If not, on to the next iteration.
*/
    if (!genCut) continue ;

#   if CGL_DEBUG > 0
    if (violation <= 0) {
      std::cout
        << "  NOT VIOLATED (strengthenCoeff): need " << needEffectiveness
        << ", violation " << violation ;
      if (!aipNonZero && aipPrimeNonZero)
	std::cout << ", new aipPrime" ;
      if (info->strengthenRow)
	std::cout << ", strengthenRow" ;
      std::cout << "." << std::endl ;
    }
#   endif
/*
  Generate the coefficients. Copy all except a<ip>. Add a<ip>' at the end, if
  it's nonzero. (a<ip>' can be zero and still have a perfectly valid cut.)
*/
    int n = 0 ;
    for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
      int k = column[kk] ;
      double aik = rowElements[kk] ;
      if (k != jProbe) {
	index[n] = k ;
	element[n++] = aik ;
      }
    }
    if (aipPrimeNonZero) {
      index[n] = jProbe ;
      element[n++] = aipPrime ;
    }
/*
  Fill in a cut structure with the cut.
*/
    OsiRowCut rc ;
    rc.setLb(blowiPrime) ;
    rc.setUb(biPrime) ;
    rc.setEffectiveness(effectiveness) ;
    rc.setRow(n,index,element,false) ;
#   if CGL_DEBUG > 1
    printf("Strengthen Cut:\n") ;
    CglProbingDebug::dump_row(irow,rc.lb(),rc.ub(),nan(""),nan(""),
        &si,true,true,4,n,index,element,1.0e-10,colLower,colUpper) ;
    printf("Original Row:\n") ;
    int k = rowStart[irow] ;
    CglProbingDebug::dump_row(irow,rowLower[irow],rowUpper[irow],
        minR[irow],maxR[irow],&si,true,true,4,rowLen[irow],
	&column[k],&rowElements[k],1.0e-10,colLower,colUpper) ;
#   endif
#   if CGL_DEBUG > 0
    if (debugger) assert(!debugger->invalidCut(rc)); 
#   endif
/*
  If we're in preprocessing, we might try to simply replace the existing
  constraint (justReplace = true). Otherwise, drop the cut into the cut set.

  realRows comes in as a parameter. This is the translation array created if
  we modified the constraint system during preprocessing in gutsOfGenerateCuts.

  Effectiveness, if we're strengthening in place, seems to be absolute size of
  coefficients; smaller is better. (Why? Proxy for fewer coefficients?)
*/
      int realRow = irow ;
      if (realRows)
	realRow = realRows[realRow] ;
      assert(realRow >= 0) ;
      rc.setWhichRow(realRow) ;
      if (!justReplace) {
#       if CGL_DEBUG > 1
	  std::cout
	    << "    strengthen coeff cut on real row "
	    << realRow << " added to cut collection." << std::endl ;
	  newCuts++ ;
#       endif
	rowCut.addCutIfNotDuplicate(rc) ;
      } else {
	double effectiveness = 0.0 ;
	for (int i = 0 ; i < n ; i++)
	  effectiveness += CoinAbs(element[i]) ;
	if (!info->strengthenRow[realRow] ||
	    info->strengthenRow[realRow]->effectiveness() > effectiveness) {
	  delete info->strengthenRow[realRow] ;
	  rc.setEffectiveness(effectiveness) ;
#         if CGL_DEBUG > 1
	  std::cout
	    << "    strengthen coeff on real row " << realRow
	    << " will replace row, eff " << effectiveness << "."
	    << std::endl ;
	  newCuts++ ;
#         endif
	  info->strengthenRow[realRow] = rc.clone() ;
	} else {
#         if CGL_DEBUG > 1
	    std::cout
	      << "    strengthen coeff cut on real row " << realRow
	      << " rejected; eff " << effectiveness << " >  eff "
	      << info->strengthenRow[realRow]->effectiveness()
	      << " of incumbent." << std::endl ;
#         endif
	}
      }
    }

# if CGL_DEBUG > 1
  int cutsAtEnd = rowCut.numberCuts() ;
  int outplaceCuts = cutsAtEnd-cutsAtStart ;
  int inplaceCuts = newCuts-outplaceCuts ;
  std::cout << "Leaving strengthenCoeff" ;
  if (newCuts > 0)
    std::cout << ", " << newCuts << " cuts, " << inplaceCuts << " in place" ;
  std::cout << "." << std::endl ;
# endif

  return ;
}




// =========================================================


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
      This is six specialised replications of the same functionality: walk
      the column of a variable from stackC, and for each affected row, try to
      tighten the bounds of other variables in the row. If we successfully
      tighten bounds on a variable, walk the column of that variable,
      tighten the bounds of affected rows, and place the variable on stackC.

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

  minR, maxR, markR  Initialised externally by calcRowBounds. For row i,
  		     minR[i] = LB(i), maxR[i] = UB(i) (row lhs lower and
		     upper bounds, respectively).  markR is set to -1 if
		     there's at least one finite lhs bound, -2 if we shouldn't
		     process this row (no finite bounds, or too many
		     coefficients).

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
bool CglProbing::probe (const OsiSolverInterface &si, 
		        OsiCuts &cs, 
		        double *const colLower, double *const colUpper, 
		        const CoinPackedMatrix *const rowCopy,
		        const CoinPackedMatrix *const columnCopy,
		        const CoinBigIndex *const rowStartPos,
		        const int *const realRows, 
		        const double *const rowLower,
		        const double *const rowUpper,
		        const char *const intVar,
		        double *const minR, double *const maxR, 
		        int *const markR, 
                        CglTreeInfo *const info,
		        bool useObj, bool useCutoff, double cutoff) const

{
# if CGL_DEBUG > 0
  std::cout << "Entering CglProbing::probe." << std::endl ;
  int numberRowCutsBefore = cs.sizeRowCuts() ;
  int numberColCutsBefore = cs.sizeColCuts() ;
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

  If what we're doing is strengthening rows in place (option 0x40 indicates
  we're in preprocessing, a good place to do this), we can generate at most
  one cut per original row (nRows) and we need the vector strengthenRow
  to pass the cut back to the caller.

  TODO: I'm asking myself `Why do we need a cut container at all, if we're
  	in `just replace' mode?'. The cuts are kept in strengthenRow. But
	at the end we copy them over into cs anyhow. This can probably be
	rationalised.  -- lh, 110209 --

        Maybe not. Now that I better understand what CglPreProcessing is doing
	with the result, there's some possibility of returning implication
	cuts in rowCut while returning coefficient strengthening cuts in
	strengthenRow. If we have 0x40 and no strengthenRow, we're basically
	limited to column cuts.   -- lh, 110305 --
*/
  bool justReplace =
      ((info->options&0x40) != 0) && (info->strengthenRow != NULL) ;
  int cutCapacity = 0 ;
  if (justReplace) {
    cutCapacity = nRows ;
  } else {
    int nRowsSafe = CoinMin(nRows,si.getNumRows()) ;
    cutCapacity = info->inTree ? nRowsSafe/3 : nRowsSafe*10 ;
    if (!info->inTree && info->pass == 0) cutCapacity *= 10 ;
  }
  CglProbingRowCut rowCut(cutCapacity,!info->inTree) ;


  // Unpack matrices
  const int *column = rowCopy->getIndices() ;
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *rowLength = rowCopy->getVectorLengths() ;
  const double *rowElements = rowCopy->getElements() ;

  const int *row = columnCopy->getIndices() ;
  const CoinBigIndex *columnStart = columnCopy->getVectorStarts() ;
  const int *columnLength = columnCopy->getVectorLengths() ;
  const double *columnElements = columnCopy->getElements() ;

/*
  Grab the column solution and reduced costs and groom them.
*/
  CoinMemcpyN(si.getReducedCost(),nCols,djs) ;
  CoinMemcpyN(si.getColSolution(),nCols,colsol) ;
  double direction = si.getObjSense() ;
  groomSoln(direction,primalTolerance_,djs,colLower,colsol,colUpper,columnGap,
	    rowCopy,rowStartPos,largestNegativeInRow,largestPositiveInRow) ;
/*
  Scan the variables, noting the ones that are fixed or have a bound too large
  to be useful.
*/
  for (int i = 0 ; i < nCols ; i++) {
    if (colUpper[i]-colLower[i] < 1.0e-8) {
      markC[i] = tightenLower|tightenUpper ;
    } else {
      markC[i] = 0 ;
      if (colUpper[i] > 1.0e10)
	markC[i] |= infUpper ;
      if (colLower[i] < -1.0e10)
	markC[i] |= infLower ;
    }
  }
/*
  jjf: If we are going to replace coefficient then we don't need to be
       effective

  Seems like this comment does not agree with CglProbingRowCut.addCuts(),
  which checks effectiveness before adding a cut to the OsiCuts collection or
  entering it into the strengthenRow array.

  But it does agree with the way coefficient strengthening cuts are handled
  down in the end-of-iteration processing for the down/up/down probe loop.

  The idea of strengthenRow is that a cut that's simply strengthening an
  existing row i is entered in slot i of strengthenRow. It's left to the
  client to deal with it on return. I need to get a better idea of how
  justReplace is handled, here and in the client.  -- lh, 101125 --

*/
  //double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3 ;
  double needEffectiveness = info->strengthenRow ? 1.0e-8 : 1.0e-3 ;
  if (justReplace && ((info->pass&1) != 0))
    needEffectiveness = -1.0e10 ;

  double tolerance = 1.0e1*primalTolerance_ ;

  /* for both way coding */
  int nstackC0 = -1 ;
  int nstackR,nstackC ;

# if CGL_DEBUG > 0
  std::cout
    << "    prior to probe loop, justReplace "
    << ((justReplace)?"true":"false") << ", required effectiveness "
    << needEffectiveness << "." << std::endl ;
# endif

/*
  PASS LOOP: HEAD

  Major block #2: Main pass loop.

  anyColumnCuts is set only in the case that we've
  fixed a variable by probing (i.e., one of the up or down probe resulted
  in infeasibility) and that probe entailed column cuts. Once set, it is
  never rescinded. In the reworked code, it's set as the return value of
  monotoneActions().
*/
  bool anyColumnCuts = false ;
  int ninfeas = 0 ;
  int rowCuts = rowCuts_ ;
  double disaggEffectiveness = 1.0e-3 ;
  int iPass = 0 ;
  int nfixed = -1 ;
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_ ;

  while (iPass < maxPass && nfixed) {
    iPass++ ;
    nfixed = 0 ;
#   if CGL_DEBUG > 0
    std::cout
      << "  probe loop: starting pass " << iPass << "." << std::endl ;
#   endif
/*
  We went through a fair bit of trouble in gutsOfGenerateCut to determine
  the set of variables to be probed and loaded the indices into lookedAt_;
  numberThisTime_ reflects that.

  The net effect here is that the root gets special treatment
  (maxProbeRoot_) and the first pass at the root gets extra special treatment
  (numberThisTime_).

  See CbcCutGenerator::generateCuts for the origin of magic number 123. 
  Comments there indicate it's intended to take effect deep in the tree, but
  the code here will also work at the root. This appeared in r629.

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
    int iLook ;
    for (iLook = 0 ; iLook < numberThisTime_ ; iLook++) {
/*
  Too successful? Consider bailing out.

  If we're generating row cuts, but haven't fixed any variables or we're in
  the tree, break from probing.

  JustFix < 0 says this is the first pass at the root; for any other
  iteration of the pass loop, it'll be initialised to 0.  If it is the first
  pass at the root, turn off row cuts, keep on fixing variables, and stop
  at the end of this pass. Otherwise, if we haven't fixed any variables, break.

  TODO: Otherwise keep going with no changes? That doesn't seem right. The
  	logic here does not cover all situations. This bit of code appeared
	at r629.   -- lh, 101126 --
*/
      if (rowCut.outOfSpace() || leftTotalStack <= 0) {
	if (!justFix && (!nfixed || info->inTree)) {
#         if CGL_DEBUG > 0
	  if (!info->inTree)
	    std::cout
	      << "    Bailing (limit A) pass " << iPass
	      << ", maxProbe " << maxProbe << "." << std::endl ;
#         endif	  
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
#           if CGL_DEBUG > 0
	    std::cout
	      << "    Bailing (limit B) pass " << iPass
	      << ", maxProbe " << maxProbe << "." << std::endl ;
#           endif	  
	    break ;
	  }
	}
      }
/*
  Have a look at the candidate probing variable. It must be integer, but
  not necessarily binary. Failure is an algorithm error.

  It's possible that previous activity may have improved the bounds or
  even fixed the variable. We're not interested if the value in the primal
  solution is no longer within bounds -- this variable has already been
  swept up in bounds propagation. There's no need to try and use it to
  initiate propagation.

  The best way to understand the following blocks of code is to draw a
  picture, marking out +/- tolerance to either side of l<j> or u<j> and
  look at the values for down and up.
*/
      int j = lookedAt_[iLook] ;
      assert(intVar[j]) ;
      double xj = colsol[j] ;
      double down = floor(xj+tolerance) ;
      double up = ceil(xj-tolerance) ;
      const double uj = colUpper[j] ;
      const double lj = colLower[j] ;
      const bool binaryProbeVar = (lj == 0.0 && uj == 1.0) ;
      const bool probeIsDichotomy = ((uj-lj) < 1.5) ;
      const double gap = columnGap[j] ;
#     if CGL_DEBUG > 1
      std::cout
	<< "  probe: looking at " << " x(" << j << ") = " << xj
	<< ", l = " << lj << ", u = " << uj << "." << std::endl ;
      if (gap < 1.0e-8)
        assert(markC[j] == (tightenUpper|tightenLower)) ;
      if ((xj > uj+tolerance) || (xj < lj-tolerance))
        std::cout << "    discard; out-of-bound." << std::endl ;
      else if (columnGap[j] < 1.0e-8)
        std::cout << "    discard: fixed." << std::endl ;
#     endif
      if (gap < 1.0e-8 || (xj > uj+tolerance) || (xj < lj-tolerance))
        continue ;
/*
  `Normalize' variables that are near their bounds or interior integral; we
  can spot them easily because down == up. We want to make sure we have a
  spread of 1 between down and up, within the current bounds.
*/
      if (down == up) {
        if (down == lj) {
	  // at l<j>
	  up = down+1 ;
	} else if (down == uj) {
	  // at u<j>
	  down = up-1 ;
	} else {
	  // interior integer; arbitrarily go up
	  up = down+1 ;
	}
      }
#     if CGL_DEBUG > 0
      if (down < lj || up > uj || down == up) {
	std::cout
	  << "ERROR! Up and down probe values out of bounds or equal!"
	  << std::endl
	  << "    x<" << j << "> " << xj << ", l " << lj
	  << ", u " << uj << ", down " << down << ", up " << up
	  << "." << std::endl ;
      }
#     endif
      assert(down >= lj && up <= uj && ((up-down) == 1)) ;

#     if CGL_DEBUG > 0
      const bool downIsLower = (CoinAbs(down-lj) < 1.0e-7) ;
      const bool upIsUpper = (CoinAbs(up-uj) < 1.0e-7) ;
#     endif
/*
  Set up to probe each variable down (1), then up (2).

  The notion is that we'll set a bit (1, 2) in the result record if we're
  feasible for a particular way. Way is defined as:
    1: we're doing a down probe, u<j> reduced to floor(x*<j>)
    2: we're doing a up probe, l<j> increased to ceil(x*<j>)
  As defined below, movement (new bound - old bound) is negative for a down
  probe, positive for an up probe.
*/
      unsigned int iWay ;
      unsigned int way[] = { probeDown, probeUp } ;
      unsigned int feasValue[] =
          { downIterFeas, upIterFeas } ;
      unsigned int feasRecord = 0 ;
      bool notFeasible = false ;
      int istackC = 0 ;
/*
  PROBE LOOP: HEAD

  Open a loop to probe up and down for the current variable. Get things
  started by placing the variable on the propagation queue (stackC).

  As with the previous loops (variables in lookedAt_, passes) this loop
  extends to the end of major block #2).
*/
      for (iWay = downIter ; iWay < oneIterTooMany ; iWay++) {
        stackC[0] = j ;
        markC[j] = way[iWay] ;
/*
  Calculate movement given the current direction. There are various
  circumstances in which we need to know that the probe has tightened by
  one from the existing bound.
*/
        double solMovement ;
        double movement ;
	bool probeDistOne = true ;

#	if CGL_DEBUG > 0
	assert(lj == colLower[j]) ;
	assert(uj == colUpper[j]) ;
	assert(downIsLower == (down == lj)) ;
	assert(upIsUpper == (up == uj)) ;
#	endif

        if (way[iWay] == probeDown) {
          movement = down-uj ;
          solMovement = down-xj ;
          assert(movement < -0.99999) ;
	  if (movement < -1.5)
	    probeDistOne = false ;
        } else {
          movement = up-lj ;
          solMovement = up-xj ;
          assert(movement > 0.99999) ;
	  if (movement > 1.5)
	    probeDistOne = false ;
        }
/*
  About those `various circumstances'. Disaggregation cuts are mathematically
  correct only across the dichotomy used to derive the cut. So we're in
  trouble if there's any polytope outside those points.

  Implication cuts are also only valid across a dichotomy. The way the
  current method is written, the row cuts are only valid for binary variables.
  The variable controls row cuts, distinct from consensus bound changes.

  Coefficient strengthening is a bit more subtle. The cut is always valid,
  but if we're replacing (tilting) an existing constraint, there cannot be
  any polytope on the side that's getting looser!
*/
	const bool doDisaggCuts =
	    (((rowCuts&0x04) == 0) && ((rowCuts&0x01) != 0) &&
	     !justReplace && probeIsDichotomy) ;
	const bool doCoeffStrengthen =
	    (((rowCuts&0x04) == 0) && ((rowCuts&0x02) != 0) &&
	     (!justReplace || (justReplace && probeDistOne))) ;
	const bool doSingletons = true ;
	const bool doImplicationCuts =
	    (((rowCuts&0x04) == 0) && !justReplace && binaryProbeVar) ;

#       if CGL_DEBUG > 1
	if (iWay == downIter)
	  std::cout
	    << "    down probe, old u " << colUpper[j]
	    << ", new u " << colUpper[j]+movement ;
	else
	  std::cout
	    << "    up probe, old l " << colLower[j]
	    << ", new l " << colLower[j]+movement ;
	std::cout
	    << ", move " << movement << ", soln " << colsol[j]
	    << "(" << xj << "), move " << solMovement << "." << std::endl ;
#       endif
/*
  Recall that we adjusted the reduced costs (djs) and current objective
  (current) to the minimisation sign convention. objVal will accumulate
  objective change due to forced bound changes.
*/
        double objVal = si.getObjValue() ;
        if (solMovement*djs[j] > 0.0)
          objVal += solMovement*djs[j] ;
        nstackC = 1 ;
        nstackR = 0 ;
        saveL[0] = colLower[j] ;
        saveU[0] = colUpper[j] ;
        assert (saveU[0] > saveL[0]) ;
        notFeasible = false ;
/*
  Set the bounds and gap for the probe. Clean up the markC flags.
*/
        if (movement < 0.0) {
          colUpper[j] = down ;
        } else {
          colLower[j] = up ;
        }
	columnGap[j] = colUpper[j]-colLower[j]-primalTolerance_ ;
        if (CoinAbs(colUpper[j]-colLower[j]) < 1.0e-6) {
          markC[j] = tightenUpper|tightenLower ;
	} else {
	  markC[j] &= ~(infUpper|infLower) ;
	  if (colUpper[j] > 1.0e10)
	    markC[j] |= infUpper ;
	  if (colLower[j] < -1.0e10)
	    markC[j] |= infLower ;
	}
        istackC = 0 ;
/*
  Update row bounds to reflect the change in variable bound.
*/
	if (!updateRowBounds(j,movement,columnStart,
			     columnLength,row,columnElements,
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
*/
        while (istackC < nstackC && nstackC < maxStack && !notFeasible) {
	  leftTotalStack-- ;
          int jway ;
          int jcol = stackC[istackC] ;
          jway = markC[jcol] ;
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
	  const CoinBigIndex jjstart = columnStart[jcol] ;
	  const CoinBigIndex jjend = jjstart+columnLength[jcol] ;
          for (CoinBigIndex jj = jjstart ; jj < jjend ; jj++) {
            if (notFeasible)
              break ;
            int irow = row[jj] ;
	    // Row already processed with these bounds.
	    if (markR[irow]/10000 > 0) continue ;

	    CoinBigIndex iistart = rowStart[irow] ;
	    const CoinBigIndex iistartPos = rowStartPos[irow] ;
	    const CoinBigIndex iiend = iistart+rowLength[irow] ;
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
	    double rowUp = rowUpper[irow] ;
	    double rowUp2 = 0.0 ;
	    bool doRowUpN ;
	    bool doRowUpP ;
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
*/
#           if CGL_DEBUG > 0
	    if (markR[irow]/10000 > 1)
	      std::cout
		<< "REDUNDANT: processing " << irow
		<< " for the " << markR[irow]/10000 << " time." << std::endl ;
#           endif
/*
  PROBE LOOP: WALK ROW

  Given the above analysis, work on the columns with negative coefficients.

  doRowUpN: a<ij> < 0, minR can violate rowUp => increase l<j> (more negative)
  doRowLoN: a<ij> < 0, maxR can violate rowLo => decrease u<j> (more positive)
*/
	    if (doRowUpN && doRowLoN) {
	      //doRowUpN=doRowLoN=false ;
	      // Start neg values loop
	      for (CoinBigIndex ii = iistart ; ii < iistartPos ; ii++) {
		const int kcol = column[ii] ;
		int markIt = markC[kcol] ;
		// Skip columns with fixed variables.
		if ((markIt&(tightenLower|tightenUpper)) !=
		                    (tightenLower|tightenUpper)) {
		  double value2 = rowElements[ii] ;
#                 if CGL_DEBUG > 0
		  if (colUpper[kcol] <= 1e10)
		    assert ((markIt&infUpper) == 0) ;
		  else
		    assert ((markIt&infUpper) != 0) ;
		  if (colLower[kcol] >= -1e10)
		    assert ((markIt&infLower) == 0) ;
		  else
		    assert ((markIt&infLower) != 0) ;
		  assert (value2 < 0.0) ;
#                 endif
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " l " << saveL[nstackC-1]
			<< " -> " << newLower << "; i " << irow << std::endl ;
#		      endif
/*
  Propagate the column bound change.

  TODO: Notice that we're not setting a variable here when we discover
        infeasibility; rather, we're setting the column bound to a ridiculous
	value which will be discovered farther down.  -- lh, 101127 --
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " u " << saveU[nstackC-1]
			<< " -> " << newUpper << "; i " << irow << std::endl ;
#		      endif
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
		      jj = jjend ;
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
	      for (CoinBigIndex ii = iistart ; ii < iistartPos ; ii++) {
		int kcol =column[ii] ;
		int markIt=markC[kcol] ;
		if ((markIt&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[ii] ;
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " l " << saveL[nstackC-1]
			<< " -> " << newLower << "; i " << irow << std::endl ;
#		      endif
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
		      jj = jjend ;
		      istackC = nstackC+1 ;
		      break ;
		    }
		  }
		}
	      } // end big loop iistart->rPos
	    } else if (doRowLoN) {
/*
  And yet again, for the case where we can only work against the lower bound.
*/
	      // Start neg values loop
	      for (CoinBigIndex ii = iistart ; ii < iistartPos ; ii++) {
		int kcol =column[ii] ;
		if ((markC[kcol]&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double moveDown=0.0 ;
		  double newUpper=-1.0 ;
		  double value2=rowElements[ii] ;
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " u " << saveU[nstackC-1]
			<< " -> " << newUpper << "; i " << irow << std::endl ;
#		      endif
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
		      jj = jjend ;
		      istackC = nstackC+1 ;
		      break ;
		    }
		  }
		}
	      }  // end big loop iistart->rPos
	    }
/*
  We've finished working on the negative coefficients of the row. Advance
  iistart to cover the positive coefficients and repeat the previous 500 lines.
*/
	    iistart = iistartPos ;
	    if (doRowUpP&&doRowLoP) {
	      //doRowUpP=doRowLoP=false ;
	      // Start pos values loop
	      for (CoinBigIndex ii = iistart ; ii < iiend ; ii++) {
		int kcol=column[ii] ;
		int markIt=markC[kcol] ;
		if ((markIt&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[ii] ;
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " l " << saveL[nstackC-1]
			<< " -> " << newLower << "; i " << irow << std::endl ;
#		      endif
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " u " << saveU[nstackC-1]
			<< " -> " << newUpper << "; i " << irow << std::endl ;
#		      endif
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
		      jj = jjend ;
		      istackC=nstackC+1 ;
		      break ;
		    }
		  }
		}
	      } // end big loop rPos->iiend
	    } else if (doRowUpP) {
	      // Start pos values loop
	      for (CoinBigIndex ii = iistart ; ii < iiend ; ii++) {
		int kcol =column[ii] ;
		int markIt=markC[kcol] ;
		if ((markIt&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[ii] ;
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " u " << saveU[nstackC-1]
			<< " -> " << newUpper << "; i " << irow << std::endl ;
#		      endif
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
		      jj = jjend ;
		      istackC=nstackC+1 ;
		      break ;
		    }
		  }
		}
	      } // end big loop rPos->iiend
	    } else if (doRowLoP) {
	      // Start pos values loop
	      for (CoinBigIndex ii = iistart ; ii < iiend ; ii++) {
		int kcol =column[ii] ;
		if ((markC[kcol]&(tightenLower|tightenUpper))!= (tightenLower|tightenUpper)) {
		  double value2=rowElements[ii] ;
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
#		      if CGL_DEBUG > 2
		      std::cout
			<< "        " << nstackC-1
			<< " col " << kcol << " " << colsol[kcol]
			<< " l " << saveL[nstackC-1]
			<< " -> " << newLower << "; i " << irow << std::endl ;
#		      endif
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
		      jj = jjend ;
		      istackC = nstackC+1 ;
		      break ;
		    }
		  }
		}
	      }      // end big loop rPos->iiend
	    }    // end processing of positive coefficients of row.
          }   // end loop to walk the column of a variable popped off stackC
          istackC++ ;
        }  // end stackC processing loop
#       if CGL_DEBUG > 1
	std::cout
	  << "    " << ((notFeasible)?"infeasible":"feasible")
	  << ", stacked " << nstackC-1 << " vars." << std::endl ;
# 	endif
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
*/
	if (notFeasible || (useCutoff && (objVal > cutoff))) {
#	  if CGL_DEBUG > 1
	  std::cout
	    << "    " << ((way[iWay] == probeDown)?"down":"up")
	    << " probe for x<" << j << "> infeasible on " ;
	  if (notFeasible && (useCutoff && (objVal > cutoff)))
	    std::cout << "primal bounds and objective cutoff" ;
	  else if (notFeasible)
	    std::cout << "primal bounds" ;
	  else
	    std::cout << "objective cutoff" ;
	  std::cout << "." << std::endl ;
#	  endif
	  notFeasible = true ;
	  if (iWay == upIter && feasRecord == 0) {
	    ninfeas = 1 ;
	    j = nCols-1 ;
	    iLook = numberThisTime_ ;
	    iPass = maxPass ;
	    break ;
	  }
	} else {
	  feasRecord |= feasValue[iWay] ; 
	}
/*
  PROBE LOOP: DOWN PROBE

  If this is the down probe, generate cuts if we're feasible and then set up
  to iterate for the up probe.
*/
	if (iWay == downIter) {
          if (!notFeasible) {
	    if (doSingletons)
	      singletonRows(j,primalTolerance_,useObj,si,rowCopy,markC,
			    nstackC,stackC,saveL,saveU,
			    colUpper,colLower,columnGap,
			    nstackR,stackR,rowUpper,rowLower,maxR,minR) ;
	    if (doDisaggCuts)
	      disaggCuts(nstackC,probeDown,primalTolerance_,
			 disaggEffectiveness,si,rowCut,stackC,colsol,
			 colUpper,colLower,saveU,saveL,index,element) ;
	    if (doCoeffStrengthen && !anyColumnCuts)
	      strengthenCoeff(j,probeDown,primalTolerance_,justReplace,
			      needEffectiveness,si,rowCut,
			      rowCopy,colUpper,colLower,colsol,nstackR,stackR,
			      rowUpper,rowLower,maxR,minR,realRows,
			      element,index,info) ;
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

          restoreRowColBounds(primalTolerance_,nstackC,stackC,markC,
	  		      saveL,saveU,colLower,colUpper,columnGap,
			      nstackR,stackR,markR,
			      saveMin,saveMax,minR,maxR) ;
	  continue ;
	}
/*
  PROBE LOOP: MONOTONE

  We've done a down probe and an up probe. Outright infeasibility was
  handled above, as was iteration after the down probe. We can only reach
  here after the up probe, with at least one direction feasible. Check
  first for monotonicity.
  * If we are feasible up but not down probe, we're monotone up.
  * If we were feasible down but not up, we're monotone down. This requires a
    bit more work, as we need to back out the up probe state and reinstall the
    down probe state.
  Either way, we will capture the results and move on to the next variable.
*/
        assert(iWay == upIter &&
	       ((feasRecord&(downIterFeas|upIterFeas)) != 0)) ;

        if (feasRecord != (downIterFeas|upIterFeas)) {
	  bool retVal ;
          nfixed++ ;
	  if (feasRecord == downIterFeas) {
	    restoreRowColBounds(primalTolerance_,nstackC,stackC,markC,
	    			saveL,saveU,colLower,colUpper,columnGap,
				nstackR,stackR,markR,
				saveMin,saveMax,minR,maxR) ;
	    restoreRowColBounds(primalTolerance_,nstackC0,stackC0,markC,
	    			lo0,up0,colLower,colUpper,columnGap,
				0,0,0,0,0,0,0) ;
	    retVal = monotoneActions(primalTolerance_,si,cs,
				     nstackC0,stackC0,intVar,
				     colLower,colsol,colUpper,
				     index,element) ;
	  } else {
	    retVal = monotoneActions(primalTolerance_,si,cs,
	  			     nstackC,stackC,intVar,
				     colLower,colsol,colUpper,
				     index,element) ;
	  }
	  if (retVal) anyColumnCuts = true ;
	  clearStacks(primalTolerance_,nstackC,stackC,markC,colLower,colUpper,
	  	      nstackR,stackR,markR) ;
	  break ;
        }
/*
  PROBE LOOP: DOWN AND UP FEASIBLE

  To reach here, this must be the up probe and both up and down probes are
  feasible. Generate row singletons and disaggregation cuts (just as for a
  feasible down probe), then call implicationCuts to look for cuts based on
  both probes (binary probe only).

  Then restore the bounds, because we can't keep these.
*/
	assert(iWay == upIter && (feasRecord == (downIterFeas|upIterFeas))) ;

	if (doSingletons)
	  singletonRows(j,primalTolerance_,useObj,si,rowCopy,markC,
			nstackC,stackC,saveL,saveU,
			colUpper,colLower,columnGap,
			nstackR,stackR,rowUpper,rowLower,maxR,minR) ;
	if (doDisaggCuts)
	  disaggCuts(nstackC,probeUp,primalTolerance_,
		     disaggEffectiveness,si,
		     rowCut,stackC,colsol,colUpper,colLower,saveU,saveL,
		     index,element) ;
	if (doCoeffStrengthen && !anyColumnCuts)
	  strengthenCoeff(j,probeUp,primalTolerance_,justReplace,
			  needEffectiveness,si,rowCut,
			  rowCopy,colUpper,colLower,colsol,nstackR,stackR,
			  rowUpper,rowLower,maxR,minR,realRows,
			  element,index,info) ;
	if (probeIsDichotomy) {
	  int jProbe = j ;
	  if (!doImplicationCuts) jProbe = -1 ;
	  implicationCuts(jProbe,primalTolerance_,nCols,cs,si,intVar,colsol,
			  nstackC,stackC,markC,colUpper,colLower,
			  saveU,saveL,nstackC0,stackC0,up0,lo0,
			  element,index) ;
        }
	restoreRowColBounds(primalTolerance_,nstackC,stackC,markC,
			    saveL,saveU,colLower,colUpper,columnGap,
			    nstackR,stackR,markR,
			    saveMin,saveMax,minR,maxR) ;

      }     // PROBE LOOP: END
    }    // LOOKEDAT LOOP: END
#   if CGL_DEBUG > 0
    std::cout
      << "  probe: ending pass " << iPass
      << ", fixed " << nfixed << " vars." << std::endl ;
#   endif
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
  replace' mode, they will be carried back to the caller in info.strengthenRow.
*/
  if (!ninfeas) {
    if (!justReplace) {
      // INPLACE rowCut.addCuts(cs,info->strengthenRow,info->pass) ;
      rowCut.addCuts(cs,0,info->pass) ;
    } else {
      assert(info->strengthenRow != 0) ;
      for (int i = 0 ; i < nRows ; i++) {
	int realRow = (realRows)?realRows[i]:i ;
	assert(realRow >= 0 && realRow < si.getNumRows()) ;
	OsiRowCut *cut = info->strengthenRow[realRow] ;
	if (cut) {
#         if CGL_DEBUG > 1
	  std::cout
	    << "    row " << realRow << " (local " << i
	    << ") strengthened, effectiveness " << cut->effectiveness()
	    << "." << std::endl ;
#         endif
	  // INPLACE cs.insert(cut) ;
	}
      }
    }
  }
# if CGL_DEBUG > 0
  int numberRowCutsAfter = cs.sizeRowCuts() ;
  int numberColCutsAfter = cs.sizeColCuts() ;
  int inPlaceCuts = 0 ;
  if (justReplace) {
    for (int i = 0 ; i < nRows ; i++)
      if (info->strengthenRow[i] != 0) inPlaceCuts++ ;
  }

  std::cout
    << "Leaving CglProbing::probe, ninfeas " << ninfeas << ", "
    << (numberRowCutsAfter-numberRowCutsBefore) << " row cuts, "
    << inPlaceCuts << " inplace, "
    << (numberColCutsAfter-numberColCutsBefore) << " col cuts."
    << std::endl ;
 
  if ((numberRowCutsAfter-numberRowCutsBefore) > 0) {
    std::cout << "  row cuts:" << std::endl ;
    for (int k = numberRowCutsBefore ; k < numberRowCutsAfter ; k++) {
      OsiRowCut thisCut = cs.rowCut(k) ;
      thisCut.print() ;
    }
  }
  if ((numberColCutsAfter-numberColCutsBefore) > 0) {
    std::cout << "  col cuts:" << std::endl ;
    for (int k = numberColCutsBefore ; k < numberColCutsAfter ; k++) {
      OsiColCut thisCut = cs.colCut(k) ;
      thisCut.print() ;
    }
  }
# endif

  return (((ninfeas)?false:true)) ;
}
