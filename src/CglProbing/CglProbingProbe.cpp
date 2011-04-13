
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
#include "CglPhic.hpp"
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

/*
  Calculate the change in objective due to bound changes that move a variable
  past the current solution (column cuts). Return true if column cuts are
  present, as this affects the generation of other cuts. Note that the probe
  variable (newBnds[0]) doesn't count in this determination.
*/
bool calcObjChg (double &z, double feasTol,
		 int bndLen, const CglPhic::CglPhicBndPair *const &newBnds,
		 const double *const colsol,
		 const double *const djs)
{
  bool colCuts = false ;
  double deltaz = 0 ;
  const int &p = newBnds[0].ndx_ ;
  if (newBnds[0].lb_-colsol[p] > feasTol) {
    deltaz += (newBnds[0].lb_-colsol[p])*djs[p] ;
  } else if (colsol[p]-newBnds[0].ub_ > feasTol) {
    deltaz -= (colsol[p]-newBnds[0].ub_)*djs[p] ;
  }

  for (int k = 1 ; k < bndLen ; k++) {
    const CglPhic::CglPhicBndPair &chgRec = newBnds[k] ;
    const int &j = chgRec.ndx_ ;
    const double &lj = chgRec.lb_ ;
    const double &uj = chgRec.ub_ ;
    const double &solj = colsol[j] ;
    if (lj-solj > feasTol) {
      colCuts = true ;
      deltaz += (lj-solj)*djs[j] ;
    } else if (solj-uj > feasTol) {
      colCuts = true ;
      deltaz -= (solj-uj)*djs[j] ;
    }
  }
  z += deltaz ;
  return (colCuts) ;
}

/*
  jjf: clean up djs and solution

  Do some cleanup of reduced costs.
*/
void groomSoln (double direction, double primalTolerance, int nCols,
		double *const djs,
		const double *const colLower,
		double *const colsol,
		const double *const colUpper)
{
/*
  Convert reduced costs to minimisation and clean up any that are on the wrong
  side of zero.
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
  }
  return ;
}

}  // end file-local namespace



/*
  This method looks for implication cuts: implications that can be generated
  given a binary probing variable x<p>. For an arbitrary binary variable x<j>,
  we're looking for things like the following:
    1) x<p> -> 0 ==> x<j> -> 0 and x<p> -> 1 ==> x<j> -> 0
       In this case, we can fix x<j> to 0. Similarly if x<j> is driven to 1.
    2) x<p> -> 0 ==> x<j> -> 0 and x<p> -> 1 ==> x<j> -> 1
       In this case, we can write x<j> = x<p>. If the sense is reversed, we
       get x<j> = (1-x<p>).
  Form 2) are actual cuts; form 1) just amounts to a bound change. Bound
  changes are installed in the CglPhic object. The return value will be true
  if we identify column cuts (this will require recalculation of row lhs
  bounds).

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

  1) is applied to all variables (integer and continuous) in the change
  record and consensus bounds are installed for all. Only integer variables
  qualify for inclusion in a column cut.

  The method uses the probe index, p, to decide if it should generate form
  2) cuts. If p >= 0, it's assumed to be the index of a binary variable
  and form 2) cuts are attempted.

  It's assumed that at the time this method is called, the up probe state has
  been reverted. The CglPhic object should contain pre-probe bounds. Consensus
  bound improvements are written directly to colLower, colUpper, which implies
  that the row lhs bounds in the CglPhic object will need to be recalculated
  before any further propagation activity.

  Clearly, if the only entry in either the up or down change record is the
  probe variable, we aren't going to see any consensus.

  TODO: Check that the above comment is correct --- do improvements to
	continuous bounds really propagate back to the caller?
	-- lh, 110201 --
*/

int CglProbing::implicationCuts (int p, int n,
		     OsiCuts &cs, const OsiSolverInterface &si,
		     CglPhic &phic,
		     const char *const intVar, const double *const colsol,
		     int numDnProbeBnds,
		     const CglPhic::CglPhicBndPair *const dnProbeBnds,
		     int numUpProbeBnds,
		     const CglPhic::CglPhicBndPair *const upProbeBnds) const
{
  double *colLower ;
  double *colUpper ;
  phic.getColBnds(colLower,colUpper) ;

# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout << "    "
      << "start implicationCuts, probe x<" << p << ">, " << numDnProbeBnds
      << " down records, " << numUpProbeBnds << " up records." << std::endl ;
  }
  int origRowCutCnt = cs.sizeRowCuts() ;
  int origColCutCnt = cs.sizeColCuts() ;
  char typLett[3] = {'c', 'b', 'g'} ;
/*
  If we have a binary probe variable, it should not be fixed. Also check that
  the first change record really is the probe variable.
*/
  if (p >= 0) {
    const CglPhic::CglPhicBndPair &dnPBnd_p = dnProbeBnds[0] ;
    const CglPhic::CglPhicBndPair &upPBnd_p = upProbeBnds[0] ;
    assert(dnPBnd_p.ndx_ == p && upPBnd_p.ndx_ == p) ;
    assert(colLower[p] == 0.0 && colUpper[p] == 1.0) ;
  }
# endif
/*
  If there's only one variable in either the up or down records, it's the
  probe variable and we aren't going to see any consensus bounds. Bail out now
  while it's still uncomplicated.
*/
  if (numDnProbeBnds == 1 || numUpProbeBnds == 1) {
#   if CGL_DEBUG > 0
    if (verbosity_ >= 3)
      std::cout << "    " << "end implicationCuts; trivial." << std::endl ;
#   endif
    return (0) ;
  }
/*
  Construct a cross-reference so that we can correlate entries in the down and
  up probe records. Then, as we scan the up probe records, we can quickly
  locate the record for x<j> (if it exists) in the down probe records. Don't
  process the probe variable.
*/
  int *up2down = new int [n] ;
  CoinZeroN(up2down,n) ;
  for (int kd = 1 ; kd < numDnProbeBnds ; kd++) {
    const CglPhic::CglPhicBndPair &dpRec_j = dnProbeBnds[kd] ;
    int j = dpRec_j.ndx_ ;
    up2down[j] = kd ;
  }
/*
  See if we have bounds improvement. Given a lower bound ld<j> from the
  down probe, a bound lu<j> from the up probe, and the original bound lo<j>,
  we can surely impose the lesser of ld<j> and lu<j> if this is more than the
  original lo<j>.  In math, max(min(ld<j>,lu<j>),lo<j>). Use a conservative
  tolerance. Improved bounds are written directly into colLower, colUpper.
  Improved bounds on an integer variable rate column cut. Make a note if the
  improved bound actually cuts off the current solution.
*/
  int nTot = 0 ;
  int nInt = 0 ;
  int cutsOffSoln = 0 ;
  const double minChange = 1.0e-4 ;
  int implIndices[2] ;
  double implCoeffs[2] ;
  CoinPackedVector consensus_lbs ;
  CoinPackedVector consensus_ubs ;
/*
  Scan the up records and look for variables that had bound improvement on
  both the up and down probe. If the variable was changed by only one probe,
  move on. Don't process the probe variable.
*/
  for (int ku = 1 ; ku < numUpProbeBnds ; ku++) {
    const CglPhic::CglPhicBndPair &upRec_j = upProbeBnds[ku] ;
    const int &j = upRec_j.ndx_ ;
    const int &kd = up2down[j] ;
    if (kd == 0) continue ;
    const CglPhic::CglPhicBndPair &dnRec_j = dnProbeBnds[kd] ;

    const double &ldj = dnRec_j.lb_ ;
    const double &udj = dnRec_j.ub_ ;
    const double &luj = upRec_j.lb_ ;
    const double &uuj = upRec_j.ub_ ;
    double &lj = colLower[j] ;
    double &uj = colUpper[j] ;
    const double &xj = colsol[j] ;
    const int typej = intVar[j] ;
/*
  Obtain the consensus lower bound. Revised bounds are written directly back
  to colLower, colUpper (via references lj, uj).
*/
    const double bestProbelj = CoinMin(ldj,luj) ;
    if (bestProbelj > lj+minChange) {
#     if CGL_DEBUG > 0
      if (verbosity_ >= 4) {
	std::cout
	  << "      consensus lb for (" << typLett[typej]
	  << ") x<" << j << "> improved from "
	  << lj << " to " << bestProbelj << "." << std::endl ;
      }
#     endif
      lj = bestProbelj ;
      nTot++ ;
      if (typej) {
        consensus_lbs.insert(j,bestProbelj) ;
	nInt++ ;
	if (xj < bestProbelj-primalTolerance_)
	  cutsOffSoln++ ;
      }
    } 
/*
  Obtain the consensus upper bound.
*/
    const double bestProbeuj = CoinMax(udj,uuj) ;
    if (bestProbeuj < uj-minChange) {
#     if CGL_DEBUG > 0
      if (verbosity_ >= 4) {
	std::cout
	  << "      consensus ub for (" << typLett[typej]
	  << ") x<" << j << "> improved from "
	  << uj << " to " << bestProbeuj << "." << std::endl ;
      }
#     endif
      uj = bestProbeuj ;
      nTot++ ;
      if (typej) {
        consensus_ubs.insert(j,bestProbeuj) ;
	nInt++ ;
	if (xj > bestProbeuj+primalTolerance_)
	  cutsOffSoln++ ;
      }
    } 
/*
  Now look for honest implication cuts (form 2 in the introductory comments).
  Recall that p >= 0 is taken as an indication that the probe variable
  is binary and these cuts should be attempted.  It's possible (miracles
  happen) that the above code has managed to fix x<j> by pulling together
  the bounds from the down and up probe but we should not lose feasibility
  here. In the case that x<j> is fixed, move on.
*/
    if (p >= 0) {
      assert(uj >= lj) ;
      if (CoinAbs(uj-lj) < primalTolerance_) continue ;
/*
  Check for x<j> = (u<j>-l<j>)x<p> + l<j>. In words, check that the down probe
  forced x<j> to l<j> and the up probe forced x<j> to u<j>.
*/
      if (udj < lj+primalTolerance_ && luj > uj-primalTolerance_) {
	OsiRowCut rc ;
	rc.setLb(lj) ;
	rc.setUb(lj) ;
	rc.setEffectiveness(1.0e-5) ;
	implIndices[0] = p ;
	implIndices[1] = j ;
	implCoeffs[0] = -(uj-lj) ;
	implCoeffs[1] = 1.0 ;
	rc.setRow(2,implIndices,implCoeffs,false) ;
#       if CGL_DEBUG > 0
        if (verbosity_ >= 4) {
	  std::cout << "      implication: " ;
	  rc.print() ;
	}
#       endif
	cs.insert(rc) ;
      } else if (uuj < lj+primalTolerance_ && ldj > uj-primalTolerance_) {
/*
  Repeat for x<j> = -(u<j>-l<j>)x<p> + u<j>
*/
	OsiRowCut rc ;
	rc.setLb(uj) ;
	rc.setUb(uj) ;
	rc.setEffectiveness(1.0e-5) ;
	implIndices[0] = p ;
	implIndices[1] = j ;
	implCoeffs[0] = uj-lj ;
	implCoeffs[1] = 1.0 ;
	rc.setRow(2,implIndices,implCoeffs,false) ;
#       if CGL_DEBUG > 0
	if (verbosity_ >= 4) {
	  std::cout << "      inverse implication: " ;
	  rc.print() ;
	}
#       endif
	cs.insert(rc) ;
      } 
    }
  }
  delete[] up2down ;
/*
  If we have column cuts, load them into the cut collection passed in as a
  parameter.
*/
  if (nInt) {
    OsiColCut cc ;
    if (consensus_lbs.getNumElements() > 0) {
      cc.setLbs(consensus_lbs) ;
    }
    if (consensus_ubs.getNumElements() > 0) {
      cc.setUbs(consensus_ubs) ;
    }
    if (cutsOffSoln) {
      cc.setEffectiveness(100.0) ;
    } else {
      cc.setEffectiveness(1.0e-5) ;
    }
#   if CGL_DEBUG > 0
    if (verbosity_ >= 3) {
      std::cout
	<< "      tightened "
	<< consensus_lbs.getNumElements() << " integer lb, "
	<< consensus_ubs.getNumElements()  << " integer ub" ;
      if (cutsOffSoln)
	std::cout << "; " << cutsOffSoln << " cuts." << std::endl ;
      else
	std::cout << "; no cut." << std::endl ;
    }
    if (paranoia_ > 0)
      assert(CglProbingDebug::checkBounds(si,cc,verbosity_) == 0) ;
#   endif
    cs.insert(cc) ;
  }
# if CGL_DEBUG > 1
  if (verbosity_ >= 3) {
    std::cout << "    "
      << "end implicationCuts, improved " << nTot << " bounds, "
      << cs.sizeRowCuts()-origRowCutCnt << " row cuts, "
      << cs.sizeColCuts()-origColCutCnt << " col cuts."
      << std::endl ;
  }
# endif

  if (nTot)
    return (true) ;
  else
    return (false) ;
}


/*
  Run through the set of variables with tightened bounds and create
  disaggregation cuts. See the typeset documentation for the full
  derivation. Note that the constraints generated here are a specialised
  form and require that the change in bound (u<p> to u'<p>, or l<p> to
  l'<p>) is distance 1.

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

  These cuts do not cast the probe value in concrete, so we can apply them
  whether or not the probe cuts off portions of the current solution.
*/
void CglProbing::disaggCuts (
		 int p, unsigned int probeDir, double disaggEffectiveness,
		 const OsiSolverInterface &si,
		 const double *const colsol,
		 int bndLen,
		 const CglPhic::CglPhicBndPair *const newBnds,
		 const CglPhic::CglPhicBndPair *const oldBnds,
		 CglProbingRowCut &rowCut,
		 int *const index, double *const element) const
{ 

# if CGL_DEBUG > 0
  if (paranoia_)
    assert((probeDir == probeDown) || (probeDir == probeUp)) ;
  if (verbosity_ >= 3) {
    std::cout << "    "
      << "start disaggCuts, "
      << ((probeDir == probeDown)?"down":"up") << " probe on "
      << "x<" << p << ">." << std::endl ;
  }

  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  int cutsAtStart = rowCut.numberCuts() ;
# endif


/*
  The change record for the probe variable should be the first entry in the
  vector.
*/
  const CglPhic::CglPhicBndPair &newBnd_p = newBnds[0] ;
  assert(p == newBnd_p.ndx_) ;

  int plusMinus = 0 ;
  double luPrime_p = 0.0 ;
  const double &x_p = colsol[p] ;
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
    luPrime_p = newBnd_p.ub_ ;
    plusMinus = -1 ;
  } else {
    luPrime_p = newBnd_p.lb_ ;
    plusMinus = 1 ;
  }
  deltaProbe_p = x_p-luPrime_p ;

# if CGL_DEBUG > 0
  if (verbosity_ >= 4) {
    const CglPhic::CglPhicBndPair &oldBnd_p = oldBnds[0] ;
    if (probeDir == probeDown) {
      std::cout
	<< "  u<" << p << "> = " << oldBnd_p.ub_
	<< " reduced to " << luPrime_p << "." << std::endl ;
    } else {
      std::cout
	<< "  l<" << p << "> = " << oldBnd_p.lb_
	<< " increased to " << luPrime_p << "." << std::endl ;
    }
    std::cout
      << "  x<" << p << "> = " << x_p << "; deltaProbe = "
      << deltaProbe_p << std::endl ;
  }
# endif
/*
  Change record 0 is the probing variable. Start with change record 1 for
  target variables x<t>.
*/
  for (int k = 1 ; k < bndLen ; k++) {
    const CglPhic::CglPhicBndPair &newBnd_t = newBnds[k] ;
    const CglPhic::CglPhicBndPair &oldBnd_t = oldBnds[k] ;
    int t = newBnd_t.ndx_ ;
    double u_t = oldBnd_t.ub_ ;
    double l_t = oldBnd_t.lb_ ;
    double x_t = colsol[t] ;
    double uPrime_t = newBnd_t.ub_ ;
    double lPrime_t = newBnd_t.lb_ ;
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
#     if CGL_DEBUG > 0
      if (verbosity_ >= 5) {
	std::cout
	  << "    Generating upper disaggregation cut for x<"
	  << t << "> (" << k << ")." << std::endl ;
	std::cout
	  << "    u<" << t << "> = " << u_t
	  << " decreased to " << uPrime_t << "." << std::endl ;
	std::cout
	  << "    x<" << t << "> = " << x_t
	  << " decreased to " << (uPrime_t-deltaProbe_p*a_p)
	  << "." << std::endl ;
	std::cout
	  << "    x<" << p << "> = " << x_p
	  << " changed to " << deltaCut_p+luPrime_p
	  << "; effectiveness = " << effectiveness
	  << "; required " << disaggEffectiveness << "." << std::endl ;
      }
#     endif
      assert(effectiveness+primalTolerance_ > 0) ;
      if (effectiveness > disaggEffectiveness) {
	OsiRowCut rc ;
	rc.setEffectiveness(effectiveness) ;
	rc.setLb(-DBL_MAX) ;
	rc.setUb(uPrime_t+a_p*luPrime_p) ;
	index[0] = t ;
	element[0] = 1.0 ;
	index[1] = p ;
	element[1] = a_p ;
	rc.setRow(2,index,element,false) ;
#	if CGL_DEBUG > 0
	if (debugger) assert(!debugger->invalidCut(rc)); 
	if (verbosity_ >= 4)
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
      if (verbosity_ >= 5) {
	std::cout
	  << "    Generating lower disaggregation cut for x<"
	  << t << "> (" << k << ")." << std::endl ;
	std::cout
	  << "    l<" << t << "> = " << l_t
	  << " increased to " << lPrime_t << "." << std::endl ;
	std::cout
	  << "    x<" << t << "> = " << x_t
	  << " increased to " << (lPrime_t-deltaProbe_p*a_p)
	  << "." << std::endl ;
	std::cout
	  << "    x<" << p << "> = " << x_p
	  << " changed to " << deltaCut_p+luPrime_p
	  << "; effectiveness = " << effectiveness
	  << "; required " << disaggEffectiveness << "." << std::endl ;
      }
#     endif
      assert(effectiveness+primalTolerance_ > 0) ;
      if (effectiveness > disaggEffectiveness) {
	OsiRowCut rc ;
	rc.setEffectiveness(effectiveness) ;
	rc.setLb(lPrime_t+a_p*luPrime_p) ;
	rc.setUb(DBL_MAX) ;
	index[0] = t ;
	element[0] = 1.0 ;
	index[1] = p ;
	element[1] = a_p ;
	rc.setRow(2,index,element,false) ;
#	if CGL_DEBUG > 1
	if (debugger) assert(!debugger->invalidCut(rc)); 
	if (verbosity_ >= 4)
	  std::cout << "    Adding lower disaggregation cut." << std::endl ;
#	endif
	rowCut.addCutIfNotDuplicate(rc) ;
      }
    }
  }

# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout << "    " << "end disaggCuts" ;
    if (rowCut.numberCuts() > cutsAtStart)
      std::cout << ", " << rowCut.numberCuts()-cutsAtStart << " cuts" ;
    std::cout << "." << std::endl ;
  }
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

void CglProbing::strengthenCoeff (
		      int p, unsigned int probeDir,
		      bool justReplace,
		      double needEffectiveness,
		      const OsiSolverInterface &si,
		      CglProbingRowCut &rowCut,
		      const CglPhic &phic,
		      const double *const colsol,
		      int numVarChgs,
		      const CglPhic::CglPhicBndPair *const varBndChgs,
		      int numLhsChgs,
		      const CglPhic::CglPhicBndPair *const lhsBndChgs,
		      const int *const realRows,
		      double *const element,
		      int *const index,
		      CglTreeInfo *const info) const
{
# if CGL_DEBUG > 0
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
# endif
# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout << "    "
      << "start strengthenCoeff, probing "
      << "x<" << p << "> "
      << ((probeDir == probeDown)?"down":"up")
      << "; " << numVarChgs << " vars, "
      << numLhsChgs << " cons." << std::endl ;
  }
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
  Unpack.
*/
  const CoinPackedMatrix *rowCopy = 0 ;
  const CoinPackedMatrix *colCopyUnused = 0 ;
  const double *rowUpper = 0 ;
  const double *rowLower = 0 ;
  phic.getSystem(rowCopy,colCopyUnused,rowLower,rowUpper) ;
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts() ;
  const int *rowLen = rowCopy->getVectorLengths() ;
  const int *colIndices = rowCopy->getIndices() ;
  const double *coeffs = rowCopy->getElements() ;

  double *colLower = 0 ;
  double *colUpper = 0 ;
  phic.getColBnds(colLower,colUpper) ;
/*
  Open up an outer loop to walk stackR and look for interesting constraints.
*/
  for (int k = 0 ; k < numLhsChgs ; k++) {
    const CglPhic::CglPhicBndPair &lhsChgRec = lhsBndChgs[k] ;
    const int &i = lhsChgRec.ndx_ ;
    const double &bi = rowUpper[i] ;
    const double &blowi = rowLower[i] ;
    const double &Li = lhsChgRec.lb_ ;
    const double &Ui = lhsChgRec.ub_ ;
/*
  We can't get anywhere unless probing has made one end or the other of the
  constraint redundant (but not too redundant --- the gap value will become
  the new coefficient a<ip>', so it must be reasonable).

  Range constraints pose the danger that we'll want to assign two different
  values to a<ip>. If we're in `justReplace' mode, this is a nonstarter.
*/
    const double uGap = bi-Ui ;
    const double lGap = blowi-Li ;
    bool useableUGap = ((uGap > primalTolerance_) && (uGap < bigCoeff)) ;
    bool useableLGap = ((-lGap > primalTolerance_) && (-lGap < bigCoeff)) ;
#   if CGL_DEBUG > 0
    if (verbosity_ >= 5)
      std::cout
        << "          " << "r(" << i << ") uGap " << uGap << " "
	<< ((useableUGap)?"useable":"unuseable") << "; lGap " << lGap << " "
	<< ((useableUGap)?"useable":"unuseable") << "." << std::endl ;
#   endif
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
    const CoinBigIndex kkstart = rowStart[i] ;
    const CoinBigIndex kkend = kkstart+rowLen[i] ;
    for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
      int k = colIndices[kk] ;
      double aik = coeffs[kk] ;
      if (k == p) {
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
      deltaMul = colUpper[p]+1.0 ; 
    } else {
      if (aip >= 0) {
        if (!useableLGap) continue ;
        delta = lGap ;
	revisebi = false ;
      } else {
        if (!useableUGap) continue ;
        delta = uGap ;
      }
      deltaMul = colLower[p]-1.0 ; 
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
    double xp = colsol[p] ;
    sum += aipPrime*xp ;
    double biPrime = DBL_MAX ;
    double blowiPrime = -DBL_MAX ;
    double violation = 0.0 ;
    bool genCut = false ;
    if (revisebi) {
      biPrime = bi+deltaMul*delta ;
      violation = sum-biPrime ;
    } else {
      blowiPrime = blowi+deltaMul*delta ;
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
    if (violation <= 0 && paranoia_ > 0) {
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
      int k = colIndices[kk] ;
      double aik = coeffs[kk] ;
      if (k != p) {
	index[n] = k ;
	element[n++] = aik ;
      }
    }
    if (aipPrimeNonZero) {
      index[n] = p ;
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
#   if CGL_DEBUG > 0
    if (verbosity_ >= 4) {
      std::cout << "    " << "Strengthen Cut: " << std::endl ;
      CglProbingDebug::dump_row(i,rc.lb(),rc.ub(),nan(""),nan(""),
	  &si,true,true,4,n,index,element,1.0e-10,colLower,colUpper) ;
      std::cout << "    " << "Original Row:" << std::endl ;
      CoinBigIndex jj = rowStart[i] ;
      CglProbingDebug::dump_row(i,blowi,bi,Li,Ui,&si,true,true,4,rowLen[i],
	  &colIndices[jj],&coeffs[jj],1.0e-10,colLower,colUpper) ;
    }
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
      int realRow = i ;
      if (realRows)
	realRow = realRows[realRow] ;
      assert(realRow >= 0) ;
      rc.setWhichRow(realRow) ;
      if (!justReplace) {
#       if CGL_DEBUG > 0
	if (verbosity_ >= 2) {
	  std::cout
	    << "    strengthen coeff cut on real row "
	    << realRow << " added to cut collection." << std::endl ;
	}
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
#         if CGL_DEBUG > 0
	  if (verbosity_ >= 2) {
	    std::cout
	      << "    strengthen coeff on real row " << realRow
	      << " will replace row, eff " << effectiveness << "."
	      << std::endl ;
	  }
	  newCuts++ ;
#         endif
	  info->strengthenRow[realRow] = rc.clone() ;
	} else {
#         if CGL_DEBUG > 0
	  if (verbosity_ >= 2) {
	    std::cout
	      << "    strengthen coeff cut on real row " << realRow
	      << " rejected; eff " << effectiveness << " >  eff "
	      << info->strengthenRow[realRow]->effectiveness()
	      << " of incumbent." << std::endl ;
	  }
#         endif
	}
      }
    }

# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    int cutsAtEnd = rowCut.numberCuts() ;
    int outplaceCuts = cutsAtEnd-cutsAtStart ;
    int inplaceCuts = newCuts-outplaceCuts ;
    std::cout << "    " << "end strengthenCoeff" ;
    if (newCuts > 0)
      std::cout << ", " << newCuts << " cuts, " << inplaceCuts << " in place" ;
    std::cout << "." << std::endl ;
  }
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
bool CglProbing::monotoneActions (
		      const OsiSolverInterface &si,
		      OsiCuts &cs,
		      const char *const intVar,
		      int bndLen,
		      const CglPhic::CglPhicBndPair *const varBndChgs,
		      const CglPhic::CglPhicBndPair *const origBnds,
		      const double *const colsol,
		      int *const index, double *const element) const
{
  OsiColCut cc ;
  int anyColCuts = 0 ;
  int nTot = 0 ;
  int nFix = 0 ;
# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    const CglPhic::CglPhicBndPair &newBnd_p = varBndChgs[0] ;
    const int &p = newBnd_p.ndx_ ;
    std::cout
      << "    " << "start monotoneActions; "
      << "x<" << p << "> monotone ["
      << newBnd_p.lb_ << ", " << newBnd_p.ub_ 
      << "] from " << colsol[p] << "." << std::endl ;
  }
# endif
/*
  Stash the improved column bounds, u<j> first, then l<j>. If any of them
  actually cut off the current solution, indicate that we have cuts.
*/
  bool ifCut = false ;
  for (int k = 0 ; k < bndLen ; k++) {
    const CglPhic::CglPhicBndPair &newBnd_j = varBndChgs[k] ;
    const CglPhic::CglPhicBndPair &origBnd_j = origBnds[k] ;
    const int &j = newBnd_j.ndx_ ;
    const double &ou_j = origBnd_j.ub_ ;
    const double &nu_j = newBnd_j.ub_ ;
    if (intVar[j]) {
      if (nu_j < ou_j-1.0e-4) {
	element[nFix] = nu_j ;
	index[nFix++] = j ;
	if (colsol[j] > nu_j+primalTolerance_) {
	  ifCut = true ;
	  anyColCuts++ ;
	} else {
	  ifCut = false ;
	}
#       if CGL_DEBUG > 0
	if (verbosity_ >= 4) {
	  std::cout
	    << "    " << " x<" << j << "> ub " << ou_j << " ==> " << nu_j ;
	  if (ifCut) std::cout << " (cut)" ;
	  std::cout << "." << std::endl ;
	}
#       endif
      }
    }
  }
  if (nFix) {
    nTot = nFix ;
    cc.setUbs(nFix,index,element) ;
    nFix = 0 ;
  }
  for (int k = 0 ; k < bndLen ; k++) {
    const CglPhic::CglPhicBndPair &newBnd_j = varBndChgs[k] ;
    const CglPhic::CglPhicBndPair &origBnd_j = origBnds[k] ;
    const int &j = newBnd_j.ndx_ ;
    const double &ol_j = origBnd_j.lb_ ;
    const double &nl_j = newBnd_j.lb_ ;
    if (intVar[j]) {
      if (nl_j > ol_j+1.0e-4) {
	element[nFix] = nl_j ;
	index[nFix++] = j ;
	if (colsol[j] < nl_j-primalTolerance_) {
	  ifCut = true ;
	  anyColCuts++ ;
	} else {
	  ifCut = false ;
	}
#       if CGL_DEBUG > 0
	if (verbosity_ >= 4) {
	  std::cout
	    << "    " << " x<" << j << "> lb " << ol_j << " ==> " << nl_j ;
	  if (ifCut) std::cout << " (cut)" ;
	  std::cout << "." << std::endl ;
	}
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
    if (paranoia_ > 0)
      assert(CglProbingDebug::checkBounds(si,cc,verbosity_) == 0) ;
#   endif
    cs.insert(cc) ;
  }

# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout << "    " << "end monotoneActions" ;
    if (nTot > 0)
      std::cout << ", " << nTot << " monotone" ;
    if (anyColCuts > 0)
      std::cout << ", " << anyColCuts << " column cuts" ;
    std::cout << "." << std::endl ;
  }
# endif
  return (anyColCuts > 0) ;
}



/*
  This method examines rows whose lhs bounds have changed, looking
  for redundant rows that consist solely of singleton variables (i.e.,
  variables that appear in just this constraint). It then looks to see if
  any of the variables in the row can be fixed at bound.

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

void CglProbing::singletonRows (int jProbe, bool useObj,
		    const OsiSolverInterface &si, CglPhic &phic, int numChgs,
		    const CglPhic::CglPhicBndPair *const lhsBndChgs) const
{
# if CGL_DEBUG > 0
  if (verbosity_ >= 3)
    std::cout << "    " << "start singletonRows." << std::endl ;
  int numSingletons = 0 ;
# endif
/*
  Unpack a few vectors from the propagator and the row-major matrix.
*/
  const CoinPackedMatrix *rowCopy ;
  const CoinPackedMatrix *colCopyUnused ;
  const double *rowLower ;
  const double *rowUpper ;
  phic.getSystem(rowCopy,colCopyUnused,rowLower,rowUpper) ;

  const int nRows = rowCopy->getNumRows() ;
  const CoinBigIndex *rowStarts = rowCopy->getVectorStarts() ;
  const int *rowLens = rowCopy->getVectorLengths() ;
  const int *colIndices = rowCopy->getIndices() ;
  const double *coeffs = rowCopy->getElements() ;

  double *colLower ;
  double *colUpper ;
  phic.getColBnds(colLower,colUpper) ;
/*
  `Singleton' must be based on the column length in the original system.
*/
  const double *c = si.getObjCoefficients() ;
  const int *columnLengths = si.getMatrixByCol()->getVectorLengths() ;
  const double objSense = si.getObjSense() ;
/*
  Open a loop to work through the rows in lhsBndChgs. Don't process the
  objective (if present, it's the last row in the constraint system).
*/
  for (int k = 0 ; k < numChgs ; k++) {
    const CglPhic::CglPhicBndPair &chgRec = lhsBndChgs[k] ;
    const int &i = chgRec.ndx_ ;
    const double &Li = chgRec.lb_ ;
    const double &Ui = chgRec.ub_ ;
    if (useObj && i == (nRows-1)) continue ;
/*
  Check the gaps. If the constraint is potentially tight in both directions,
  there's nothing more to do here.
*/
    const double uGap = rowUpper[i]-Ui ;
    const double lGap = Li-rowLower[i] ;
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
    const CoinBigIndex kkstart = rowStarts[i] ;
    const CoinBigIndex kkend = kkstart+rowLens[i] ;
    for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
      int j = colIndices[kk] ;
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
#   if CGL_DEBUG > 0
    if (verbosity_ >= 5) {
      std::cout << "    " << "singleton row: " << std::endl ;
      CglProbingDebug::dump_row(i,rowLower[i],rowUpper[i],
	  phic.getRowLhsLB(i),phic.getRowLhsUB(i),
	  &si,true,true,4,rowLens[i],&colIndices[rowStarts[i]],
	  &coeffs[rowStarts[i]],primalTolerance_,colLower,colUpper) ;
    }
#   endif
/*
  If we've passed the tests, we can look for variables suitable to drive to
  bound. Work against U<i> first. We're looking for variables with a<ij>
  > 0 that will be naturally driven to u<j>, or variables with a<ij> <
  0 that will be naturally driven to l<j>.

  Call chgColBnd to do the work. By definition, this change isn't going
  to propagate, but bookkeeping will be updated so that revert will work
  properly. Note that the change will be reflected back in colLower and
  colUpper because these are loaned to the propagator.
*/
    if (uGap > primalTolerance_) {
      for (CoinBigIndex kk = kkstart ; kk < kkend ; kk++) {
	int j = colIndices[kk] ;
	const double lj = colLower[j] ;
	const double uj = colUpper[j] ;
	if (uj > lj) {
	  double value = coeffs[kk] ;
	  if (objSense*c[j]*value < 0.0) {
	    bool feas = true ;
	    if (value > 0.0) {
	      phic.chgColBnd(j,'l',uj,feas) ;
	    } else {
	      phic.chgColBnd(j,'u',lj,feas) ;
	    }
#	    if CGL_DEBUG > 0
	    numSingletons++ ;
	    if (verbosity_ >= 4) {
	      std::cout
		<< "    row " << i << " uGap " << rowUpper[i] << "-"
		<< Ui << " = " << uGap
		<< " " << " x(" << j
		<< ") c = " << c[j] << ", a = " << value ;
	      if (colLower[j] == uj)
		std::cout << " l " << lj << " -> " << uj ;
	      else
		std::cout << " u " << uj << " -> " << lj ;
	      std::cout << std::endl ;
	    }
#	    endif
	    assert(feas) ;
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
	int j = colIndices[kk] ;
	const double lj = colLower[j] ;
	const double uj = colUpper[j] ;
	if (uj > lj) {
	  double value = coeffs[kk] ;
	  if (objSense*c[j]*value > 0.0)
	  { 
	    bool feas = true ;
	    if (value < 0.0) {
	      phic.chgColBnd(j,'l',uj,feas) ;
	    } else {
	      phic.chgColBnd(j,'u',lj,feas) ;
	    }
#	    if CGL_DEBUG > 0
	    numSingletons++ ;
	    if (verbosity_ >= 4) {
	      std::cout
		<< "    row " << i << " lGap " << Li << "-"
		<< rowLower[i] << " = " << lGap
		<< " " << " x(" << j
		<< ") c = " << c[j] << ", a = " << value ;
	      if (colLower[j] == uj)
		std::cout << " l " << lj << " -> " << uj ;
	      else
		std::cout << " u " << uj << " -> " << lj ;
	      std::cout << std::endl ;
	    }
#	    endif
	    assert(feas) ;
	  }
	}
      }
    }
  }
# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout << "    " << "end singletonRows" ;
    if (numSingletons)
      std::cout << ", found " << numSingletons << " singletons" ;
    std::cout << "." << std::endl ;
  }
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


  TODO: The goal is to completely eliminate stackR, minR, maxR, markR. They
  	should be hidden away in CglPhic. Similarly for rowStartPos,
	largestPositiveInRow, largestNegativeInRow.    -- lh, 110405 --

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
			CglPhic &phic,
		        const int *const realRows,
                        CglTreeInfo *const info,
		        bool useObj, bool useCutoff, double cutoff) const

{
# if CGL_DEBUG > 0
  int numberRowCutsBefore = cs.sizeRowCuts() ;
  int numberColCutsBefore = cs.sizeColCuts() ;
  if (verbosity_ >= 3) {
    std::cout << "  " << "start CglProbing::probe." << std::endl ;
  }
# endif
/*
  Parameter replacement. Unpack from the CglPhic object various things that
  used to come in as individual parameters.
*/
  double *colLower = 0 ;
  double *colUpper = 0 ;
  phic.getColBnds(colLower,colUpper) ;
  const CoinPackedMatrix *rowCopy = 0 ;
  const CoinPackedMatrix *colCopyUnused = 0 ;
  const double *rowLower = 0 ;
  const double *rowUpper = 0 ;
  phic.getSystem(rowCopy,colCopyUnused,rowLower,rowUpper) ;
  const char *intVar = 0 ;
  phic.getColType(intVar) ;

/*
  PREPARATION

  All the code down through `PASS LOOP: HEAD' is preparation. Do all of the
  setup that will not change over the nested loops that do the work of probing.
  Note that rowCopy may have considerably fewer rows than the original system.
*/

  int nRows = rowCopy->getNumRows() ;
  int nCols = rowCopy->getNumCols() ;
  int maxStack = (info->inTree)?maxStack_:maxStackRoot_ ;

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

/*
  TODO: Remove useless arrays largestPositiveInRow, largestNegativeInRow,
  	columnGap
*/
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
  Used for grabbing bound changes.
*/
  CoinPackedVector varlbChgs,varubChgs ;
  CoinPackedVector lhsLBChgs,lhsUBChgs ;
  int numVarBndChgs = 0 ;
  CglPhic::CglPhicBndPair *varBndChgs = 0 ;
  CglPhic::CglPhicBndPair *origBnds = 0 ;
  int numRowLhsBndChgs = 0 ;
  CglPhic::CglPhicBndPair *rowLhsBndChgs = 0 ;

  int numVarBndChgsDown = 0 ;
  CglPhic::CglPhicBndPair *varBndChgsDown = 0 ;
  CglPhic::CglPhicBndPair *origBndsDown = 0 ;

  bool recalcRowLhsBnds = false ;
  bool overallFeasible = true ;

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
/*
  Grab the column solution and reduced costs and groom them.
*/
  CoinMemcpyN(si.getReducedCost(),nCols,djs) ;
  CoinMemcpyN(si.getColSolution(),nCols,colsol) ;
  double direction = si.getObjSense() ;
# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    std::cout
      << "      " << "grooming solution and reduced costs." << std::endl ;
  }
# endif
  groomSoln(direction,primalTolerance_,nCols,djs,colLower,colsol,colUpper) ;
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

# if CGL_DEBUG > 0
  if (verbosity_ > 4) {
    std::cout
      << "    prior to probe loop, justReplace "
      << ((justReplace)?"true":"false") << ", required effectiveness "
      << needEffectiveness << "." << std::endl ;
  }
# endif

/*
  PASS LOOP: HEAD

  Major block #2: Main pass loop.
*/
  int rowCuts = rowCuts_ ;
  double disaggEffectiveness = 1.0e-3 ;
  int iPass = 0 ;
  int nfixed = -1 ;
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_ ;

  while (iPass < maxPass && nfixed) {
    iPass++ ;
    nfixed = 0 ;
#   if CGL_DEBUG > 0
    if (verbosity_ > 4) {
    std::cout
      << "  probe loop: starting pass " << iPass << "." << std::endl ;
    }
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

  TODO: Figure out something to replace leftTotalStack as a limit on effort.
  	-- lh, 110405 --
*/
      if (rowCut.outOfSpace()) {
	if (!justFix && (!nfixed || info->inTree)) {
#         if CGL_DEBUG > 0
	  if (verbosity_ > 2 && !info->inTree)
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
	    if (verbosity_ > 2)
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
      const double gap = uj-lj ;
#     if CGL_DEBUG > 0
      if (verbosity_ > 2) {
	std::cout
	  << "  probe: looking at " << " x<" << j << "> = " << xj
	  << ", l = " << lj << ", u = " << uj << "." << std::endl ;
	if ((xj > uj+tolerance) || (xj < lj-tolerance))
	  std::cout << "    discard; out-of-bound." << std::endl ;
	else if (gap < 1.0e-8)
	  std::cout << "    discard: fixed." << std::endl ;
      }
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
      if (verbosity_ > 0 && (down < lj || up > uj || down == up)) {
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

  If we need to recalculate the row lhs bounds, we can't put it off any
  longer.
*/
      unsigned int iWay ;
      unsigned int way[] = { probeDown, probeUp } ;
      unsigned int feasValue[] =
          { downIterFeas, upIterFeas } ;
      unsigned int feasRecord = 0 ;
      bool probeFeasible = true ;
      if (recalcRowLhsBnds) {
        phic.initLhsBnds() ;
	recalcRowLhsBnds = false ;
      }
/*
  PROBE LOOP: HEAD

  Open a loop to probe up and down for the current variable.

  As with the previous loops (variables in lookedAt_, passes) this loop
  extends to the end of major block #2).
*/
      for (iWay = downIter ; iWay < oneIterTooMany ; iWay++) {
/*
  Calculate movement given the current direction. There are various
  circumstances in which we need to know that the probe has tightened by
  one from the existing bound.
*/
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
          assert(movement < -0.99999) ;
	  if (movement < -1.5)
	    probeDistOne = false ;
        } else {
          movement = up-lj ;
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
	double solMovement ;
	if (verbosity_ >= 3) {
	  if (iWay == downIter) {
	    solMovement = down-xj ;
	    std::cout
	      << "  down probe, old u " << colUpper[j]
	      << ", new u " << colUpper[j]+movement ;
	  } else {
	    solMovement = xj-up ;
	    std::cout
	      << "  up probe, old l " << colLower[j]
	      << ", new l " << colLower[j]+movement ;
	  }
	  std::cout
	      << ", move " << movement << ", soln " << colsol[j]
	      << "(" << xj << "), move " << solMovement << "." << std::endl ;
	}
#       endif
/*
  Now that we know what we're doing, push the probe's bound change and
  propagate. If we're still feasible at the end of propagation, acquire
  the variable bound changes and adjust the objective.
*/
	if (iWay == downIter)
	  phic.chgColBnd(j,'u',down,probeFeasible) ;
	else
	  phic.chgColBnd(j,'l',up,probeFeasible) ;
        phic.propagate(probeFeasible) ;
/*
  Whether we're feasible or not, we don't need these any longer.
*/
	delete[] varBndChgs ;
	varBndChgs = 0 ;
	delete[] origBnds ;
	origBnds = 0 ;
	delete[] rowLhsBndChgs ;
	rowLhsBndChgs = 0 ;
/*
  If we're primal feasible, we need the variable bound change records. We may
  also need the row lhs bound changes (depending on the cuts we can generate).
  
  With the bound changes, calculate the change in objective. The calculation
  is conservative: we consider only column cuts (i.e., the new bound cuts
  off the existing solution) and adjust by that amount.
*/
	double objVal = 0.0 ;
	bool anyColCuts = false ;
        if (probeFeasible)
	{ objVal = si.getObjValue() ;
	  numVarBndChgs = phic.getColBndChgs(&varBndChgs,&origBnds) ;
	  anyColCuts = calcObjChg(objVal,primalTolerance_,
				  numVarBndChgs,varBndChgs,colsol,djs) ;
	  if (doSingletons || (doCoeffStrengthen && !anyColCuts))
	    numRowLhsBndChgs = phic.getRowLhsBndChgs(&rowLhsBndChgs,0) ;
#         if CGL_DEBUG > 0
	  if (verbosity_ > 2) {
	    std::cout << "    " << "probe feasible" ;
	    if (anyColCuts) std::cout << "; cuts off solution" ;
	    std::cout << "." << std::endl ;
	  }
#         endif
	}
/*
  PROBE LOOP: INFEASIBILITY

  Feasibility check. Primal feasibility is recorded in probeFeasible; we still
  need to test against the objective bound. If we have shown up and down
  infeasibility, we're infeasible, period. Break from the probe loop and
  terminate the lookedAt and pass loops by forcing their iteration counts past
  the maximum.
*/
	if (!probeFeasible || (useCutoff && (objVal > cutoff))) {
#	  if CGL_DEBUG > 0
	  if (verbosity_ > 2) {
	    std::cout
	      << "    " << ((way[iWay] == probeDown)?"down":"up")
	      << " probe for x<" << j << "> infeasible on " ;
	    if (!probeFeasible && (useCutoff && (objVal > cutoff)))
	      std::cout << "primal bounds and objective cutoff" ;
	    else if (!probeFeasible)
	      std::cout << "primal bounds" ;
	    else
	      std::cout
	        << "objective cutoff, z = " << objVal << ", cutoff "
		<< cutoff  ;
	    std::cout << "." << std::endl ;
	  }
#	  endif
	  probeFeasible = false ;
	  if (iWay == upIter && feasRecord == 0) {
	    overallFeasible = false ;
	    iLook = numberThisTime_ ;
	    iPass = maxPass ;
	    break ;
	  }
	} else {
	  feasRecord |= feasValue[iWay] ; 
	}
/*
  PROBE LOOP: DOWN PROBE

  End of the down probe iteration. If we're feasible, generate cuts and
  save the bounds changes for use after the up probe.  Then reset for the
  up probe and iterate.
*/
	if (iWay == downIter) {
          if (probeFeasible) {
	    if (doSingletons)
	      singletonRows(j,useObj,si,phic,numRowLhsBndChgs,rowLhsBndChgs) ;
	    if (doDisaggCuts)
	      disaggCuts(j,probeDown,disaggEffectiveness,
	      		 si,colsol,numVarBndChgs,varBndChgs,origBnds,
			 rowCut,index,element) ;
	    if (doCoeffStrengthen && !anyColCuts)
	      strengthenCoeff(j,probeDown,justReplace,
			      needEffectiveness,si,rowCut,phic,
			      colsol,numVarBndChgs,varBndChgs,
			      numRowLhsBndChgs,rowLhsBndChgs,
			      realRows,element,index,info) ;
	    numVarBndChgsDown = numVarBndChgs ;
	    delete[] varBndChgsDown ;
	    varBndChgsDown = varBndChgs ;
	    varBndChgs = 0 ;
	    delete[] origBndsDown ;
	    origBndsDown = origBnds ;
	    origBnds = 0 ;
	  }
	  phic.revert(true,true) ;
	  continue ;
	}
/*
  PROBE LOOP: MONOTONE

  We've done a down probe and an up probe. Outright infeasibility was
  handled above, as was iteration after the down probe. We can only reach
  here after the up probe, with at least one direction feasible. Check
  first for monotonicity.
  * If we are feasible up but not down, we're monotone up. Easy --- we simply
    record the resulting column bounds as column cuts and keep the current
    state.
  * If we were feasible down but not up, we're monotone down. This requires
    a bit more work, as we need to back out the up probe state and reinstall
    the down probe state. We only care about column bounds in monotoneActions,
    but we need to recalculate row lhs bounds if we have more variables to
    probe.
  Either way, we capture the results and move on to the next variable.
*/
        assert(iWay == upIter &&
	       ((feasRecord&(downIterFeas|upIterFeas)) != 0)) ;

        if (feasRecord != (downIterFeas|upIterFeas)) {
          nfixed++ ;
	  if (feasRecord == downIterFeas) {
	    phic.revert(true,false) ;
	    phic.editColBnds(numVarBndChgsDown,varBndChgsDown) ;
	    monotoneActions(si,cs,intVar,numVarBndChgsDown,varBndChgsDown,
			    origBndsDown,colsol,index,element) ;
	    recalcRowLhsBnds = true ;
	  } else {
	    monotoneActions(si,cs,intVar,numVarBndChgs,varBndChgs,
			    origBnds,colsol,index,element) ;
	  }
	  phic.clearPropagation() ;
	  break ;
        }
/*
  PROBE LOOP: DOWN AND UP FEASIBLE

  To reach here, this must be the up probe and both up and down probes are
  feasible. Generate row singletons and disaggregation cuts (just as for a
  feasible down probe). Next, revert to the bounds prior to the up probe.
  Finally, call implicationCuts to look for cuts based on both probes
  (binary probe variable only). If we manage to improve bounds
  (implicationCuts returns true) we'll need to recalculate row lhs bounds.
*/
	assert(iWay == upIter && (feasRecord == (downIterFeas|upIterFeas))) ;

	if (doSingletons)
	  singletonRows(j,useObj,si,phic,numRowLhsBndChgs,rowLhsBndChgs) ;
	if (doDisaggCuts)
	  disaggCuts(j,probeUp,disaggEffectiveness,
		     si,colsol,numVarBndChgs,varBndChgs,origBnds,
		     rowCut,index,element) ;
	if (doCoeffStrengthen && !anyColCuts)
	  strengthenCoeff(j,probeUp,justReplace,
			  needEffectiveness,si,rowCut,phic,
			  colsol,numVarBndChgs,varBndChgs,
			  numRowLhsBndChgs,rowLhsBndChgs,
			  realRows,element,index,info) ;
	phic.revert(true,true) ;
	if (probeIsDichotomy) {
	  int jProbe = j ;
	  bool foundColCuts = false ;
	  if (!doImplicationCuts) jProbe = -1 ;
	  foundColCuts = implicationCuts(jProbe,nCols,
	  				 cs,si,phic,intVar,colsol,
					 numVarBndChgsDown,varBndChgsDown,
					 numVarBndChgs,varBndChgs) ;
	  if (foundColCuts) recalcRowLhsBnds = true ;
        }
      }     // PROBE LOOP: END
    }    // LOOKEDAT LOOP: END
#   if CGL_DEBUG > 0
    if (verbosity_ > 3)
      std::cout
	<< "  probe: ending pass " << iPass
	<< ", fixed " << nfixed << " variables." << std::endl ;
#   endif
  }   // PASS LOOP: END
/*
  Free up the space we've been using.

  If this was allocated as one big block, it suffices to delete colsol.
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

  delete[] varBndChgs ;
  delete[] origBnds ;
  delete[] rowLhsBndChgs ;
  delete[] varBndChgsDown ;
  delete[] origBndsDown ;

/*
  jjf:  Add in row cuts

  Transfer row cuts from the local container into something that can escape
  this method. If we're not in `just replace' mode, simply copy them over to
  the container (cs) passed as a parameter using addCuts. If we're in `just
  replace' mode, they will be carried back to the caller in info.strengthenRow.
*/
  if (overallFeasible) {
    if (!justReplace) {
      rowCut.addCuts(cs,0,info->pass) ;
    } else {
      assert(info->strengthenRow != 0) ;
      for (int i = 0 ; i < nRows ; i++) {
	int realRow = (realRows)?realRows[i]:i ;
	assert(realRow >= 0 && realRow < si.getNumRows()) ;
	OsiRowCut *cut = info->strengthenRow[realRow] ;
	if (cut) {
#         if CGL_DEBUG > 0
	  if (verbosity_ > 2) {
	    std::cout
	      << "    row " << realRow << " (local " << i
	      << ") strengthened, effectiveness " << cut->effectiveness()
	      << "." << std::endl ;
	  }
#         endif
	}
      }
    }
  }
# if CGL_DEBUG > 0
  if (verbosity_ >= 3) {
    int numberRowCutsAfter = cs.sizeRowCuts() ;
    int numberColCutsAfter = cs.sizeColCuts() ;
    int inPlaceCuts = 0 ;
    if (justReplace) {
      for (int i = 0 ; i < nRows ; i++)
	if (info->strengthenRow[i] != 0) inPlaceCuts++ ;
    }
    std::cout << "  "
      << "end probe, "
      << ((overallFeasible)?"feasible":"infeasible")
      << ", " << (numberRowCutsAfter-numberRowCutsBefore) << " row cuts, "
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
  }
# endif

  return (overallFeasible) ;
}
