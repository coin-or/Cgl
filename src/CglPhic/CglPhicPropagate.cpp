/*
  Copyright (C) 2011, Lou Hafer
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file contains the methods that deal with bound propagation proper.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CglConfig.h"

#include "CoinFinite.hpp"
#include "CoinHelperFunctions.hpp"
#include "CglPhic.hpp"

namespace {

/*
  A pair of rounding methods that will do 'safe' rounding to the specified
  tolerance. For example, say we're rounding a value that's used as an upper
  bound. round_ub(1.15,.1) = 1.2.
*/
inline double roundU (double val, double tol) {
  double result = val-fmod(val,tol) ;
  if (val > 0.0) result += tol ;
  return (result) ;
}

inline double roundL (double val, double tol) {
  double result = val-fmod(val,tol) ;
  if (val < 0.0) result -= tol ;
  return (result) ;
}

}    // end file-local namespace

/*
  Add a constraint to the heap of constraints waiting to be processed. 
*/
void CglPhic::addToPending (int i, double ddelta, double dmetric)
{
  const float delta = static_cast<float>(ddelta) ;
  const float metric = static_cast<float>(dmetric) ;
/*
  Grab the metric for the candidate at the top of the heap. We'll use this to
  decide if we need to rebuild.
*/
  int bestPending = -1 ;
  float bestMetric = static_cast<float>(-infty_) ;
  if (!pending_.empty()) {
    bestPending = *pending_.begin() ;
    CglPhicCand &bestCand = candInfo_[bestPending] ;
    bestMetric = bestCand.metric_ ;
  }
/*
  Is this constraint already in the heap? If so, simply add to the existing
  entry and indicate we need to rebuild the heap.
*/
  CglPhicCand &candi = candInfo_[i] ;
  bool alreadyPendingi = candi.isPending_ ;

  if (alreadyPendingi) {
    candi.delta_ += delta ;
    if (candi.metric_ < infty_) candi.accum_ += metric ;
    if (candi.accum_ > disturbTol_) {
      candi.metric_ += candi.accum_ ;
      candi.accum_ = 0 ;
      if (candi.metric_ > bestMetric) rebuildHeap_ = true ;
    }
  } else if (i != inProcess_) {
/*
  Not in the heap, and not currently in process. Initialise the entry in
  candInfo_, add a pointer to the end of pending_ and set rebuildHeap_.
*/
    candi.delta_ = delta ;
    candi.accum_ = metric ;
    if (candi.accum_ > disturbTol_) {
      candi.metric_ = candi.accum_ ;
      candi.accum_ = 0 ;
      candi.isPending_ = true ;
      pending_.push_back(i) ;
      if (candi.metric_ > bestMetric) rebuildHeap_ = true ;
    }
  }

  if (verbosity_ >= 4 && candi.metric_ > disturbTol_) {
    std::cout << "          " ;
    if (alreadyPendingi)
      std::cout << "updated" ;
    else
      std::cout << "queued" ;
    std:: cout
      << " r(" << i << "), delta " << candi.delta_ << " (" << delta
      << "), metric " << candi.metric_ << " (" << metric << "); pending set "
      << pending_.size() << " entries." << std::endl ;
    if (rebuildHeap_ == true) {
      std::cout << "          " ;
      std::cout
        << "Requesting heap rebuild; best " << bestMetric << "." << std::endl ;
    }
  }
}

/*
  Remove the constraint with the largest metric from the pile of constraints
  waiting to be processed.

  This initial implementation is straightforward: If rebuildHeap_ is true,
  rebuild before popping the first entry.
*/
int CglPhic::getFromPending ()
{
  typedef CglPhicCandHeap::iterator HeapIter ;
  HeapIter heapStart = pending_.begin() ;
  HeapIter heapEnd = pending_.end() ;
/*
  Do we need to rebuild the heap?
*/
  if (rebuildHeap_) {
    if (verbosity_ >= 4) {
      std::cout << std::endl ;
      std::cout
        << "      rebuilding heap, " << pending_.size() << " entries ... " ;
    }
    std::make_heap(heapStart,heapEnd,heapCmpObj_) ;
    rebuildHeap_ = false ;
  }
/*
  Pop the best element: read the value of the entry at the top of the heap,
  then remove the element from the heap and rebuild the heap (pop_heap), and
  finally remove the element from pending_ (pop_back).
*/
  int bestPending = -1 ;
  if (!pending_.empty()) {
    bestPending = *heapStart ;
    std::pop_heap(heapStart,heapEnd,heapCmpObj_) ;
    pending_.pop_back() ;
  }

  if (verbosity_ >= 5) {
    std::cout << std::endl ;
    std::cout << "          " ;
    if (bestPending >= 0) {
      CglPhicCand &candi = candInfo_[bestPending] ;
      std::cout
	<< "selected r(" << bestPending << "), delta "
	<< candi.delta_ << ", metric " << candi.metric_ ;
    } else {
      std::cout << "pending set empty." ;
    }
    std::cout << std::endl ;
  }
  return (bestPending) ;
}


/*
  Process a change in a bound for a variable. We need to walk the column,
  updating L<i> or U<i> as necessary. For a <= constraint (finite b<i>),
  only L<i> matters. For a >= constraint (finite blow<i>), only U<i> matters.
  Range and equality constraints need both.

  For variable x<j>, bnd should be 'l' or 'u', and nbndj the new bound value.
  Loss of feasibility is reported back through feas.
*/
void CglPhic::chgColBnd (int j, char bnd, double nbndj, bool &feas)
{
  assert(bnd == 'l' || bnd == 'u') ;
  assert(0 <= j && j < n_) ;

  feas = true ;

  // Create separate deltal and deltau for readability.
  const bool deltal = (bnd == 'l') ;
  const bool deltau = !deltal ;
/*
  In general, the change L(i) or U(i) is to adjust bnd_ (the finite portion
  of the row bound) by subtracting the contribution of the old variable
  bound and adding the contribution of the new variable bound.

  When we make the transition from a single infinite contribution to
  completely finite, the finite portion of the row bound is already correct
  except for the component from the variable contributing the infinity
  (because we recalculated it when infCnt_ dropped from 2 to 1). So the
  correct delta is simply the new finite bound value (and may well be zero).

  No need to go further if the bounds are equal.
*/
  const double obndj = (deltal)?colL_[j]:colU_[j] ;
  const bool toFinite = (deltal)?(obndj <= -infty_):(obndj >= infty_) ;
  double delta = 0.0 ;
  if (toFinite) {
    delta = nbndj ;
  } else {
    delta = nbndj-obndj ;
    if (CoinAbs(delta) < zeroTol_) delta = 0.0 ;
  }
  assert((deltal && (delta > 0)) || (deltau && (delta < 0)) || toFinite) ;
/*
  Create/update a change record and change the bound in colLower or colUpper.
*/
  recordVarBndChg(j,bnd,nbndj) ;
/*
  See if we're propagating this variable type.
*/
  if (((1<<intVar_[j])&propType_) == 0) {
    if (verbosity_ >= 5) {
      char vartypelet[3] = {'c','b','g'} ;
      std::cout
        << "          suppressing propagation of x(" << j << ") (" 
	<< vartypelet[intVar_[j]] << "); propagating" ;
	if ((CglPhic::PropCon&propType_) != 0) std::cout << " (c)" ;
	if ((CglPhic::PropGenInt&propType_) != 0) std::cout << " (g)" ;
	if ((CglPhic::PropBinary&propType_) != 0) std::cout << " (b)" ;
	std::cout << "." << std::endl ;
    }
    return ;
  }
    
/*
  Set up a loop to walk the column
*/
  const CoinBigIndex &colStart = cm_colStarts_[j] ;
  const CoinBigIndex colEnd = colStart+cm_colLens_[j] ;
  for (CoinBigIndex ii = colStart ; ii < colEnd ; ii++) {
    int i = cm_rowIndices_[ii] ;
    if (i == inProcess_) continue ;
    double aij = cm_coeffs_[ii] ;
    if (CoinAbs(aij) < zeroTol_) continue ;
    const double &blowi = rhsL_[i] ;
    const double &bi = rhsU_[i] ;
    CglPhicLhsBnd newLi = lhsL_[i] ;
    CglPhicLhsBnd newUi = lhsU_[i] ;
/*
  Which lhs bound will change?
  * U<i> changes if (a<ij> > 0 && delta(u<j>)) || (a<ij> < 0 && delta(l<j>))
  * L<i> changes if (a<ij> > 0 && delta(l<j>)) || (a<ij> < 0 && delta(u<j>))
  Again, create separate changeLi and changeUi for readability.
*/
    const bool changeLi = (((aij > 0) && deltal) || ((aij < 0) && deltau)) ;
    const bool changeUi = !changeLi ;
/*
  And what sort of change are we talking about? The possibilities are:

  * A reduction in the number of infinite contributions, where the end
    result is still two or more.
  * A reduction in the number of infinite contributions from two to one
    (full recalculation required to identify the single remaining infinite
    contribution).
  * A reduction in the number of infinite contributions from one to zero
    (incremental update, because the change to a single infinity triggered a
    full bound recalculation). Note that infCnt_ contains -(ndx+1) for the
    offending variable.
  * A simple change in a finite bound (incremental update unless full
    recalculation is required because we've reached the revision limit).

  A >= constraint (finite blow<i>) uses only U<i> in propagation, while a <=
  constraint (finite b<i>) uses only L<i>. Range and equality use both. But
  there are other uses where the bound not used in propagation is still
  important. Both are calculated here, but the constraint will not be queued
  for propagation unless the relevant bound is updated.

  Note that at most one bound revision block will execute.
*/
    bool fullRecalc = false ;
    double deltaLi = 0.0 ;
    double deltaUi = 0.0 ;
    bool updatedLi = false ;
    bool updatedUi = false ;
    if (changeUi) {
      if (!toFinite) {
	if (newUi.infCnt_ == 0) {
	  if (newUi.revs_ < revLimit_) {
	    deltaUi = aij*delta ;
	    if (CoinAbs(deltaUi) > zeroTol_) {
	      newUi.bnd_ += deltaUi ;
	      newUi.revs_++ ;
	      updatedUi = true ;
	    }
	  } else {
	    fullRecalc = true ;
	    updatedUi = true ;
	  }
	}
      } else {
        assert(newUi.infCnt_ != 0) ;
        if (newUi.infCnt_ > 2) {
	  newUi.infCnt_-- ;
	} else if (newUi.infCnt_ == 2) {
	  fullRecalc = true ;
	} else {
	  if (paranoia_ >= 3) {
	    if (!(0 <= (-newUi.infCnt_)-1) && ((-newUi.infCnt_)-1 < m_)) {
	      std::cout
	        << "INV FAIL: (U) single infinity fails 0 <= "
		<< (-newUi.infCnt_)-1 << " < " << m_ << std::endl ;
	      assert(false) ;
	    }
	  }
	  deltaUi = aij*delta ;
	  if (CoinAbs(deltaUi) < zeroTol_) deltaUi = 0 ;
	  newUi.infCnt_ = 0 ;
	  newUi.bnd_ += deltaUi ;
	  newUi.revs_++ ;
	}
	updatedUi = true ;
      }
    }
    if (changeLi) {
      if (!toFinite) {
	if (newLi.infCnt_ == 0) {
	  if (newLi.revs_ < revLimit_) {
	    deltaLi = aij*delta ;
	    if (CoinAbs(deltaLi) > zeroTol_) {
	      newLi.bnd_ += deltaLi ;
	      newLi.revs_++ ;
	      updatedLi = true ;
	    }
	  } else {
	    fullRecalc = true ;
	    updatedLi = true ;
	  }
	}
      } else {
        assert(newLi.infCnt_ != 0) ;
        if (newLi.infCnt_ > 2) {
	  newLi.infCnt_-- ;
	} else if (newLi.infCnt_ == 2) {
	  fullRecalc = true ;
	} else {
	  if (paranoia_ >= 3) {
	    if (!(0 <= (-newLi.infCnt_)-1) && ((-newLi.infCnt_)-1 < m_)) {
	      std::cout
	        << "INV FAIL: (L) single infinity fails 0 <= "
		<< (-newLi.infCnt_)-1 << " < " << m_ << std::endl ;
	      assert(false) ;
	    }
	  }
	  deltaLi = aij*delta ;
	  if (CoinAbs(deltaLi) < zeroTol_) deltaLi = 0 ;
	  newLi.infCnt_ = 0 ;
	  newLi.bnd_ += deltaLi ;
	  newLi.revs_++ ;
	}
	updatedLi = true ;
      }
    }
    if (!(updatedLi || updatedUi)) continue ;


    if (paranoia_ >= 3) {
      if (updatedLi && updatedUi) {
	std::cout
	  << "INVARIANT FAIL: updating both Li and Ui!" << std::endl ;
	assert(false) ;
      }
      if (updatedLi && !(toFinite || (deltaLi >= 0.0))) {
	std::cout
	  << "INVARIANT FAIL: update Li, aij " << aij << ", delta " << delta
	  << ", deltaLi " << deltaLi << std::endl ;
	assert(false) ;
      }
      if (updatedUi && !(toFinite || (deltaUi <= 0.0))) {
	std::cout
	  << "INVARIANT FAIL: update Ui, aij " << aij << ", delta " << delta
	  << ", deltaUi " << deltaUi << std::endl ;
	assert(false) ;
      }
    }
/*
  Make a record of the change and change the bound in lhsL or lhsU. If we're
  actually infeasible, we do a bit of extra work by doing this prior to the
  feasibility check, but it makes the whole update process much easier if we
  can bury a full recalculation within the call to recordLhsBndChg.
*/
    if (updatedLi) {
      recordLhsBndChg(i,fullRecalc,'L',newLi) ;
    } else {
      recordLhsBndChg(i,fullRecalc,'U',newUi) ;
    }
/*
  Are we still feasible? Apparent loss of feasibility requires a more careful
  check. Recalculate the row bound (unless we just did that), then do a
  toleranced check.
*/
    const CglPhicLhsBnd &Li = lhsL_[i] ;
    const CglPhicLhsBnd &Ui = lhsU_[i] ;
    if ((Ui.infCnt_ == 0 && Ui.bnd_ < blowi) ||
        (Li.infCnt_ == 0 && Li.bnd_ > bi)) {
      double infLow = 0.0 ;
      double infUp  = 0.0 ;
      if (Ui.infCnt_ == 0 && Ui.bnd_ < blowi)
        infLow = blowi-Ui.bnd_ ;
      if (Li.infCnt_ == 0 && Li.bnd_ > bi)
        infUp = Li.bnd_-bi ;
      if (verbosity_ >= 3) {
        std::cout
	  << "            " << "apparent loss of feasibility, r("
	  << i << ")" ;
	if (infLow != 0.0)
	  std::cout
	    << "; U " << Ui.bnd_ << " < blow " << blowi << " by " << infLow ;
	if (infUp != 0.0)
	  std::cout
	    << "; L " << Li.bnd_ << " > b " << bi << " by " << infUp ;
	std::cout << "; rechecking ..." ;
      }
      if (!fullRecalc) calcLhsBnds(i) ;
      if (paranoia_ >= 3) {
	if (Ui.infCnt_ != 0 || Li.infCnt_ != 0) {
	  std::cout
	    << std::endl
	    << "ERROR: infinity count changed for row bound!" << std::endl ;
	  assert (false) ;
	}
      }
      if (Li.bnd_-bi > feasTol_) {
        feas = false ;
	if (verbosity_ >= 3)
	  std::cout
	    << std::endl
	    << "      " << "r(" << i << ") infeasible; " << Li
	    << " > b<" << i << "> = " << bi << ", infeas "
	    << Li.bnd_-bi << ", tol " << feasTol_ << "." << std::endl ;
      } else if (blowi-Ui.bnd_ > feasTol_) {
        feas = false ;
	if (verbosity_ >= 3)
	  std::cout
	    << std::endl
	    << "      " << "r(" << i << ") infeasible; " << Ui
	    << " < blow<" << i << "> = " << blowi << ", infeas "
	    << blowi-Ui.bnd_ << ", tol " << feasTol_ << "." << std::endl ;
      } else {
        if (verbosity_ >= 3)
	  std::cout << " false alarm." << std::endl ;
	if (paranoia_ >= 3 && (infLow > feasTol_ || infUp > feasTol_)) {
	  std::cout
	    << "ERROR: unacceptable error in lhs bounds." << std::endl ;
	  assert(false) ;
	}
      }
      if (!feas) return ;
    }
/*
  Is this change worth propagating? It has to be useful and large enough to
  matter.
*/
    const double &posGapi = info_[i].posGap_ ;
    const double &negGapi = info_[i].negGap_ ;
    double metric = 0.0 ;
    if ((blowi > -infty_) && updatedUi) {
      if (Ui.infCnt_ == 0) {
	if  (toFinite) {
	  addToPending(i,infty_,infty_) ;
	} else {
	  double gapLE = blowi-Ui.bnd_ ;
	  if (posGapi > gapLE || -negGapi < gapLE) {
	    metric = -deltaUi/info_[i].l1norm_ ;
	    if (metric > rowPropTol_) addToPending(i,deltaUi,metric) ;
	  }
	}
      }
    } else {
      if ((bi < infty_) && Li.infCnt_ == 0) {
	if  (toFinite) {
	  addToPending(i,infty_,infty_) ;
	} else {
	  double gapGE = bi-Li.bnd_ ;
	  if (-posGapi > gapGE || negGapi < gapGE) {
	    metric = deltaLi/info_[i].l1norm_ ;
	    if (metric > rowPropTol_) addToPending(i,deltaLi,metric) ;
	  }
	}
      }
    }
  }
}


/*
  Helper function to tighten an upper bound. The calculation is one of
    x<t> <= ( b<i> - (L<i>-a<it>l<t>))/a<it>	finite b<i>, a<it> > 0
    x<t> <= (blow<i>-(U<i>-a<it>l<t>))/a<it>	finite blow<i>, a<it> < 0

  In use, gap should be one of b<i>-L<i> or blow<i>-U<i> and the actual
  calculation will be x<t> <= gap/a<it>+l<t>.

  When the constraint bound (L<i>, U<i>) contains a single infinity
  contributed by x<t>, there's no need for the correction term a<it>l<t>
  (Li.bnd_ and Ui.bnd_ are calculated without the contribution of l<t>).
  Just pass a zero for the correction term (corr).

  Returns true if the variable remains feasible, false otherwise.
*/
bool CglPhic::tightenVu (int t, double gap, double ait, double corr)
{
  bool feas ;
/*
  Calculate a clean value.
*/
  double n_ut = gap/ait+corr ;
  if (CoinAbs(n_ut) < zeroTol_) n_ut = 0 ;
  n_ut = roundU(n_ut,10*zeroTol_) ;
/*
  If this is an integer, we may be able to tighten further.
*/
  if (intVar_[t]) n_ut = floor(n_ut+feasTol_) ;
/*
  Did we make any progress? If so, check to see whether we're still feasible.
  If we're feasible, propagate the change.
*/
  const double &o_ut = colU_[t] ;
  const double &o_lt = colL_[t] ;
  feas = true ;
  const double diff = o_ut-n_ut ;
  if (diff > zeroTol_) {
    if ((o_lt-n_ut) > feasTol_) {
      feas = false ;
      if (verbosity_ >= 1)
	std::cout
	  << "      " << "x(" << t << ") infeasible; "
	  << n_ut << " = u(" << t << ") < l(" << t << ") = " << o_lt
	  << ", infeas = " << o_lt-n_ut << ", tol = " << feasTol_
	  << "." << std::endl ;
    } else if (diff > colPropTol_) {
      chgColBnd(t,'u',n_ut,feas) ;
    }
  }
  return (feas) ;
}
/*
  Helper function to tighten a lower bound. The calculation is one of
    x<t> >= ( b<i> - (L<i>-a<it>u<t>))/a<it>	finite b<i>, a<it> < 0
    x<t> >= (blow<i>-(U<i>-a<it>u<t>))/a<it>	finite blow<i>, a<it> > 0

  In use, gap should be one of b<i>-L<i> or blow<i>-U<i> and the actual
  calculation will be x<t> <= gap/a<it>+u<t>.

  When the constraint bound (L<i>, U<i>) contains a single infinity
  contributed by x<t>, there's no need for the correction term a<it>u<t>
  (Li.bnd_ and Ui.bnd_ are calculated without the contribution of u<t>).
  Just pass a zero for the correction term (corr).

  Returns true if the variable remains feasible, false otherwise.
*/
bool CglPhic::tightenVl (int t, double gap, double ait, double corr)
{
  bool feas ;
/*
  Calculate a clean value.
*/
  double n_lt = gap/ait+corr ;
  if (CoinAbs(n_lt) < zeroTol_) n_lt = 0 ;
  n_lt = roundL(n_lt,10*zeroTol_) ;
/*
  If this is an integer, we may be able to tighten further.
*/
  if (intVar_[t]) n_lt = ceil(n_lt-feasTol_) ;
/*
  Did we make any progress? If so, check to see whether we're still feasible.
  If we're feasible, propagate the change.
*/
  const double &o_ut = colU_[t] ;
  const double &o_lt = colL_[t] ;
  feas = true ;
  const double diff = n_lt-o_lt ;
  if (diff > zeroTol_) {
    if ((n_lt-o_ut) > feasTol_) {
      feas = false ;
      if (verbosity_ >= 1)
	std::cout
	  << "      " << "x(" << t << ") infeasible; "
	  << n_lt << " = l(" << t << ") > u(" << t << ") = " << o_ut
	  << ", infeas = " << n_lt-o_ut << ", tol = " << feasTol_
	  << "." << std::endl ;
    } else if (diff > colPropTol_) {
      chgColBnd(t,'l',n_lt,feas) ;
    }
  }
  return (feas) ;
}



/*
  This method attempts to propagate a change in row lhs bound(s) to tighten
  bounds on variables. This is really an exercise in case analysis, and the
  trick here is to try to leave as little decision-making as possible in the
  inner loop that walks the constraint. Fortunately, the inner loop is fairly
  compact so we can replicate it.

  The case analysis for what bound we can tighten looks like this:
  * finiteUGap && a<it> > 0  =>  u<t>
  * finiteUGap && a<it> < 0  =>  l<t>
  * finiteLGap && a<it> > 0  =>  l<t>
  * finiteLGap && a<it> < 0  =>  u<t>
  Once we've nailed down the case, call the appropriate tighten method, with
  appropriate parameters. If the bound is tightened, chgColBnd will be
  invoked to propagate the change.

  There's one additional complication: It may be that L<i> or U<i> contains a
  single infinity, contributed by x<p>. Case analysis says that even in this
  case we might be able to improve the finite bound of x<p>, but the loop is
  sufficiently specialised that it's better to hold it separate.

  The end result is that we have four replications of the inner loop, and the
  only decision we have to make in the loop is based on a<it>.
*/
void CglPhic::processRow (int i, bool &feas)
{
  if (verbosity_ >= 5) {
    std::cout << "          " << "processing r(" << i << ")" << std::endl ;
  }
  feas = true ;

  const CglPhicLhsBnd &Li = lhsL_[i] ;
  const CglPhicLhsBnd &Ui = lhsU_[i] ;

  assert(Li.infCnt_ <= 0 || Ui.infCnt_ <= 0) ;

  const double &blowi = rhsL_[i] ;
  const double &bi = rhsU_[i] ;
  const bool finiteUGap = ((Li.infCnt_ == 0) && (bi < infty_))?true:false ;
  const double upperGap = bi-Li.bnd_ ;
  const bool finiteLGap = ((Ui.infCnt_ == 0) && (blowi > -infty_))?true:false ;
  const double lowerGap = blowi-Ui.bnd_ ;
/*
  We now have some idea of what we can do, but the final decision can't be
  made until we see each coefficient a<ij> as we walk the row.
*/
  const CoinBigIndex &rowStart = rm_rowStarts_[i] ;
  const CoinBigIndex rowEnd = rowStart+rm_rowLens_[i] ;
/*
  If the upper gap is finite, life is easy. But if Li contains a contribution
  of infinity from a single variable x<p>, L<i\t> will be finite only for t =
  p. We can try to improve the other (finite) bound on x<p>, but that's it.
  It's the sign of a<ip> that tells us whether l<p> or u<p> is infinite.
*/
  if (finiteUGap) {
    for (CoinBigIndex tt = rowStart ; tt < rowEnd && feas ; tt++) {
      const int &t = rm_colIndices_[tt] ;
      const double &ait = rm_coeffs_[tt] ;
      if (CoinAbs(ait) < zeroTol_) continue ;
      if (ait > 0) {
	feas = tightenVu(t,upperGap,ait,colL_[t]) ;
      } else {
	feas = tightenVl(t,upperGap,ait,colU_[t]) ;
      }
    }
  } else if (bi < infty_) {
    int p = (-Li.infCnt_)-1 ;
    for (CoinBigIndex tt = rowStart ; tt < rowEnd && feas ; tt++) {
      const int &t = rm_colIndices_[tt] ;
      if (t == p) {
	const double &aip = rm_coeffs_[tt] ;
	if (aip > zeroTol_)
	  feas = tightenVl(t,upperGap,aip,0) ;
	else if (aip < zeroTol_)
	  feas = tightenVu(t,upperGap,aip,0) ;
	break ;
      }
    }
  }
  if (!feas) return ;
/*
  If we're still feasible, repeat, using the lower gap.
*/
  if (finiteLGap) {
    for (CoinBigIndex tt = rowStart ; tt < rowEnd && feas ; tt++) {
      const int &t = rm_colIndices_[tt] ;
      const double &ait = rm_coeffs_[tt] ;
      if (CoinAbs(ait) < zeroTol_) continue ;
      if (ait > 0) {
	feas = tightenVl(t,lowerGap,ait,colU_[t]) ;
      } else {
	feas = tightenVu(t,lowerGap,ait,colL_[t]) ;
      }
    }
  } else if (blowi > -infty_) {
    int p = (-Ui.infCnt_)-1 ;
    for (CoinBigIndex tt = rowStart ; tt < rowEnd && feas ; tt++) {
      const int &t = rm_colIndices_[tt] ;
      if (t == p) {
	const double &aip = rm_coeffs_[tt] ;
	if (aip > zeroTol_)
	  feas = tightenVu(t,lowerGap,aip,0) ;
	else
	  feas = tightenVl(t,lowerGap,aip,0) ;
	break ;
      }
    }
  }
}



/*
  This method propagates bound changes over the set of pending constraints.
  The pending set must be loaded with an initial set of constraints by one
  or more direct calls to addToPending; one or more calls to chgColBnd;
  or a call to tightenAbInitio (which will also make the call to propagate).

  The main loop removes a constraint from the pending set and processes
  it to see if it's possible to tighten bounds on any of the variables
  involved in the constraint. When a bound is tightened, chgColBnd is
  called to walk the column and tighten bounds on constraints entangled
  with the variable.  This, in turn, will add constraints to the pending
  set for further processing.
*/
void CglPhic::propagate (bool &feas)
{
  if (verbosity_ >= 3) {
    if (verbosity_ >= 4) std::cout << "   " ;
    std::cout << " propagating ... " ;
  }
/*
  Open a loop to remove a constraint from the pending set and process it for
  changes to variable bounds.
*/
  feas = true ;
  while (!pending_.empty() && feas) {
    int i = getFromPending() ;
    inProcess_ = i ;
    candInfo_[i].isPending_ = false ;
    processRow(i,feas) ;
  }
  return ;
}


/*
  Set up for ab initio bound tightening. The constraint bound arrays are
  scanned and all constraints where bound propagation might be productive are
  placed on the propagation stack. Constraints with one infinite contribution
  to the row bound are given priority, since processing these constraints will
  result in a variable bound going from infinite to finite. After that, the
  metric is calculated using looseness and gap. See the typeset documentation
  for details.
*/
void CglPhic::tightenAbInitio (bool &feas)
{
  if (verbosity_ >= 3) {
    std::cout << "    " << "Ab initio bound tightening ..." ;
    if (verbosity_ >= 4) std::cout << std::endl ;
  }
/*
  Open a loop to step through the rows.
*/
  for (int i = 0 ; i < m_ ; i++) {
    const double &normi = info_[i].l1norm_ ;
    const double &posGapi = info_[i].posGap_ ;
    const double &negGapi = info_[i].negGap_ ;
    const double &bi = rhsU_[i] ;
    const double &blowi = rhsL_[i] ;
    double gapLE = 0.0 ;
    double gapGE = 0.0 ;
    double metricLE = 0.0 ;
    double metricGE = 0.0 ;
    bool propagate = false ;
    if (verbosity_ >= 5) {
      std::cout
	<< "          " << "considering " << blowi << " < "
	<< lhsL_[i] << " <= r(" << i << ") <= " << lhsU_[i] << " < "
	<< bi << std::endl ;
    }
/*
  The top-level decision is based whether or not the rhs upper bound b<i>
  and lhs lower bound L(i) are finite. After that, we have two cases:
    * U(i) is not finite: In this case we should be able to change infinite
      bounds to finite bounds. Queue this constraint with high priority.
    * U(i) is finite: In this case we can queue based on the gap, which will
      be finite. If (bi-Li)-posGap < 0 we will be able to tighten the upper
      bound on some variable with a positive coefficient. If (bi-Li)+negGap
      > 0, we will be able to tighten the lower bound on some variable with
      a negative coefficient.
*/
    if (bi < infty_ && lhsL_[i].infCnt_ == 0) {
      if (lhsU_[i].infCnt_ != 0) {
        gapLE = infty_ ;
	metricLE = infty_ ;
	propagate = true ;
      } else {
	gapLE = bi-lhsL_[i].bnd_ ;
	if (posGapi > gapLE || -negGapi < gapLE) {
	  propagate = true ;
	  metricLE = CoinMax((posGapi-gapLE),(gapLE+negGapi)) ;
	  metricLE /= normi ;
	}
      }
      if (verbosity_ >= 4 && propagate) {
	std::cout << "        queuing r(" << i << ") on ugap, " ;
	if (lhsU_[i].infCnt_ != 0)
	  std::cout
	    << "reducing " << ((lhsU_[i].infCnt_ < 0)?1:lhsU_[i].infCnt_)
	    << " infinities" ;
	else
	  std::cout
	    << "gap " << metricLE*(-normi) << ", metric " << metricLE ;
	std::cout << std::endl ;
      }
    }
/*
  For the >= side, it's (blowi-Ui) and the metric should be > 0 to propagate.
*/
    if (blowi > -infty_ && lhsU_[i].infCnt_ == 0) {
      if (lhsL_[i].infCnt_ != 0) {
        gapGE = infty_ ;
	metricGE = infty_ ;
	propagate = true ;
      } else {
	gapGE = blowi-lhsU_[i].bnd_ ;
	if (-posGapi > gapGE || negGapi < gapGE) {
	  propagate = true ;
	  metricGE = CoinMax(-(gapGE+posGapi),(gapGE-negGapi)) ;
	  metricGE /= normi ;
	}
      }
      if (verbosity_ >= 4 && propagate) {
	std::cout << "        queuing r(" << i << ") on lgap, " ;
	if (lhsL_[i].infCnt_ != 0)
	  std::cout
	    << "reducing " << ((lhsL_[i].infCnt_ < 0)?1:lhsL_[i].infCnt_)
	    << " infinities" ;
	else
	  std::cout
	    << "gap " << metricGE*normi << ", metric " << metricGE ;
	std::cout << std::endl ;
      }
    }
/*
  If this constraint should be processed, add it to the pile.
*/
    if (propagate) {
      const double metric = CoinMax(metricLE,metricGE) ;
      addToPending(i,0,metric) ;
    }
  }
/*
  Ok, we're loaded and ready to go ...
*/
  feas = false ;
  propagate(feas) ;
/*
  Report the results.
*/
  if (verbosity_ >= 3) {
    if (verbosity_ >= 4) std::cout << std::endl << "    " ;
    std::cout
      << ((feas)?"feasible":"infeasible") << "; tightened bounds on "
      << numVarBndChgs_ << " variables." << std::endl ;
  }
  if (verbosity_ >= 2 && numVarBndChgs_ > 0) {
    std::cout << "    Variable bound changes:" << std::endl ;
    for (int k = 0 ; k < numVarBndChgs_ ; k++)
      std::cout << "      " << varBndChgs_[k] << std::endl ;
  }
}
