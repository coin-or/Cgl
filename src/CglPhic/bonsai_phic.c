/*
  This file is a portion of the bonsaiG MILP code.

        Copyright (C) 2004 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation; either version 2 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc., 59
  Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
  This file contains the routines which implement the arc consistency
  algorithms used in bonsai.
*/

#include <math.h>
#include "bonsai.h"
#include "milp.h"

static char UNUSED sccsid[] = "@(#)consistency.c	3.7	11/11/04" ;

/*
  As an overview, the arc consistency algorithm propagates the effects of
  changes in the upper and lower bounds of variables. It maintains full
  k-consistency (i.e., it propagates changes until convergence).

  In operation, the package is notified of changes in variable bounds by
  calls to arc_chgvarbnd. It does some preliminary processing and then, if
  warranted, stacks the change for further propagation. The routine arc_phic
  processes the stacked changes, propagating them (and likely stacking new
  changes) until the system stablilises.

  Initialisation consists of data structure allocation (arc_init) and initial
  calculation of the constraint (row) lhs bounds (arc_initrowbnds). Where the
  initial lhs bounds can be used to derive variable bounds, the constraints
  are stacked for propagation, and arc_phic should be run.

  The package works directly on the variable and constraint bounds vectors
  associated with the constraint system. The vectors are modified, and when
  a variable is fixed the upper and lower variable bounds are set equal.

  As presently written, the package is initialised by arc_init to work on a
  specific constraint system. To allow the package to work on multiple systems
  with the constraint system passed as a parameter, the following modifications
  would be needed:
    * the bound change records should be allocated dynamically in arc_phic,
      instead of attached to the constraint system in arc_init.
    * externally visible routines would need to add consys as a parameter, and
      use it to set phicsys.
    * the fixed variable vector would need to come in as a parameter to
      arc_phic, which would set fxvars.
  While not strictly necessary:
    * a paranoid check should be added to make sure that constraint systems
      don't get interleaved (i.e., constraints are pushed onto the stack from
      one system and then the call to arc_phic specifies another).
*/

/*
  A lot of words about the mathematics used in the constraint propagation
  routines.

  The constraints are linear, of the form
    SUM<A+>[a<i,j>x<j>] + SUM<A->[a<i,j>x<j>] =  b<i>
  where
    A+ is the set of indices j such that a<i,j> is positive, and
    A- is the set of indices j such that a<i,j> is negative.
  The "=" can be either "<=", "=", or ">=".

  Let UB(x<j>) be the upper bound on x<j>, and LB(x<j>) be the lower bound on
  x<j>. Then we can calculate upper and lower bounds on the value of the
  left-hand side (lhs) of the constraint as
    UB(i) = SUM<A+>[a<i,j>UB(x<j>)] + SUM<A->[a<i,j>LB(x<j>)]
    LB(i) = SUM<A+>[a<i,j>LB(x<j>)] + SUM<A->[a<i,j>UB(x<j>)].

  Let x<t> be the current target variable (i.e., the one whose bound we're
  working on), with ub(x<t>) its upper bound and lb(x<t>) its lower bound.
  Suppose, for now, that we're looking at a >= constraint and that a<i,t>
  is positive. An isolation of x<t> (i.e., solving the constraint for x<t>)
  is
     x<t> >= (1/a<i,t>)[b<i> - (SUM<A+\t>[a<i,j>x<j>] + SUM<A->[a<i,j>x<j>])]
  This gives a lower bound on x<t>. To make that lower bound as small as
  possible, the two SUMs should be maximised. But then the two sums can be
  replaced, as
    SUM<A+\t>[a<i,j>UB(x<j>)]+SUM<A->[a<i,j>LB(x<j>)] = UB(i)-a<i,t>UB(x<t>)
  so what we actually end up calculating is
    LB(x<t>) = (1/a<i,t>)[b<i> - [UB(i)-a<i,t>UB(x<t>)]]
  If the new value for LB(x<t>) is larger than the old one, we tighten the
  bound and update LB(i).

  A similar line of reasoning can be followed if a<i,t> is negative. The
  calculation becomes
    UB(x<t>) = (1/a<i,t>)[b<i> - [UB(i)-a<i,t>LB(x<t>)]].
  
  Note that in both cases, we use UB(i); LB(i) is, in fact, not needed for a
  >= constraint.
  
  For <= constraints, the same reasoning applies, with appropriate changes in
  sign and the use of lower and upper bounds (left as an exercise for the
  maintainer, or read the supplementary documentation). Equalities can
  strengthen both bounds of a variable.
*/

/*
  A few words about accuracy and rounding. There doesn't seem to be any way
  to check the accuracy of the arc consistency calculation, nor any way to
  restore accuracy (analogous to refactoring the basis in LP). The code tries
  to cope with this by snapping to the nearest integer based on a tolerance
  (typically main_tols->arcsnap, scaled by the maximum of the target value
  and the currrent delta).  Otherwise it will round up (for upper bounds) or
  down (for lower bounds) to the nearest multiple of the rounding tolerance.

  This rounding is performed before applying floor/ceiling operations to
  integer variables, to avoid doing something like taking n-epsilon to n-1
  for epsilon smaller than the snap tolerance.

  It's not a perfect system, but then it's not a perfect world.
*/

extern double rint(double x) ;

#define RND_TOL main_tols->arcrnd
#define ROUND_UB(zz_val,zz_tol) \
  ((fabs((zz_val)-rint(zz_val)) < zz_tol)?\
     rint(zz_val):\
     ((zz_val) < 0.0)?\
	 ((zz_val)-fmod((zz_val),RND_TOL)):\
	 ((zz_val)-fmod((zz_val),RND_TOL)+RND_TOL))

#define ROUND_LB(zz_val,zz_tol) \
  ((fabs((zz_val)-rint(zz_val)) <= zz_tol)?\
     rint(zz_val):\
     ((zz_val) < 0.0)?\
	 ((zz_val)-fmod((zz_val),RND_TOL)-RND_TOL):\
	 ((zz_val)-fmod((zz_val),RND_TOL)))



/*
  The main data structure for the arc consistency algorithm is a stack of
  constraints waiting to be evaluated to see if we can use the bounds on the
  constraint to improve the bounds on the variables in the constraint. The
  arrays propstack and on_propstack handle this. propstack is worked as a
  stack, with the top pointed to by propstack_top.
  
  on_propstack provides an efficient check to see if a constraint is already
  on the stack.  For efficiency, last_popped tracks the constraint currently
  being used (since we don't want to push it back on the stack).
  
  The function propstack_cmp is used to sort the propagation stack (largest
  changes to the top) prior to popping a constraint for processing. need_sort
  is used to pass an indication from pushcon to phic when the stack needs to
  be sorted.

  The stack is handled pre-increment on push, post-decrement on pop.
  Allocated in arc_init, expanded if necessary in pushcon. on_propstack is
  attached to the constraint system and is grown as needed when the constraint
  system increases.

  Field		Description
  -----		-----------
  cndx		constraint index
  delta		accumulated change to UB(i) and LB(i) while the constraint
		has been waiting for propagation
  pct		sort key for the stack entries, delta/(UB(i)-LB(i))

  Both delta and pct are abused at times, particularly during initialisation
  of row bounds and, more generally, any time a bound changes from infinite
  to finite.
*/

typedef struct { int cndx ;
		 float delta ;
		 float pct ; } propstack_struct ;

static propstack_struct *propstack ;
static int propstack_top,propstack_sze,last_popped ;
static bool *on_propstack,need_sort ;

static int propstack_cmp (const void *v_prop1, const void *v_prop2)

{ propstack_struct *prop1,*prop2 ;

  prop1 = (propstack_struct *) v_prop1 ;
  prop2 = (propstack_struct *) v_prop2 ;

  if (prop1->pct < prop2->pct)
    return (-1) ;
  else
  if (prop1->pct > prop2->pct)
    return (1) ;
  else
    return (0) ; }

/*
  We also keep a few vectors and vars global to the propagation routines.
    * phicsys is the constraint system we're working on. Set from the
      system passed as a parameter to arc_init.
    * fxvars is used to return information about variables fixed in a given
      round of propagation. fxvars_sze is the capacity of fxvars. fxvars is
      allocated in arc_init, expanded if necessary in typ[L|E]_revise.
*/

static consys_struct *phicsys ;
static fxvar_struct *fxvars ;
static int fxvars_sze ;

#ifndef NDEBUG
/*
  If debug printing is enabled, we need some additional vectors. These are
  used to record changes to variable and constraint bounds. They are attached
  to the constraint system in arc_init.
*/

typedef struct { bool empty ;
		 int cnt ;
		 double olb ;
		 double oub ;
		 double nlb ;
		 double nub ; } vbndchg_struct ;

typedef struct { bool empty ;
		 int cnt ;
		 conbnd_struct olb ;
		 conbnd_struct oub ;
		 conbnd_struct nlb ;
		 conbnd_struct nub ; } cbndchg_struct ;

static vbndchg_struct *varbnd_chgs ;
static cbndchg_struct *conbnd_chgs ;
#endif



static bool pushcon (int cndx, double delta, double pct)

/*
  This routine pushes a constraint onto the constraint propagation stack,
  provided
    1) it isn't already there, and
    2) if the push request is for the constraint currently being evaluated,
       the constraint is an equality or range constraint.
  If the constraint is already stacked, the delta value is increased by the
  amount of this change in the bounds.

  All >= constraints should have been converted to <= on input, so if we
  see one, it's an error.

  Parameters:
    cndx:	the index of a constraint
    delta:	the absolute change in the constraint's bound
    pct:	the percentage change in the constraint's bound

  Returns: TRUE if the constraint is processed without error, FALSE if there's
	   a problem.
*/


{ int ndx ;
  contyp_enum contyp ;
  propstack_struct *tmpstk ;
  char *rtnnme = "pushcon" ;

# ifdef PARANOIA
/*
  Check that the constraint index is legitimate, that it's an equality,
  less-than-equal, or range constraint, and that the percent change is large
  enough for propagation.
*/
  if (cndx < 1 || cndx > phicsys->concnt)
  { errmsg(102,rtnnme,phicsys->nme,"constraint",cndx,1,phicsys->concnt) ;
    return (FALSE) ; }
  contyp = phicsys->ctyp[cndx] ;
  if (!(contyp == contypLE || contyp == contypEQ || contyp == contypRNG))
  { errmsg(453,rtnnme,phicsys->nme,consys_nme(phicsys,'c',cndx,FALSE,NULL),
	   cndx,consys_prtcontyp(contyp)) ;
    return (FALSE) ; }
# else
  contyp = phicsys->ctyp[cndx] ;
# endif

  if (pct < main_tols->arcprop)
  {
# ifndef NDEBUG
    if (main_opts->print.arcstk >= 1)
      warn(454,rtnnme,phicsys->nme,pct*100,"constraint",
	   consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx,main_tols->arcprop) ;
# endif
    return (TRUE) ; }
/*
  If the constraint isn't on the stack, attempt a push, provided that
  conditions 1) and 2) are met. If we have to expand the stack, grow by
  10% (and force at least one additional entry).

  Currently, we fail hard if the stack can't grow, but this might not be
  necessary. What would be the consequences of incomplete propagation for
  bonsai?
*/
  if (on_propstack[cndx] == FALSE &&
      (cndx != last_popped || contyp == contypEQ || contyp == contypRNG))
  { if (propstack_top >= propstack_sze-1)
    { ndx = 1+propstack_sze*1.1 ;
      tmpstk = (propstack_struct *)
			REALLOC(propstack,ndx*sizeof(propstack_struct)) ;
      if (tmpstk == NULL)
      { errmsg(455,rtnnme,propstack_sze,ndx) ;
	return (FALSE) ; }
#     ifndef NDEBUG
      if (main_opts->print.arcstk >= 1)
      { outfmt(logchn,gtxecho,
	       "\n\tpropstack size boosted from %d to %d entries.",
	       propstack_sze,ndx) ; }
#     endif
      propstack = tmpstk ;
      propstack_sze = ndx ; }
#   ifndef NDEBUG
    if (main_opts->print.arcstk >= 1)
    { outfmt(logchn,gtxecho,
	     "\n\tpushing constraint %s (%d), delta %.3g (%.3g%%).",
	     consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx,delta,pct) ; }
#   endif
    propstack[++propstack_top].cndx = cndx ;
    propstack[propstack_top].delta = delta ;
    propstack[propstack_top].pct = pct ;
    if (propstack_top > 0)
      if (propstack[propstack_top-1].pct >= pct) need_sort = TRUE ;
    on_propstack[cndx] = TRUE ;
    return (TRUE) ; }
/*
  If the constraint is already on the stack, but not the active constraint,
  update its delta and percent values. It's a fatal error to fail to find
  the constraint.
*/
  else
  if (on_propstack[cndx] == TRUE)
  { for (ndx = 0 ; ndx <= propstack_top ; ndx++)
    { if (propstack[ndx].cndx == cndx)
      { propstack[ndx].delta += delta ;
	propstack[ndx].pct += pct ;
	if (propstack_top > ndx)
	  if (propstack[ndx+1].pct < propstack[ndx].pct)
	    need_sort = TRUE ;
#	ifndef NDEBUG
	if (main_opts->print.arcstk >= 2)
	{ outfmt(logchn,gtxecho,"\n\tnew delta %g (%g%%), ",
		 propstack[ndx].delta,propstack[ndx].pct) ;
	  outfmt(logchn,gtxecho,"constraint %s (%d), pos'n %d.",
		 consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx,ndx) ; }
#	endif
	return (TRUE) ; } }
    errmsg(456,rtnnme,consys_nme(phicsys,'c',cndx,TRUE,NULL),cndx) ;
    return (FALSE) ; }
/*
  The final case: the constraint is the active constraint but not an equality
  or range constraint.
*/
  else
  if (cndx == last_popped)
  {
#   ifndef NDEBUG
    if (main_opts->print.arcstk >= 2)
      outfmt(logchn,gtxecho,
	     "\n\trequest to restack active %s constraint %s (%d) ignored.",
	     consys_prtcontyp(contyp),consys_nme(phicsys,'c',cndx,FALSE,NULL),
	     cndx) ;
#   endif
    return (TRUE) ; }
/*
  If we get here, we're in deep trouble.
*/
  errmsg(1,rtnnme,__LINE__) ;
  return (FALSE) ; }



static int popcon (void)

/*
  This routine pops a constraint from the constraint propagation stack.

  Parameters: none

  Returns: the constraint index from the top of the stack, -1 if the
	   stack is empty, -2 if there's an error.
*/

{ int cndx ;
  char *rtnnme = "popcon" ;

  if (propstack_top < 0)
  { last_popped = -1 ;
    return (-1) ; }

  cndx = propstack[propstack_top--].cndx ;
# ifdef PARANOIA
  if (cndx < 1 || cndx > phicsys->concnt)
  { errmsg(102,rtnnme,phicsys->nme,"constraint",cndx,1,phicsys->concnt) ;
    return (-2) ; }
# endif
  last_popped = cndx ;
  on_propstack[cndx] = FALSE ;

# ifndef NDEBUG
  if (main_opts->print.arcstk >= 1)
    outfmt(logchn,gtxecho,
	   "\n    popping %s constraint %s (%d), delta %g (%g%%).",
	   consys_prtcontyp(phicsys->ctyp[cndx]),
	   consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx,
	   propstack[propstack_top+1].delta,propstack[propstack_top+1].pct) ;
# endif

  return (cndx) ; }





static bool calc_rowbnd (int cndx, conbnd_struct *clb, conbnd_struct *cub,
			 double *p_delta, double *pct)

/*
  This routine calculates the upper and/or lower bound for the lhs of a
  constraint, returning them in the parameters. For a >= constraint, only
  UB(i) is calculated; for a <= constraint, only LB(i) is calculated.

  The parameter pct, if supplied, is used to indicate if the constraint
  should be pushed on the propagation stack for further processing. The value
  returned is meant to be handed to pushcon as the sort key. If the lhs bound
  is finite for only one variable, pct is set to infinity, as propagating the
  lhs bound will result in changing an infinite variable bound to a finite
  bound. If the lhs bound is finite for all variables in the constraint, pct
  is set to a value reflecting the tightness of the constraint. First the
  absolute looseness delta = min(|rhs-LB|,|UB-rhs|) is calculated. This is
  scaled to 1/min(delta/fabs(rhs),delta/1-norm(a<i>)). If the constraint is
  violated, pct is set to infinity.

  Parameters:
    cndx:	the index of the constraint
    clb:	(o) the lower bound for the left-hand side of the row (not
		used for a >= constraint)
    cub:	(o) the upper bound for the left-hand side of the row (not
		used for a <= constraint)
    p_delta:	(o) tightness of the constraint; min(|rhs-LB|,|UB-rhs|) for
		equalities; appropriate difference for <= and >= constraints.
    pct:	(o) propagation control, as outlined above

  Returns: TRUE if the calculation proceeds without error, FALSE otherwise.
*/

{ int pkndx,vndx,inflb,infub,inflndx,infundx ;
  double *vlb,*vub,lbnd,ubnd,bi,aij,delta,ainorm ;
  contyp_enum ctyp ;
  bool propagate ;
  pkvec_struct *pkrow ;
  pkcoeff_struct *coeff ;
  char *rtnnme = "calc_rowbnd" ;

  inflndx = -1 ;
  infundx = -1 ;
  bi = 0.0 ;

# ifdef PARANOIA
  if (cndx < 1 || cndx > phicsys->concnt)
  { errmsg(102,rtnnme,phicsys->nme,"constraint",cndx,1,phicsys->concnt) ;
    return (FALSE) ; }
# endif
  ctyp = phicsys->ctyp[cndx] ;
# ifdef PARANOIA
  if (!VALID_CONTYPE(ctyp))
  { errmsg(4,rtnnme,"constraint type",consys_prtcontyp(ctyp)) ;
    return (FALSE) ; }
  if (ctyp != contypGE && clb == NULL)
  { errmsg(2,rtnnme,"clb") ;
    return (FALSE) ; }
  if (ctyp != contypLE && cub == NULL)
  { errmsg(2,rtnnme,"cub") ;
    return (FALSE) ; }
# endif
  vlb = phicsys->vlb ;
  vub = phicsys->vub ;
# ifdef PARANOIA
  if (vlb == NULL)
  { errmsg(101,rtnnme,phicsys->nme,consys_assocnme(NULL,CONSYS_VLB)) ;
    return (FALSE) ; }
  if (vub == NULL)
  { errmsg(101,rtnnme,phicsys->nme,consys_assocnme(NULL,CONSYS_VUB)) ;
    return (FALSE) ; }
# endif
/*
  Initialise the bounds to 0 (we'll accumulate them), and hope they'll be
  finite. Then get the coefficients for the constraint.
*/
  ubnd = 0.0 ;
  infub = 0 ;
  lbnd = 0.0 ;
  inflb = 0 ;
  pkrow = NULL ;
  if (consys_getrow_pk(phicsys,cndx,&pkrow) == FALSE)
  { errmsg(112,rtnnme,phicsys->nme,"retrieve","constraint",
	   consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
    if (pkrow != NULL) pkvec_free(pkrow) ;
    return (FALSE) ; }
/*
  Open up a loop to step through the variables. For each variable, check the
  sign of the coefficient. If the coefficient is positive, then the
  variable's upper bound contributes to the upper bound of the constraint,
  and the variable's lower bound contributes to the lower bound of the
  constraint. If the coefficient is negative, the contributions are
  reversed.  If a variable's bound is not finite (i.e., +infinity or
  -infinity), then it contributes an infinity to the bound.
*/
  ainorm = 0 ;
  for (pkndx = 0, coeff = pkrow->coeffs ; pkndx < pkrow->cnt ; coeff++, pkndx++)
  { aij = coeff->val ;
    vndx = coeff->ndx ;
#   ifdef PARANOIA
    if (vndx < 1 || vndx > phicsys->archvcnt)
    { errmsg(102,rtnnme,phicsys->nme,"variable",vndx,1,phicsys->archvcnt) ;
      return (FALSE) ; }
    if (fabs(aij) < main_lptols->zero && flgon(phicsys->opts,CONSYS_WRNZERO))
      warn(450,rtnnme,cndx,vndx,consys_nme(phicsys,'c',cndx,FALSE,NULL),
	   consys_nme(phicsys,'v',vndx,FALSE,NULL)) ;
#   endif
    ainorm += fabs(aij) ;
    if (aij > 0)
    { if (ctyp != contypLE)
      { if (vub[vndx] < main_lptols->inf)
	  ubnd += aij*vub[vndx] ;
	else
	{ infub++ ;
	  infundx = vndx ; } }
      if (ctyp != contypGE)
      { if (vlb[vndx] > -main_lptols->inf)
	  lbnd += aij*vlb[vndx] ;
	else
	{ inflb++ ;
	  inflndx = vndx ; } } }
    else
    { if (ctyp != contypLE)
      { if (vlb[vndx] > -main_lptols->inf)
	  ubnd += aij*vlb[vndx] ;
	else
	{ infub++ ;
	  infundx = vndx ; } }
      if (ctyp != contypGE)
      { if (vub[vndx] < main_lptols->inf)
	  lbnd += aij*vub[vndx] ;
	else
	{ inflb++ ;
	  inflndx = vndx ; } } }

#   ifndef NDEBUG
    if (main_opts->print.arcinit >= 2)
    { outfmt(logchn,gtxecho,"\n\ta<%d,%d> = %g, ",cndx,vndx,aij) ;
      outfmt(logchn,gtxecho," %g <= %s <= %g, ",vlb[vndx],
	     consys_nme(phicsys,'v',vndx,FALSE,NULL),vub[vndx]) ;
      if (ctyp != contypGE)
      { if (inflb > 0) outfmt(logchn,gtxecho,"%d*%g+",inflb,-main_lptols->inf) ;
	outfmt(logchn,gtxecho,"%g",lbnd) ; }
      else
      { outfmt(logchn,gtxecho,"n/a") ; }
      outfmt(logchn,gtxecho," <= %s <= ",
	     consys_nme(phicsys,'c',cndx,FALSE,NULL)) ;
      if (ctyp != contypLE)
      { if (infub > 0) outfmt(logchn,gtxecho,"%d*%g+",infub,main_lptols->inf) ;
	outfmt(logchn,gtxecho,"%g.",ubnd) ; }
      else
      { outfmt(logchn,gtxecho,"n/a.") ; } }
#   endif

    }
/*
  We've processed the variables. Clean up and recode the bounds.
*/
  if (ctyp != contypLE)
  { if (infub == 1)
      cub->inf = -infundx ;
    else
      cub->inf = infub ;
    setcleanzero(ubnd,main_lptols->zero) ;
    cub->bnd = ubnd ;
    cub->revs = 0 ; }
  if (ctyp != contypGE)
  { if (inflb == 1)
      clb->inf = -inflndx ;
    else
      clb->inf = inflb ;
    setcleanzero(lbnd,main_lptols->zero) ;
    clb->bnd = lbnd ;
    clb->revs = 0 ; }

# ifndef NDEBUG
  if (main_opts->print.arcinit >= 2)
  { if (ctyp != contypLE && infub == 1)
      outfmt(logchn,gtxecho,"\n\tUB(%s) finite for %s (%d) only.",
	     consys_nme(phicsys,'c',cndx,FALSE,NULL),
	     consys_nme(phicsys,'v',-cub->inf,FALSE,NULL),-cub->inf) ;
    if (ctyp != contypGE && inflb == 1)
      outfmt(logchn,gtxecho,"\n\tLB(%s) finite for %s (%d) only.",
	     consys_nme(phicsys,'c',cndx,FALSE,NULL),
	     consys_nme(phicsys,'v',-clb->inf,FALSE,NULL),-clb->inf) ; }
# endif

/*
  Now figure out whether further propagation is necessary, and how urgent it
  is to get to this constraint (see comments in header). If the calculated
  value of pct (b-LB, or UB-blow) is negative, the constraint is violated.
*/
  if (pct != NULL)
  { propagate = TRUE ;
    ctyp = phicsys->ctyp[cndx] ;
    if (ctyp == contypLE)
    { if (inflb == 1)
      { delta = main_lptols->inf ; }
      else
      if (inflb == 0)
      { bi = phicsys->rhs[cndx] ;
	delta = bi-lbnd ; }
      else
      { propagate = FALSE ; } }
    else
    if (ctyp == contypGE)
    { if (infub == 1)
      { delta = main_lptols->inf ; }
      else
      if (infub == 0)
      { bi = phicsys->rhs[cndx] ;
	delta = ubnd-bi ; }
      else
      { propagate = FALSE ; } }
    else
    if (ctyp == contypEQ)
    { if (infub == 1 || inflb == 1)
      { delta = main_lptols->inf ; }
      else
      if (infub == 0 || inflb == 0)
      { bi = phicsys->rhs[cndx] ;
	if (infub == 0 && inflb == 0)
	  delta = minn(bi-lbnd,ubnd-bi) ;
	else
	if (inflb == 0)
	  delta = bi-lbnd ;
	else
	  delta = ubnd-bi ; }
      else
      { propagate = FALSE ; } }
    else	/* contypRNG */
    { if (infub == 1 || inflb == 1)
      { delta = main_lptols->inf ; }
      else
      if (infub == 0 || inflb == 0)
      { bi = phicsys->rhs[cndx] ;
	if (infub == 0 && inflb == 0)
	{ if (bi-lbnd < ubnd-phicsys->rhslow[cndx])
	  { delta = bi-lbnd ; }
	  else
	  { bi = phicsys->rhslow[cndx] ;
	    delta = ubnd-bi ; } }
	else
	if (inflb == 0)
	{ delta = bi-lbnd ; }
	else
	{ bi = phicsys->rhslow[cndx] ;
	  delta = ubnd-bi ; } }
      else
      { propagate = FALSE ; } }
    if (propagate == TRUE)
    { if (delta < 0 || delta >= main_lptols->inf)
	*pct = main_lptols->inf ;
      else
      { *pct = minn(delta/(1+fabs(bi)),delta/ainorm) ;
	*pct = 1/(1+*pct) ; }
      if (p_delta != NULL) *p_delta = delta ; }
    else
    { *pct = 0 ; } }

  pkvec_free(pkrow) ;

  return (TRUE) ; }



bool arc_initrowbnds (void)

/*
  This routine works through the constraints, calculating the initial upper and
  lower bounds for each constraint. Each constraint which ends up with a finite
  upper or lower bound is pushed onto the constraint propagation stack for the
  initial round of constraint propagation.

  Parameters: none

  Returns: TRUE if the calculation proceeds without error, FALSE otherwise.
*/

{ int cndx,pushcnt ;
  conbnd_struct conub,conlb ;
  double delta,pct ;
  char *rtnnme = "arc_initrowbnds" ;

  pushcnt = 0 ;

# ifndef NDEBUG
  if (main_opts->print.arcinit >= 1)
    outfmt(logchn,gtxecho,"\n  initialising constraint bounds.") ;
# endif
/*
  Open up a loop to step through the constraints (rows). calc_rowbnd does the
  actual work of calculating the bound. For a <= constraint, only LB(i) is
  calculated; for a >= constraint, only UB(i).
*/
  for (cndx = 1 ; cndx <= phicsys->concnt ; cndx++)
  {
#   ifndef NDEBUG
    if (main_opts->print.arcinit >= 1)
    { outfmt(logchn,gtxecho,"\n    %s (%d): ",
	     consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ; }
#   endif

    if (calc_rowbnd(cndx,&conlb,&conub,&delta,&pct) == FALSE)
    { errmsg(451,rtnnme,phicsys->nme,
	     consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
      return (FALSE) ; }
/*
  Set the bounds in the data structure. The paranoia code here inserts
  {inf,bnd} = {1,SNaN} constraint bounds entries that should not be accessed
  if the arc consistency code is working properly. This makes both fields
  invalid and should be sufficient to catch any access.
*/
    if (phicsys->ctyp[cndx] == contypLE)
    { conub.revs = 0 ;
      conub.inf = 1 ;
      conub.bnd = signaling_nan(0) ; }
    if (phicsys->ctyp[cndx] == contypGE)
    { conlb.revs = 0 ;
      conlb.inf = 1 ;
      conlb.bnd = signaling_nan(0) ; }
    phicsys->cub[cndx] = conub ;
    phicsys->clb[cndx] = conlb ;

#   ifndef NDEBUG
    if (main_opts->print.arcinit >= 1)
    { if (main_opts->print.arcinit >= 2)
	outfmt(logchn,gtxecho,"\n\tfinal bounds: ") ;
      if (phicsys->ctyp[cndx] != contypGE)
      { outfmt(logchn,gtxecho,"%s = %s",consys_conbndval(&phicsys->clb[cndx]),
	       consys_conbndnme('L',cndx,&phicsys->clb[cndx])) ; }
      else
      { outfmt(logchn,gtxecho,"n/a = LB(%d) ",cndx) ; }
      if (phicsys->ctyp[cndx] != contypLE)
      { outfmt(logchn,gtxecho,"%s = %s",
	       consys_conbndnme('U',cndx,&phicsys->cub[cndx]),
	       consys_conbndval(&phicsys->cub[cndx])) ; }
      else
      { outfmt(logchn,gtxecho,"UB(%d) = n/a",cndx) ; } }
#   endif
/*
  If calc_rowbnd says push, do it. Either one or both of the bounds is
  finite, or the constraint is an explicit upper/lower bound constraint.
*/
    if (pct > 0)
    { if (pushcon(cndx,delta,pct) == FALSE)
      { errmsg(452,rtnnme,phicsys->nme,
	       consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
	return (FALSE) ; }
      pushcnt++ ; } }

# ifndef NDEBUG
/*
  If debug printing is enabled, we need the `off' row bounds also, to avoid
  accessing unitialised bound entries.
*/
  if (arc_rowoffbnds(phicsys) == FALSE)
  { errmsg(853,rtnnme) ;
    return (FALSE) ; }
  if (main_opts->print.arcinit >= 1)
  { outfmt(logchn,gtxecho,
	   "\n    %d constraints pushed for arc consistency processing.",
	   pushcnt) ; }
# endif

  return (TRUE) ; }



bool arc_init (consys_struct *consys)

/*
  This routine sets the constraint system to be processed, attaches vectors,
  allocates space and does other such initialisations for the constraint
  propagation machinery.

  Parameters:
    consys:     the constraint system

  Returns: TRUE if initialisation completes without error, FALSE otherwise
*/

{ int ndx ;
  char *rtnnme = "arc_init" ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (main_opts->phic == 0)
  { errmsg(467,rtnnme) ;
    return (FALSE) ; }
# endif
/*
  Set the constraint system.
*/
  phicsys = consys ;
/*
  Establish the constraint propagation stack, and an index vector that records
  whether a constraint is present on the stack.
*/
  propstack_sze = phicsys->concnt/4 ;
  propstack =
	(propstack_struct *) MALLOC(propstack_sze*sizeof(propstack_struct)) ;
  propstack_top = -1 ;
  last_popped = -1 ;
  on_propstack = NULL ;
  if (consys_attach(phicsys,CONSYS_COL,
		    sizeof(bool),(void **) &on_propstack) == FALSE)
  { errmsg(100,rtnnme,phicsys->nme,"propagation stack presence vector") ;
    return (FALSE) ; }
/*
  Create a fixed variable vector that will be used to return the indices of
  any variables fixed during a round of constraint propagation.
*/
  fxvars_sze = phicsys->varcnt/4 ;
  fxvars = (fxvar_struct *) CALLOC(fxvars_sze,sizeof(fxvar_struct)) ;
/*
  Attach arrays to hold the row bounds, if they're not already present. These
  are associated vectors.
*/
# ifdef PARANOIA
  if (flgon(phicsys->parts,CONSYS_CUB) && phicsys->cub == NULL)
  { errmsg(113,rtnnme,phicsys->nme,consys_assocnme(NULL,CONSYS_CUB)) ;
    return (FALSE) ; }
# endif
  if (flgoff(phicsys->parts,CONSYS_CUB))
  { phicsys->cub = NULL ;
    if (consys_attach(phicsys,CONSYS_COL,
		      sizeof(conbnd_struct),(void **) &phicsys->cub) == FALSE)
    { errmsg(100,rtnnme,phicsys->nme,consys_assocnme(NULL,CONSYS_CUB)) ;
      return (FALSE) ; }
    setflg(phicsys->parts,CONSYS_CUB) ; }
# ifdef PARANOIA
  if (flgon(phicsys->parts,CONSYS_CLB) && phicsys->cub == NULL)
  { errmsg(113,rtnnme,phicsys->nme,consys_assocnme(NULL,CONSYS_CLB)) ;
    return (FALSE) ; }
# endif
  if (flgoff(phicsys->parts,CONSYS_CLB))
  { phicsys->clb = NULL ;
    if (consys_attach(phicsys,CONSYS_COL,
		      sizeof(conbnd_struct),(void **) &phicsys->clb) == FALSE)
    { errmsg(100,rtnnme,phicsys->nme,consys_assocnme(NULL,CONSYS_CLB)) ;
      return (FALSE) ; }
    setflg(phicsys->parts,CONSYS_CLB) ; }

# ifndef NDEBUG
/*
  Attach arrays that will record changes to variable and constraint bounds
  for use in debugging prints.
*/
  varbnd_chgs = (vbndchg_struct *)
			MALLOC((phicsys->colsze+1)*sizeof(vbndchg_struct)) ;
  if (consys_attach(phicsys,CONSYS_ROW,
		    sizeof(vbndchg_struct),(void **) &varbnd_chgs) == FALSE)
  { errmsg(100,rtnnme,phicsys->nme,"variable bound change record") ;
    return (FALSE) ; }
  for (ndx = 0 ; ndx <= phicsys->colsze ; ndx++) varbnd_chgs[ndx].empty = TRUE ;
  conbnd_chgs = (cbndchg_struct *)
			MALLOC((phicsys->rowsze+1)*sizeof(cbndchg_struct)) ;
  if (consys_attach(phicsys,CONSYS_COL,
		    sizeof(cbndchg_struct),(void **) &conbnd_chgs) == FALSE)
  { errmsg(100,rtnnme,phicsys->nme,"constraint bound change record") ;
    return (FALSE) ; }
  for (ndx = 0 ; ndx <= phicsys->rowsze ; ndx++) conbnd_chgs[ndx].empty = TRUE ;
# endif

  return (TRUE) ; }



void arc_free (void)

/*
  This routine frees the vectors used by the arc consistency package.

  Parameters: none
  Returns: undefined
*/

{ if (propstack != NULL) FREE(propstack) ;
  if (on_propstack != NULL) FREE(on_propstack) ;
  if (fxvars != NULL) FREE(fxvars) ;
# ifndef NDEBUG
  if (varbnd_chgs != NULL) FREE(varbnd_chgs) ;
  if (conbnd_chgs != NULL) FREE(conbnd_chgs) ;
# endif
  return ; }


void arc_clear (void)

/*
  This routine cleans off the constraint propagation stack. It is intended to
  be used to clean up when constraints have been pushed for propagation but
  propagation never takes place (typically because infeasibility is detected
  somewhere else). It's really just a streamlined version of popcon.

  We can only fail here if we're paranoid, hence the void return.

  It might be faster to implement this by setting propstack_top = -1 and using
  memset to clear on_propstack. The code below is betting that there won't be
  many constraints on propstack.

  Parameters: none

  Returns: undefined
*/

{ int cndx ;
  char *rtnnme = "arc_clear" ;

  last_popped = -1 ;

  while (propstack_top >= 0)
  { cndx = propstack[propstack_top--].cndx ;
#   ifdef PARANOIA
    if (cndx < 1 || cndx > phicsys->concnt)
    { errmsg(102,rtnnme,phicsys->nme,"constraint",cndx,1,phicsys->concnt) ;
      return ; }
#   endif
    on_propstack[cndx] = FALSE ;
#   ifndef NDEBUG
    if (main_opts->print.arcstk >= 1)
      outfmt(logchn,gtxecho,
	     "\n    clearing %s constraint %s (%d).",
	     consys_prtcontyp(phicsys->ctyp[cndx]),
	     consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
#   endif
  }

  return ; }



static bool recalc_rowbnd (int i, char cbnd, conbnd_struct *cbi)

/*
  This routine recalculates a row bound and does some paranoid checks.  It's
  split out of arc_chgvarbnd to keep that code from becoming completely
  unreadable. There's only about 10 lines of nonparanoid code.

  The paranoid part looks for differences in the value of the iteratively
  updated bound and the recalculated bound. Differences in the finite part of
  the bound rate a warning. Differences in the infinity part of the bound are
  considered fatal.

  Parameters:
    i:		constraint index
    cbnd:	the bound we're updating; 'L' for LB(i), 'U' for UB(i)
    cbi:	the new bound

  Returns: TRUE if the recalculation goes ok, FALSE otherwise
*/

{ char *rtnnme = "recalc_rowbnd" ;

# ifdef PARANOIA
  conbnd_struct alt_cbi ;
  char buf1[128],buf2[128] ;
# endif

  if (cbnd == 'L')
  {
#   ifdef PARANOIA
    alt_cbi = phicsys->cub[i] ;
#   endif
    if (calc_rowbnd(i,cbi,&phicsys->cub[i],NULL,NULL) == FALSE)
    { errmsg(451,rtnnme,phicsys->nme,consys_nme(phicsys,'c',i,FALSE,NULL),i) ;
      return (FALSE) ; }
#   ifdef PARANOIA
    if (alt_cbi.inf != phicsys->cub[i].inf)
    { (void) outfxd(buf1,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('U',i,&alt_cbi),
		    consys_conbndval(&alt_cbi)) ;
      (void) outfxd(buf2,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('U',i,&phicsys->cub[i]),
		    consys_conbndval(&phicsys->cub[i])) ;
      errmsg(464,rtnnme,phicsys->nme,buf1,buf2) ;
      return (FALSE) ; }
    if (alt_cbi.inf == 0 && !atbnd(alt_cbi.bnd,phicsys->cub[i].bnd))
    { (void) outfxd(buf1,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('U',i,&alt_cbi),
		    consys_conbndval(&alt_cbi)) ;
      (void) outfxd(buf2,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('U',i,&phicsys->cub[i]),
		    consys_conbndval(&phicsys->cub[i])) ;
      warn(465,rtnnme,phicsys->nme,buf1,buf2,
	   fabs(alt_cbi.bnd-phicsys->cub[i].bnd)) ; }
#   endif
  }
  else
  {
#   ifdef PARANOIA
    alt_cbi = phicsys->clb[i] ;
#   endif
    if (calc_rowbnd(i,&phicsys->clb[i],cbi,NULL,NULL) == FALSE)
    { errmsg(451,rtnnme,phicsys->nme,consys_nme(phicsys,'c',i,FALSE,NULL),i) ;
      return (FALSE) ; }
#   ifdef PARANOIA
    if (alt_cbi.inf != phicsys->clb[i].inf)
    { (void) outfxd(buf1,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('L',i,&alt_cbi),
		    consys_conbndval(&alt_cbi)) ;
      (void) outfxd(buf2,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('L',i,&phicsys->clb[i]),
		    consys_conbndval(&phicsys->clb[i])) ;
      errmsg(464,rtnnme,phicsys->nme,buf1,buf2) ;
      return (FALSE) ; }
    if (alt_cbi.inf == 0 && !atbnd(alt_cbi.bnd,phicsys->clb[i].bnd))
    { (void) outfxd(buf1,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('L',i,&alt_cbi),
		    consys_conbndval(&alt_cbi)) ;
      (void) outfxd(buf2,-(sizeof(buf1)-1),'l',"%s = %s",
		    consys_conbndnme('L',i,&phicsys->clb[i]),
		    consys_conbndval(&phicsys->clb[i])) ;
      warn(465,rtnnme,phicsys->nme,buf1,buf2,
	   fabs(alt_cbi.bnd-phicsys->clb[i].bnd)) ; }
#   endif
  }

  return (TRUE) ; }



bool arc_calcrowbnds (consys_struct *sys, int cndx,
		      conbnd_struct *clb, conbnd_struct *cub)

/*
  This routine is made available for clients outside of consistency.c that
  want to calculate upper and lower bounds on the lhs of a constraint (most
  often for installing a cutting plane). It calculates the upper and lower
  bound for the lhs of a constraint, returning them in the parameters. Both
  LB(i) and UB(i) are always calculated, but returned only if the parameter
  is non-NULL.

  This routine is a sort of stripped-down version of calc_rowbnd. It also takes
  an explicit constraint system parameter. (At some point in the future, this
  routine and calc_rowbnd should merge.)

  Parameters:
    cndx:	the index of the constraint
    clb:	(o) if non-NULL, used to return LB(i)
    cub:	(o) if non-NULL, used to return UB(i)

  Returns: TRUE if the calculation proceeds without error, FALSE otherwise.
*/

{ int pkndx,vndx,inflb,infub,inflndx,infundx ;
  double *vlb,*vub,lbnd,ubnd,aij ;
  pkvec_struct *pkrow ;
  pkcoeff_struct *coeff ;

  char *rtnnme = "arc_calcrowbnds" ;

  inflndx = -1 ;
  infundx = -1 ;

# ifdef PARANOIA
  if (sys == NULL)
  { errmsg(2,rtnnme,"constraint system") ;
    return (FALSE) ; }
  if (cndx < 1 || cndx > sys->concnt)
  { errmsg(102,rtnnme,sys->nme,"constraint",cndx,1,sys->concnt) ;
    return (FALSE) ; }
# endif

  vlb = sys->vlb ;
  vub = sys->vub ;

# ifdef PARANOIA
  if (vlb == NULL)
  { errmsg(101,rtnnme,sys->nme,consys_assocnme(NULL,CONSYS_VLB)) ;
    return (FALSE) ; }
  if (vub == NULL)
  { errmsg(101,rtnnme,sys->nme,consys_assocnme(NULL,CONSYS_VUB)) ;
    return (FALSE) ; }
# endif
/*
  Initialise the bounds to 0 (we'll accumulate them), and hope they'll be
  finite. Then get the coefficients for the constraint.
*/
  ubnd = 0.0 ;
  infub = 0 ;
  lbnd = 0.0 ;
  inflb = 0 ;
  pkrow = NULL ;
  if (consys_getrow_pk(sys,cndx,&pkrow) == FALSE)
  { errmsg(112,rtnnme,sys->nme,"retrieve","constraint",
	   consys_nme(sys,'c',cndx,FALSE,NULL),cndx) ;
    if (pkrow != NULL) pkvec_free(pkrow) ;
    return (FALSE) ; }
/*
  Open up a loop to step through the variables. For each variable, check the
  sign of the coefficient. If the coefficient is positive, then the
  variable's upper bound contributes to the upper bound of the constraint,
  and the variable's lower bound contributes to the lower bound of the
  constraint. If the coefficient is negative, the contributions are
  reversed.  If a variable's bound is not finite (i.e., +infinity or
  -infinity), then it contributes an infinity to the bound.
*/
  for (pkndx = 0, coeff = pkrow->coeffs ; pkndx < pkrow->cnt ; coeff++, pkndx++)
  { aij = coeff->val ;
    vndx = coeff->ndx ;
#   ifdef PARANOIA
    if (vndx < 1 || vndx > sys->archvcnt)
    { errmsg(102,rtnnme,sys->nme,"variable",vndx,1,sys->archvcnt) ;
      return (FALSE) ; }
    if (fabs(aij) < main_lptols->zero && flgon(sys->opts,CONSYS_WRNZERO))
      warn(450,rtnnme,cndx,vndx,consys_nme(sys,'c',cndx,FALSE,NULL),
	   consys_nme(sys,'v',vndx,FALSE,NULL)) ;
#   endif
    if (aij > 0)
    { if (vub[vndx] < main_lptols->inf)
	ubnd += aij*vub[vndx] ;
      else
      { infub++ ;
	infundx = vndx ; }
      if (vlb[vndx] > -main_lptols->inf)
	lbnd += aij*vlb[vndx] ;
      else
      { inflb++ ;
	inflndx = vndx ; } }
    else
    { if (vlb[vndx] > -main_lptols->inf)
	ubnd += aij*vlb[vndx] ;
      else
      { infub++ ;
	infundx = vndx ; }
      if (vub[vndx] < main_lptols->inf)
	lbnd += aij*vub[vndx] ;
      else
      { inflb++ ;
	inflndx = vndx ; } } }
  pkvec_free(pkrow) ;
/*
  We've processed the variables. Clean up and recode the bounds.
*/
  if (cub != NULL)
  { if (infub == 1)
      cub->inf = -infundx ;
    else
      cub->inf = infub ;
    setcleanzero(ubnd,main_lptols->zero) ;
    cub->bnd = ubnd ;
    cub->revs = 0 ; }
  if (clb != NULL)
  { if (inflb == 1)
      clb->inf = -inflndx ;
    else
      clb->inf = inflb ;
    setcleanzero(lbnd,main_lptols->zero) ;
    clb->bnd = lbnd ;
    clb->revs = 0 ; }

  return (TRUE) ; }




bool arc_rowoffbnds (consys_struct *consys)

/*
  This routine calculates the `off' bound for <= and >= constraints, adding
  them to the clb and cub arrays of consys. The reason for this routine is
  that when propagating bounds changes, only UB(i) is used for >= constraints,
  and only LB(i) for <= constraints. Hence there's no utility in incrementally
  updating these bounds during bound propagation. But, we do need these values
  when we're considering the suitability of a constraint as a branching
  hyperplane. So, calculate them in one swoop, at the point where we need
  them.

  Unlike other routines in consistency.c, arc_rowoffbnds already takes a
  constraint system as a parameter --- makes more sense in the context of
  its use.

  Parameters:
    consys	constraint system to be processed.

  Returns: TRUE if the calculation proceeds without error, FALSE otherwise.
*/

{ int i,pkndx,j,inflb,infub,inflndx,infundx ;
  double *vlb,*vub,lbnd,ubnd,aij ;
  contyp_enum ctypi ;
  pkvec_struct *pkrow ;
  pkcoeff_struct *coeff ;
  conbnd_struct *cubi,*clbi ;
  char *rtnnme = "arc_rowoffbnds" ;

  inflb = -1 ;
  infub = -1 ;
  inflndx = -1 ;
  infundx = -1 ;
  lbnd = 0.0 ;
  ubnd = 0.0 ;
  clbi = NULL ;
  cubi = NULL ;

# ifdef PARANOIA
  if (consys == NULL)
  { errmsg(2,rtnnme,"consys") ;
    return (FALSE) ; }
  if (consys->cub == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CUB)) ;
    return (FALSE) ; }
  if (consys->clb == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_CLB)) ;
    return (FALSE) ; }
# endif

  vlb = consys->vlb ;
  vub = consys->vub ;

# ifdef PARANOIA
  if (vlb == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VLB)) ;
    return (FALSE) ; }
  if (vub == NULL)
  { errmsg(101,rtnnme,consys->nme,consys_assocnme(NULL,CONSYS_VUB)) ;
    return (FALSE) ; }
# endif
# ifndef NDEBUG
  if (main_opts->print.conclass >= 4)
    outfmt(logchn,gtxecho,"\n      calculating constraint off bounds.") ;
# endif
/*
  The outer loop steps through the constraints of the system, looking for
  <= and >= constraints.
*/
  pkrow = NULL ;
  for (i = 1 ; i <= consys->concnt ; i++)
  { ctypi = consys->ctyp[i] ;
#   ifdef PARANOIA
    if (!VALID_CONTYPE(ctypi))
    { errmsg(4,rtnnme,"constraint type",consys_prtcontyp(ctypi)) ;
      return (FALSE) ; }
#   endif
#   ifndef NDEBUG
    if (main_opts->print.conclass >= 4)
    { outfmt(logchn,gtxecho,"\n      %s (%d) (%s): ",
	     consys_nme(consys,'c',i,FALSE,NULL),i,
	     consys_prtcontyp(ctypi)) ; }
#   endif
/*
  If we're looking at a >= or <= constraint, initialise the appropriate bound
  to 0. Then get the coefficients for the constraint.
*/
    switch (ctypi)
    { case contypLE:
      { cubi = &consys->cub[i] ;
	ubnd = 0.0 ;
	infub = 0 ;
	break ; }
      case contypGE:
      { clbi = &consys->cub[i] ;
	lbnd = 0.0 ;
	inflb = 0 ;
	break ; }
      default:
      { continue ; } }
/*
  Get the coefficients for the constraint.
*/
    if (consys_getrow_pk(consys,i,&pkrow) == FALSE)
    { errmsg(112,rtnnme,consys->nme,"retrieve","constraint",
	     consys_nme(consys,'c',i,FALSE,NULL),i) ;
      if (pkrow != NULL) pkvec_free(pkrow) ;
      return (FALSE) ; }
/*
  Open up a loop to step through the variables. For each variable, check the
  sign of the coefficient. If the coefficient is positive, then the
  variable's upper bound contributes to the upper bound of the constraint,
  and the variable's lower bound contributes to the lower bound of the
  constraint. If the coefficient is negative, the contributions are
  reversed.  If a variable's bound is not finite (i.e., +infinity or
  -infinity), then it contributes an infinity to the bound.
*/
    for (pkndx = 0, coeff = pkrow->coeffs ;
	 pkndx < pkrow->cnt ;
	 coeff++, pkndx++)
    { aij = coeff->val ;
      j = coeff->ndx ;
#     ifdef PARANOIA
      if (j < 1 || j > consys->archvcnt)
      { errmsg(102,rtnnme,consys->nme,"variable",j,1,consys->archvcnt) ;
	return (FALSE) ; }
      if (fabs(aij) < main_lptols->zero && flgon(consys->opts,CONSYS_WRNZERO))
	warn(450,rtnnme,i,j,consys_nme(consys,'c',i,FALSE,NULL),
	     consys_nme(consys,'v',j,FALSE,NULL)) ;
#     endif
      if (aij > 0)
      { if (ctypi == contypLE)
	{ if (vub[j] < main_lptols->inf)
	    ubnd += aij*vub[j] ;
	  else
	  { infub++ ;
	    infundx = j ; } }
	else
	{ if (vlb[j] > -main_lptols->inf)
	    lbnd += aij*vlb[j] ;
	  else
	  { inflb++ ;
	    inflndx = j ; } } }
      else
      { if (ctypi == contypLE)
	{ if (vlb[j] > -main_lptols->inf)
	    ubnd += aij*vlb[j] ;
	  else
	  { infub++ ;
	    infundx = j ; } }
	else
	{ if (vub[j] < main_lptols->inf)
	    lbnd += aij*vub[j] ;
	  else
	  { inflb++ ;
	    inflndx = j ; } } }

#     ifndef NDEBUG
      if (main_opts->print.conclass >= 5)
      { outfmt(logchn,gtxecho,"\n\ta<%d,%d> = %g, ",i,j,aij) ;
	outfmt(logchn,gtxecho," %g <= %s <= %g, ",vlb[j],
	       consys_nme(consys,'v',j,FALSE,NULL),vub[j]) ;
	if (ctypi == contypGE)
	{ if (inflb > 0)
	    outfmt(logchn,gtxecho,"%d*%g+",inflb,-main_lptols->inf) ;
	  outfmt(logchn,gtxecho,"%g",lbnd) ; }
	else
	{ outfmt(logchn,gtxecho,"n/a") ; }
	outfmt(logchn,gtxecho," <= %s <= ",
	       consys_nme(consys,'c',i,FALSE,NULL)) ;
	if (ctypi == contypLE)
	{ if (infub > 0)
	    outfmt(logchn,gtxecho,"%d*%g+",infub,main_lptols->inf) ;
	  outfmt(logchn,gtxecho,"%g.",ubnd) ; }
	else
	{ outfmt(logchn,gtxecho,"n/a.") ; } }
#     endif

    }
/*
  We've processed the variables. Clean up and recode the bound.
*/
  if (ctypi == contypLE)
  { if (infub == 1)
      cubi->inf = -infundx ;
    else
      cubi->inf = infub ;
    setcleanzero(ubnd,main_lptols->zero) ;
    cubi->bnd = ubnd ;
    cubi->revs = 0 ; }
  else
  { if (inflb == 1)
      clbi->inf = -inflndx ;
    else
      clbi->inf = inflb ;
    setcleanzero(lbnd,main_lptols->zero) ;
    clbi->bnd = lbnd ;
    clbi->revs = 0 ; }

# ifndef NDEBUG
  if (main_opts->print.conclass >= 5)
  { if (ctypi == contypLE)
    { if (infub == 1)
      { outfmt(logchn,gtxecho,"\n\tUB(%s) finite for %s (%d) only.",
	       consys_nme(consys,'c',i,FALSE,NULL),
	       consys_nme(consys,'v',-consys->cub[i].inf,FALSE,NULL),
	       -consys->cub[i].inf) ; } }
    else
    { if (inflb == 1)
      { outfmt(logchn,gtxecho,"\n\tLB(%s) finite for %s (%d) only.",
	       consys_nme(consys,'c',i,FALSE,NULL),
	       consys_nme(consys,'v',-consys->clb[i].inf,FALSE,NULL),
	       -consys->clb[i].inf) ; } } }
  if (main_opts->print.conclass >= 4)
  { if (main_opts->print.conclass >= 5)
      outfmt(logchn,gtxecho,"\n      final bound: ") ;
    if (consys->ctyp[i] == contypGE)
    { outfmt(logchn,gtxecho,"%s = %s",consys_conbndnme('L',i,&consys->clb[i]),
	     consys_conbndval(&consys->clb[i])) ; }
    else
    { outfmt(logchn,gtxecho,"%s = %s",
	     consys_conbndnme('U',i,&consys->cub[i]),
	     consys_conbndval(&consys->cub[i])) ; } }
# endif
  }

/*
  That's it --- all constraints have been processed.
*/
  if (pkrow != NULL) pkvec_free(pkrow) ;

  return (TRUE) ; }



bool arc_chgvarbnd (int t, char bnd, double ovbt, double nvbt, bool *p_feas)

/*
  This routine is responsible for propagating the change in the upper or lower
  bound of a variable. It also checks the change against the lp solution to
  see if it's likely to be significant.

  The basic criteria for a change to be significant is that it makes the
  current optimum infeasible. (If it doesn't do that, then the current optimum
  is still the optimum, eh?) In practical terms, this comes down to
    * A change to a bound of a nonbasic variable. (We're by definition
      infeasible, but whether it's serious requires recalculating the basic
      variables, to see if any have become infeasible due to the change. We'll
      run the lp code to check.)
    * A change to a bound which makes the current value of a basic variable
      infeasible. This'll definitely require at least one pivot to fix.
  
  To propagate the change, the routine retrieves the column for the variable,
  then steps through it. For each constraint involving the variable, the row
  type is recovered. Based on the row type and the sign of the variable, the
  code determines if the bound should be corrected and makes the correction.
  The bound needs to be corrected according to the following logic:

  * A >= constraint uses its upper bound in typGE_revise; the upper bound has
    changed if a<i,t> > 0 and ub(x<t>) has changed, or if a<i,t> < 0 and
    lb(x<t>) has changed.
  * A <= constraint uses its lower bound in typLE_revise; the lower bound has
    changed if a<i,t> > 0 and lb(x<t>) has changed, or if a<i,t> < 0 and
    ub(x<t>) has changed.
  * Equality and range constraints are viewed as a <=, >= pair.

  Propagation of a variable bound change is based on the percentage change
  it will have relative to max(|rhs-bound|,|bound|). The first term considers
  the effect of the change in use in an isolation, and the second attempts
  to deal with the problem of large but nearly equal rhs and bound values.
  When a change in variable bound results in the constraint bound changing from
  infinite to one of the B(i) or B(i\t) forms, propagation is forced with
  high priority.

  Whether or not we propagate a change, it is installed; the small changes
  are accumulated. Eventually the row bounds will be used, and all these small
  tightenings will be put to use.

  If the bound is actually corrected, the constraint is stacked for further
  propagation. Bounds are recalculated from scratch with a frequency determined
  by mipopts.arcrecalc.

  Parameters:
    t		the index of the target variable x<t>
    bnd		the bound which has changed; 'l' for lower, 'u' for upper
    ovbt	the old value of the bound
    nvbt	the new value of the bound
    p_feas	(o) set to FALSE if feasibility is lost, TRUE otherwise

  Returns: TRUE if the calculation proceeds without error, FALSE otherwise.
*/

{ double ait,delta,pct,tol,rhs,ainorm,valt ;
  conbnd_struct cbi,*ocbi ;
  int pkndx,i ;
  bool feas ;
  flags statt ;
  contyp_enum ctyp ;
  char cbnd ;
  pkvec_struct *pkcol ;
  char *rtnnme = "arc_chgvarbnd" ;

# ifndef NDEBUG
  cbndchg_struct *chg ;
# endif

# ifdef PARANOIA
  if (main_opts->phic == 0)
  { errmsg(467,rtnnme) ;
    return (FALSE) ; }
  if (t < 1 || t > phicsys->varcnt)
  { errmsg(102,rtnnme,phicsys->nme,"variable",t,1,phicsys->varcnt) ;
    return (FALSE) ; }
  if (bnd != 'l' && bnd != 'u')
  { errmsg(3,rtnnme,"bound type",bnd) ;
    return (FALSE) ; }
  if (bnd == 'l')
  { if (!abovebnd(nvbt,ovbt))
    { errmsg(458,rtnnme,phicsys->nme,"l<i>",
	     consys_nme(phicsys,'v',t,FALSE,NULL),t,nvbt,ovbt,ovbt-nvbt) ;
    return (FALSE) ; } }
  else
  { if (!belowbnd(nvbt,ovbt))
    { errmsg(458,rtnnme,phicsys->nme,"u<i>",
	     consys_nme(phicsys,'v',t,FALSE,NULL),t,nvbt,ovbt,nvbt-ovbt) ;
    return (FALSE) ; } }
  if (p_feas == NULL)
  { errmsg(2,rtnnme,"&feas") ;
    return (FALSE) ; }
# endif
# ifndef NDEBUG
  if (main_opts->print.arcprop >= 5)
  { outfmt(logchn,gtxecho,"\n\t    propagating %cb(%s) from %g to %g ...",
	   bnd,consys_nme(phicsys,'v',t,FALSE,NULL),ovbt,nvbt) ; }
#   endif
/*
  Check against the current solution. We want to err on the side of caution.
  Running the occasional unnecessary lp is better than missing a necessary one.
  For basic variables, it's a straightforward matter of checking if the bound
  change makes the current value infeasible.
  For nonbasic variables, if we're changing the bound that the variable is
  at, there's no real need to check the value (if the change wasn't nontrivial,
  we wouldn't be propagating it). It's remotely possible a nonbasic free or
  superbasic will show up here. We should never see NBFX, or anything else.
*/
  if (main_relax->chklpsol == TRUE && main_relax->needslp == FALSE)
  { statt = main_relax->lp->status[t] ;
    i = (int) statt ;
    if (i < 0)
    { i = -i ;
      valt = main_relax->lp->x[i] ;
      if ((bnd == 'l' && nvbt > valt) || (bnd == 'u' && nvbt < valt))
      {
#	ifndef NDEBUG
	if (main_opts->print.arcprop >= 5)
	{ outfmt(logchn,gtxecho," (nontrivial basic)") ; }
#	endif
	main_relax->needslp = TRUE ; } }
    else
    { switch (statt)
      { case vstatNBLB:
	{ if (bnd == 'l')
	  {
#	    ifndef NDEBUG
	    if (main_opts->print.arcprop >= 5)
	    { outfmt(logchn,gtxecho," (nontrivial %s)",dy_prtvstat(statt)) ; }
#	    endif
	    main_relax->needslp = TRUE ; }
	  break ; }
	case vstatNBUB:
	{ if (bnd == 'u')
	  {
#	    ifndef NDEBUG
	    if (main_opts->print.arcprop >= 5)
	    { outfmt(logchn,gtxecho," (nontrivial %s)",dy_prtvstat(statt)) ; }
#	    endif
	    main_relax->needslp = TRUE ; }
	  break ; }
	case vstatNBFR:
	case vstatSB:
	{ 
#	  ifndef NDEBUG
	  if (main_opts->print.arcprop >= 5)
	  { outfmt(logchn,gtxecho," (nontrivial %s)",dy_prtvstat(statt)) ; }
#	  endif
	  main_relax->needslp = TRUE ;
	  break ; }
	default:
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (FALSE) ; } } } } 
/*
  Acquire the packed column for the variable.
*/
  pkcol = NULL ;
  if (consys_getcol_pk(phicsys,t,&pkcol) == FALSE)
  { errmsg(122,rtnnme,phicsys->nme,
	   "column",consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
    if (pkcol != NULL) pkvec_free(pkcol) ;
    return (FALSE) ; }
/*
  Fire up a loop to step down the packed column. 
*/
  feas = TRUE ;
  for (pkndx = 0 ; pkndx < pkcol->cnt && feas == TRUE ; pkndx++)
  { i = pkcol->coeffs[pkndx].ndx ;
#   ifdef PARANOIA
    if (i < 1 || i > phicsys->concnt)
    { errmsg(102,rtnnme,phicsys->nme,"constraint",i,1,phicsys->concnt) ;
      return (FALSE) ; }
#   endif
    ait = pkcol->coeffs[pkndx].val ;
#   ifndef NDEBUG
    if (fabs(ait) < main_lptols->zero && flgon(phicsys->opts,CONSYS_WRNZERO))
      warn(450,rtnnme,i,t,consys_nme(phicsys,'c',i,FALSE,NULL),
	   consys_nme(phicsys,'v',t,FALSE,NULL)) ;
#   endif
    if (fabs(ait) < main_lptols->zero) continue ;
    ctyp = phicsys->ctyp[i] ;
#   ifndef NDEBUG
    if (main_opts->print.arcprop >= 6)
    { outfmt(logchn,gtxecho,"\n\t    %s constraint %s: ",
	     consys_prtcontyp(ctyp),consys_nme(phicsys,'c',i,FALSE,NULL)) ; }
#   endif
/*
  The first question to answer is whether the change in x<t>'s bound will
  affect the relevant bound of the constraint. For contypGE we're interested
  only in UB(i); for contypLE we're interested only in LB(i). For contypEQ
  and contypRNG, we're interested in both bounds, but only one will change in
  response to a change in one of a variable's bounds. Hence we'll execute at
  most one of the cases below.

  If the change in x<t>'s bound is from infinite to finite, we decrease the
  inf count and put ait*nvbt into the bnd portion (bnd did not contain any
  contribution from x<t> prior to the change to a finite bound). There are
  special cases to handle the change to 1 and 0 infinite bounds (the first
  requires we hunt down the remaining infinite variable with calc_rowbnd, the
  second is just a matter of encoding).

  If the change in x<t>'s bound is from one finite value to another, then we
  correct bnd by -ait*ovbt+ait*nvbt = ait(nvbt-ovbt).

  All of this is irrelevant if the revs count says it's time to recalculate
  the bound from scratch.
*/
    if (ctyp != contypGE &&
	 ((ait > 0 && bnd == 'l') || (ait < 0 && bnd == 'u')))
    { ocbi = &phicsys->clb[i] ;
      cbnd = 'L' ; }
    else
    if (ctyp != contypLE &&
	 ((ait > 0 && bnd == 'u') || (ait < 0 && bnd == 'l')))
    { ocbi = &phicsys->cub[i] ;
      cbnd = 'U' ; }
    else
    { 
#     ifndef NDEBUG
      if (main_opts->print.arcprop >= 6)
      { outfmt(logchn,gtxecho,"n/a.") ; }
#     endif
      continue ; }

    if ((fabs(ovbt) < main_lptols->inf || ocbi->inf != 2) &&
	ocbi->revs < main_opts->arcrecalc)
    { if (fabs(ovbt) < main_lptols->inf)
      { delta = ait*(nvbt-ovbt) ;
	cbi.inf = ocbi->inf ; }
      else
      { if (ocbi->inf > 2)
	{ cbi.inf = ocbi->inf-1 ; }
	else
	{ cbi.inf = 0 ; }
	delta = ait*nvbt ; }
      cbi.bnd = ocbi->bnd+delta ;
      setcleanzero(cbi.bnd,main_lptols->zero) ;
      cbi.revs = ocbi->revs+1 ; }
    else
    { if (recalc_rowbnd(i,cbnd,&cbi) == FALSE)
      { if (pkcol != NULL) pkvec_free(pkcol) ;
	return (FALSE) ; } }
/*
  The bound is considered to have tightened if either of the infinite or
  finite contributions has tightened.  The change needs to be propagated if
  the infinity contribution has dropped to 1 or 0, or if it has remained at 1
  or 0 and the finite contribution has tightened.
*/
    if (ocbi->inf != cbi.inf || !atbnd(ocbi->bnd,cbi.bnd))
    { 
#     ifndef NDEBUG
      if (main_opts->print.arcprop >= 4)
      { chg = &conbnd_chgs[i] ;
	if (chg->empty == TRUE)
	{ chg->olb = phicsys->clb[i] ;
	  chg->nlb = phicsys->clb[i] ;
	    chg->oub = phicsys->cub[i] ;
	    chg->nub = phicsys->cub[i] ;
	  chg->empty = FALSE ;
	  chg->cnt = 0 ; }
	chg->cnt++ ;
	if (cbnd == 'L')
	  chg->nlb = cbi ;
	else
	  chg->nub = cbi ;
	if (main_opts->print.arcprop >= 6)
	{ if (ocbi->inf != cbi.inf)
	    delta = main_lptols->inf ;
	  else
	    delta = fabs(ocbi->bnd-cbi.bnd) ;
	  outfmt(logchn,gtxecho,"%s = %s tightened by %g to ",
		 consys_conbndnme(cbnd,i,ocbi),consys_conbndval(ocbi),delta) ;
	  outfmt(logchn,gtxecho,"%s = %s.",
		 consys_conbndnme(cbnd,i,&cbi),consys_conbndval(&cbi)) ; } }
#     endif
      if (cbi.inf > 1)
      { *ocbi = cbi ;
	continue ; }
/*
  Check for loss of feasibility, if the bound is strictly finite. The first
  level tests with the values immediately available. If it looks like we've
  lost feasibility, recalculate the bound from scratch, retrieve the norm
  of the coefficients, and do a proper scaled test. If we lose feasibility,
  there's no need to proceed further. Clear off the stack and we're out of
  here.

  If we pass the feasibility test, check for the case where cbi.bnd is right
  at the rhs value. If it is, set it to the rhs value. This removes any error
  below the zero tolerance threshold.
*/
      if (ctyp == contypRNG && cbnd == 'U')
	rhs = phicsys->rhslow[i] ;
      else
	rhs = phicsys->rhs[i] ;
      ainorm = consys_1normrow(phicsys,i) ;
      if (cbi.inf == 0)
      { if ((cbnd == 'L' && abovebnd(cbi.bnd,rhs)) ||
	    (cbnd == 'U' && belowbnd(cbi.bnd,rhs)))
	{ if (recalc_rowbnd(i,cbnd,&cbi) == FALSE)
	  { if (pkcol != NULL) pkvec_free(pkcol) ;
	    return (FALSE) ; }
	  tol = main_lptols->zero*(1+ainorm) ;
	  if ((cbnd == 'L' && cbi.bnd > rhs+tol) ||
	      (cbnd == 'U' && cbi.bnd < rhs-tol))
	  { feas = FALSE ;
	    arc_clear() ;
#	    ifndef NDEBUG
	    if (main_opts->print.arcprop >= 1)
	    { outfmt(logchn,gtxecho,
		     "\n\t    feasibility lost while propagating %cb(%s).",
		     bnd,consys_nme(phicsys,'v',t,FALSE,NULL)) ;
	      outfmt(logchn,gtxecho,"\n\t    %s constraint %s: %s = %s %c ",
		     consys_prtcontyp(ctyp),
		     consys_nme(phicsys,'c',i,FALSE,NULL),
		     consys_conbndval(&cbi),consys_conbndnme(cbnd,i,&cbi),
		     (cbnd == 'L')?'>':'<') ;
	      outfmt(logchn,gtxecho,"b%s<%d> = %g by %g, tol = %g.",
		     (ctyp == contypRNG && cbnd == 'U')?"low":"",i,rhs,
		     fabs(rhs-cbi.bnd),tol) ; }
#	    endif
	    *ocbi = cbi ;
	    break ; } }
	  else
	  if (atbnd(cbi.bnd,rhs)) cbi.bnd = rhs ; }
      if (ocbi->inf != cbi.inf)
      { delta = main_lptols->inf ;
	pct = main_lptols->inf ; }
      else
      { delta = fabs(ocbi->bnd-cbi.bnd) ;
	pct = delta/ainorm ; }
      *ocbi = cbi ; }
    else
    { 
#     ifndef NDEBUG
      if (main_opts->print.arcprop >= 6)
      { outfmt(logchn,gtxecho,"no change.") ; }
#     endif
      continue ; }
/*
  If we're here, we're considering propagating the bound.
  Stack the constraint for further propagation, if the change is
  big enough to justify it.
*/
    if (pct >= main_tols->arcprop)
    { 
#     ifndef NDEBUG
      if (main_opts->print.arcprop >= 6)
      { outfmt(logchn,gtxecho," pushed at %g/%g = %g%%.",delta,ainorm,pct) ; }
#     endif
      if (pushcon(i,delta,pct) == FALSE)
      { errmsg(452,rtnnme,phicsys->nme,consys_nme(phicsys,'c',i,FALSE,NULL),i) ;
	if (pkcol != NULL) pkvec_free(pkcol) ;
	return (FALSE) ; } }
#   ifndef NDEBUG
    else
    { if (main_opts->print.arcprop >= 4)
	warn(454,rtnnme,phicsys->nme,pct*100,"constraint",
	     consys_nme(phicsys,'c',i,FALSE,NULL),i,main_tols->arcprop) ; }
#   endif
  }
/*
  That's the end of the column. Return success.
*/
  if (pkcol != NULL) pkvec_free(pkcol) ;
  *p_feas = feas ;

  return (TRUE) ; }



bool arc_chgrhs (int i, char cbnd, double orhs, double nrhs, bool recalc_bnds)

/*
  This routine propagates a change in b<i> or blow<i> by pushing the
  constraint onto the propagation stack. If recalc_bnds is true, LB(i) and
  UB(i) are recalculated (this is necessary when converting a <= or >=
  constraint to a <=> or = constraint).

  Note that we might be changing either or rhs or rhslow, depending on the
  bound and constraint type. The caller deals with this; here the task is
  more straightforward. Note also that we're not actually changing LB(i)
  or UB(i) --- these are calculated based on upper and lower bounds on
  the variables in the lhs of the constraint. The bounds on variables
  will change because the value of the isolation will change when b<i> or
  blow<i> changes. See the notes at the top of the file for the isolation
  expressions.

  Parameters:
    i:		index of the constraint being altered
    cbnd:	'L' if we're changing the lower bound on LHS(i), 'U' if we're
		changing the upper bound.
    orhs:	old rhs value
    nrhs:	new rhs value
    recalc_bnds: TRUE to force recalculation of LB(i) and UB(i).

  Returns: TRUE if the propagation is completed without error, FALSE otherwise.
*/

{ double deltai,ainorm,pct ;
  conbnd_struct *ocbi ;
  char lhsbnd ;

  char *rtnnme = "arc_chgrhs" ;

# ifdef PARANOIA
  if (main_opts->phic == 0)
  { errmsg(467,rtnnme) ;
    return (FALSE) ; }
  if (i < 1 || i > phicsys->concnt)
  { errmsg(102,rtnnme,phicsys->nme,"constraint",i,1,phicsys->varcnt) ;
    return (FALSE) ; }
  if (cbnd != 'L' && cbnd != 'U')
  { errmsg(3,rtnnme,"bound type",cbnd) ;
    return (FALSE) ; }
  if (cbnd == 'L')
  { if (!abovebnd(nrhs,orhs))
    { errmsg(458,rtnnme,phicsys->nme,"blow<i>",
	     consys_nme(phicsys,'c',i,FALSE,NULL),i,nrhs,orhs,orhs-nrhs) ;
    return (FALSE) ; } }
  else
  { if (!belowbnd(nrhs,orhs))
    { errmsg(458,rtnnme,phicsys->nme,"b<i>",
	     consys_nme(phicsys,'c',i,FALSE,NULL),i,nrhs,orhs,nrhs-orhs) ;
    return (FALSE) ; } }
# endif
# ifndef NDEBUG
  if (main_opts->print.arcprop >= 4)
  { outfmt(logchn,gtxecho,"\n\t    propagating %s<%s> from %g to %g ...",
	   (cbnd == 'L')?"blow":"b",consys_nme(phicsys,'c',i,FALSE,NULL),
	   orhs,nrhs) ; }
#   endif
/*
  If we're supposed to recalculate LB<i> and UB<i>, do it now.
*/
  if (recalc_bnds == TRUE)
  { if (calc_rowbnd(i,&phicsys->clb[i],&phicsys->cub[i],NULL,NULL) == FALSE)
    { errmsg(451,rtnnme,consys_nme(phicsys,'c',i,FALSE,NULL),i) ;
      return (FALSE) ; } }
/*
  Haul out the appropriate lhs bound (LB<i> for b<i>, UB<i> for blow<i>) and
  check that it's finite. We can only propagate changes to a finite bound.
*/
# ifndef NDEBUG
  if (main_opts->print.arcprop >= 6)
  { outfmt(logchn,gtxecho,"\n\t    %s constraint %s: ",
	   consys_prtcontyp(phicsys->ctyp[i]),
	   consys_nme(phicsys,'c',i,FALSE,NULL)) ; }
# endif
  if (cbnd == 'U')
  { ocbi = &phicsys->clb[i] ;
    lhsbnd = 'L' ; }
  { ocbi = &phicsys->cub[i] ;
    lhsbnd = 'U' ; }
  if (ocbi->inf > 0)
  {
#   ifndef NDEBUG
    if (main_opts->print.arcprop >= 6)
    { outfmt(logchn,gtxecho,"%s = %s; not propagated.",
	     consys_conbndnme(lhsbnd,i,ocbi),consys_conbndval(ocbi)) ; }
#   endif
    return (TRUE) ; }
/*
  Calculate delta<i> and 1-norm(a<i>).
*/
  deltai = fabs(nrhs-orhs) ;
  ainorm = consys_1normrow(phicsys,i) ;
  pct = deltai/ainorm ;
  setcleanzero(pct,main_lptols->zero) ;
/*
  If we're here, we're considering propagating the bound.
  Stack the constraint for further propagation, if the change is
  big enough to justify it.
*/
  if (pct >= main_tols->arcprop)
  { 
#   ifndef NDEBUG
    if (main_opts->print.arcprop >= 6)
    { outfmt(logchn,gtxecho,"pushed at %g/%g = %g%%.",deltai,ainorm,pct) ; }
#   endif
    if (pushcon(i,deltai,pct) == FALSE)
    { errmsg(452,rtnnme,phicsys->nme,
	     consys_nme(phicsys,'c',i,FALSE,NULL),i) ;
      return (FALSE) ; } }
# ifndef NDEBUG
  else
  { if (main_opts->print.arcprop >= 4)
      warn(454,rtnnme,phicsys->nme,pct*100,"constraint",
	   consys_nme(phicsys,'c',i,FALSE,NULL),i,main_tols->arcprop) ; }
# endif

  return (TRUE) ; }



static bool typLE_revise (int i, int *fxcnt, bool *p_feas)

/*
  This routine handles less than or equal inequalities. We're here under one of
  two circumstances: LB(i) is finite, and all variables in the constraint
  are to be evaluated, or only LB(i\t) is finite, hence only x<t> is to be
  evaluated.

  The routine retrieves the constraint coefficients and steps through them
  evaluating the bound on each variable x<t>. When a bound is tightened, the
  variable bounds are written and the change is propagated to the bounds of
  other constraints which use the variable. When a variable is fixed, it is
  recorded in fxvars (static global) and fxcnt is incremented.

  If feasibility is lost, feas is set to FALSE. This can happen if the domain
  of one of the variables in constraint i becomes empty, or if arc_chgvarbnd
  detects that some constraint has become infeasible while it's propagating
  a change in the bound of a variable.

  See the opening comments of the file for a longer explanation of the math.
  Basically, we're processing constraint i, calculating isolations of variable
  x<t>, where t ranges over all variables in the constraint.

  Parameters:
    i		the constraint (row) index
    fxcnt	(i/o) the number of variables fixed during this run of the
		arc consistency algorithm
    p_feas	(o) set to FALSE if feasibility is lost during bounds
		revision, TRUE otherwise

  Returns: TRUE if the individual variables of the constraint are successfully
	   processed, FALSE otherwise.
*/


{ double vlbt,vubt,bi,ait,vlbt_n,vubt_n,rngt,tol ;
  int pkndx,t ;
  conbnd_struct *clbi ;
  pkvec_struct *pkrow ;
  bool feas ;
  char *rtnnme = "typLE_revise" ;

# ifndef NDEBUG
  vbndchg_struct *chg ;
# endif

# ifdef PARANOIA
  if (i < 1 || i > phicsys->concnt)
  { errmsg(102,rtnnme,phicsys->nme,"constraint",i,1,phicsys->concnt) ;
    return (FALSE) ; }
  if (fxcnt == NULL)
  { errmsg(2,rtnnme,"&fxcnt") ;
    return (FALSE) ; }
  if (p_feas == NULL)
  { errmsg(2,rtnnme,"&feas") ;
    return (FALSE) ; }
  if (phicsys->clb[i].inf > 0)
  { errmsg(457,rtnnme,phicsys->nme,consys_nme(phicsys,'c',i,FALSE,NULL),
	   consys_conbndnme('L',i,&phicsys->clb[i]),
	   consys_conbndval(&phicsys->clb[i])) ;
    return (FALSE) ; }
  if (phicsys->clb[i].inf < 0 && -phicsys->clb[i].inf > phicsys->varcnt)
  { errmsg(102,rtnnme,phicsys->nme,"LB(i\t) tgt ndx",-phicsys->clb[i].inf,
	   1,phicsys->concnt) ;
    return (FALSE) ; }
# endif

/*
  Fetch the coefficients of the constraint a<i>, the row lower bound LB(i),
  and the right-hand side value b<i>.
*/
  pkrow = NULL ;
  if (consys_getrow_pk(phicsys,i,&pkrow) == FALSE)
  { errmsg(122,rtnnme,phicsys->nme,
	   "constraint",consys_nme(phicsys,'c',i,FALSE,NULL),i) ;
    if (pkrow != NULL) pkvec_free(pkrow) ;
    return (FALSE) ; }
  clbi = &phicsys->clb[i] ;
  bi = phicsys->rhs[i] ;
/*
  Step through the row.  If we're here with LB(i\t), we're looking
  specifically for a<i,t>, and nothing else. If we're here with LB(i), we're
  looking at each variable x<t>.

  We really shouldn't see a<i,t> == 0, but it's not precisely wrong. Warn the
  user, but don't abort.
*/
  feas = TRUE ;
  for (pkndx = 0 ; pkndx < pkrow->cnt && feas == TRUE ; pkndx++)
  { t = pkrow->coeffs[pkndx].ndx ;
    if (clbi->inf < 0 && t != -clbi->inf) continue ;
    ait = pkrow->coeffs[pkndx].val ;
#   ifndef NDEBUG
    if (fabs(ait) < main_lptols->zero && flgon(phicsys->opts,CONSYS_WRNZERO))
      warn(450,rtnnme,i,t,consys_nme(phicsys,'c',i,FALSE,NULL),
	   consys_nme(phicsys,'v',t,FALSE,NULL)) ;
#   endif
    if (fabs(ait) < main_lptols->zero) continue ;
    vlbt = phicsys->vlb[t] ;
    vubt = phicsys->vub[t] ;
    rngt = vubt-vlbt ;
/*
  Do the calculation to see if x<t>'s upper bound can be tightened. If we're
  working with LB(i\t), clbi->bnd is correct as is. If we're dealing with
  LB(i), we have to correct for x<t>'s contribution. Then round the new value
  according to the rounding tolerance and integrality.

  If there's a change, we have several possibilities:
    * vubt_n >= vubt. Nice try, but no cigar.
    * vubt_n is a bit tighter. The domain has shrunk, but that's all.
    * vubt_n = vlbt (within tolerance). We have a fixed variable.
    * vubt_n < vlbt. x<t>'s domain is empty and we have infeasibility.
  In the middle two cases, we'll want to propagate the bound change. The
  first case requires no action; continue to the next variable. The final
  case means we're done.
*/
    if (ait > 0)
    { if (clbi->inf < 0)
	vubt_n = (bi-clbi->bnd)/ait ;
      else
	vubt_n = (bi-(clbi->bnd-ait*vlbt))/ait ;
      setcleanzero(vubt_n,main_lptols->zero) ;
      vubt_n = ROUND_UB(vubt_n,main_tols->arcsnap) ;
      if (INT_VARTYPE(phicsys->vtyp[t])) vubt_n = floor(vubt_n) ;

      if (belowbnd(vubt_n,vubt))
      { if (rngt < main_lptols->inf)
	  tol = main_lptols->zero*(1+maxx(fabs(vlbt),rngt)) ;
	else
	  tol = main_lptols->zero*(1+fabs(vlbt)) ;
	if (vlbt-vubt_n > tol)
	{ feas = FALSE ;
#	  ifndef NDEBUG
	  if (main_opts->print.arcprop >= 1)
	  { outfmt(logchn,gtxecho,"\n\tempty domain for %s (%d); ",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    outfmt(logchn,gtxecho,"%g = lb > ub = %g, lb-ub = %g.",
		   vlbt,vubt_n,vlbt-vubt_n) ; }
#	  endif
	}
	else
	{ if (withintol(vlbt,vubt_n,tol))
	  { vubt_n = vlbt ;
	    if ((*fxcnt)+1 > fxvars_sze)
	    { fxvars_sze = maxx(main_sys->varcnt,.1*(fxvars_sze+10)) ;
	      fxvars = (fxvar_struct *)
		       REALLOC(fxvars,fxvars_sze*sizeof(fxvar_struct)) ; }
	    fxvars[*fxcnt].ndx = t ;
	    fxvars[(*fxcnt)++].val = vubt_n ;
#	    ifndef NDEBUG
	    if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\t%s (%d) fixed at %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),t,vlbt) ; }
#	    endif
	  }
#	  ifndef NDEBUG
	  else
	  { if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\tub(%s) tightened by %g from %g to %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),vubt-vubt_n,
		     vubt,vubt_n) ; } }
	  if (main_opts->print.arcprop >= 1)
	  { chg = &varbnd_chgs[t] ;
	    if (chg->empty == TRUE)
	    { chg->olb = vlbt ;
	      chg->oub = vubt ;
	      chg->nlb = vlbt ;
	      chg->nub = vubt ;
	      chg->empty = FALSE ;
	      chg->cnt = 0 ; }
	    chg->cnt++ ;
	    chg->nub = vubt_n ; }
#         endif
	  phicsys->vub[t] = vubt_n ;
	  setflg(main_lp->ctlopts,lpctlUBNDCHG) ;
	  if (arc_chgvarbnd(t,'u',vubt,vubt_n,&feas) == FALSE)
	  { errmsg(459,rtnnme,phicsys->nme,"ub<j>",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    if (pkrow != NULL) pkvec_free(pkrow) ;
	    return (FALSE) ; } } } }
/*
  The same flow, but in this case we're working on x<t>'s lower bound.
*/
    else
    { if (clbi->inf < 0)
	vlbt_n = (bi-clbi->bnd)/ait ;
      else
	vlbt_n = (bi-(clbi->bnd-ait*vubt))/ait ;
      setcleanzero(vlbt_n,main_lptols->zero) ;
      vlbt_n = ROUND_LB(vlbt_n,main_tols->arcsnap) ;
      if (INT_VARTYPE(phicsys->vtyp[t])) vlbt_n = ceil(vlbt_n) ;

      if (abovebnd(vlbt_n,vlbt))
      { if (rngt < main_lptols->inf)
	  tol = main_lptols->zero*(1+maxx(fabs(vubt),rngt)) ;
	else
	  tol = main_lptols->zero*(1+fabs(vubt)) ;
	if (vlbt_n-vubt > tol)
	{ feas = FALSE ;
#	  ifndef NDEBUG
	  if (main_opts->print.arcprop >= 1)
	  { outfmt(logchn,gtxecho,"\n\tempty domain for %s (%d); ",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    outfmt(logchn,gtxecho,"%g = lb > ub = %g, lb-ub = %g.",
		   vlbt_n,vubt,vlbt_n-vubt) ; }
#	  endif
	}
	else
	{ if (withintol(vlbt_n,vubt,tol))
	  { vlbt_n = vubt ;
	    if ((*fxcnt)+1 > fxvars_sze)
	    { fxvars_sze = maxx(main_sys->varcnt,.1*(fxvars_sze+10)) ;
	      fxvars = (fxvar_struct *)
		       REALLOC(fxvars,fxvars_sze*sizeof(fxvar_struct)) ; }
	    fxvars[*fxcnt].ndx = t ;
	    fxvars[(*fxcnt)++].val = vlbt_n ;
#	    ifndef NDEBUG
	    if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\t%s (%d) fixed at %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),t,vubt) ; }
#	    endif
	  }
#	  ifndef NDEBUG
	  else
	  { if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\tlb(%s) tightened by %g from %g to %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),vlbt_n-vlbt,
		     vlbt,vlbt_n) ; } }
	  if (main_opts->print.arcprop >= 1)
	  { chg = &varbnd_chgs[t] ;
	    if (chg->empty == TRUE)
	    { chg->olb = vlbt ;
	      chg->oub = vubt ;
	      chg->nlb = vlbt ;
	      chg->nub = vubt ;
	      chg->empty = FALSE ;
	      chg->cnt = 0 ; }
	    chg->cnt++ ;
	    chg->nlb = vlbt_n ; }
#         endif
	  phicsys->vlb[t] = vlbt_n ;
	  setflg(main_lp->ctlopts,lpctlLBNDCHG) ;
	  if (arc_chgvarbnd(t,'l',vlbt,vlbt_n,&feas) == FALSE)
	  { errmsg(459,rtnnme,phicsys->nme,"lb<j>",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    if (pkrow != NULL) pkvec_free(pkrow) ;
	    return (FALSE) ; } } } } }
/*
  We're done. Return success.
*/
  if (pkrow != NULL) pkvec_free(pkrow) ;
  *p_feas = feas ;

  return (TRUE) ; }



static bool typGE_revise (int i, int *fxcnt, bool *p_feas)

/*
  This routine handles greater than or equal inequalities. The flow of the
  algorithm is identical to typLE_revise, but the number of detail changes
  in the math makes it preferable to duplicate the routine. The code is much
  more readable this way.

  Parameters:
    i		the constraint (row) index
    fxcnt	(i/o) the number of variables fixed during this run of the
		arc consistency algorithm
    p_feas	(o) set to FALSE if feasibility is lost during bounds
		revision, TRUE otherwise

  Returns: TRUE if the individual variables of the constraint are successfully
	   processed, FALSE otherwise.
*/


{ double vlbt,vubt,bi,ait,vlbt_n,vubt_n,rngt,tol ;
  int pkndx,t ;
  conbnd_struct *cubi ;
  pkvec_struct *pkrow ;
  bool feas ;
  char *rtnnme = "typGE_revise" ;

# ifndef NDEBUG
  vbndchg_struct *chg ;
# endif

# ifdef PARANOIA
  if (i < 1 || i > phicsys->concnt)
  { errmsg(102,rtnnme,phicsys->nme,"constraint",i,1,phicsys->concnt) ;
    return (FALSE) ; }
  if (fxcnt == NULL)
  { errmsg(2,rtnnme,"&fxcnt") ;
    return (FALSE) ; }
  if (p_feas == NULL)
  { errmsg(2,rtnnme,"&feas") ;
    return (FALSE) ; }
  if (phicsys->cub[i].inf > 0)
  { errmsg(457,rtnnme,phicsys->nme,consys_nme(phicsys,'c',i,FALSE,NULL),
	   consys_conbndnme('U',i,&phicsys->cub[i]),
	   consys_conbndval(&phicsys->cub[i])) ;
    return (FALSE) ; }
  if (phicsys->cub[i].inf < 0 && -phicsys->cub[i].inf > phicsys->varcnt)
  { errmsg(102,rtnnme,phicsys->nme,"LB(i\t) tgt ndx",-phicsys->cub[i].inf,
	   1,phicsys->concnt) ;
    return (FALSE) ; }
# endif

/*
  Fetch the coefficients of the constraint a<i>, the row upper bound UB(i),
  and the right-hand-side b<i>.  One quirk here --- if we're processing a
  range constraint, we need to acquire blow<i> instead of b<i>.
*/
  pkrow = NULL ;
  if (consys_getrow_pk(phicsys,i,&pkrow) == FALSE)
  { errmsg(122,rtnnme,phicsys->nme,
	   "constraint",consys_nme(phicsys,'c',i,FALSE,NULL),i) ;
    if (pkrow != NULL) pkvec_free(pkrow) ;
    return (FALSE) ; }
  if (phicsys->ctyp[i] == contypRNG)
    bi = phicsys->rhslow[i] ;
  else
    bi = phicsys->rhs[i] ;
  cubi = &phicsys->cub[i] ;
/*
  Step through the row.  If we're here with UB(i\t), we're looking
  specifically for a<i,t>, and nothing else. If we're here with UB(i), we're
  looking at each variable x<t>.

  We really shouldn't see a<i,t> == 0, but it's not precisely wrong. Warn the
  user, but don't abort.
*/
  feas = TRUE ;
  for (pkndx = 0 ; pkndx < pkrow->cnt && feas == TRUE ; pkndx++)
  { t = pkrow->coeffs[pkndx].ndx ;
    if (cubi->inf < 0 && t != -cubi->inf) continue ;
    ait = pkrow->coeffs[pkndx].val ;
#   ifndef NDEBUG
    if (fabs(ait) < main_lptols->zero && flgon(phicsys->opts,CONSYS_WRNZERO))
      warn(450,rtnnme,i,t,consys_nme(phicsys,'c',i,FALSE,NULL),
	   consys_nme(phicsys,'v',t,FALSE,NULL)) ;
#   endif
    if (fabs(ait) < main_lptols->zero) continue ;
    vlbt = phicsys->vlb[t] ;
    vubt = phicsys->vub[t] ;
    rngt = vubt-vlbt ;
/*
  Do the calculation to see if x<t>'s lower bound can be tightened. If we're
  working with UB(i\t), cubi->bnd is correct as is. If we're dealing with
  UB(i), we have to correct for x<t>'s contribution. Then round the new value
  according to the rounding tolerance and integrality.

  If there's a change, we have several possibilities:
    * vlbt_n > vubt, and x<t>'s domain is empty. We have infeasibility.
    * vlbt_n = vubt (within tolerance). We have a fixed variable.
    * vlbt_n is a bit tighter, but that's all.
  In the latter two cases, we'll want to propagate the bound change.
*/
    if (ait > 0)
    { if (cubi->inf < 0)
	vlbt_n = (bi-cubi->bnd)/ait ;
      else
	vlbt_n = (bi-(cubi->bnd-ait*vubt))/ait ;
      setcleanzero(vlbt_n,main_lptols->zero) ;
      vlbt_n = ROUND_LB(vlbt_n,main_tols->arcsnap) ;
      if (INT_VARTYPE(phicsys->vtyp[t])) vlbt_n = ceil(vlbt_n) ;

      if (abovebnd(vlbt_n,vlbt))
      { if (rngt < main_lptols->inf)
	  tol = main_lptols->zero*(1+maxx(fabs(vubt),rngt)) ;
	else
	  tol = main_lptols->zero*(1+fabs(vubt)) ;
	if (vlbt_n-vubt > tol)
	{ feas = FALSE ;
#	  ifndef NDEBUG
	  if (main_opts->print.arcprop >= 1)
	  { outfmt(logchn,gtxecho,"\n\tempty domain for %s (%d); ",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    outfmt(logchn,gtxecho,"%g = lb > ub = %g, lb-ub = %g.",
		   vlbt_n,vubt,vlbt_n-vubt) ; }
#	  endif
	}
	else
	{ if (withintol(vlbt_n,vubt,tol))
	  { vlbt_n = vubt ;
	    if ((*fxcnt)+1 > fxvars_sze)
	    { fxvars_sze = maxx(main_sys->varcnt,.1*(fxvars_sze+10)) ;
	      fxvars = (fxvar_struct *)
		       REALLOC(fxvars,fxvars_sze*sizeof(fxvar_struct)) ; }
	    fxvars[*fxcnt].ndx = t ;
	    fxvars[(*fxcnt)++].val = vlbt_n ;
#	    ifndef NDEBUG
	    if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\t%s (%d) fixed at %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),t,vubt) ; }
#	    endif
	  }
#	  ifndef NDEBUG
	  else
	  { if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\tlb(%s) tightened by %g from %g to %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),vlbt_n-vlbt,
		     vlbt,vlbt_n) ; } }
	  if (main_opts->print.arcprop >= 1)
	  { chg = &varbnd_chgs[t] ;
	    if (chg->empty == TRUE)
	    { chg->olb = vlbt ;
	      chg->oub = vubt ;
	      chg->nlb = vlbt ;
	      chg->nub = vubt ;
	      chg->empty = FALSE ;
	      chg->cnt = 0 ; }
	    chg->cnt++ ;
	    chg->nlb = vlbt_n ; }
#         endif
	  phicsys->vlb[t] = vlbt_n ;
	  setflg(main_lp->ctlopts,lpctlLBNDCHG) ;
	  if (arc_chgvarbnd(t,'l',vlbt,vlbt_n,&feas) == FALSE)
	  { errmsg(459,rtnnme,phicsys->nme,"lb<j>",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    if (pkrow != NULL) pkvec_free(pkrow) ;
	    return (FALSE) ; } } } }
    else
    { if (cubi->inf < 0)
	vubt_n = (bi-cubi->bnd)/ait ;
      else
	vubt_n = (bi-(cubi->bnd-ait*vlbt))/ait ;
      setcleanzero(vubt_n,main_lptols->zero) ;
      vubt_n = ROUND_UB(vubt_n,main_tols->arcsnap) ;
      if (INT_VARTYPE(phicsys->vtyp[t])) vubt_n = floor(vubt_n) ;

      if (belowbnd(vubt_n,vubt))
      { if (rngt < main_lptols->inf)
	  tol = main_lptols->zero*(1+maxx(fabs(vlbt),rngt)) ;
	else
	  tol = main_lptols->zero*(1+fabs(vlbt)) ;
	if (vlbt-vubt_n > tol)
	{ feas = FALSE ;
#	  ifndef NDEBUG
	  if (main_opts->print.arcprop >= 1)
	  { outfmt(logchn,gtxecho,"\n\tempty domain for %s (%d); ",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    outfmt(logchn,gtxecho,"%g = lb > ub = %g, lb-ub = %g.",
		   vlbt,vubt_n,vlbt-vubt_n) ; }
#	  endif
	}
	else
	{ if (withintol(vlbt,vubt_n,tol))
	  { vubt_n = vlbt ;
	    if ((*fxcnt)+1 > fxvars_sze)
	    { fxvars_sze = maxx(main_sys->varcnt,.1*(fxvars_sze+10)) ;
	      fxvars = (fxvar_struct *)
		       REALLOC(fxvars,fxvars_sze*sizeof(fxvar_struct)) ; }
	    fxvars[*fxcnt].ndx = t ;
	    fxvars[(*fxcnt)++].val = vubt_n ;
#	    ifndef NDEBUG
	    if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\t%s (%d) fixed at %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),t,vlbt) ; }
#	    endif
	  }
#	  ifndef NDEBUG
	  else
	  { if (main_opts->print.arcprop >= 4)
	    { outfmt(logchn,gtxecho,
		     "\n\tub(%s) tightened by %g from %g to %g.",
		     consys_nme(phicsys,'v',t,FALSE,NULL),vubt-vubt_n,
		     vubt,vubt_n) ; } }
	  if (main_opts->print.arcprop >= 1)
	  { chg = &varbnd_chgs[t] ;
	    if (chg->empty == TRUE)
	    { chg->olb = vlbt ;
	      chg->oub = vubt ;
	      chg->nlb = vlbt ;
	      chg->nub = vubt ;
	      chg->empty = FALSE ;
	      chg->cnt = 0 ; }
	    chg->cnt++ ;
	    chg->nub = vubt_n ; }
#         endif
	  phicsys->vub[t] = vubt_n ;
	  setflg(main_lp->ctlopts,lpctlUBNDCHG) ;
	  if (arc_chgvarbnd(t,'u',vubt,vubt_n,&feas) == FALSE)
	  { errmsg(459,rtnnme,phicsys->nme,"ub<j>",
		   consys_nme(phicsys,'v',t,FALSE,NULL),t) ;
	    if (pkrow != NULL) pkvec_free(pkrow) ;
	    return (FALSE) ; } } } } }
/*
  We're done. Return success.
*/
  if (pkrow != NULL) pkvec_free(pkrow) ;
  *p_feas = feas ;

  return (TRUE) ; }




bool arc_phic (bool *p_feas, int *p_fxcnt, fxvar_struct **p_fxvars)

/*
  This is the controlling routine for a round of constraint propagation. It
  initialises the data structures, then processes the constraints on the
  propagation stack.  Constraint propagation ends when the stack is empty, or
  when the constraint system loses feasibility because some variable's domain
  is empty.  The bound vectors associated with phicsys are updated to reflect
  new bounds.

  Parameters:
    p_feas:	(o) TRUE if the system is still feasible, FALSE if it's not.
    p_fxcnt:	(o) the number of variables fixed on this round
    p_fxvars:	(o) vector of {index,value} pairs for variables fixed on
		this round

  Returns: TRUE if the propagation proceeds without error, FALSE otherwise.
*/

{ int cndx ;
  bool success = FALSE ;
  char *rtnnme = "phic" ;

/* system quicker sort routine */

  void qsort(void *base, size_t nel, size_t width,
	     int (*compar) (const void *, const void *)) ;

# ifndef NDEBUG
  int ndx ;
  vbndchg_struct *vbndchg ;
  cbndchg_struct *cbndchg ;
# endif

# ifdef PARANOIA
  if (main_opts->phic == 0)
  { errmsg(467,rtnnme) ;
    return (FALSE) ; }
  if (p_fxcnt == NULL)
  { errmsg(2,rtnnme,"&fxcnt") ;
    return (FALSE) ; }
  if (p_fxvars == NULL)
  { errmsg(2,rtnnme,"&fxvars") ;
    return (FALSE) ; }
# endif
# ifndef NDEBUG
  if (main_opts->print.arcprop >= 1)
    outfmt(logchn,gtxecho,"\n(%s) propagating bounds for system %s.",
	   rtnnme,phicsys->nme) ;
# endif
/*
  Initialisation. We'll assume the problem will remain feasible, and that
  we won't fix any variables.
*/
  *p_feas = TRUE ;
  *p_fxcnt = 0 ;

# ifdef NDEBUG
/*
  Sort the entries on the propagation stack -- this turns out to be important
  to improve the rate of convergence.
*/
  qsort(propstack,propstack_top+1,sizeof(propstack_struct),propstack_cmp) ;
# else
/*
  If the debug prints are on, then we have much more work to do. We have to
  initialise the arrays that record which variables and constraints actually
  change bounds. We may also need to print the constraint stack before and
  after the initial sort.
*/
  if (main_opts->print.arcprop >= 1)
  { for (cndx = 1, vbndchg = &varbnd_chgs[1] ;
	 cndx <= phicsys->varcnt ;
	 cndx++, vbndchg++)
    { vbndchg->empty = TRUE ;
      vbndchg->cnt = 0 ; }
    
    if (main_opts->print.arcprop >= 3)
    { for (cndx = 1, cbndchg = &conbnd_chgs[1] ;
	   cndx <= phicsys->concnt ;
	   cndx++, cbndchg++)
      { cbndchg->empty = TRUE ;
	cbndchg->cnt = 0 ; } } }

  if (main_opts->print.arcstk >= 3)
  { outfmt(logchn,gtxecho,"\n    stack prior to sort:") ;
    for (ndx = 0 ; ndx <= propstack_top ; ndx++)
    { if (ndx%3 == 0)
	outfmt(logchn,gtxecho,"\n\t") ;
      else
	outchr(logchn,gtxecho,',') ;
      outfmt(logchn,gtxecho,"[%s (%d),%.3g (%.3g%%)]",
	     consys_nme(phicsys,'c',propstack[ndx].cndx,FALSE,NULL),
	     propstack[ndx].cndx,propstack[ndx].delta,propstack[ndx].pct) ; } }

  qsort((void *) propstack,propstack_top+1,
	sizeof(propstack_struct),propstack_cmp) ;

  if (main_opts->print.arcstk >= 2)
  { outfmt(logchn,gtxecho,
	   "\n    sorted %d constraints on the propagation stack.",
	   propstack_top+1) ;
    if (main_opts->print.arcstk >= 3)
    { outfmt(logchn,gtxecho,"\n    stack after sort:",rtnnme) ;
      for (ndx = 0 ; ndx <= propstack_top ; ndx++)
      { if (ndx%3 == 0)
	  outfmt(logchn,gtxecho,"\n\t") ;
	else
	  outchr(logchn,gtxecho,',') ;
	outfmt(logchn,gtxecho,"[%s (%d),%.3g (%.3g%%)]",
	       consys_nme(phicsys,'c',propstack[ndx].cndx,FALSE,NULL),
	       propstack[ndx].cndx,
	       propstack[ndx].delta,propstack[ndx].pct) ; } } }
# endif
/*
  Start up a loop to pop a constraint off of the stack and process it. The
  correct processing routine is selected based on the type of constraint.
  The processing routines will, in general, push new constraints onto the
  stack for processing. Equalities and range constraints are handled by
  considering them equivalent to a <= and a >= constraint. The only real loss
  of efficiency is that the constraint coefficients will be fetched twice.
  (But we do have to be careful of the case where only one of the bounds is
  finite, and watch for loss of feasibility.)

  We shouldn't see contypGE --- they've been converted to contypLE in mpsin.
  Nor should we see contypNB --- mpsin ignores them.
*/

  for (cndx = popcon(), need_sort = FALSE ;
       cndx > 0 ;
       cndx = popcon(), need_sort = FALSE)
  {
#   ifndef NDEBUG
    if (main_opts->print.arcprop >= 3)
    { outfmt(logchn,gtxecho,"\n    processing %s constraint %s (%d).",
	     consys_prtcontyp(phicsys->ctyp[cndx]),
	     consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ; }
#   endif
    switch (phicsys->ctyp[cndx])
    { case contypLE:
      { success = typLE_revise(cndx,p_fxcnt,p_feas) ;
	break ; }
      case contypEQ:
      case contypRNG:
      { 
#	ifdef PARANOIA
	if (phicsys->clb[cndx].inf > 0 && phicsys->cub[cndx].inf > 0)
	{ errmsg(466,rtnnme,phicsys->nme,
		 consys_prtcontyp(phicsys->ctyp[cndx]),
		 consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
	  success = FALSE ;
	  break ; }
#	endif
	if (phicsys->clb[cndx].inf <= 0)
	{ success = typLE_revise(cndx,p_fxcnt,p_feas) ;
	  if (*p_feas == FALSE || success == FALSE) break ; }
	if (phicsys->cub[cndx].inf <= 0)
	{ success = typGE_revise(cndx,p_fxcnt,p_feas) ;
	  if (success == FALSE) break ; }
	break ; }
      case contypGE:
      case contypNB:
      { errmsg(463,rtnnme,consys_prtcontyp(phicsys->ctyp[cndx]),
	       consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
	success = FALSE ;
	break ; }
      default:
      { errmsg(461,rtnnme,(int) phicsys->ctyp[cndx],
	       consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
	success = FALSE ;
	break ; } }
    if (success == FALSE)
    { errmsg(462,rtnnme,consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx) ;
      return (FALSE) ; }
    if (*p_feas == FALSE)
    {
#     ifndef NDEBUG
      if (main_opts->print.arcprop >= 1)
        outfmt(logchn,gtxecho,
	     "\nT%d:#%d: (%s) feasibility lost, %d constraints unprocessed.\n",
	       main_status->tid,main_status->subprob,rtnnme,propstack_top+1) ;
#     endif
      arc_clear() ;
      break ; }
    if (need_sort == TRUE)
    {
#     ifndef NDEBUG
      if (main_opts->print.arcstk >= 3)
      { outfmt(logchn,gtxecho,"\n    stack prior to sort:") ;
	for (ndx = 0 ; ndx <= propstack_top ; ndx++)
	{ if (ndx%3 == 0)
	    outfmt(logchn,gtxecho,"\n\t") ;
	  else
	    outchr(logchn,gtxecho,',') ;
	  outfmt(logchn,gtxecho,"[%s (%d),%.3g (%.3g%%)]",
		 consys_nme(phicsys,'c',propstack[ndx].cndx,FALSE,NULL),
		 propstack[ndx].cndx,
		 propstack[ndx].delta,propstack[ndx].pct) ; } }
#     endif
      qsort(propstack,propstack_top+1,sizeof(propstack_struct),propstack_cmp) ;
#     ifndef NDEBUG
      if (main_opts->print.arcstk >= 2)
      { outfmt(logchn,gtxecho,
	       "\n    sorted %d constraints on the propagation stack.",
	       propstack_top+1) ;
	if (main_opts->print.arcstk >= 3)
	{ outfmt(logchn,gtxecho,"\n    stack after sort:",rtnnme) ;
	  for (ndx = 0 ; ndx <= propstack_top ; ndx++)
	  { if (ndx%3 == 0)
	      outfmt(logchn,gtxecho,"\n\t") ;
	    else
	      outchr(logchn,gtxecho,',') ;
	    outfmt(logchn,gtxecho,"[%s (%d),%.3g (%.3g%%)]",
	           consys_nme(phicsys,'c',propstack[ndx].cndx,FALSE,NULL),
		   propstack[ndx].cndx,
		   propstack[ndx].delta,propstack[ndx].pct) ; } } }
#     endif
    } }

# ifndef NDEBUG
  if (main_opts->print.arcprop >= 1)
  { outfmt(logchn,gtxecho,"\n    constraint propagation completed.") ;
    if (main_opts->print.arcprop >= 2)
    { bool first ;
    
      first = TRUE ;
      for (cndx = 1, vbndchg = &varbnd_chgs[1] ;
	   cndx <= phicsys->varcnt ;
	   cndx++, vbndchg++)
	if (vbndchg->empty == FALSE)
	{ if (first == TRUE)
	  { outfmt(logchn,gtxecho,"\n    updated variable bounds:") ;
	    first = FALSE ; }
	  outfmt(logchn,gtxecho,
		 "\n\t%s (%d) tightened %d times from [%g, %g] to [%g, %g].",
		 consys_nme(phicsys,'v',cndx,FALSE,NULL),cndx,vbndchg->cnt,
		 vbndchg->olb,vbndchg->oub,vbndchg->nlb,vbndchg->nub) ; }
      if (first == TRUE)
	outfmt(logchn,gtxecho,"\n    no updated variable bounds.") ;
      first = TRUE ;
      if (main_opts->print.arcprop >= 4)
      { for (cndx = 1, cbndchg = &conbnd_chgs[1] ;
	     cndx <= phicsys->concnt ;
	     cndx++, cbndchg++)
	  if (cbndchg->empty == FALSE)
	  { if (first == TRUE)
	    { outfmt(logchn,gtxecho,"\n    updated constraint bounds:") ;
	      first = FALSE ; }
	    outfmt(logchn,gtxecho,
		   "\n\t%s (%d) tightened %d times from [%s = %s, ",
		   consys_nme(phicsys,'c',cndx,FALSE,NULL),cndx,cbndchg->cnt,
		   consys_conbndnme('L',cndx,&cbndchg->olb),
		   consys_conbndval(&cbndchg->olb)) ;
	    outfmt(logchn,gtxecho,"%s = %s] to ",
		   consys_conbndnme('U',cndx,&cbndchg->oub),
		   consys_conbndval(&cbndchg->oub)) ;
	    outfmt(logchn,gtxecho,"[%s = %s, ",
		   consys_conbndnme('L',cndx,&cbndchg->nlb),
		   consys_conbndval(&cbndchg->nlb)) ;
	    outfmt(logchn,gtxecho,"%s = %s].",
		   consys_conbndnme('U',cndx,&cbndchg->nub),
		   consys_conbndval(&cbndchg->nub)) ; }
	if (first == TRUE)
	  outfmt(logchn,gtxecho,"\n    no updated constraint bounds.") ; } } }
# endif
/*
  Set the return values and return to the caller.
*/
  *p_fxvars = fxvars ;

  return (TRUE) ; }
