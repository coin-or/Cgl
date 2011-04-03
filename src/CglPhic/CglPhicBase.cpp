/*
  Copyright (C) 2011, Lou Hafer
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
/*
  This file contains class boilerplate, setup, initialisation, and reporting
  methods.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"

#include "CglPhic.hpp"

namespace {

/*! \brief Default zero tolerance

  This value is used for coefficients, bounds, variable values, etc.
*/
const double dfltZeroTol = 1.0e-12 ;

/* \brief Default feasibility tolerance

  This value is used when comparing constraint and variable bounds in
  feasibility tests.
*/
const double dfltFeasTol = 1.0e-8 ;

// Default value for infinity
const double dfltInfinity = COIN_DBL_MAX ;

/*! \brief Default propagation tolerance

  This value is used to determine if the change to a lhs bound is worth
  propagating.
*/
const double dfltPropTol = 1.0e-8 ;

/*! \brief Default revision limit for constraint lhs bounds

  The value comes from previous experience in bonsaiG.
*/
const int dfltRevLimit = 10 ;

}    // end file-local namespace


/*
  Default constructor
*/
CglPhic::CglPhic ()
  : zeroTol_(dfltZeroTol),
    feasTol_(dfltFeasTol),
    propTol_(dfltPropTol),
    infty_(dfltInfinity),
    revLimit_(dfltRevLimit),
    verbosity_(0),
    paranoia_(0),
    m_(0),
    n_(0),
    ourRowMtx_(false),
    rowMtx_(0),
    rm_rowStarts_(0),
    rm_rowLens_(0),
    rm_colIndices_(0),
    rm_coeffs_(0),
    ourColMtx_(0),
    colMtx_(0),
    cm_colStarts_(0),
    cm_colLens_(0),
    cm_rowIndices_(0),
    cm_coeffs_(0),
    rhsL_(0),
    rhsU_(0),
    ourColL_(false),
    colL_(0),
    ourColU_(false),
    colU_(0),
    intVar_(0),
    lhsL_(0),
    lhsU_(0),
    info_(0),
    szePending_(0),
    numPending_(0),
    pending_(0),
    isPending_(0),
    inProcess_(0),
    szeVarBndChgs_(0),
    numVarBndChgs_(0),
    varBndChgs_(0),
    hasChanged_(0)
{ }

/*
  Constructor with constraint system.
*/
CglPhic::CglPhic (const CoinPackedMatrix *const rowMtx,
		  const CoinPackedMatrix *const colMtx,
		  const double *const rhsLower, const double *const rhsUpper)
  : zeroTol_(dfltZeroTol),
    feasTol_(dfltFeasTol),
    propTol_(dfltPropTol),
    infty_(dfltInfinity),
    revLimit_(dfltRevLimit),
    verbosity_(0),
    paranoia_(0),
    m_(0),
    n_(0),
    ourRowMtx_(false),
    rowMtx_(0),
    rm_rowStarts_(0),
    rm_rowLens_(0),
    rm_colIndices_(0),
    rm_coeffs_(0),
    ourColMtx_(0),
    colMtx_(0),
    cm_colStarts_(0),
    cm_colLens_(0),
    cm_rowIndices_(0),
    cm_coeffs_(0),
    rhsL_(0),
    rhsU_(0),
    ourColL_(false),
    colL_(0),
    ourColU_(false),
    colU_(0),
    intVar_(0),
    lhsL_(0),
    lhsU_(0),
    info_(0),
    szePending_(0),
    numPending_(0),
    pending_(0),
    isPending_(0),
    inProcess_(0),
    szeVarBndChgs_(0),
    numVarBndChgs_(0),
    varBndChgs_(0),
    hasChanged_(0)
{
  loanSystem(rowMtx,colMtx,rhsLower,rhsUpper) ;
}
  

/*
  Destructor
*/
CglPhic::~CglPhic ()
{
 if (ourRowMtx_) delete rowMtx_ ;
 if (ourColMtx_) delete colMtx_ ;
 delete[] lhsL_ ;
 delete[] lhsU_ ;
 delete[] info_ ;
 if (ourColL_) delete[] colL_ ;
 if (ourColU_) delete[] colU_ ;
 delete[] pending_ ;
 delete[] isPending_ ;
 delete[] varBndChgs_ ;
 delete[] hasChanged_ ;
}



/*
  Install constraint system on loan from the client.
*/
bool CglPhic::loanSystem (const CoinPackedMatrix *const rowMtx,
				 const CoinPackedMatrix *const colMtx,
				 const double *const rhsLower,
				 const double *const rhsUpper)
{
  assert(rhsLower != 0 && rhsUpper != 0) ;
  assert(rowMtx != 0 || colMtx != 0) ;
/*
  Install the lower and upper rhs row bounds.
*/
  rhsL_ = rhsLower ;
  rhsU_ = rhsUpper ;
/*
  Install row- and column-major copies of the constraint matrix.
*/
  if (rowMtx != 0 && colMtx != 0) {
    rowMtx_ = rowMtx ;
    ourRowMtx_ = false ;
    colMtx_ = colMtx ;
    ourColMtx_ = false ;
  } else if (colMtx != 0) {
    colMtx_ = colMtx ;
    ourColMtx_ = false ;
    CoinPackedMatrix *tmp = new CoinPackedMatrix() ;
    tmp->reverseOrderedCopyOf(*colMtx_) ;
    rowMtx_ = tmp ;
    ourRowMtx_ = true ;
  } else if (rowMtx != 0) {
    rowMtx_ = rowMtx ;
    ourRowMtx_ = false ;
    CoinPackedMatrix *tmp = new CoinPackedMatrix() ;
    tmp->reverseOrderedCopyOf(*colMtx_) ;
    colMtx_ = tmp ;
    ourColMtx_ = true ;
  }
/*
  If the lhs bound arrays no longer fit, remove them.
*/
  if (colMtx_->getNumRows() > m_) {
    delete[] lhsL_ ;
    lhsL_ = 0 ;
    delete[] lhsU_ ;
    lhsU_ = 0 ;
  }
  m_ = colMtx_->getNumRows() ;
  n_ = colMtx_->getNumCols() ;
/*
  Unpack the matrix structure vectors --- we'll use these a lot, no sense
  doing it every time.
*/
  rm_rowStarts_ = rowMtx_->getVectorStarts() ;
  rm_rowLens_ = rowMtx_->getVectorLengths() ;
  rm_colIndices_ = rowMtx_->getIndices() ;
  rm_coeffs_ = rowMtx_->getElements() ;
  cm_colStarts_ = colMtx_->getVectorStarts() ;
  cm_colLens_ = colMtx_->getVectorLengths() ;
  cm_rowIndices_ = colMtx_->getIndices() ;
  cm_coeffs_ = colMtx_->getElements() ;

  return (true) ;
}


/*
  Install column bounds on loan from client.
*/
void CglPhic::loanColBnds (double *const colLower,
				  double *const colUpper)
{
  assert(colLower && colUpper) ;
  colL_ = colLower ;
  ourColL_ = false ;
  colU_ = colUpper ;
  ourColU_ = false ;
}


/*
  Copy column bounds provided by client.
*/
void CglPhic::setColBnds (const double *const colLower,
			  const double *const colUpper)
{
  assert(colLower && colUpper) ;
  CoinMemcpyN(colLower,n_,colL_) ;
  ourColL_ = true ;
  CoinMemcpyN(colUpper,n_,colU_) ;
  ourColU_ = true ;
}

/*
  Return copies of the row lhs bounds arrays.

  Takes a little bit of work because the internal structure is a bit more
  complicated.
*/
void CglPhic::getRowLhsBnds (double *&lhsLower, double *&lhsUpper) const
{
  assert(lhsL_ && lhsU_) ;
  if (!lhsLower) lhsLower = new double [m_] ;
  for (int i = 0 ; i < m_ ; i++) {
    if (lhsL_[i].infCnt_ == 0)
      lhsLower[i] = lhsL_[i].bnd_ ;
    else
      lhsLower[i] = -infty_ ;
  }
  if (!lhsUpper) lhsUpper = new double [m_] ;
  for (int i = 0 ; i < m_ ; i++) {
    if (lhsU_[i].infCnt_ == 0)
      lhsUpper[i] = lhsU_[i].bnd_ ;
    else
      lhsUpper[i] = infty_ ;
  }
}

/*
  Calculate upper and lower lhs bounds for a given row, and the row measures
  (gaps and norm). The gaps are calculated over variables with both bounds
  finite. An infinite gap causes problems elsewhere, and there are explicit
  mechanisms for converting infinite bounds to finite bounds where possible.

  For a <= constraint, we will only ever use L(i) during propagation.
  Similarly, for a >= constraint, only U(i). The coefficient norm won't
  change, and the gaps have limited utility. Arguably this method should
  pay attention. The question is whether the additional tests to avoid
  calculation would disrupt the execution pipeline and slow down execution
  more than the additional work.  -- lh, 110326 --
*/
void CglPhic::calcLhsBnds (int i)
{
  assert(0 <= i && i < m_) ;
  assert(colL_ && colU_ && rowMtx_) ;
  assert(lhsL_ && lhsU_ && info_) ;

  double l1norm = 0.0 ;
  double posGap = 0.0 ;
  double negGap = 0.0 ;
  int infU = 0 ;
  int ndxU = -1 ;
  double bndU = 0.0 ;
  int infL = 0 ;
  int ndxL = -1 ;
  double bndL = 0.0 ;
/*
  Walk the row, accumulating the bound values.
*/
  const CoinBigIndex firstjj = rm_rowStarts_[i] ;
  const CoinBigIndex lastjj = firstjj+rm_rowLens_[i] ;
  for (CoinBigIndex jj = firstjj ; jj < lastjj ; jj++) {
    const int j = rm_colIndices_[jj] ;
    const double aij = rm_coeffs_[jj] ;
    const double lj = colL_[j] ;
    const double uj = colU_[j] ;
    if (aij > zeroTol_) {
      l1norm += aij ;
      if (uj >= infty_) {
        infU++ ;
	ndxU = j ;
      } else {
	bndU += aij*uj ;
	if (lj > -infty_) posGap = CoinMax(posGap,aij*(uj-lj)) ;
      }
      if (lj <= -infty_) {
        infL++ ;
	ndxL = j ;
      } else {
	bndL += aij*lj ;
	if (uj < infty_) posGap = CoinMax(posGap,aij*(uj-lj)) ;
      }
    } else if (aij < -zeroTol_) {
      l1norm -= aij ;
      if (uj >= infty_) {
        infL++ ;
	ndxL = j ;
      } else {
	bndL += aij*uj ;
	if (lj > -infty_) negGap = CoinMin(negGap,aij*(uj-lj)) ;
      }
      if (lj <= -infty_) {
        infU++ ;
	ndxU = j ;
      } else {
	bndU += aij*lj ;
	if (uj < infty_) negGap = CoinMin(negGap,aij*(uj-lj)) ;
      }
    }
  }
/*
  Clean up and fill in the bound and info entries.
*/
  CglPhicConBnd &Li = lhsL_[i] ;
  CglPhicConBnd &Ui = lhsU_[i] ;
  CglPhicConInfo &infoi = info_[i] ;

  Li.bnd_ = bndL ;
  if (infL == 1)
    Li.infCnt_ = -(ndxL+1) ;
  else
    Li.infCnt_ = infL ;
  Li.revs_ = 0 ;

  Ui.bnd_ = bndU ;
  if (infU == 1)
    Ui.infCnt_ = -(ndxU+1) ;
  else
    Ui.infCnt_ = infU ;
  Ui.revs_ = 0 ;

  infoi.l1norm_ = l1norm ;
  infoi.posGap_ = posGap ;
  infoi.negGap_ = negGap ;
/*
  And a bit of debug printing
*/
  if (verbosity_ >= 4) {
    std::cout
      << "        " << "init " << rhsL_[i] << " < " << Li << " <= r(" << i
      << ") <= " << Ui << " < " << rhsU_[i] << ", l1 " << l1norm << ", pGap "
      << posGap << ", nGap " << negGap << std::endl ;
  }
}

/*
  Calculate upper and lower row lhs bounds.
*/
void CglPhic::initLhsBnds ()
{
  if (verbosity_ >= 1)
    std::cout << "  " << "Initialising row info and lhs bounds ... " ;
  assert(colL_ && colU_ && rowMtx_) ;
/*
  Allocate fresh bound and info arrays if we need them.
*/
  if (!lhsL_) lhsL_ = new CglPhicConBnd [m_] ;
  if (!lhsU_) lhsU_ = new CglPhicConBnd [m_] ;
  if (!info_) info_ = new CglPhicConInfo [m_] ;
/*
  Step through the major vectors and calculate measures and bounds.
*/
  if (verbosity_ >= 5) std::cout << std::endl ;
  for (int i = 0 ; i < m_ ; i++) {
    calcLhsBnds(i) ;
  }
  if (verbosity_ >= 1) {
    if (verbosity_ >= 5) std::cout << "    " ;
    std::cout << "done." << std::endl ;
  }
}

/*
  Initialise the propagation data structures
*/
void CglPhic::initPropagation ()
{
  assert(colMtx_) ;

/*
  Allocate and/or reset the data structures for the set of pending
  constraints: the set itself, and a cross reference. Reset is simply a matter
  of setting numPending back to zero.
*/
  int newSize = CoinMin(m_/4+10,m_) ;
  if (newSize > szePending_) {
    szePending_ = newSize ;
    delete[] pending_ ;
    pending_ = new CglPhicCand [szePending_] ;
  }
  if (!isPending_) isPending_ = new int [m_] ;
  CoinZeroN(isPending_,m_) ;
  numPending_ = 0 ;
  inProcess_ = -1 ;
/*
  And again, for the data structures that record variable bound changes.
*/
  newSize = CoinMin(n_/4+10,n_) ;
  if (newSize > szeVarBndChgs_) {
    szeVarBndChgs_ = newSize ;
    delete[] varBndChgs_ ;
    varBndChgs_ = new CglPhicBndChg [szeVarBndChgs_] ;
  }
  if (!hasChanged_) hasChanged_ = new int [n_] ;
  CoinZeroN(hasChanged_,n_) ;
  numVarBndChgs_ = 0 ;
}


/*
  Report out column bound changes as a pair of CoinPackedVectors.
*/
void CglPhic::getVarBoundChanges (CoinPackedVector &lbs,
				  CoinPackedVector &ubs)
{
  for (int ndx = 0 ; ndx < numVarBndChgs_ ; ndx++) {
    const CglPhicBndChg &chgRec = varBndChgs_[ndx] ;
    if (chgRec.revl_ > 0) lbs.insert(chgRec.ndx_,chgRec.nl_) ;
    if (chgRec.revu_ > 0) ubs.insert(chgRec.ndx_,chgRec.nu_) ;
  }
}


/*
  Print method for a constraint bound.
*/
std::ostream& operator<< (std::ostream &strm, CglPhic::CglPhicConBnd lhs)
{
  strm << "(" ;
  if (lhs.infCnt_ < 0)
    strm << "x(" << (-lhs.infCnt_)-1 << ")" ;
  else if (lhs.infCnt_ == 0 || lhs.infCnt_ >= 2)
    strm << lhs.infCnt_ ;
  else
    strm << "invalid!" ;
  strm << "," << lhs.bnd_ << ")" ;

  return (strm) ;
}

/*
  Print method for a variable bound change record.
*/
std::ostream& operator<< (std::ostream &strm, CglPhic::CglPhicBndChg chg)
{
  char typlet[3] = {'c','b','g'} ;
  strm
    << "x(" << chg.ndx_ << ") " << typlet[chg.type_]
    << " [" << chg.ol_ << "," << chg.ou_ << "] --#"
    << chg.revl_ << "," << chg.revu_
    << "#-> [" << chg.nl_ << "," << chg.nu_ << "]" ;

  return (strm) ;
}
    
