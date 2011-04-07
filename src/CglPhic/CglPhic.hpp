/*
  Copyright (C) 2011, Lou Hafer
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

#ifndef CglPhic_H
#define CglPhic_H

#include "CoinPackedMatrix.hpp"

/*! \file CglPhic.hpp

  This file contains classes used to implement bound propagation. Changes
  in the bounds of variables are propagated using knowledge of integrality.

  In constraint propagation jargon, we're doing partial hyperarc integer
  consistency (PHIC) over linear constraints. The algorithm enforces
  k-consistency if allowed to propagate to quiescence.
*/

/*! \brief An implementation of bound propagation.

  This class implements bound propagation (partial hyperarc integer
  consistency, PHIC, in constraint propagation jargon). Given a set of
  changes to bounds on variables, it will propagate the changes, taking
  integrality into account. The algorithm enforces k-consistency if allowed to
  propagate to quiescence.

  The constraint system data structures are loaned to the propagator. They
  are not modified, and ownership remains with the client.

  To use this class, follow this sequence of actions:
  - Create the CglPhic object.
  - Supply pointers to (const) row- and column-major copies of the constraint
    matrix and upper and lower row bounds. If only one matrix is provided, the
    other will be synthesized, but it's more efficient to provide pointers.
  - Supply pointers to arrays for the upper and lower bounds of variables.
*/
class CglPhic {

public:

  /*! \name Constructors and destructor */
  //@{

  /// Default constructor
  CglPhic() ;

  /*! \brief Constructor with constraint system

    Works as #loanSystem.
  */
  CglPhic(const CoinPackedMatrix *const rowMtx,
	  const CoinPackedMatrix *const colMtx,
	  const double *const rhsLower, const double *const rhsUpper) ;


  /// Destructor
  virtual ~CglPhic() ;

  //@}

  /*! \name Set and get methods */
  //@{

  /*! \brief Loan constraint system to the propagator

    You must supply at least one matrix and both rhs vectors. If only one
    matrix is provided, the other will be synthesized.  Responsibility
    for destruction remains with the client.  The client cannot obtain
    ownership of a matrix synthesized by the propagator.  Really, it's to
    your advantage to provide both if you have them available.
  */
  bool loanSystem(const CoinPackedMatrix *const rowMtx,
  		  const CoinPackedMatrix *const colMtx,
		  const double *const rhsLower, const double *const rhsUpper) ;

  /// Return const pointers to the constraint system
  inline void getSystem(const CoinPackedMatrix *&rowMtx,
  			const CoinPackedMatrix *&colMtx,
			const double *&rhsLower,
			const double *&rhsUpper) const
  { rowMtx = rowMtx_ ;
    colMtx = colMtx_ ;
    rhsLower = rhsL_ ;
    rhsUpper = rhsU_ ;
  }

  /*! \brief Loan variable bounds arrays to the propagator

    These will be modified with tightened bounds. Handy if you already have
    bounds arrays suitable for use.
  */
  void loanColBnds(double *const colLower, double *const colUpper) ;

  /*! \brief Set variable bounds in the propagator

    The bounds arrays are copied for internal use.
  */
  void setColBnds(const double *const colLower, const double *const colUpper) ;

  /*! \brief Return pointers to the modified variable bounds arrays

    If you loaned variable bounds arrays to the propagator with #loanColBnds,
    you'll get back the same arrays, so there's no need to ask.
  */
  inline void getColBnds(double *&colLower, double *&colUpper) const
  { colLower = colL_ ;
    colUpper = colU_ ; }

  /*! \brief Retrieve variable bound changes as packed vectors

    The vectors are suitable as inputs for #editColBnds.
  */
  void getColBndChgs(CoinPackedVector &lbs, CoinPackedVector &ubs) ;

  /*! \brief Apply a set of changes to the variable bounds arrays

    Changes are supplied as a pair of packed vectors (viz. #getVarBndChanges)
  */
  void editColBnds(const CoinPackedVector *const lbs,
  		   const CoinPackedVector *const ubs) ;

  /*! \brief Paired bound structure

    It's often handy to get both both bounds together for a given row or
    column, even when only one has changed.
  */
  struct CglPhicBndPair {
    /// index
    int ndx_ ;
    /// changes (bit set): 0x1 if lower bound is changed, 0x2 for upper
    int changed_ ;
    /// lower bound
    double lb_ ;
    /// upper bound
    double ub_ ;
  } ;

  /*! \brief Retrieve a vector of bound changes for variables where one or
  	     both bounds have changed.

    The pointers \p newBnds and \p oldBnds will point to correlated vectors of
    bound pairs with the tightened and original bounds, respectively. The
    return value of the method is the length of each vector. If \p newBnds or
    \p oldBnds is 0, the corresponding vector is not returned.
    
    By default, all changed bounds are returned. The boolean values \p
    binVars, \p intVars, and \p conVars can be used to control inclusion
    for binary, general integer, and continuous variables, respectively.
  */
  int getColBndChgs(CglPhicBndPair **newBnds, CglPhicBndPair **oldBnds,
  		    bool binVars = true, bool intVars = true,
		    bool conVars = true) ;

  /*! \brief Apply a set of changes to the variable bounds arrays

    Changes are supplied as a vector of CglPhicBndPair entries of length
    \p len.
  */
  void editColBnds(int len, const CglPhicBndPair *const newBnds) ;

  /*! \brief Return the row lhs lower bound for the specified row.

    Returns the finite lhs bound L(i) if it exists, -infinity otherwise.
    If the internal lhs bound arrays are not initialised, returns -infinity.
  */
  inline double getRowLhsLB (int i) const
  { double retval = -infty_ ;
    if (lhsL_ && lhsL_[i].infCnt_ == 0) retval = lhsL_[i].bnd_ ;
    return (retval) ;
  }
      
  /*! \brief Return the row lhs upper bound for the specified row.

    Returns the finite lhs bound U(i) if it exists, infinity otherwise.
    If the internal lhs bound arrays are not initialised, returns infinity.
  */
  inline double getRowLhsUB (int i) const
  { double retval = infty_ ;
    if (lhsU_ && lhsU_[i].infCnt_ == 0) retval = lhsU_[i].bnd_ ;
    return (retval) ;
  }

  /*! \brief Return copies of the constraint lhs bounds arrays
  
    The information returned is incomplete, in the sense that it does not
    have separate finite and infinite contributions (cf. #CglPhicLhsBnd).
    If either \p lhsLower or \p lhsUpper is 0, an array of length #m_ will be
    supplied.
  */
  void getRowLhsBnds(double *&lhsLower, double *&lhsUpper) const ;

  /*! \brief Retrieve variable bound changes as packed vectors

    As with #getRowLhsBnds, there's a loss of information.
  */
  void getRowLhsBndChgs (CoinPackedVector &Lbs, CoinPackedVector &Ubs) ;

  /*! \brief Retrieve a vector of row lhs bound changes for constraints where
  	     one or both bounds have changed.

    The pointers \p newBnds and \p oldBnds will point to correlated vectors of
    bound pairs with the tightened and original bounds, respectively. The
    return value of the method is the length of each vector. If \p newBnds or
    \p oldBnds is 0, the corresponding vector is not returned.
  */
  int getRowLhsBndChgs(CglPhicBndPair **newBnds, CglPhicBndPair **oldBnds) ;


  /*! \brief Loan variable type array to the propagator

    This array will not be changed.
  */
  inline void loanColType(const char *const colType) { intVar_ = colType ; }

  /// Return a pointer to the variable type array
  inline void getColType(const char *&colType) { colType = intVar_ ; }

  //@}

  /*! \name Propagation methods */
  //@{

  /// Calculate lhs bounds for the specified row
  void calcLhsBnds(int i) ;

  /// Initialise row lhs bounds
  void initLhsBnds() ;

  /// Initialise propagation data structures
  void initPropagation() ;

  /*! \brief Clear propagation data structures

    This clears the pending set and bound change records.
  */
  void clearPropagation() ;

  /// Attempt to tighten the upper bound on a variable
  bool tightenVu(int t, double gap, double ait, double corr) ;

  /// Attempt to tighten the lower bound on a variable
  bool tightenVl(int t, double gap, double ait, double corr) ;

  /// Process a single constraint to propagate changes in the lhs bounds
  void processRow (int i, bool &feas) ;

  /*! \brief Propagate the effect of change in a variable bound

    This method walks the column for the specified variable x_j looking for
    changes to row lhs bounds L_i and U_i. If it finds significant changes
    the constraint is added to the pending set for further processing. The
    parameter \p feas will be set on return to indicate feasibility.
  */
  void chgColBnd(int j, char bnd, double nbndj, bool &feas) ;

  /*! \brief Propagate bound changes

    This method is responsible for processing the set of pending constraints
    to propagate bound changes. The parameter \p feas will be set on return to
    indicate feasibility.
  */
  void propagate (bool &feas) ;

  /*! \brief Revert bound changes

    This method will revert the current set of bound changes. The two
    parameters control whether variable bounds or constraint lhs bounds (or
    both) are reverted.
  */
  void revert (bool revertColBnds, bool revertRowBnds) ;

  /*! \brief Tighten up the existing system

    This method attempts to tighten constraint and variable bounds in the
    system as it exists. It will select a promising subset of constraints (if
    such exist) and attempt to generate and propagate bound changes. The
    parameter \p feas will be set on return to indicate feasibility.
  */
  void tightenAbInitio(bool &feas) ;

  /// Add a constraint to the set of constraints waiting to be processed.
  void addToPending(int i, double delta, double metric) ;

  /*! \brief Retrieve the constraint with the best metric from the set of
             constraints waiting to be processed.
  */
  int getFromPending() ;

  //@}

  /*! \name Control parameter methods */

  //@{

  /// Set the left-hand-side bound revision limit
  inline void setLhsRevisionLimit(int revLimit) { revLimit_ = revLimit ; }
  /// Get the left-hand-side bound revision limit
  inline int getLhsRevisionLimit() const { return (revLimit_) ; }

  /// Set the zero tolerance
  inline void setZeroTol(double tol) { zeroTol_ = tol ; }
  /// Get the zero tolerance
  inline double getZeroTol() const { return (zeroTol_) ; }

  /// Set the feasibility tolerance
  inline void setFeasTol(double tol) { feasTol_ = tol ; }
  /// Get the feasibility tolerance
  inline double getFeasTol() const { return (feasTol_) ; }

  /// Set the propagation tolerance
  inline void setPropTol(double tol) { propTol_ = tol ; }
  /// Get the propagation tolerance
  inline double getPropTol() const { return (propTol_) ; }

  /// Set the value of infinity
  inline void setInfinity(double val) { infty_ = val ; }
  /// Get the value of infinity
  inline double getInfinity() const { return (infty_) ; }

  /// Set the print level (#verbosity_)
  inline void setVerbosity(int val) { verbosity_ = val ; }
  /// Get the print level
  inline double getVerbosity() const { return (verbosity_) ; }

  /// Set the level of paranoia (#paranoia_)
  inline void setParanoia(int val) { paranoia_ = val ; }
  /// Get the vlevel of paranoia
  inline double getParanoia() const { return (paranoia_) ; }
  //@}

private:

  /*! \name Control parameters */
  //@{
  /*! \brief Zero tolerance
  
    This value is used to determine if coefficients, bounds, variable values,
    etc., are zero.
  */
  double zeroTol_ ;
  /*! \brief Feasibility tolerance

    This value is used to determine if a change in a constraint or variable
    bound has resulted in loss of feasibility.
  */
  double feasTol_ ;
  /*! \brief Propagation tolerance

    This value is used to determine if a change to a lhs bound is worth
    propagating.
  */
  double propTol_ ;
  /*! \brief Infinity

    The value used as infinity during propagation.
  */
  double infty_ ;
  /*! \brief Revision limit for constraint lhs bounds

    Finite lower and upper bounds on the constraint lhs are recalculated from
    scratch after this many revisions to remove accumulated numerical error.
  */
  int revLimit_ ;
  /*! \brief verbosity

    Controls the amount of output:
    -  0: catatonic
    -  1: summary messages
    -  2: print a list of changed variable bounds
    -  3: print revisions to bounds on variables
    -  4: print revisions to bounds on rows
    -  5: print constraints queued for processing
    -  6: add additional processing traces
    Anything above 1 is likely too much unless you're debugging.
  */
  int verbosity_ ;
  /*! \brief paranoia

    Controls the level of paranoia:
    -  0: all clients are careful, conscientious, and perfect
    -  1: everyone makes the occasional mistake
    -  2: trust no one
    -  3: trust no one, not even myself!

    Defaults to 0; 1 is useful for application development. Higher values are
    useful for debugging CglPhic. 
  */
  int paranoia_ ;
  //@}

  /*! \name Constraint system */
  //@{
  /// Number of constraints (rows)
  int m_ ;
  /// Number of variables (columns)
  int n_ ;
  /// Ownership of row-major copy
  bool ourRowMtx_ ;
  /// Row-major constraint matrix
  const CoinPackedMatrix *rowMtx_ ;
  /// Row-major row start vector
  const CoinBigIndex *rm_rowStarts_ ;
  /// Row-major row length vector
  const int *rm_rowLens_ ;
  /// Row-major column indices
  const int *rm_colIndices_ ;
  /// Row-major coefficients
  const double *rm_coeffs_ ;
  /// Ownership of column-major copy
  bool ourColMtx_ ;
  /// Column-major constraint matrix
  const CoinPackedMatrix *colMtx_ ;
  /// Column-major column start vector
  const CoinBigIndex *cm_colStarts_ ;
  /// Column-major column length vector
  const int *cm_colLens_ ;
  /// Column-major row indices
  const int *cm_rowIndices_ ;
  /// Column-major coefficients
  const double *cm_coeffs_ ;
  /// Row rhs lower bounds array (always owned by client)
  const double *rhsL_ ;
  /// Row rhs upper bounds array (always owned by client)
  const double *rhsU_ ;
  /// Ownership of variable lower bounds array
  bool ourColL_ ;
  /// Column (variable) lower bounds array
  double *colL_ ;
  /// Ownership of variable upper bounds array
  bool ourColU_ ;
  /// Column (variable) upper bounds array
  double *colU_ ;
  /*! \brief Variable type

    Coded as:
    - 0: continuous ('c')
    - 1: binary ('b')
    - 2: general integer ('g')
    The character is shown in various print statements.
  */
  const char *intVar_ ;
  //@}

  /*! \name Propagation data structures */
  //@{

  /*! \brief Constraint lhs bound
  
    A count of the number of infinite variable bounds contributing to the
    lhs bound makes it easier to detect when a bound can be changed
    from infinite to finite. The l1norm puts propagation candidates on an
    equal footing.

    Lhs bounds are updated when variable bounds change, and numerical error
    accumulates. The bound must be periodically recalculated to minimise this
    error. Revs_ keeps track.
  */
  struct CglPhicLhsBnd {
    /// number of times this bound has been revised
    int revs_ ;
    /* \brief number of infinities contributing to the bound

      The value 1 is invalid. If there's a single variable x_j with an
      infinite bound, this field contains -(j+1).
    */
    int infCnt_ ;
    /// finite component of the bound
    double bnd_ ;
  } ;

  /// Row lhs lower bounds array
  CglPhicLhsBnd *lhsL_ ;
  /// Row lhs upper bounds array
  CglPhicLhsBnd *lhsU_ ;

  /*! \brief Constraint characteristics

    Various measures that are useful when determining whether to process a
    constraint. The two gap values are used for ab initio setup only.
  */
  struct CglPhicConInfo {
    /// L1-norm of coefficients
    double l1norm_ ;
    /// max{k | a_ik > 0}(a_ik(u_k-l_k))  (over finite bounds only)
    double posGap_ ;
    /// min{k | a_ik < 0}(a_ik(u_k-l_k))  (over finite bounds only)
    double negGap_ ;
  } ;

  /// Constraint characteristics array
  CglPhicConInfo *info_ ;

  /*! \brief Propagation candidate
  
    It's important to the efficiency of the algorithm that large changes in
    lhs bounds be propagated first. Delta_ and metric_ are used for this
    purpose.
  */
  struct CglPhicCand {
    /// constraint index
    int ndx_ ;
    /// accumulated bound change
    double delta_ ;
    /// percentage bound change (sort key)
    double metric_ ;
  } ;
  
  /// Capacity of #pending_
  int szePending_ ;
  /// The number of constraints waiting for bound propagation
  int numPending_ ;
  /// The set of constraints waiting for bound propagation.
  CglPhicCand *pending_ ;
  /*! \brief Cross-reference from constraint index to pending entry.

    Contains (index+1) of the entry in #pending_, or zero if the constraint is
    not pending.
  */
  int *isPending_ ;
  /// The constraint currently in process
  int inProcess_ ;

  /*! \brief Constraint lhs bound change record

    For a given constraint, tracks changes to the lhs bounds. Created when a
    bound is first tightened.
  */
  struct CglPhicLhsBndChg {
    /// Constraint index
    int ndx_ ;
    /// Number of revisions to lower bound
    int revL_ ;
    /// Original lower bound
    CglPhicLhsBnd oL_ ;
    /// New lower bound
    CglPhicLhsBnd nL_ ;
    /// Number of revisions to upper bound
    int revU_ ;
    /// Original lower bound
    CglPhicLhsBnd oU_ ;
    /// New lower bound
    CglPhicLhsBnd nU_ ;
  } ;

  /// Capacity of #lhsBndChgs_
  int szeLhsBndChgs_ ;
  /// Number of constraint lhs change records
  int numLhsBndChgs_ ;
  /// Constraint lhs bound change records
  CglPhicLhsBndChg *lhsBndChgs_ ;
  /*! \brief Cross-reference from constraint index to change record

    Contains (index+1) of the entry in #lhsBndChgs_, or zero if there's no
    change record.
  */
  int *lhsHasChanged_ ;


  /*! \brief Variable bound change record

    For a given variable, tracks changes to the bounds. Created when a bound
    is first tightened.
  */
  struct CglPhicVarBndChg {
    /// Variable index
    int ndx_ ;
    /// Variable type
    int type_ ;
    /// Number of revisions to lower bound
    int revl_ ;
    /// Original lower bound
    double ol_ ;
    /// New lower bound
    double nl_ ;
    /// Number of revisions to upper bound
    int revu_ ;
    /// Original lower bound
    double ou_ ;
    /// New lower bound
    double nu_ ;
  } ;

  /// Capacity of #varBndChgs_
  int szeVarBndChgs_ ;
  /// Number of variable bound change records
  int numVarBndChgs_ ;
  /// Variable bound change records
  CglPhicVarBndChg *varBndChgs_ ;
  /*! \brief Cross-reference from variable index to change record

    Contains (index+1) of the entry in #varBndChgs_, or zero if there's no
    change record.
  */
  int *varHasChanged_ ;

  //@}

  friend std::ostream& operator<< (std::ostream &out,
				   CglPhic::CglPhicLhsBnd lhs) ;

  friend std::ostream& operator<< (std::ostream &out,
  				   CglPhic::CglPhicLhsBndChg chg) ;

  friend std::ostream& operator<< (std::ostream &out,
  				   CglPhic::CglPhicVarBndChg chg) ;

} ;

/*! \brief Print a constraint bound

  Usually, "(infCnt,bnd)". When there is a single infinity remaining, the
  form is "(x(j),bnd)".
*/
std::ostream& operator<< (std::ostream &out, CglPhic::CglPhicLhsBnd lhs) ;
/*! \brief Print a bound change record

  The form is "x(j) [oldlj,olduj] --#lj revs,uj revs#-> [newlj,newuj]" where
  revs are the number of times the bound was revised.
*/
std::ostream& operator<< (std::ostream &out, CglPhic::CglPhicVarBndChg chg) ;

#endif
