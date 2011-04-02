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

  /*! \brief Loan variable bounds arrays to the propagator

    These will be changed; handy if you just need the updated bounds arrays.
  */
  void loanColBnds(double *const colLower, double *const colUpper) ;

  /*! \brief Set variable bounds in the propagator

    The bounds arrays are copied for internal use.
  */
  void setColBnds(const double *const colLower, const double *const colUpper) ;

  /*! \brief Return pointers to the modified variable bounds arrays

    If you loaned variable bounds arrays to the propagator, you'll get back
    the same arrays.
  */
  inline void getColBnds(double *&colLower, double *&colUpper) const
  { colLower = colL_ ;
    colUpper = colU_ ; }

  /*! \brief Return copies of the constraint lhs bounds arrays
  
    If either \p lhsLower or \p lhsUpper is 0, an array of length #m_ will be
    supplied.
  */
  void getRowLhsBnds(double *&lhsLower, double *&lhsUpper) const ;

  /*! \brief Loan variable type array to the propagator

    This array will not be changed.
  */
  inline void loanColType(const char *const colType) { intVar_ = colType ; }


  /// Retrieve variable bound changes as packed vectors
  void getVarBoundChanges(CoinPackedVector &lbs, CoinPackedVector &ubs) ;

  //@}

  /*! \name Propagation methods */
  //@{

  /// Calculate lhs bounds for the specified row
  void calcLhsBnds(int i) ;

  /// Initialise row lhs bounds
  void initLhsBnds() ;

  /// Initialise propagation data structures
  void initPropagation() ;

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
  void changeVarBound(int j, char bnd, double nbndj, bool &feas) ;

  /*! \brief Propagate bound changes

    This method is responsible for processing the set of pending constraints
    to propagate bound changes. The parameter \p feas will be set on return to
    indicate feasibility.
  */
  void propagate (bool &feas) ;

  /*! \brief Tighten up the existing system

    This method attempts to tighten constraint and variable bounds in the
    system as it exists. It will select a promising subset of constraints (if
    such exist) and attempt to generate and propagate bound changes.
  */
  void tightenAbInitio() ;

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
    -  2: add changes to bounds on variables
    -  3: add changes to bounds on rows
    -  4: add constraints queued for processing
    -  5: add processing traces
    Anything above 1 is likely too much unless you're debugging.
  */
  int verbosity_ ;
  /*! \brief paranoia

    Controls the level of paranoia:
    -  0: all clients are careful, conscientious, and perfect
    -  1: everyone makes the occasional mistake
    -  2: trust no one
    -  3: trust no one, not even myself
    Defaults to 0 and that's where you should leave it unless you're
    debugging.
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
    - 0: continuous
    - 1: binary
    - 2: general integer
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
  struct CglPhicConBnd {
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
  CglPhicConBnd *lhsL_ ;
  /// Row lhs upper bounds array
  CglPhicConBnd *lhsU_ ;

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

  /*! \brief Bound change record

    For a given variable, tracks changes to the bounds. Created when a bound
    if first tightened.
  */
  struct CglPhicBndChg {
    /// Variable index
    int ndx_ ;
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
  /// Number of change records
  int numVarBndChgs_ ;
  /// Variable bound change records
  CglPhicBndChg *varBndChgs_ ;
  /*! \brief Cross-reference from variable index to change record

    Contains (index+1) of the entry in #varBndChgs_, or zero if there's no
    change record.
  */
  int *hasChanged_ ;

  //@}

  friend std::ostream& operator<< (std::ostream &out,
				   CglPhic::CglPhicConBnd lhs) ;

} ;

/*! \brief Print a constraint bound

  Usually, "(infCnt,bnd)". When there is a single infinity remaining, the
  form is "(x(j),bnd)".
*/
std::ostream& operator<< (std::ostream &out, CglPhic::CglPhicConBnd lhs) ;

#endif
