// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglProbing_H
#define CglProbing_H

#include <string>

#include "CglCutGenerator.hpp"

/*! \brief Probing Cut Generator

  CglProbing borrows a subset of probing ideas put into OSL about ten years
  ago and extends them.  The only known documentation is a copy of a talk
  handout - Robin Lougee and Lou Hafer have copies.

  For selected integer variables (probe variables, see #setMode), the effect
  of changing bounds to force them up or down to the next integer is
  investigated.  Changing a bound on the probe variable may in turn affect
  other variables (continuous as well as integer).  There are various
  possible results:
  <ol>
    <li>
    It is shown that the problem is infeasible or can be fathomed by bound
    in one or both directions. If one direction remains viable (monotone),
    we can generate one or more column cuts and continue probing.
    </li>
    <li>
    It can happen that x -> 0 implies y -> 1, or any of three similar
    implications. This can be captured in a disaggregation cut.
    </li>
    <li>
    If both ways are feasible, it can happen that x -> 0 implies y -> 1
    *and* x -> 1 implies y -> 1 (or similar; again, a column cut).  More
    common is that x -> 0 implies y -> 1 and x -> 1 implies y -> 0 (or
    similar). This can be captured in an implication cut.
    </li>
    <li>
    x -> 1 causes a constraint to become slack by some amount c.  We can
    tighten the constraint ax + .... <= b (where a may be zero) to
    (a+c)x + .... <= b. This is coefficient strengthening. During
    preprocessing, it can be done in place, altering the coefficient in the
    existing constraint.
    </li>
  </ol>

  \todo
  The reader will note that when an implication is discovered we could
  directly substitute x for y, or vice versa, during preprocessing. This is
  not done at present as there is no mechanism for returning the information.

  The generated cuts are inserted into an OsiCut object. Both row cuts and
  column cuts may be returned.

  In all modes, user specified limits (number of variables probed, types of
  cuts generated, etc.) are (mostly) observed.
*/
class CglProbing : public CglCutGenerator {
   friend int CglProbingUnitTest(const OsiSolverInterface *siP,
				 const std::string mpsDir) ;
 
public:

  /*! \name Generate Cuts */
  //@{
  /*! \brief Generate probing/disaggregation cuts and tighter bounds

    Cuts are generated for the constraint system held in \p si and returned in
    \p cs.
  */
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
			    const CglTreeInfo info = CglTreeInfo()) const ;
			    
  /*! \brief Generate probing/disaggregation cuts and tighter bounds

    Cuts are generated for the constraint system held in \p si and returned in
    \p cs. In addition, the client will have access to the tightened variable
    bound arrays through #tightLower and #tightUpper.

    Return value is 0 for feasible, 1 for infeasible.
  */
  int generateCutsAndModify(const OsiSolverInterface &si, OsiCuts &cs, 
			    CglTreeInfo *info) ;
  //@}

  /**@name Get tighter column bounds */
  //@{
  /// Lower
  const double *tightLower() const ;
  /// Upper
  const double *tightUpper() const ;
  //@}

  /*! \name Behavioural Options
  
    This group of methods controls the type of cuts generated and other
    aspects of overall behaviour.
  */
  //@{
  /*! \brief Set cut generation mode

    The mode options are:

    - 1: Look at unsatisfied integer variables, using current bounds.
	 Probing will be done on all looked at.
    - 2: Look at all integer variables, using current bounds.
	 Probing will be done on all.

    See #mode_ for additional detail.
  */
  void setMode(int mode) ;
  /// Get cut generation mode
  int getMode() const ;

  /*! \brief Specify whether to generate row cuts

    When row cuts are disabled, CglProbing is simply fixing variables based on
    probing. The parameter is interpreted as single-bit flags:

    - 0x0: do not generate row cuts
    - 0x1: generate disaggregation row cuts
    - 0x2: generate coefficient (strengthening) row cuts
    - 0x4: generate column cuts only

    A negative value limits row cuts to the root; only variable fixing will be
    performed in the search tree.
  */
  void setRowCuts(int type) ;
  /// Get the row cut generation setting
  int rowCuts() const ;

  /*! \brief Inquire if row cuts can be generated other than at the root node

     Returns true if cut generation may produce row cuts in the tree, as
     opposed to only at the root node. For CglProbing, this is equivalent to
     #rowCuts() > 0.
     Provided so that a client can know if the constraint matrix will change
     in a search tree.
  */
  virtual bool mayGenerateRowCutsInTree() const ;

  /*! \brief Specify whether to use the objective as a constraint

      -  1: do add the objective as a constraint while probing
      -  0: don't add the objective as a constraint while probing
      - -1: don't use the objective as a constraint; further, don't compare
	    estimated degradation of the objective against a cutoff value
   
    Normally probing compares the estimated objective after forcing integrality
    to a cutoff value, when a cutoff is available. A value of -1 prevents this.
  */
  void setUsingObjective(int yesNo) ;
  /// Get use of objective
  int getUsingObjective() const ;

  //@}

  /*! \name Generator Limits
  
    This group of methods sets limits on cut generator activity.
  */
  //@{
  /// Set maximum number of passes per node
  void setMaxPass(int value) ;
  /// Get maximum number of passes per node
  int getMaxPass() const ;
  /// Set maximum number of passes per node  (root node)
  void setMaxPassRoot(int value) ;
  /// Get maximum number of passes per node (root node)
  int getMaxPassRoot() const ;
  /*! \brief Set log level
  
    - 0: none
    - 1 - a bit
    - 2 - more details
  */
  void setLogLevel(int value) ;
  /// Get log level
  int getLogLevel() const ;
  /// Set maximum number of unsatisfied variables to look at
  void setMaxProbe(int value) ;
  /// Get maximum number of unsatisfied variables to look at
  int getMaxProbe() const ;
  /// Set maximum number of variables to look at in one probe
  /// Set maximum number of unsatisfied variables to look at (root node)
  void setMaxProbeRoot(int value) ;
  /// Get maximum number of unsatisfied variables to look at (root node)
  int getMaxProbeRoot() const ;
  void setMaxLook(int value) ;
  /// Get maximum number of variables to look at in one probe
  int getMaxLook() const ;
  /// Set maximum number of variables to look at in one probe (root node)
  void setMaxLookRoot(int value) ;
  /// Get maximum number of variables to look at in one probe (root node)
  int getMaxLookRoot() const ;
  /// Set maximum number of elements in a row for it to be considered
  void setMaxElements(int value) ;
  /// Get maximum number of elements in a row for it to be considered
  int getMaxElements() const ;
  /// Set maximum number of elements in row for it to be considered (root node)
  void setMaxElementsRoot(int value) ;
  /// Get maximum number of elements in row for it to be considered (root node)
  int getMaxElementsRoot() const ;
  /// Set primal feasibility tolerance
  void setPrimalTolerance(double tol) { primalTolerance_ = tol ; }
  /// Get primal feasibility tolerance
  double getPrimalTolerance() { return (primalTolerance_) ; }
  //@}

  /**@name Get information back from probing */
  //@{
  /// Number looked at this time
  inline int numberThisTime() const
  { return numberThisTime_; }
  /// Which ones looked at this time
  inline const int *lookedAt() const
  { return lookedAt_; }
  //@}

  /**@name Mark which continuous variables are to be tightened */
  //@{
  /// Mark variables to be tightened
  void tightenThese(const OsiSolverInterface & solver, int number, const int * which) ;
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglProbing() ;
 
  /// Copy constructor 
  CglProbing(const CglProbing &) ;

  /// Clone
  virtual CglCutGenerator *clone() const ;

  /// Assignment operator 
  CglProbing &operator=(const CglProbing& rhs) ;
  
  /// Destructor 
  virtual ~CglProbing () ;

  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp) ;
  //@}
      
private:
  
 // Private member methods
  /**@name probe */
  //@{
  /*! \brief Probe and generate cuts

    This method does the heavy lifting. Given a set of variables to probe, do
    the probing and generate cuts. Returns false if the problem is infeasible.
  */
  bool probe(const OsiSolverInterface &si, 
	     OsiCuts &cs, 
	     double *const colLower, double *const colUpper,
	     const CoinPackedMatrix *const rowCopy,
	     const CoinPackedMatrix *const columnCopy,
	     const CoinBigIndex *const rowStartPos,
	     const int *const realRow, const double *const rowLower,
	     const double *const rowUpper,
	     const char *const intVar,
	     double *const minR, double *const maxR, int *const markR, 
	     CglTreeInfo *const info,
	     bool useObj, bool useCutoff, double cutoff) const ;

  /*! \brief Core method for overall control of probing.

    Common worker for generateCuts and generateCutsAndModify. Does common
    setup and calls #probe.  Returns false if the problem is infeasible.
  */
  bool gutsOfGenerateCuts(const OsiSolverInterface &si, 
			  OsiCuts &cs,
			  double *rowLower, double *rowUpper,
			  double *colLower, double *colUpper,
                          CglTreeInfo *info) const ;

  /*! \brief Calculate constraint left-hand-side bounds

    Calculate lower and upper bounds on the constraint lhs based on variable
    lower and upper bounds. The new bounds are loaded into minR and maxR.
    The array markR will be set as:
    - -1: constraint is useful for further propagation (at least one lhs
          bound finite)
    - -2: constraint is not useful (both lhs bounds infinite).
    In practice, finite translates as `reasonable magnitude'.
  */
  void calcRowBounds(const double *const colLower,
  		     const double *const colUpper,
		     const int *const column, const double *const rowElements, 
		     const CoinBigIndex *const rowStart,
		     const int *const rowLength,
		     const double *const rowLower,
		     const double *const rowUpper, 
		     double *const minR, double *const maxR, int *const markR,
		     int nRows) const ;
  //@}

  // Private member data

  /**@name Private member data */
  //@{

  /*! \brief Arbitrary bound to impose on rows and columns.

    This really shouldn't be necessary, but it's wired in and will require
    a fair bit of work to remove. As the code is broken up, it needs to be
    available across multiple files, so hide it as a private constant.
    -- lh, 101124 --
  */
  static const double CGL_REASONABLE_INTEGER_BOUND ;

  /*! \brief Arbitrarily large bound to impose on rows.

    This one is completely bogus. It gets installed as a lower or upper
    row bound in tighten2 if the existing row bound is larger than 1e20,
    or if a column bound exceeds 1e12.
    
    This bound leaks out of CglProbing as bogus `improved' bounds on
    variables, and that's a problem.
  */
  static const double CGL_BOGUS_1E60_BOUND ;

  /*
    These are vestigial at this point. I'm keeping them around on the premise
    that after generateCutsAndModify it will be useful to set these to point
    to the working colLower and colUpper. Not yet sure I want to keep them as
    mutable.  -- lh, 110402 --
  */ 
  /// Lower bounds on columns
  mutable double * colLower_ ;
  /// Upper bounds on columns
  mutable double * colUpper_ ;

  /*
    Not entirely sure why these need to be mutable. Check. -- lh, 101007 --
    Possibly because they're modified to reflect the size of the working copy?
    -- lh, 110402 --
  */
  /// Number of rows in snapshot (or when cliqueRow stuff computed)
  mutable int numberRows_ ;
  /// Number of columns in problem ( must == current)
  mutable int numberColumns_ ;

  /// Tolerance to see if infeasible
  double primalTolerance_ ;

  /*! \brief Probing mode

    Determines which variables will be probed. Coded as:
    - 1: just unsatisfied integer variables
    - 2: all integer variables
  */
  int mode_ ;

  /*! \brief Row cuts flags

    Broadly, the possibilities are:
    - generate row cuts that will be returned in an OsiCuts collection ;
    - generate row cuts that can replace (strengthen) rows in the original
      constraint system, returned through a CglTreeInfo structure ;
    - generate column cuts
    If the CglTreeInfo structure contains a strengthenRow vector, and you
    specify 0x4 here, no row cuts are returned in the cut collection, but rows
    are strengthened if possible.

    - 0x0: do not generate row cuts
    - 0x1: generate disaggregation row cuts
    - 0x2: generate coefficient (strengthening) row cuts
    - 0x4: generate column cuts only

    A negative value -n is interpreted as +n but only fixes variables unless
    at the root node
  */
  mutable int rowCuts_ ;

  /// Maximum number of passes to do in probing
  int maxPass_ ;
  /// Log level - 0 none, 1 - a bit, 2 - more details
  int logLevel_ ;
  /// Maximum number of unsatisfied variables to probe
  int maxProbe_ ;
  /// Maximum number of variables to look at in one probe
  int maxStack_ ;
  /// Maximum number of elements in row for scan
  int maxElements_ ;
  /// Maximum number of passes to do in probing at root
  int maxPassRoot_ ;
  /// Maximum number of unsatisfied variables to probe at root
  int maxProbeRoot_ ;
  /// Maximum number of variables to look at in one probe at root
  int maxStackRoot_ ;
  /// Maximum number of elements in row for scan at root
  int maxElementsRoot_ ;
  /// Whether to include the objective as a constraint
  int usingObjective_ ;
  /// Number of integer variables
  int numberIntegers_ ;
  /// Number of 0-1 integer variables
  int number01Integers_ ;
  /// Number looked at this time
  mutable int numberThisTime_ ;
  /*! \brief Total number of times called */
  mutable int totalTimesCalled_ ;
  /// Which ones looked at this time
  mutable int * lookedAt_ ;

  //@}
} ;

//#############################################################################
/** A function that tests the methods in the CglProbing class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
int CglProbingUnitTest(const OsiSolverInterface *siP,
		       const std::string mpsDir) ;

#endif
