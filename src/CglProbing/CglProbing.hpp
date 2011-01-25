// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglProbing_H
#define CglProbing_H

#include <string>

#include "CglCutGenerator.hpp"

  /** Only useful type of disaggregation is most normal
      For now just done for 0-1 variables
      Can be used for building cliques

      ?!?? What is this?
   */
  typedef struct {
    //unsigned int zeroOne:1; // nonzero if affected variable is 0-1
    //unsigned int whenAtUB:1; // nonzero if fixing happens when this variable at 1
    //unsigned int affectedToUB:1; // nonzero if affected variable fixed to UB
    //unsigned int affected:29; // If 0-1 then 0-1 sequence, otherwise true
    unsigned int affected;
  } disaggregationAction;

/*! \brief Probing Cut Generator

  This is a simplification of probing ideas put into OSL about ten years
  ago.  The only known documentation is a copy of a talk handout - we think
  Robin Lougee has a copy!

  For selected integer variables (e.g. unsatisfied ones) the effect of
  setting them up or down is investigated.  Setting a variable up may in turn
  set other variables (continuous as well as integer).  There are various
  possible results:

  1) It is shown that the problem is infeasible or can be fathomed by bound
     in one or both directions. If one direction remains viable, we can
     generate a column cut (and continue probing).

  2) If both ways are feasible, it can happen that x -> 0 implies y -> 1
     *and* x -> 1 implies y -> 1 (again a column cut).  More common
     is that x -> 0 implies y -> 1 and x -> 1 implies y -> 0 and we can
     substitute for y, which might lead later to more powerful cuts.

     Sadly, this is not done in this code as there is no mechanism for
     returning the information. This is a big TODO.

  3) x -> 1 causes a constraint to become slack by some amount c.  We can
     tighten the constraint ax + .... <= b (where a may be zero) to
     (a+c)x + .... <= b.  If this cut is violated then it is generated.

  4) Similarly, we can generate implied disaggregation cuts

  Note the following differences from cuts in OSL:

  a) OSL had structures intended to make this faster.
  b) The "chaining" in 2) was implemented.
  c) Row cuts modified the original constraint rather than adding a cut.
  b) This code can cope with general integer variables.

  The generated cuts are inserted into an OsiCut object. Both row cuts and
  column cuts may be returned.

  If a "snapshot" of the matrix exists then this will be used.  Presumably
  this will give global cuts and will be faster.  No check is done to see if
  cuts will be global. Otherwise use the current matrix held in the solver
  interface.

  The mode options are:

  - 0: Only unsatisfied integer variables will be looked at.  If no
       information exists for that variable then probing will be done so as a
       by-product you "may" get a fixing or infeasibility.  This will be fast
       and is only available if a snapshot exists (otherwise as 1).  The
       bounds in the snapshot are the ones used.

  - 1: Look at unsatisfied integer variables, using current bounds.  Probing
       will be done on all looked at.

  - 2: Look at all integer variables, using current bounds.  Probing will be
       done on all
*/
class CglProbing : public CglCutGenerator {
   friend int CglProbingUnitTest(const OsiSolverInterface *siP,
				 const std::string mpsDir) ;
 
public:

  /*! \name Generate Cuts */
  //@{
  /*! \brief Generate probing/disaggregation cuts

    Cuts are generated for the constraint system held in \p si unless a
    snapshot matrix exists, in which case the snapshot is used.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo()) const;
			    
  /*! \brief Generate probing/disaggregation cuts and tighter bounds

    New relaxed row bounds and tightened column bounds are generated.
    The number of infeasibilities is returned. 
  */
  int generateCutsAndModify( const OsiSolverInterface & si, OsiCuts & cs, 
			     CglTreeInfo * info);
  //@}

  /*! \name Snapshot and Clique Utilities */
  //@{
  /*! \brief Create a copy of the constraint matrix which will be used for
	     probing.

      This is to speed up the process and to give global cuts. Column bounds
      are tightened on the copy and remembered from one call to the next, as
      are the results of probing.  Returns 1 if infeasible, otherwise 0.

      If given, \p possible specifies which rows should be used; set entry i
      to 1 to retain row i in the copy, 0 to drop it. In addition, if \p
      possible is provided, rows that are determined to be ineffective will
      be dropped.  If \p withObjective is true, the objective will be added
      as an additional constraint.
  */
  int snapshot ( const OsiSolverInterface & si,
		  char * possible=NULL,
                  bool withObjective=true);

  /*! \brief Refresh an existing snapshot from the solver.
  
    Yes, the name is totally misleading.
  */
  virtual void refreshSolver(OsiSolverInterface * solver);

  /// Delete the snapshot
  void deleteSnapshot ( );

  /*! \brief Create cliques for use by probing.

      Only cliques >= minimumSize and < maximumSize are created.  The method
      will also try to extend cliques as a result of probing (root node
      only).

      Returns the number of cliques found.
  */
  int createCliques( OsiSolverInterface & si, 
		     int minimumSize=2, int maximumSize=100);

  /// Delete all clique information
  void deleteCliques();
  //@}

  /**@name Get tighter column bounds */
  //@{
  /// Lower
  const double * tightLower() const;
  /// Upper
  const double * tightUpper() const;
  /// Array which says tighten continuous
  const char * tightenBounds() const
  { return tightenBounds_;}
  //@}

  /**@name Get possible freed up row bounds - only valid after mode==3 */
  //@{
  /// Lower
  const double * relaxedRowLower() const;
  /// Upper
  const double * relaxedRowUpper() const;
  //@}

  /*! \name Behavioural Options
  
    This group of methods controls the type of cuts generated and other
    aspects of overall behaviour.
  */
  //@{
  /*! \brief Set cut generation mode

    The mode options are:

    - 0: Only unsatisfied integer variables will be looked at.  If no
	 information exists for that variable then probing will be done so as
	 a by-product you "may" get a fixing or infeasibility.  This will be
	 fast and is only available if a snapshot exists (otherwise as 1).
	 The bounds in the snapshot are the ones used.

    - 1: Look at unsatisfied integer variables, using current bounds.
	 Probing will be done on all looked at.

    - 2: Look at all integer variables, using current bounds.
	 Probing will be done on all.

    - 3: As mode 2, but will also return tightened row bounds when called
	 through #generateCutsAndModify.
    
    See #mode_ for additional detail.
  */
  void setMode(int mode);
  /// Get cut generation mode
  int getMode() const;

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
  void setRowCuts(int type);
  /// Get the row cut generation setting
  int rowCuts() const;

  /*! \brief Inquire if row cuts can be generated other than at the root node

     Returns true if cut generation may produce row cuts in the tree, as
     opposed to only at the root node. For CglProbing, this is equivalent to
     #rowCuts() > 0.
     Provided so that a client can know if the constraint matrix will change
     in a search tree.
  */
  virtual bool mayGenerateRowCutsInTree() const;

  /*! \brief Specify whether to use the objective as a constraint

      -  1: do add the objective as a constraint while probing
      -  0: don't add the objective as a constraint while probing
      - -1: don't use the objective as a constraint; further, don't compare
	    estimated degradation of the objective against a cutoff value
   
    Normally probing compares the estimated objective after forcing integrality
    to a cutoff value, when a cutoff is available. A value of -1 prevents this.
  */
  void setUsingObjective(int yesNo);
  /// Get use of objective
  int getUsingObjective() const;

  //@}

  /*! \name Generator Limits
  
    This group of methods sets limits on cut generator activity.
  */
  //@{
  /// Set maximum number of passes per node
  void setMaxPass(int value);
  /// Get maximum number of passes per node
  int getMaxPass() const;
  /// Set maximum number of passes per node  (root node)
  void setMaxPassRoot(int value);
  /// Get maximum number of passes per node (root node)
  int getMaxPassRoot() const;
  /*! \brief Set log level
  
    - 0: none
    - 1 - a bit
    - 2 - more details
  */
  void setLogLevel(int value);
  /// Get log level
  int getLogLevel() const;
  /// Set maximum number of unsatisfied variables to look at
  void setMaxProbe(int value);
  /// Get maximum number of unsatisfied variables to look at
  int getMaxProbe() const;
  /// Set maximum number of variables to look at in one probe
  /// Set maximum number of unsatisfied variables to look at (root node)
  void setMaxProbeRoot(int value);
  /// Get maximum number of unsatisfied variables to look at (root node)
  int getMaxProbeRoot() const;
  void setMaxLook(int value);
  /// Get maximum number of variables to look at in one probe
  int getMaxLook() const;
  /// Set maximum number of variables to look at in one probe (root node)
  void setMaxLookRoot(int value);
  /// Get maximum number of variables to look at in one probe (root node)
  int getMaxLookRoot() const;
  /// Set maximum number of elements in a row for it to be considered
  void setMaxElements(int value);
  /// Get maximum number of elements in a row for it to be considered
  int getMaxElements() const;
  /// Set maximum number of elements in row for it to be considered (root node)
  void setMaxElementsRoot(int value);
  /// Get maximum number of elements in row for it to be considered (root node)
  int getMaxElementsRoot() const;
  /// Set primal feasibility tolerance
  void setPrimalTolerance(double tol) { primalTolerance_ = tol ; }
  /// Get primal feasibility tolerance
  double getPrimalTolerance() { return (primalTolerance_) ; }
  //@}

  /**@name Get information back from probing */
  //@{
  /// Number looked at this time
  inline int numberThisTime() const
  { return numberThisTime_;}
  /// Which ones looked at this time
  inline const int * lookedAt() const
  { return lookedAt_;}
  //@}

  /**@name Mark which continuous variables are to be tightened */
  //@{
  /// Mark variables to be tightened
  void tightenThese(const OsiSolverInterface & solver, int number, const int * which);
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglProbing ();
 
  /// Copy constructor 
  CglProbing (
    const CglProbing &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglProbing &
    operator=(
    const CglProbing& rhs);
  
  /// Destructor 
  virtual
    ~CglProbing ();

  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  //@}
      
private:
  
 // Private member methods
  /**@name probe */
  //@{
  /// Does probing and adding cuts (without cliques and mode_!=0)
  int probe(const OsiSolverInterface &si, 
	    OsiCuts &cs, 
	    double *colLower, double *colUpper, CoinPackedMatrix *rowCopy,
	    CoinPackedMatrix *columnCopy, const CoinBigIndex *rowStartPos,
	    const int *realRow, const double *rowLower, const double *rowUpper,
	    const char *intVar, double *minR, double *maxR, int *markR, 
	    CglTreeInfo *info) const;
  /// Does probing and adding cuts (with cliques)
  int probeCliques( const OsiSolverInterface & si, 
	     const OsiRowCutDebugger * debugger, 
	     OsiCuts & cs, 
	     double * colLower, double * colUpper, CoinPackedMatrix *rowCopy,
		    CoinPackedMatrix *columnCopy, const int * realRow,
	     double * rowLower, double * rowUpper,
	     char * intVar, double * minR, double * maxR, int * markR, 
             CglTreeInfo * info) const;
  /// Does probing and adding cuts for clique slacks
  int probeSlacks( const OsiSolverInterface & si, 
                    const OsiRowCutDebugger * debugger, 
                    OsiCuts & cs, 
                    double * colLower, double * colUpper, CoinPackedMatrix *rowCopy,
		   CoinPackedMatrix *columnCopy,
                    double * rowLower, double * rowUpper,
                    char * intVar, double * minR, double * maxR,int * markR,
                     CglTreeInfo * info) const;
  /** Does most of work of generateCuts 
      Returns number of infeasibilities */
  int gutsOfGenerateCuts( const OsiSolverInterface & si, 
			  OsiCuts & cs,
			  double * rowLower, double * rowUpper,
			  double * colLower, double * colUpper,
                           CglTreeInfo * info) const;
  /// Sets up clique information for each row
  void setupRowCliqueInformation(const OsiSolverInterface & si);
  /*! \brief Tighten column bounds
  
     Use bound propagation to tighten column bounds. Can declare infeasibility
     and  may also declare rows to be redundant
  */
  int tighten(double *colLower, double * colUpper,
              const int *column, const double *rowElements, 
              const CoinBigIndex *rowStart,const CoinBigIndex * rowStartPos,
	      const int * rowLength,
              double *rowLower, double *rowUpper, 
              int nRows,int nCols,char * intVar,int maxpass,
              double tolerance) const;
  /*! \brief Tighten column bounds using clique information
  
     Use bound propagation augmented with clique information to tighten column
     bounds. Can declare infeasibility and  may also declare rows to be
     redundant
  */
  int tightenClique(double *colLower, double * colUpper,
              const int *column, const double *rowElements, 
              const CoinBigIndex *rowStart,const CoinBigIndex * rowStartPos,
	      const int * rowLength,
              double *rowLower, double *rowUpper, 
              int nRows,int nCols,char * intVar,int maxpass,
              double tolerance) const;

  /*! \brief Calculate constraint left-hand-side bounds

    Calculate lower and upper bounds on the constraint lhs based on variable
    lower and upper bounds. The new bounds are loaded into minR and maxR.
    The array markR will be set as:
    - -1: constraint is useful for further propagation (at least one lhs
          bound finite)
    - -2: constraint is not useful (both lhs bounds infinite).
    In practice, finite translates as `reasonable magnitude'.
  */
  void calcRowBounds(double *colLower, double *colUpper,
		     const int *column, const double *rowElements, 
		     const CoinBigIndex *rowStart, const int *rowLength,
		     double *rowLower, double *rowUpper, 
		     double *minR, double *maxR, int *markR,
		     int nRows) const;
  //@}

  // Private member data

  /**@name Private member data */
  //@{

  struct disaggregation_struct_tag ;
  friend struct CglProbing::disaggregation_struct_tag ;

  /*! \brief Row copy
  
    This is created only as part of a snapshot.

    \todo The more I browse in CglProbing, the more it seems to me that
	  it'd be a good idea to create a snapshot object. Then we'd have
	  one mutable pointer here and could simplify parameter passing
	  in several places in the code.  -- lh, 100917 --
  */
  CoinPackedMatrix * rowCopy_;
  /*! \brief Column copy
  
    This is created only as part of a snapshot. More striking, as far as I can
    tell it's never used! (CglProbing indeed uses a column copy, but it's
    created and destroyed in gutsOfGenerateCut.) (lh, 101021)
  */
  CoinPackedMatrix * columnCopy_;

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
    These don't have to be mutable because the only time probing tightens
    row bounds is when called as generateCutsAndModify, a non-const method.
  */
  /// Lower bounds on rows
  double * rowLower_;
  /// Upper bounds on rows
  double * rowUpper_;

  /*
    These have to be mutable because generateCuts is a const method, and
    tightening bounds is one of the main activities of generateCuts.
  */ 
  /// Lower bounds on columns
  mutable double * colLower_;
  /// Upper bounds on columns
  mutable double * colUpper_;

  /*
    Not entirely sure why these need to be mutable. Check. -- lh, 101007 --
  */
  /// Number of rows in snapshot (or when cliqueRow stuff computed)
  mutable int numberRows_;
  /// Number of columns in problem ( must == current)
  mutable int numberColumns_;

  /// Tolerance to see if infeasible
  double primalTolerance_;

  /*! \brief cut generation mode

    The basic mode is the first hex digit. Higher bits are variations.
  
    - 0x00: lazy using snapshot
    - 0x01: just unsatisfied integer variables
    - 0x02: all integer variables
    - 0x03: as 0x02, but when called as generateCutsAndModify will
	    return tightened row bounds as well as column bounds.
    - 0x10: set if want to extend cliques at root node
  */
  int mode_;

  /*! \brief Row cuts flags

    - 0x0: do not generate row cuts
    - 0x1: generate disaggregation row cuts
    - 0x2: generate coefficient (strengthening) row cuts
    - 0x4: generate column cuts only

    A negative value -n is interpreted as +n but only fixes variables unless
    at the root node
  */
  mutable int rowCuts_;

  /// Maximum number of passes to do in probing
  int maxPass_;
  /// Log level - 0 none, 1 - a bit, 2 - more details
  int logLevel_;
  /// Maximum number of unsatisfied variables to probe
  int maxProbe_;
  /// Maximum number of variables to look at in one probe
  int maxStack_;
  /// Maximum number of elements in row for scan
  int maxElements_;
  /// Maximum number of passes to do in probing at root
  int maxPassRoot_;
  /// Maximum number of unsatisfied variables to probe at root
  int maxProbeRoot_;
  /// Maximum number of variables to look at in one probe at root
  int maxStackRoot_;
  /// Maximum number of elements in row for scan at root
  int maxElementsRoot_;
  /// Whether to include the objective as a constraint
  int usingObjective_;
  /// Number of integer variables
  int numberIntegers_;
  /// Number of 0-1 integer variables
  int number01Integers_;
  /// Number looked at this time
  mutable int numberThisTime_;
  /*! \brief Total number of times called */
  mutable int totalTimesCalled_;
  /// Which ones looked at this time
  mutable int * lookedAt_;

  /// For disaggregation cuts and for building cliques
  typedef struct disaggregation_struct_tag {
    /// Integer variable column number
    int sequence;
    /// length of new value
    int length;
    /// columns whose bounds will be changed
    disaggregationAction *index;
  } disaggregation;

  /*! \brief Disaggregation cuts

    Initialised when probing starts, or when a snapshot is created.
  */
  disaggregation * cutVector_;

  /************************  Cliques  ************************/
  /*
    Cliques detected by CglProbing are always a single row and have a
    very stylised form. See the description with cliqueEntry in
    CglTreeInfo.hpp.
  */

  /// Number of cliques
  int numberCliques_;

  /*! \brief Clique type

    Currently a single flag indicating if the relevant row is an equality.

    \todo See if there's any indication that this will become more complex
    in the future. If so, it's good to leave it as a defined type. If this
    is all there'll ever be, make it a boolean. Or a char, which would give
    7 bits for future expansion. Almost certainly a single-bit field is extra
    work.
    -- lh, 100924 --
  */
  typedef struct {
    unsigned int equality:1; //  nonzero if clique row is an equality
  } cliqueType;

  /*! \brief Clique type
  
    See the description with the cliqueType structure.
  */
  cliqueType *cliqueType_;

  /*! \brief Start of each clique

    Each entry points to the start of the block of entries for this clique
    in #cliqueEntry_.
  */
  int *cliqueStart_;

  /*! \brief Clique membership

    Contains one block of entries for each clique, specifying the variables
    in the clique. See cliqueEntry for a description of the information in each
    entry.
  */
  cliqueEntry *cliqueEntry_;

  /*! \brief Start of strong-1 block for a variable

    One entry per variable. An entry contains the index of the strong-1
    subblock for this variable in #whichClique_, or -1 if the variable is
    not in a clique.
  */
  int * oneFixStart_;

  /*! \brief Start of strong-0 block for a variable

    One entry per variable. An entry contains the index of the strong-0
    subblock for this variable in #whichClique_, or -1 if the variable is
    not in a clique.
  */
  int * zeroFixStart_;

  /*! \brief End of block for a variable

    One entry per variable. An entry contains the index one past the end of
    the block for a variable in #whichClique_. Note that it's not necessarily
    true that oneFixStart_[j] is one past the end of the block for x<j-1>
    because there's no guarantee x<j> is in any clique; oneFixStart_[j] could
    well be -1. So we need an array to specify the end of each block.
  */
  int * endFixStart_;

  /*! \brief Cross-reference variables back to cliques

    Divided into blocks for each variable; each block is divided into a
    strong-1 subblock (see #oneFixStart_) and a strong-0 subblock (see
    #zeroFixStart_). Each entry is the index of a clique. See also
    #endFixStart_.
  */
  int * whichClique_;

  /*! \brief Preprocessed clique information for use when calculating row
	     bounds.
  
    See the explanation with #setupRowCliqueInformation. For each row, there's
    a block of entries, in the order in which columns are referenced in the
    row. Each entry specifies the strong direction for the variable, and a
    locally-valid clique id that identifies variables in the row that are part
    of the same clique.

    When is this valid? That's a bit unclear. Created by
    #setupRowCliqueInformation at the end of #generateCutsAndModify, so only
    useable on subsequent passes. Seems to be a bit of a work in progress.
  */
  cliqueEntry * cliqueRow_;
  /// cliqueRow_ starts for each row
  int * cliqueRowStart_;
  /// If not null and [i] !=0 then also tighten even if continuous
  char * tightenBounds_;
  //@}
};
inline int affectedInDisaggregation(const disaggregationAction & dis)
{ return dis.affected&0x1fffffff;}
inline void setAffectedInDisaggregation(disaggregationAction & dis,
					   int affected)
{ dis.affected = affected|(dis.affected&0xe0000000);}
#ifdef NDEBUG
inline bool zeroOneInDisaggregation(const disaggregationAction & )
{ return true;}
#else
inline bool zeroOneInDisaggregation(const disaggregationAction & dis)
//{ return (dis.affected&0x80000000)!=0;}
{ assert ((dis.affected&0x80000000)!=0); return true;}
#endif
inline void setZeroOneInDisaggregation(disaggregationAction & dis,bool zeroOne)
{ dis.affected = (zeroOne ? 0x80000000 : 0)|(dis.affected&0x7fffffff);}
inline bool whenAtUBInDisaggregation(const disaggregationAction & dis)
{ return (dis.affected&0x40000000)!=0;}
inline void setWhenAtUBInDisaggregation(disaggregationAction & dis,bool whenAtUB)
{ dis.affected = (whenAtUB ? 0x40000000 : 0)|(dis.affected&0xbfffffff);}
inline bool affectedToUBInDisaggregation(const disaggregationAction & dis)
{ return (dis.affected&0x20000000)!=0;}
inline void setAffectedToUBInDisaggregation(disaggregationAction & dis,bool affectedToUB)
{ dis.affected = (affectedToUB ? 0x20000000 : 0)|(dis.affected&0xdfffffff);}

//#############################################################################
/** A function that tests the methods in the CglProbing class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
int CglProbingUnitTest(const OsiSolverInterface *siP,
		       const std::string mpsDir);

#endif
