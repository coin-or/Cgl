// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglDuplicateRow_H
#define CglDuplicateRow_H

#include <string>

#include "CglCutGenerator.hpp"
class CglStored;

/*! \brief DuplicateRow Cut Generator Class

  This is a very simple minded idea but I (JJF) am using it in a project
  so thought I might as well add it.  It should really be called before
  first solve and I may modify CBC to allow for that.

  This is designed for problems with few rows and many integer variables
  where the rhs are <= or == and all coefficients and rhs are small integers.

  If effective rhs is K then we can fix all variables with coefficients > K
  to their lower bounds (effective rhs just means original with variables with
  nonzero lower bounds subtracted out).

  If one row is a subset of another and the effective rhs are same we can
  fix some variables and then the two rows are identical.

  The generator marks identical rows so can be taken out in solve.
*/
class CglDuplicateRow : public CglCutGenerator {
 
public:
  
  /*! \name Generate Cuts */
  //@{
  /// Look for duplicate/dominated rows and variables that can be fixed.
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo()) const;

  /*! \brief Fix variables and find duplicate/dominated rows

    This version performs the deletions and fixings and may return stored
    cuts for dominated columns 
  */
  CglStored * outDuplicates( OsiSolverInterface * solver);

  //@}

public:

  /**@name Get problem information */
  //@{
  /*! \brief Get the duplicate row list.
  
    Entries in the list are coded as:
    - -1: still in
    - -2: out (all fixed)
    - -3: not clique (mode 0x08 only?)
    -  k >= 0: same as row k 
  */
  inline const int * duplicate() const { return duplicate_;}

  /// Size of dynamic program
  inline int sizeDynamic() const { return sizeDynamic_;}

  /// Number of rows in original problem
  inline int numberOriginalRows() const { return matrix_.getNumRows();}
  //@}

  /*! \name Control methods */
  //@{

  /// Get mode
  inline int mode() const { return mode_;}

  /*! \brief Set mode

    - 0x01: check for duplicate rows
    - 0x02: check for duplicate columns
    - 0x04: undocumented
    - 0x08: undocumented
  */
  inline void setMode(int value) { mode_=value;}

  /// Get right-hand-side magnitude limit
  inline int maximumRhs() const { return maximumRhs_;}

  /*! \brief Set right-hand-side magnitude limit

    We only check for duplicates amongst rows where the magnitude of the
    effective rhs is less than the given value.
  */
  inline void setMaximumRhs(int value) { maximumRhs_=value;}

  /// Get column size limit
  inline int maximumDominated() const { return maximumDominated_;}

  /*! \brief Set column size limit

    We only check for dominated columns amongst groups of columns whose size
    is less than the given value.
  */
  inline void setMaximumDominated(int value) { maximumDominated_=value;}

  /// Get log level
  inline int logLevel() const { return logLevel_;}

  /// Set log level
  inline void setLogLevel(int value) { logLevel_ = value;}
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglDuplicateRow ();
 
  /// Useful constructor 
  CglDuplicateRow (OsiSolverInterface * solver);
 
  /// Copy constructor 
  CglDuplicateRow (
    const CglDuplicateRow & rhs);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglDuplicateRow &
    operator=(
    const CglDuplicateRow& rhs);
  
  /// Destructor 
  virtual
    ~CglDuplicateRow ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);

  /// This can be used to refresh any information
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}
      
protected:
  

  // Protected member data

  /**@name Protected member data */
  //@{
  /// Matrix
  CoinPackedMatrix matrix_;
  /// Matrix by row
  CoinPackedMatrix matrixByRow_; 
  /// Possible rhs (if 0 then not possible)
  int * rhs_;
  /// Marks duplicate rows
  mutable int * duplicate_;
  /// To allow for <= rows
  int * lower_;
  /// Stored cuts if we found dominance cuts
  mutable CglStored * storedCuts_;
  /// Check dominated columns if less than this number of candidates
  int maximumDominated_;
  /// Check duplicates if effective rhs <= this
  int maximumRhs_;
  /// Size of dynamic program
  mutable int sizeDynamic_;
  /// 1 do rows, 2 do columns, 3 do both
  int mode_;
  /// Controls print out
  int logLevel_;
  //@}

private:

  /*! \name Cut generation worker methods */
  //@{
  /// Does work for modes 1,2
  void generateCuts12( const OsiSolverInterface & si, OsiCuts & cs,
		       const CglTreeInfo info = CglTreeInfo()) const;
  /// Does work for mode 4
  void generateCuts4( const OsiSolverInterface & si, OsiCuts & cs,
		       const CglTreeInfo info = CglTreeInfo()) const;
  /// Does work for mode 8
  void generateCuts8( const OsiSolverInterface & si, OsiCuts & cs,
		       const CglTreeInfo info = CglTreeInfo()) const;
  //@}
};
#endif
