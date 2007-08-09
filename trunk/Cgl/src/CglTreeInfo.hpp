// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef CglTreeInfo_H
#define CglTreeInfo_H

#include "OsiCuts.hpp"
#include "OsiSolverInterface.hpp"

/** Information about where the cut generator is invoked from. */

class CglTreeInfo {
public:
  /// The level of the search tree node 
  int level;
  /** How many times the cut generator was already invoked in this search tree
      node */
  int pass;
  /** The number of rows in the original formulation. Some generators may not
      want to consider already generated rows when generating new ones. */
  int formulation_rows;
  /// Set true if in tree (to avoid ambiguity at first branch)
  bool inTree;
  /** Replacement array.  Before Branch and Cut it may be beneficial to strengthen rows
      rather than adding cuts.  If this array is not NULL then the cut generator can
      place a pointer to the stronger cut in this array which is number of rows in size.

      A null (i.e. zero elements and free rhs) cut indicates that the row is useless 
      and can be removed.

      The calling function can then replace those rows.
  */
  OsiRowCut ** strengthenRow;
  /// Default constructor 
  CglTreeInfo ();
 
  /// Copy constructor 
  CglTreeInfo (
    const CglTreeInfo &);
  /// Clone
  virtual CglTreeInfo * clone() const;

  /// Assignment operator 
  CglTreeInfo &
    operator=(
    const CglTreeInfo& rhs);
  
  /// Destructor 
  virtual
    ~CglTreeInfo ();
  /// Take action if cut generator can fix a variable (toValue -1 for down, +1 for up)
  virtual void fixes(int variable, int toValue, int fixedVariable,double fixedToValue) {}
  /// Initalizes fixing arrays etc - returns true if we want to save info
  virtual bool initializeFixing() {return false;}
  
};

/** Derived class to pick up probing info. */
typedef struct {
  unsigned int oneFixed:1; //  nonzero if variable fixed to 1
  unsigned int sequence:31; //  variable (in matrix)
} fixEntry;

class CglTreeProbingInfo : public CglTreeInfo {
public:
  /// Default constructor 
  CglTreeProbingInfo ();
  /// Constructor from model
  CglTreeProbingInfo (const OsiSolverInterface * model);
 
  /// Copy constructor 
  CglTreeProbingInfo (
    const CglTreeProbingInfo &);
  /// Clone
  virtual CglTreeInfo * clone() const;

  /// Assignment operator 
  CglTreeProbingInfo &
    operator=(
    const CglTreeProbingInfo& rhs);
  
  /// Destructor 
  virtual
    ~CglTreeProbingInfo ();
  /// Take action if cut generator can fix a variable (toValue -1 for down, +1 for up)
  virtual void fixes(int variable, int toValue, int fixedVariable,double fixedToValue);
  /// Initalizes fixing arrays etc - returns true if we want to save info
  virtual bool initializeFixing() ;
  /// Entries for fixing variables
  inline fixEntry * fixEntries() const
  { return fixEntry_;}
  /// Starts of integer variable going to zero
  inline int * toZero() const
  { return toZero_;}
  /// Starts of integer variable going to one
  inline int * toOne() const
  { return toOne_;}
  /// List of 0-1 integer variables
  inline int * integerVariable() const
  { return integerVariable_;}
  /// Backward look up
  inline int * backward() const
  { return backward_;}
  /// Number of variables
  inline int numberVariables() const
  { return numberVariables_;}
  /// Number of 0-1 variables
  inline int numberIntegers() const
  { return numberIntegers_;}
protected:
  /// Entries for fixing variables
  fixEntry * fixEntry_;
  /// Starts of integer variable going to zero
  int * toZero_;
  /// Starts of integer variable going to one
  int * toOne_;
  /// List of 0-1 integer variables
  int * integerVariable_;
  /// Backward look up
  int * backward_;
  /// Number of variables
  int numberVariables_;
  /// Number of 0-1 variables
  int numberIntegers_;
  /// Maximum number in fixEntry_
  int maximumEntries_;
  /// Last integer done in pass
  int lastInteger_;
};

#endif
