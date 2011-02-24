// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglTreeProbingInfo_H
#define CglTreeProbingInfo_H

#include "CoinPackedMatrix.hpp"
#include "OsiCuts.hpp"
#include "CglTreeInfo.hpp"
#include "CglGeneralisedClique.hpp"


/*! \brief Specialisation of CglTreeInfo for the CglProbing cut generator.

  It's unclear that this is actually in use. The relevant guard symbol in
  CglProbing, PROBING100, has been set to a disabling value since 080722. In
  cbc, code relevant to CglTreeProbingInfo is also mostly disabled.
  -- lh, 101127 --

  Well, not entirely, but I'm working on it. As best I can tell, the
  functionality supported by CglTreeProbingInfo --- generalised cliques ---
  was never fully functional, and correspondence with John confirms this.
  I'm working on lopping it out.
  -- lh, 110224 --
*/
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
  OsiSolverInterface * analyze(const OsiSolverInterface & si, int createSolver=0);
  /** Take action if cut generator can fix a variable 
      (toValue -1 for down, +1 for up)
      Returns true if still room, false if not  */
  virtual bool fixes(int variable, int toValue, int fixedVariable,bool fixedToLower);
  /*! \brief Initializes fixing arrays
  
    Erases existing implication information, if any, and initialises the
    object with binary variable information.

    jjf: Returns > 0 if we want to save info, 0 if we don't, and -1 if is to be
    used (?)

    The above doesn't seem to match reality. The method returns these values:
    -  2: arrays already initialised
    -  1: arrays initialised by this call.
    -  0: unused
    - -1: unused
    - -2: `correct style'? Whatever it is, it's apparently imposed by
          executing #convert() and indicates that the object was not
	  reinitialised.
  */
  virtual int initializeFixing(const OsiSolverInterface * model) ;
  /// Fix entries in a solver using implications
  int fixColumns(OsiSolverInterface & si) const;
  /// Fix entries in a solver using implications for one variable
  int fixColumns(int iColumn, int value, OsiSolverInterface & si) const;
  /// Packs down entries
  int packDown();
  /// Generate cuts from implications
  void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
		    const CglTreeInfo info) const;
  /// Entries for fixing variables
  inline cliqueEntry * fixEntries() const
  { convert(); return fixEntry_;}
  /// Starts of integer variable going to zero
  inline int * toZero() const
  { convert(); return toZero_;}
  /// Starts of integer variable going to one
  inline int * toOne() const
  { convert(); return toOne_;}
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
private:
  /// Converts to ordered
  void convert() const;
protected:
  /// Entries for fixing variables
  mutable cliqueEntry * fixEntry_;
  /// Starts of integer variable going to zero
  mutable int * toZero_;
  /// Starts of integer variable going to one
  mutable int * toOne_;
  /// List of 0-1 integer variables
  int * integerVariable_;
  /// Backward look up
  int * backward_;
  /// Entries for fixing variable when collecting
  mutable int * fixingEntry_;
  /// Number of variables
  int numberVariables_;
  /// Number of 0-1 variables
  int numberIntegers_;
  /// Maximum number in fixEntry_
  int maximumEntries_;
  /// Number entries in fixingEntry_ (and fixEntry_) or -2 if correct style
  mutable int numberEntries_;
};

#endif
