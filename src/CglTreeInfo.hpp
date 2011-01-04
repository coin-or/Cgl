// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglTreeInfo_H
#define CglTreeInfo_H

#include "OsiCuts.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinHelperFunctions.hpp"
class CglStored;
/*! \brief Context for an invocation of a cut generator

  This class is intended to provide the cut generator with some information
  about the context in which it is invoked (<i>e.g.</i>, preprocessing,
  evaluation of the root node, or evaluation of a non-root node). It provides
  an assortment of other information of general utility.

  Where a cut generator requires additional specialised information, a derived
  class may be appropriate (<i>cf.</i> CglTreeProbingInfo).
*/

class CglTreeInfo {
public:
  /// The level of the search tree node 
  int level ;

  /*! \brief Number of times the cut generator has been invoked at this
  	     search tree node
  */
  int pass ;

  /*! \brief Number of rows in the original formulation.
  
    Some generators may not want to consider already generated rows when
    generating new ones.
  */
  int formulation_rows ;

  /*! \brief Generic options.

      -  0x01 - treat costed integers as important
      -  0x02 - switch off some stuff as variables semi-integer
      -  0x04 - set global cut flag if at root node
      -  0x08 - set global cut flag if at root node and first pass
      -  0x10 - set global cut flag and make cuts globally valid
      -  0x20 - last round of cuts did nothing - maybe be more aggressive
      -  0x40 - in preprocessing stage
      -  0x80 - looks like solution
  */
  int options ;

  /// Set true if in tree (to avoid ambiguity at first branch)
  bool inTree ;

  /*! \brief Row replacement array.
  
    Before Branch and Cut it may be beneficial to strengthen (replace) rows
    rather than adding cuts. If this array is not NULL then the cut generator
    can place a pointer to the stronger cut in this array which is number of
    rows in size.  The calling function can then replace those rows.

    A null cut (i.e. zero elements and free rhs) indicates that the row is
    useless and can be removed.
  */
  OsiRowCut ** strengthenRow;

  /// Optional pointer to thread specific random number generator
  CoinThreadRandom * randomNumberGenerator;

  /// Default constructor 
  CglTreeInfo ();
 
  /// Copy constructor 
  CglTreeInfo (const CglTreeInfo &);

  /// Clone
  virtual CglTreeInfo *clone() const;

  /// Assignment operator 
  CglTreeInfo &operator=(const CglTreeInfo& rhs);
  
  /// Destructor 
  virtual ~CglTreeInfo ();

  /*! \brief Take action if the cut generator can fix a variable

    jjf: (toValue -1 for down, +1 for up)

    It's unclear why this is here, except that it conveniently solves some
    coding issues related to nonuse of CglTreeProbingInfo. Really should be
    removed. -- lh, 101127 --
  */
  virtual bool fixes(int , int , int ,bool) {return false;}
  
  /*! \brief Initializes fixing arrays
  
    jjf: Returns > 0 if we want to save info, 0 if we don't, and -1 if is to be
    used (?)

    It's unclear why this is here, except that it conveniently solves some
    coding issues related to nonuse of CglTreeProbingInfo. Really should be
    removed. -- lh, 101127 --
  */
  virtual int initializeFixing(const OsiSolverInterface * ) {return 0;}
  
};

/*! \brief Utility structure to encode clique information.

  The only purpose is to hide the details of the encoding. The msb of the word
  is set to 1 if setting the variable to 1 fixes all other variables in the
  clique, 0 if setting the variable to 0 fixes all other variables in the
  clique. The remaining 31 bits are the index of the variable.
*/
typedef struct {
  //unsigned int oneFixed:1; //  nonzero if variable to 1 fixes all
  //unsigned int sequence:31; //  variable (in matrix) (but also see cliqueRow_)
  unsigned int fixes;
} cliqueEntry;

/*! \brief Specialisation of CglTreeInfo for the CglProbing cut generator.

  It's unclear that this is actually in use. The relevant guard symbol in
  CglProging, PROBING100, has been set to a disabling value since 080722. In
  cbc, code relevant to CglTreeProbingInfo is also mostly disabled.
  -- lh, 101127 --
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
inline int sequenceInCliqueEntry(const cliqueEntry & cEntry)
{ return cEntry.fixes&0x7fffffff;}
inline void setSequenceInCliqueEntry(cliqueEntry & cEntry,int sequence)
{ cEntry.fixes = sequence|(cEntry.fixes&0x80000000);}
inline bool oneFixesInCliqueEntry(const cliqueEntry & cEntry)
{ return (cEntry.fixes&0x80000000)!=0;}
inline void setOneFixesInCliqueEntry(cliqueEntry & cEntry,bool oneFixes)
{ cEntry.fixes = (oneFixes ? 0x80000000 : 0)|(cEntry.fixes&0x7fffffff);}

#endif
