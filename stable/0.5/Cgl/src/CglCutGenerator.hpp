// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglCutGenerator_H
#define CglCutGenerator_H

#include "OsiCuts.hpp"
#include "OsiSolverInterface.hpp"

/** Information about where the cut generator is invoked from. */

struct CglTreeInfo {
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
  CglTreeInfo() : level(-1), pass(-1), formulation_rows(-1), inTree(false),
                   strengthenRow(NULL) {}
};

//-------------------------------------------------------------------
//
// Abstract base class for generating cuts.
//
//-------------------------------------------------------------------
///
/** Cut Generator Base Class

This is an abstract base class for generating cuts.  A specific cut 
generator will inherit from this class.
*/
class CglCutGenerator  {
  
public:
    
  /**@name Generate Cuts */
  //@{
  /** Generate cuts for the model data contained in si.
  The generated cuts are inserted into and returned in the
  collection of cuts cs.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo()) const=0; 
  //@}

    
  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglCutGenerator (); 
 
  /// Copy constructor 
  CglCutGenerator ( const CglCutGenerator &);

  /// Clone
  virtual CglCutGenerator * clone() const = 0;

  /// Assignment operator 
  CglCutGenerator & operator=(const CglCutGenerator& rhs);

  /// Destructor 
  virtual ~CglCutGenerator ();

  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp) {return "";};

  /// This can be used to refresh any inforamtion
  virtual void refreshSolver(OsiSolverInterface * solver) {};
  //@}
  
  /**@name Gets and Sets */
  //@{
  /**
     Get Aggressiveness - 0 = neutral, 100 is normal root node.
     Really just a hint to cut generator
  */
  inline int getAggressiveness() const
  { return aggressive_;};

  /**
     Set Aggressiveness - 0 = neutral, 100 is normal root node.
     Really just a hint to cut generator
  */
  inline void setAggressiveness(int value)
  { aggressive_=value;};

  /**
     Returns true if may generate Row cuts in tree (rather than root node).
     Used so know if matrix will change in tree.  Really
     meant so column cut generators can still be active
     without worrying code.
     Default is true
  */
  virtual bool mayGenerateRowCutsInTree() const;
  /// Return true if needs optimal basis to do cuts
  virtual bool needsOptimalBasis() const;
  //@}
  
  // test this class
  //static void unitTest();
  
// private:
  
  /**
     Aggressiveness - 0 = neutral, 100 is normal root node.
     Really just a hint to cut generator
  */
  int aggressive_;
};


#endif
