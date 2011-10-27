// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglTreeInfo_H
#define CglTreeInfo_H

#include "CoinHelperFunctions.hpp"

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
      256 - want alternate cuts
      512 - in sub tree (i.e. parent model)
      1024 - in must call again mode or after everything mode
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

    The only reason for this method is that it sidesteps coding issues
    related to nonuse of CglTreeProbingInfo. Really should be removed.
    -- lh, 110224 --
  */
  virtual bool fixes(int, int, int, bool) { return false ; }
  
  /*! \brief Initializes fixing arrays
  
    jjf: Returns > 0 if we want to save info, 0 if we don't, and -1 if is to be
    used (?)

    The only reason for this method is that it sidesteps coding issues
    related to nonuse of CglTreeProbingInfo. Really should be removed.
    -- lh, 110224 --
  */
  virtual int initializeFixing(const OsiSolverInterface*) {return (0) ; }
  
} ;

#endif
