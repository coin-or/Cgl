
/*
  Copyright (C) 2010, International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
#ifndef CglProbingRowCut_H
#define CglProbingRowCut_H

#include <OsiCuts.hpp>

/*
  \brief Augmented row cut class to accumulate probing cuts.

  A collection of OsiRowCut objects, specialised for use in CglProbing.
  It supports fast identification of duplicate cuts and some bookkeeping
  that's useful in CglProbing.

  The value given to the constructor for \p nRows will be used to to set
  the capacity of the collection, based on a complicated calculation that
  tries to come up with something that's unlikely to overflow during the
  probe. Expansion is costly because the hash table must be expanded and
  repopulated.  The actual capacity will have no obvious relation to \p
  nRows. Check the code for the constructor if this concerns you.
*/
class CglProbingRowCut {

  public:

  /*!  \brief Default constructor

    The value given for \p nRows will be used to set #nRows_. It's also
    used in a magic calculation to set the capacity of the collection. See
    the notes for the class for details.  Setting \p initialPass to true
    will reduce the size of the hash table by half, in the context of a
    similar magic calculation.
  */
  CglProbingRowCut(int nRows, bool initialPass) ;

  /// Destructor
  ~CglProbingRowCut() ;

  inline OsiRowCut *cut (int i) const { return rowCut_[i] ; }

  /// The number of cuts in the collection.
  inline int numberCuts() const { return numberCuts_ ; }

  /// Check if space remains to add more cuts.
  inline bool outOfSpace() const { return (maxSize_ == numberCuts_) ; }

  /*! \brief Add a cut if it's not a duplicate
  
    Add a cut if it's not a duplicate of a cut already in the collection.

    Return value:
     - -1 if the cut is not added because there's no space,
     -  0 if the cut is added,
     -  1 if the cut is not added because it's a duplicate
  */
  int addCutIfNotDuplicate(OsiRowCut &cut) ;

  /*! \brief Transfer cuts into the given OsiCuts collection, updating
  	     whichRow.

    Transfer cuts from this collection into \p cs, deleting all cuts from this
    collection when finished. Nominally, up to #nRows_ cuts are transferred
    out but this limit can be exceeded if the cuts straddling the limit have
    equal effectiveness.

    During the transfer, if a cut has an associated row index, \p whichRow is
    checked. If the entry for the row has no associated cut, a pointer to the
    cut is stored. The \p iPass parameter gives some randomness in lieu of
    sorting when the number of cuts in the container is less than the transfer
    limit.
  */
  void addCuts(OsiCuts &cs, OsiRowCut **whichRow, int iPass) ;

  private:

  /*! \brief find an entry for a cut in the hash table

    If the cut is already in the hash table, returns the index of the entry.
    If the cut is not in the hash table, returns the index of an appropriate
    empty entry. If there's no space, returns -1.
  */
    int findHashTableEntry (const OsiRowCut &newCut) ;

  /*! \brief hash table entry
  
    index is the cut index, next is the next hash table entry for this
    hash value
  */

  typedef struct { int index ;
  		   int next ; } HashLink ;

  /// Fudge factor used by constructor to calculate #maxSize_
  static const int SIZE_ROW_MULT = 4 ;
  /// Another fudge factor used by constructor to calculate #maxSize_
  static const int SIZE_ROW_ADD = 2000 ;

  /// The number of cuts in the collection
  int numberCuts_ ;
  /// The maximum number of cuts that can be held in the collection.
  int maxSize_ ;
  /// The initial size of the cut storage vector
  int size_ ;
  /// The cut storage vector
  OsiRowCut **rowCut_ ;
  /// Nominal number of cuts transferred out by #addCuts
  int nRows_ ;

  /// Size of hash table
  int hashSize_ ;
  /// Hash table for fast duplicate check
  HashLink *hash_ ;
  /*! \brief Most recently filled hash table entry
  
    And the starting point to search for a new empty entry when there's a hash
    collision.
  */
  int lastHash_ ;

} ;

#endif
