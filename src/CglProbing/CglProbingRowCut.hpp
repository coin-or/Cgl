
/*
  Copyright (C) 2010, International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/
#ifndef CglProbingRowCut_H
#define CglProbingRowCut_H

/*
  TODO: What does this provide that's better / distinct from OsiCuts?
	There's a cryptic comment down in gutsOfGenerateCuts, `keep cuts out
	of cs until end so we can find duplicates quickly', but that doesn't
	quite justify a new class.  -- lh, 100923 --

  TODO: This should be a separate file. It's about 200 lines.
	-- lh, 100923 --

  TODO: This should be renamed so that it follows the usual conventions for a
	class name.  -- lh, 100923 --

  Reading a bit more, it would appear that this class has two reasons for
  existence: OsiCuts is more general, handling row and column cuts, and it
  uses OsiRowCut, instead of OsiRowCut2 used here. OsiRowCut2 adds a
  reference (simple integer index) to the original constraint. Remains to be
  seen whether this really justifies a new class or whether it'd be better to
  add the reference to the original constraint in OsiRowCut.  -- lh, 101125 --

  Separated and renamed CglProbingRowCut  -- lh, 101202 --
*/

#include <OsiCuts.hpp>

class CglProbingRowCut {

  public:

  CglProbingRowCut(int nRows, bool initialPass) ;

  ~CglProbingRowCut() ;

  inline OsiRowCut2 *cut (int i) const { return rowCut_[i] ; }

  inline int numberCuts() const { return numberCuts_ ; }

  inline bool outOfSpace() const { return (maxSize_ == numberCuts_) ; }

  // Return 0 if added, 1 if not, -1 if not added because of space
  int addCutIfNotDuplicate(OsiRowCut &cut, int whichRow = -1) ;

  void addCuts(OsiCuts &cs, OsiRowCut **whichRow, int iPass) ;

  private:

  typedef struct { int index, next ; } HashLink;

  OsiRowCut2 **rowCut_ ;
  /// Hash table
  HashLink *hash_ ;
  int size_ ;
  int maxSize_ ;
  int hashSize_ ;
  int nRows_ ;
  int numberCuts_ ;
  int lastHash_ ;

} ;

#endif
