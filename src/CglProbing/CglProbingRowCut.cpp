
/*
  Copyright (C) 2002, International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinHelperFunctions.hpp"
#include "CglProbingRowCut.hpp"

#include "CglProbingDebug.hpp"

/* \file CglProbingRowCut.cpp
   \brief Implementation of the CglProbingRowCut class.
*/

/*
  Hash functions.
  
  TODO: Why is there no CoinUtils hash package?  -- lh, 101202 --
*/
namespace {

const double multiplier[] = { 1.23456789e2, -9.87654321 } ;

/*
  Calculate a hash value for the cut, modulo the size of the hash table.
  Note that the hash value is sensitive to coefficient order.
*/
int hashCut (const OsiRowCut &x, int size) 
{
  int xN = x.row().getNumElements() ;
  double xLb = x.lb() ;
  double xUb = x.ub() ;
  const int *xIndices = x.row().getIndices() ;
  const double *xElements = x.row().getElements() ;
  double value = 1.0 ;
  if (xLb > -1.0e10)
    value += xLb*multiplier[0] ;
  if (xUb < 1.0e10)
    value += xUb*multiplier[1] ;
  // This calculation is order-sensitive
  for (int j = 0 ; j < xN ; j++) {
    int xColumn = xIndices[j] ;
    double xValue = xElements[j] ;
    int k = (j&1) ;
    value += (j+1)*multiplier[k]*(xColumn+1)*xValue ;
  }
/*
  It really would be nice if sizeof where legal in the preprocessor. But a
  good optimizing compiler shouldn't have any trouble.
*/
  unsigned int hashValue ;
  if (sizeof(double) == sizeof(int)) {
    union { double d; int i; } xx ;
    xx.d = value ;
    hashValue = xx.i ;
  } else if (sizeof(double) == 2*sizeof(int)) {
    union { double d; int i[2]; } xx ;
    xx.d = value ;
    hashValue = (xx.i[0] + xx.i[1]) ;
  } else {
    std::cout << "Hey! Time to implement another case, eh?" ;
    assert(false) ;
  }

  return (hashValue%size) ;

}

/*
  Compare two cuts to see if they're the same. Invoked when there's a hash
  collision. Coefficient order is critical --- the coefficients are compared
  in the order given.

  Arguably we could use OsiRowCut::operator==, but it's not a toleranced
  comparison, and that can be a problem. Also, we don't check efficiency here.
*/
bool same (const OsiRowCut &x, const OsiRowCut &y) 
{
  int xN = x.row().getNumElements() ;
  int yN = y.row().getNumElements() ;
  bool identical = false ;
  if (xN == yN) {
    double xLb = x.lb() ;
    double xUb = x.ub() ;
    double yLb = y.lb() ;
    double yUb = y.ub() ;
    if (fabs(xLb-yLb) < 1.0e-8 && fabs(xUb-yUb) < 1.0e-8) {
      const int * xIndices = x.row().getIndices() ;
      const double * xElements = x.row().getElements() ;
      const int * yIndices = y.row().getIndices() ;
      const double * yElements = y.row().getElements() ;
      int j ;
      for (j = 0 ; j < xN ; j++) {
	if (xIndices[j] != yIndices[j]) break ;
	if (fabs(xElements[j]-yElements[j]) > 1.0e-12) break ;
      }
      identical = (j == xN) ;
    }
  }
  return (identical) ;
}

} // end file-local namespace

/*
  Note that while nRows is used to size the collection, the connection won't
  be obvious to the casual observer. Expanding the capacity later is
  expensive, so don't scrimp on the initial allocation.
*/
CglProbingRowCut::CglProbingRowCut (int nRows, bool initialPass)
{
  // Calculate the maximum size of the collection
  if (nRows < 500) {
    maxSize_ = SIZE_ROW_MULT*nRows+SIZE_ROW_ADD ;
  } else if (nRows < 5000) {
    maxSize_ = (SIZE_ROW_MULT*nRows+SIZE_ROW_ADD)>>1 ;
  } else if (nRows < 10000) {
    maxSize_ = (SIZE_ROW_MULT*(nRows>>1)+SIZE_ROW_ADD)>>1 ;
  } else {
    maxSize_ = (SIZE_ROW_MULT*CoinMin(nRows,100000)+SIZE_ROW_ADD)>>2 ;
  }

  // Calculate the initial size of the cut collection
  size_ = (maxSize_>>3)+10 ;
  if (initialPass)
    size_ = size_>>1 ;
  // Set up the vector of cuts
  nRows_ = nRows ;
  rowCut_ = new  OsiRowCut* [size_] ;
  numberCuts_ = 0 ;

  // Set up the hash table
  if (size_ < 1000)
    hashSize_ = 4*size_ ;
  else
    hashSize_ = 2*size_ ;
  hash_ = new HashLink [hashSize_] ;
  for (int i = 0 ; i < hashSize_ ; i++) {
    hash_[i].index = -1 ;
    hash_[i].next = -1 ;
  }
  lastHash_ = -1 ;
}

/*
  Destructor. Delete any cuts, then remove the cut vector and hash table.
*/
CglProbingRowCut::~CglProbingRowCut ()
{
  for (int i = 0 ; i < numberCuts_ ; i++) delete rowCut_[i] ;
  delete [] rowCut_ ;
  delete [] hash_ ;
}

/*
  If the cut is not in the table, and there is room for the cut, returns
  the index of the appropriate empty hash table entry (index field = -1).
  If the cut is already in the table, returns the entry for the cut (index
  field >= 0).  If there's no room, returns -1
*/
int CglProbingRowCut::findHashTableEntry (const OsiRowCut &newCut)
{
/*
  Hash the cut and hope we hit an empty entry. If so, we're done.
*/
  int newNdx = hashCut(newCut,hashSize_) ;
  if (hash_[newNdx].index == -1) return (newNdx) ;
/*
  Collision. Walk the chain of entries starting from here, looking for a
  duplicate cut or the end of the chain.
*/
  int dupNdx = -1 ;
  for (int ndx = newNdx ; ndx >= 0 ; ndx = hash_[ndx].next) {
    int cutNdx = hash_[ndx].index ;
    if (same(newCut,*rowCut_[cutNdx])) {
      dupNdx = ndx ;
      break ;
    }
  }
  if (dupNdx >= 0) return (dupNdx) ;
/*
  We didn't find a duplicate, so we need to add to the chain. Start a linear
  search for an empty entry from the point where we ended the last search.
*/
  assert(lastHash_ < hashSize_) ;
  int startSearch = lastHash_ ; 
  for (lastHash_++ ;
       lastHash_ < hashSize_ && hash_[lastHash_].index >= 0 ;
       lastHash_++) ;
  if (lastHash_ >= hashSize_) {
    for (lastHash_ = 0 ;
	 lastHash_ < startSearch && hash_[lastHash_].index >= 0 ;
	 lastHash_++) ;
  }
  if (hash_[lastHash_].index >= 0) {
    std::cout
      << "  CglProbingRowCut::insertInHashTable: no room in hash table!"
      << std::endl ;
    return (-1) ;
  }
  assert(hash_[lastHash_].index == -1) ;
/*
  Insert the new link at the beginning of the chain
*/
  hash_[lastHash_].next = hash_[newNdx].next ;
  hash_[newNdx].next = lastHash_ ;

  return (lastHash_) ;
}


/*
  Add a cut only if it's not a duplicate of an existing cut.
  
  Return 0 if added, 1 if not, -1 if not added because of space
*/
int CglProbingRowCut::addCutIfNotDuplicate(OsiRowCut &cut)
{
/*
  Time to expand the cut storage vector and hash table? We have to repopulate
  the hash table, so this is a fairly expensive operation. If we're at the
  max, return -1.
*/
  if (numberCuts_ >= size_) {
    if (numberCuts_ >= maxSize_) return (-1) ;
    size_ = CoinMin(2*size_+100,maxSize_) ;
    if (size_ < 1000)
      hashSize_ = 4*size_ ;
    else
      hashSize_ = 2*size_ ;
#   if CGL_DEBUG > 0
    std::cout
      << "  addCutIfNotDuplicate: increasing size from " << numberCuts_
      << " to " << size_ << " (max " << maxSize_ << ", hash " << hashSize_
      << "." << std::endl ;
#   endif
    OsiRowCut **oldRowCut = rowCut_ ;
    rowCut_ = new  OsiRowCut* [size_] ;
    HashLink *oldHash = hash_ ;
    hash_ = new HashLink[hashSize_] ;
    for (int i = 0 ; i < hashSize_ ; i++) {
      hash_[i].index = -1 ;
      hash_[i].next = -1 ;
    }
/*
  Repopulate the hash table. A return of -1 indicates no space in the hash
  table, which really shouldn't happen here.
*/
    bool failed = false ;
    for (int i = 0 ; i < numberCuts_ ; i++) {
      rowCut_[i] = oldRowCut[i] ;
      int hashNdx = findHashTableEntry(*rowCut_[i]) ;
      if (hashNdx < 0) {
        failed = true ;
	break ;
      }
      if (hash_[hashNdx].index < 0) {
        hash_[hashNdx].index = i ;
      }
#     if CGL_DEBUG > 0
      else {
        std::cout
	  << "    addCutIfNotDuplicate: encountered duplicate cut while "
	  << "repopulating hash table after expansion." << std::endl ;
      }
#     endif
    }
    if (failed) {
      delete [] rowCut_ ;
      rowCut_ = oldRowCut ;
      delete [] hash_ ;
      hash_ = oldHash ;
      return (-1) ;
    }
    delete [] oldHash ;
    delete [] oldRowCut ;
  }
/*
  There is now space to add the cut, or we would have already returned.
*/
  CoinPackedVector &cutVec = cut.mutableRow() ;
# if 1 // CGL_DEBUG > 0
  int bad = 0 ;
  int coeffCnt = cutVec.getNumElements() ;
  const double *constCoeffs = cutVec.getElements() ;
  for (int i = 0 ; i < coeffCnt ; i++) {
    double coeff = fabs(constCoeffs[i]) ;
    if (CoinAbs(coeff) < 1.0e-12 || CoinAbs(coeff) > 1.0e12) bad++ ;
  }
  if (bad) {
    std::cout
      << "    addCutIfNotDuplicate: BAD: " << bad << " coefficients in cut."
      << std::endl ;
    assert(false) ;
  }
# endif
/*
  Find the hash table entry that matches the cut. We need to sort the indices
  to ensure a standard order, because the hash function is order-sensitive.
  The exact order isn't important.
*/
  cutVec.sortIncrIndex() ;
  int hashNdx = findHashTableEntry(cut) ;
  if (hash_[hashNdx].index < 0)
  { OsiRowCut *newCut = cut.clone() ;
    rowCut_[numberCuts_] = newCut ;
    hash_[hashNdx].index = numberCuts_ ;
    numberCuts_++ ;
  }

  return (0) ;
}

/*
  Transfer up to nRows_ cuts out of this CglProbingRowCut structure
  into the standard OsiCuts structure given as a parameter. If the total
  number of cuts is less than nRows_, take them all, otherwise, just take
  the best. nRows_ is a nominal limit which can be exceeded if the cuts
  straddling the limit all have the same effectiveness.

  If whichRows is supplied, entries which do not already associate a cut with
  the row will be filled in, if a cut is transferred that identifies itself
  with the row.

  When we're taking all cuts, scan from first or last cut on alternate passes,
  so that cuts generated later in the pass have a chance of being associated
  with a row. When we sort for the best cuts, use the best.
*/
void CglProbingRowCut::addCuts (OsiCuts &cs, OsiRowCut **whichRow, int iPass)
{
  int numberCuts = cs.sizeRowCuts() ;
  int i ;
/*
  Take all, working from first or last on alternate passes so that later cuts
  have a chance to claim places in whichRow.
*/
  if (numberCuts_ < nRows_) {
    if ((iPass&1) == 1) {
      for (i = 0 ; i < numberCuts_ ; i++) {
	cs.insert(*rowCut_[i]) ;
	if (whichRow) {
	  int iRow = rowCut_[i]->whichRow() ;
	  if (iRow >= 0 && !whichRow[iRow])
	    whichRow[iRow] = cs.rowCutPtr(numberCuts) ;
	}
	numberCuts++ ;
      }
    } else {
      for (i = numberCuts_-1 ; i >= 0 ; i--) {
	cs.insert(*rowCut_[i]) ;
	if (whichRow) {
	  int iRow = rowCut_[i]->whichRow() ;
	  if (iRow >= 0 && !whichRow[iRow])
	    whichRow[iRow] = cs.rowCutPtr(numberCuts) ;
	}
	numberCuts++ ;
      }
    }
  } else {
/*
  Transfer just the best cuts. Load effectiveness values into a vector, sort,
  select the threshold, then scan the cuts and take the ones over the
  threshold. Cuts that are associated with rows get star billing.
*/
    double *effectiveness = new double[numberCuts_] ;
    int iCut = 0 ;
    for (iCut = 0 ; iCut < numberCuts_ ; iCut++) {
      double value = -rowCut_[iCut]->effectiveness() ;
      if (whichRow) {
	int iRow = rowCut_[iCut]->whichRow() ;
	if (iRow >= 0)
	  value -= 1.0e10 ;
      }
      effectiveness[iCut] = value ;
    }
    std::sort(effectiveness,effectiveness+numberCuts_) ;
    double threshold = effectiveness[nRows_-1] ;
    for (i = 0 ; i < numberCuts_ ; i++) {
      if (rowCut_[i]->effectiveness() > threshold) {
	cs.insert(*rowCut_[i]) ;
	if (whichRow) {
	  int iRow = rowCut_[i]->whichRow() ;
	  if (iRow >= 0 && !whichRow[iRow])
	    whichRow[iRow] = cs.rowCutPtr(numberCuts) ;
	}
	numberCuts++ ;
      }
    }
    delete[] effectiveness ;
  }
  for (i = 0 ; i < numberCuts_ ; i++)
  { delete rowCut_[i] ;
    rowCut_[i] = 0 ; }
  numberCuts_ = 0 ;
}

