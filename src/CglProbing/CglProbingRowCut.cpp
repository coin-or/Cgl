
/*
  Copyright (C) 2002, International Business Machines Corporation and others.
  All Rights Reserved.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinHelperFunctions.hpp"
#include "CglProbingRowCut.hpp"

#define SIZE_ROW_MULT 4
#define SIZE_ROW_ADD 2000

/*
  Hash functions. Why is there no CoinUtils hash package? Apparently it was
  easier to just copy.  -- lh, 101202 --
*/
namespace {

double multiplier[] = { 1.23456789e2, -9.87654321 } ;

int hashCut (const OsiRowCut2 & x, int size) 
{
  int xN =x.row().getNumElements();
  double xLb = x.lb();
  double xUb = x.ub();
  const int * xIndices = x.row().getIndices();
  const double * xElements = x.row().getElements();
  unsigned int hashValue;
  double value=1.0;
  if (xLb>-1.0e10)
    value += xLb*multiplier[0];
  if (xUb<1.0e10)
    value += xUb*multiplier[1];
  for( int j=0;j<xN;j++) {
    int xColumn = xIndices[j];
    double xValue = xElements[j];
    int k=(j&1);
    value += (j+1)*multiplier[k]*(xColumn+1)*xValue;
  }
  // should be compile time but too lazy for now
  if (sizeof(value)>sizeof(hashValue)) {
    assert (sizeof(value)==2*sizeof(hashValue));
    union { double d; int i[2]; } xx;
    xx.d = value;
    hashValue = (xx.i[0] + xx.i[1]);
  } else {
    assert (sizeof(value)==sizeof(hashValue));
    union { double d; unsigned int i[2]; } xx;
    xx.d = value;
    hashValue = xx.i[0];
  }
  return hashValue%(size);
}

bool same (const OsiRowCut2 & x, const OsiRowCut2 & y) 
{
  int xN =x.row().getNumElements();
  int yN =y.row().getNumElements();
  bool identical=false;
  if (xN==yN) {
    double xLb = x.lb();
    double xUb = x.ub();
    double yLb = y.lb();
    double yUb = y.ub();
    if (fabs(xLb-yLb)<1.0e-8&&fabs(xUb-yUb)<1.0e-8) {
      const int * xIndices = x.row().getIndices();
      const double * xElements = x.row().getElements();
      const int * yIndices = y.row().getIndices();
      const double * yElements = y.row().getElements();
      int j;
      for( j=0;j<xN;j++) {
	if (xIndices[j]!=yIndices[j])
	  break;
	if (fabs(xElements[j]-yElements[j])>1.0e-12)
	  break;
      }
      identical =  (j==xN);
    }
  }
  return identical;
}

} // end file-local namespace


CglProbingRowCut::CglProbingRowCut (int nRows, bool initialPass)
{
  numberCuts_ = 0 ;
  if (nRows < 500) {
    maxSize_ = SIZE_ROW_MULT*nRows+SIZE_ROW_ADD ;
  } else if (nRows<5000) {
    maxSize_ = (SIZE_ROW_MULT*nRows+SIZE_ROW_ADD)>>1 ;
  } else if (nRows<10000) {
    maxSize_ = (SIZE_ROW_MULT*(nRows>>1)+SIZE_ROW_ADD)>>1 ;
  } else {
    maxSize_ = (SIZE_ROW_MULT*CoinMin(nRows,100000)+SIZE_ROW_ADD)>>2 ;
  }
  size_ = (maxSize_>>3)+10 ;
/*
  TODO: When this is the initial pass, cut the capacity by half? Seems
	counterintuitive. For that matter, I don't expect the constructor
	to be making this decision for me.  -- lh, 100923 --
*/
  if (initialPass)
    size_ = size_>>1 ;
  if (size_<1000)
    hashSize_=4*size_ ;
  else
    hashSize_=2*size_ ;
  nRows_ = nRows ;
  rowCut_ = new  OsiRowCut2 * [size_] ;
  hash_ = new HashLink[hashSize_] ;
  for (int i=0;i<hashSize_;i++) {
    hash_[i].index=-1 ;
    hash_[i].next=-1 ;
  }
  numberCuts_=0 ;
  lastHash_=-1 ;
}


CglProbingRowCut::~CglProbingRowCut ()
{
  for (int i=0;i<numberCuts_;i++)
    delete rowCut_[i] ;
  delete [] rowCut_ ;
  delete [] hash_ ;
}

// Return 0 if added, 1 if not, -1 if not added because of space
int CglProbingRowCut::addCutIfNotDuplicate(OsiRowCut & cut,int whichRow)
{
  if (numberCuts_==size_&&numberCuts_<maxSize_) {
    size_ = CoinMin(2*size_+100,maxSize_) ;
    if (size_<1000)
      hashSize_=4*size_ ;
    else
      hashSize_=2*size_ ;
#ifdef COIN_DEVELOP
    printf("increaing size from %d to %d (hash size %d, maxsize %d)\n",
	   numberCuts_,size_,hashSize_,maxSize_) ;
#endif
    OsiRowCut2 ** temp = new  OsiRowCut2 * [size_] ;
    delete [] hash_ ;
    hash_ = new HashLink[hashSize_] ;
    for (int i=0;i<hashSize_;i++) {
      hash_[i].index=-1 ;
      hash_[i].next=-1 ;
    }
    for (int i=0;i<numberCuts_;i++) {
      temp[i]=rowCut_[i] ;
      int ipos = hashCut(*temp[i],hashSize_) ;
      int found = -1 ;
      int jpos=ipos ;
      while ( true ) {
	int j1 = hash_[ipos].index ;
	
	if ( j1 >= 0 ) {
	  if ( !same(*temp[i],*temp[j1]) ) {
	    int k = hash_[ipos].next ;
	    if ( k != -1 )
	      ipos = k ;
	    else
	      break ;
	  } else {
	    found = j1 ;
	    break ;
	  }
	} else {
	  break ;
	}
      }
      if (found<0) {
	assert (hash_[ipos].next==-1) ;
	if (ipos==jpos) {
	  // first
	  hash_[ipos].index=i ;
	} else {
	  // find next space 
	  while ( true ) {
	    ++lastHash_ ;
	    assert (lastHash_<hashSize_) ;
	    if ( hash_[lastHash_].index == -1 ) 
	      break ;
	  }
	  hash_[ipos].next = lastHash_ ;
	  hash_[lastHash_].index = i ;
	}
      }
    }
    delete [] rowCut_ ;
    rowCut_ = temp ;
  }
  if (numberCuts_<size_) {
    double newLb = cut.lb() ;
    double newUb = cut.ub() ;
    CoinPackedVector vector = cut.row() ;
    int numberElements =vector.getNumElements() ;
    int * newIndices = vector.getIndices() ;
    double * newElements = vector.getElements() ;
    CoinSort_2(newIndices,newIndices+numberElements,newElements) ;
    int i ;
    bool bad=false ;
    for (i=0;i<numberElements;i++) {
      double value = fabs(newElements[i]) ;
      if (value<1.0e-12||value>1.0e12) 
	bad=true ;
    }
    if (bad)
      return 1 ;
    OsiRowCut2 newCut(whichRow) ;
    newCut.setLb(newLb) ;
    newCut.setUb(newUb) ;
    newCut.setRow(vector) ;
    int ipos = hashCut(newCut,hashSize_) ;
    int found = -1 ;
    int jpos=ipos ;
    while ( true ) {
      int j1 = hash_[ipos].index ;

      if ( j1 >= 0 ) {
	if ( !same(newCut,*rowCut_[j1]) ) {
	  int k = hash_[ipos].next ;
	  if ( k != -1 )
	    ipos = k ;
	  else
	    break ;
	} else {
	  found = j1 ;
	  break ;
	}
      } else {
	break ;
      }
    }
    if (found<0) {
      assert (hash_[ipos].next==-1) ;
      if (ipos==jpos) {
	// first
	hash_[ipos].index=numberCuts_ ;
      } else {
	// find next space 
	while ( true ) {
	  ++lastHash_ ;
	  assert (lastHash_<hashSize_) ;
	  if ( hash_[lastHash_].index == -1 ) 
	    break ;
	}
	hash_[ipos].next = lastHash_ ;
	hash_[lastHash_].index = numberCuts_ ;
      }
      OsiRowCut2 * newCutPtr = new OsiRowCut2(whichRow) ;
      newCutPtr->setLb(newLb) ;
      newCutPtr->setUb(newUb) ;
      newCutPtr->setRow(vector) ;
      rowCut_[numberCuts_++]=newCutPtr ;
      return 0 ;
    } else {
      return 1 ;
    }
  } else {
    return -1 ;
  }
}

void CglProbingRowCut::addCuts (OsiCuts &cs, OsiRowCut **whichRow, int iPass)
{
  int numberCuts=cs.sizeRowCuts() ;
  int i ;
  if (numberCuts_<nRows_) {
    if ((iPass&1)==1) {
      for (i=0;i<numberCuts_;i++) {
	cs.insert(*rowCut_[i]) ;
	if (whichRow) {
	  int iRow= rowCut_[i]->whichRow() ;
	  if (iRow>=0&&!whichRow[iRow])
	    whichRow[iRow]=cs.rowCutPtr(numberCuts); ;
	}
	numberCuts++ ;
      }
    } else {
      for (i=numberCuts_-1;i>=0;i--) {
	cs.insert(*rowCut_[i]) ;
	if (whichRow) {
	  int iRow= rowCut_[i]->whichRow() ;
	  if (iRow>=0&&!whichRow[iRow])
	    whichRow[iRow]=cs.rowCutPtr(numberCuts); ;
	}
	numberCuts++ ;
      }
    }
  } else {
    // just best
    double * effectiveness = new double[numberCuts_] ;
    int iCut=0 ;
    for (i=0;i<numberCuts_;i++) {
      double value = -rowCut_[i]->effectiveness() ;
      if (whichRow) {
	int iRow= rowCut_[i]->whichRow() ;
	if (iRow>=0)
	  value -= 1.0e10 ;
      }
      effectiveness[iCut++]=value ;
    }
    std::sort(effectiveness,effectiveness+numberCuts_) ;
    double threshold = -1.0e20 ;
    if (iCut>nRows_)
      threshold = effectiveness[nRows_] ;
    for ( i=0;i<numberCuts_;i++) {
      if (rowCut_[i]->effectiveness()>threshold) {
	cs.insert(*rowCut_[i]) ;
	if (whichRow) {
	  int iRow= rowCut_[i]->whichRow() ;
	  if (iRow>=0&&!whichRow[iRow])
	    whichRow[iRow]=cs.rowCutPtr(numberCuts); ;
	}
	numberCuts++ ;
      }
    }
    delete[] effectiveness ;
  }
  for (i = 0 ; i < numberCuts_ ; i++)
  { delete rowCut_[i] ;
    rowCut_[i] = 0 ; }
  numberCuts_=0 ;
}

