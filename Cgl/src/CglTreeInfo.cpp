// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CglTreeInfo.hpp"
#include "CoinHelperFunctions.hpp"

// Default constructor 
CglTreeInfo::CglTreeInfo ()
  : level(-1), pass(-1), formulation_rows(-1), inTree(false),
    strengthenRow(NULL) {}

// Copy constructor 
CglTreeInfo::CglTreeInfo (const CglTreeInfo & rhs)
  : level(rhs.level), 
    pass(rhs.pass), 
    formulation_rows(rhs.formulation_rows), 
    inTree(rhs.inTree),
    strengthenRow(rhs.strengthenRow)
{
}
// Clone
CglTreeInfo * 
CglTreeInfo::clone() const
{
  return new CglTreeInfo(*this);
}

// Assignment operator 
CglTreeInfo &
CglTreeInfo::operator=(const CglTreeInfo& rhs)
{
  if (this != &rhs) {
    //CglCutGenerator::operator=(rhs);
    level = rhs.level; 
    pass = rhs.pass; 
    formulation_rows = rhs.formulation_rows; 
    inTree = rhs.inTree;
    strengthenRow = rhs.strengthenRow;
  }
  return *this;
}
  
 // Destructor 
CglTreeInfo::~CglTreeInfo ()
{
}

// Default constructor 
CglTreeProbingInfo::CglTreeProbingInfo ()
  : CglTreeInfo(),
    fixEntry_(NULL),
    toZero_(NULL),
    toOne_(NULL),
    integerVariable_(NULL),
    backward_(NULL),
    fixingEntry_(NULL),
    numberVariables_(0),
    numberIntegers_(0),
    maximumEntries_(0),
    numberEntries_(-1)
{
}
// Constructor from model
CglTreeProbingInfo::CglTreeProbingInfo (const OsiSolverInterface * model)
  : CglTreeInfo(),
    fixEntry_(NULL),
    toZero_(NULL),
    toOne_(NULL),
    integerVariable_(NULL),
    backward_(NULL),
    fixingEntry_(NULL),
    numberVariables_(0),
    numberIntegers_(0),
    maximumEntries_(0),
    numberEntries_(-1)
{
  numberVariables_=model->getNumCols(); 
  // Too many ... but
  integerVariable_ = new int [numberVariables_];
  backward_ = new int [numberVariables_];
  int i;
  for (i=0;i<numberVariables_;i++) {
    backward_[i]=-1;
    if (model->isInteger(i)) {
      if (model->isBinary(i)) {
	backward_[i]=numberIntegers_;
	integerVariable_[numberIntegers_++]=i;
      } else {
	backward_[i]=-2;
      }
    }
  }
  // Set up to arrays
  toOne_ = new int[numberIntegers_];
  toZero_ = new int[numberIntegers_+1];
  // zero out
  CoinZeroN(toOne_,numberIntegers_);
  CoinZeroN(toZero_,numberIntegers_+1);
}

// Copy constructor 
CglTreeProbingInfo::CglTreeProbingInfo (const CglTreeProbingInfo & rhs)
  : CglTreeInfo(rhs),
    fixEntry_(NULL),
    toZero_(NULL),
    toOne_(NULL),
    integerVariable_(NULL),
    backward_(NULL),
    fixingEntry_(NULL),
    numberVariables_(rhs.numberVariables_),
    numberIntegers_(rhs.numberIntegers_),
    maximumEntries_(rhs.maximumEntries_),
    numberEntries_(rhs.numberEntries_)
{
  if (numberVariables_) {
    fixEntry_ = new fixEntry [maximumEntries_];
    memcpy(fixEntry_,rhs.fixEntry_,maximumEntries_*sizeof(fixEntry));
    if (numberEntries_<0) {
      // in order
      toZero_ = CoinCopyOfArray(rhs.toZero_,numberIntegers_+1);
      toOne_ = CoinCopyOfArray(rhs.toOne_,numberIntegers_);
    } else {
      // not in order
      fixingEntry_ = CoinCopyOfArray(rhs.fixingEntry_,maximumEntries_);
    }
    integerVariable_ = CoinCopyOfArray(rhs.integerVariable_,numberIntegers_);
    backward_ = CoinCopyOfArray(rhs.backward_,numberVariables_);
  }
}
// Clone
CglTreeInfo * 
CglTreeProbingInfo::clone() const
{
  return new CglTreeProbingInfo(*this);
}

// Assignment operator 
CglTreeProbingInfo &
CglTreeProbingInfo::operator=(const CglTreeProbingInfo& rhs)
{
  if (this != &rhs) {
    CglTreeInfo::operator=(rhs);
    delete [] fixEntry_;
    delete [] toZero_;
    delete [] toOne_;
    delete [] integerVariable_;
    delete [] backward_;
    delete [] fixingEntry_;
    numberVariables_ = rhs.numberVariables_;
    numberIntegers_ = rhs.numberIntegers_;
    maximumEntries_ = rhs.maximumEntries_;
    numberEntries_ = rhs.numberEntries_;
    if (numberVariables_) {
      fixEntry_ = new fixEntry [maximumEntries_];
      memcpy(fixEntry_,rhs.fixEntry_,maximumEntries_*sizeof(fixEntry));
      if (numberEntries_<0) {
	// in order
	toZero_ = CoinCopyOfArray(rhs.toZero_,numberIntegers_+1);
	toOne_ = CoinCopyOfArray(rhs.toOne_,numberIntegers_);
	fixingEntry_ = NULL;
      } else {
	// not in order
	fixingEntry_ = CoinCopyOfArray(rhs.fixingEntry_,maximumEntries_);
	toZero_ = NULL;
	toOne_ = NULL;
      }
      toZero_ = CoinCopyOfArray(rhs.toZero_,numberIntegers_+1);
      toOne_ = CoinCopyOfArray(rhs.toOne_,numberIntegers_);
      integerVariable_ = CoinCopyOfArray(rhs.integerVariable_,numberIntegers_);
      backward_ = CoinCopyOfArray(rhs.backward_,numberVariables_);
    } else {
      fixEntry_ = NULL;
      toZero_ = NULL;
      toOne_ = NULL;
      integerVariable_ = NULL;
      backward_ = NULL;
      fixingEntry_ = NULL;
    }
  }
  return *this;
}
  
 // Destructor 
CglTreeProbingInfo::~CglTreeProbingInfo ()
{
  delete [] fixEntry_;
  delete [] toZero_;
  delete [] toOne_;
  delete [] integerVariable_;
  delete [] backward_;
  delete [] fixingEntry_;
}
// Take action if cut generator can fix a variable (toValue -1 for down, +1 for up)
void 
CglTreeProbingInfo::fixes(int variable, int toValue, int fixedVariable,double fixedToValue)
{
  //printf("%d going to %d fixes %d at %g\n",variable,toValue,fixedVariable,fixedToValue);
  int intVariable = backward_[variable];
  if (intVariable<0) // off as no longer in order FIX
    return; // not 0-1 (well wasn't when constructor was called)
  int intFix = backward_[fixedVariable];
  if (intFix<0)
    return; // not 0-1
  int fixedTo = (int) fixedToValue;
  if (numberEntries_==maximumEntries_) {
    maximumEntries_ += 100 +maximumEntries_/2;
    fixEntry * temp1 = new fixEntry [maximumEntries_];
    memcpy(temp1,fixEntry_,numberEntries_*sizeof(fixEntry));
    delete [] fixEntry_;
    fixEntry_ = temp1;
    int * temp2 = new int [maximumEntries_];
    memcpy(temp2,fixingEntry_,numberEntries_*sizeof(int));
    delete [] fixingEntry_;
    fixingEntry_ = temp2;
  }
  fixEntry entry1;
  entry1.oneFixed=fixedTo;
  entry1.sequence=intFix;
  fixEntry_[numberEntries_] = entry1;
  assert (toValue==-1||toValue==1);
  assert (fixedTo==0||fixedTo==1);
  if (toValue<0)
    fixingEntry_[numberEntries_++] = intVariable << 1;
  else
    fixingEntry_[numberEntries_++] = (intVariable << 1) | 1;
}
// Initalizes fixing arrays etc - returns true if we want to save info
bool 
CglTreeProbingInfo::initializeFixing() 
{
  delete [] toZero_;
  delete [] toOne_;
  delete [] fixingEntry_;
  toZero_ = NULL;
  toOne_ = NULL;
  fixingEntry_ = new int[maximumEntries_];
  numberEntries_ = 0;
  return true;
}
// Converts to ordered and takes out duplicates
void 
CglTreeProbingInfo::convert() const
{
  if (numberEntries_>=0) {
    CoinSort_2( fixingEntry_, fixingEntry_+numberEntries_, fixEntry_);
    assert (!toZero_);
    toZero_ = new int [numberIntegers_+1];
    toOne_ = new int [numberIntegers_];
    toZero_[0]=0;
    int n=0;
    int put=0;
    for (int intVariable = 0;intVariable<numberIntegers_;intVariable++) {
      int last = n;
      for (;n<numberEntries_;n++) {
	int value = fixingEntry_[n];
	int iVar = value>>1;
	int way = value &1;
	if (intVariable!=iVar||way)
	  break;
      }
      if (n>last) {
	// sort
	assert (sizeof(int)==4);
	std::sort((unsigned int *) fixEntry_+last,(unsigned int *) fixEntry_+n);
	fixEntry temp2;
	temp2.sequence=numberVariables_+1;
	for (int i=last;i<n;i++) {
	  if (temp2.sequence!=fixEntry_[i].sequence||temp2.oneFixed||fixEntry_[i].oneFixed) {
	    temp2 = fixEntry_[i];
	    fixEntry_[put++]=temp2;
	  }
	}
      }
      toOne_[intVariable]=put;
      last = n;
      for (;n<numberEntries_;n++) {
	int value = fixingEntry_[n];
	int iVar = value>>1;
	if (intVariable!=iVar)
	  break;
      }
      if (n>last) {
	// sort
	assert (sizeof(int)==4);
	std::sort((unsigned int *) fixEntry_+last,(unsigned int *) fixEntry_+n);
	fixEntry temp2;
	temp2.sequence=numberVariables_+1;
	for (int i=last;i<n;i++) {
	  if (temp2.sequence!=fixEntry_[i].sequence||temp2.oneFixed||fixEntry_[i].oneFixed) {
	    temp2 = fixEntry_[i];
	    fixEntry_[put++]=temp2;
	  }
	}
	last=n;
      }
      toZero_[intVariable+1]=put;
    }
    delete [] fixingEntry_;
    fixingEntry_ = NULL;
    numberEntries_ = -1;
  }
}
