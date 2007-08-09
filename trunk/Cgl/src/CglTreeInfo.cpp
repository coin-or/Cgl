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
    numberVariables_(0),
    numberIntegers_(0),
    maximumEntries_(0),
    lastInteger_(-1)
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
    numberVariables_(0),
    numberIntegers_(0),
    maximumEntries_(0),
    lastInteger_(-1)
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
    numberVariables_(rhs.numberVariables_),
    numberIntegers_(rhs.numberIntegers_),
    maximumEntries_(rhs.maximumEntries_),
    lastInteger_(rhs.lastInteger_)
{
  if (numberVariables_) {
    fixEntry_ = new fixEntry [maximumEntries_];
    memcpy(fixEntry_,rhs.fixEntry_,maximumEntries_*sizeof(fixEntry));
    toZero_ = CoinCopyOfArray(rhs.toZero_,numberIntegers_+1);
    toOne_ = CoinCopyOfArray(rhs.toOne_,numberIntegers_);
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
    numberVariables_ = rhs.numberVariables_;
    numberIntegers_ = rhs.numberIntegers_;
    maximumEntries_ = rhs.maximumEntries_;
    lastInteger_ = rhs.lastInteger_;
    if (numberVariables_) {
      fixEntry_ = new fixEntry [maximumEntries_];
      memcpy(fixEntry_,rhs.fixEntry_,maximumEntries_*sizeof(fixEntry));
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
}
// Take action if cut generator can fix a variable (toValue -1 for down, +1 for up)
void 
CglTreeProbingInfo::fixes(int variable, int toValue, int fixedVariable,double fixedToValue)
{
  //printf("%d going to %d fixes %d at %g\n",variable,toValue,fixedVariable,fixedToValue);
  // should be more sophisticated
  int intVariable = backward_[variable];
  if (intVariable<0)
    return; // not 0-1 (well wasn't when constructor was called)
  int intFix = backward_[fixedVariable];
  if (intFix<0)
    return; // not 0-1
  int fixedTo = (int) fixedToValue;
  assert (intVariable>=lastInteger_);
  int n = toZero_[lastInteger_+1];
  while (intVariable>lastInteger_) {
    lastInteger_++;
    toOne_[lastInteger_]=n;
    toZero_[lastInteger_+1]=n;
  }
  if (n==maximumEntries_) {
    maximumEntries_ += 100 +maximumEntries_/2;
    fixEntry * temp = new fixEntry [maximumEntries_];
    memcpy(temp,fixEntry_,n*sizeof(fixEntry));
    delete [] fixEntry_;
    fixEntry_ = temp;
  }
  fixEntry entry;
  entry.oneFixed=fixedTo;
  entry.sequence=intFix;
  if (toValue==-1) {
    // to 0
    int k = toOne_[intVariable];
    if (toOne_[intVariable]<n) {
      fixEntry entry2=fixEntry_[k];
      fixEntry_[k]=entry;
      entry=entry2;
    }
    toOne_[intVariable]++;
  } else {
    assert (toValue==1);
    // to 1
  }
  fixEntry_[n++]=entry;
  toZero_[intVariable+1]=n;
}
// Initalizes fixing arrays etc - returns true if we want to save info
bool 
CglTreeProbingInfo::initializeFixing() 
{
  // zero out
  CoinZeroN(toOne_,numberIntegers_);
  CoinZeroN(toZero_,numberIntegers_+1);
  lastInteger_=0;
  return true;
}
