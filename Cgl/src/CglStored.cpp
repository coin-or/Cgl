// Copyright (C) 2005, International Business Machines
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
//#define CGL_DEBUG 2
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglStored.hpp"
#include "CglTreeInfo.hpp"
#include "CoinFinite.hpp"
//-------------------------------------------------------------------
// Generate Stored cuts
//------------------------------------------------------------------- 
void 
CglStored::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info) const
{
  // Get basic problem information
  const double * solution = si.getColSolution();
  int numberRowCuts = cuts_.sizeRowCuts();
  for (int i=0;i<numberRowCuts;i++) {
    const OsiRowCut * rowCutPointer = cuts_.rowCutPtr(i);
    double violation = rowCutPointer->violated(solution);
    if (violation>=requiredViolation_)
      cs.insert(*rowCutPointer);
  }
  if (probingInfo_) {
    int number01 = probingInfo_->numberIntegers();
    const fixEntry * entry = probingInfo_->fixEntries();
    const int * toZero = probingInfo_->toZero();
    const int * toOne = probingInfo_->toOne();
    const int * integerVariable = probingInfo_->integerVariable();
    const double * lower = si.getColLower();
    const double * upper = si.getColUpper();
    OsiRowCut cut;
    int column[2];
    double element[2];
    for (int i=0;i<number01;i++) {
      int iColumn=integerVariable[i];
      if (upper[iColumn]==lower[iColumn])
	continue;
      double value1 = solution[iColumn];
      for (int j=toZero[i];j<toOne[i];j++) {
	int jColumn=entry[j].sequence;
	if (jColumn<number01) {
	  jColumn=integerVariable[jColumn];
	  assert (jColumn>=0);
	  double value2 = solution[jColumn];
	  if (entry[j].oneFixed) {
	    double violation = 1.0-value1-value2;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %d + %d >=1\n",iColumn,jColumn);
	      cut.setLb(1.0);
	      cut.setUb(COIN_DBL_MAX);
	      column[0]=iColumn;
	      element[0]=1.0;
	      column[1]=jColumn;
	      element[1]= 1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  } else {
	    double violation = value2-value1;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %d >= %d\n",iColumn,jColumn);
	      cut.setLb(0.0);
	      cut.setUb(COIN_DBL_MAX);
	      column[0]=iColumn;
	      element[0]=1.0;
	      column[1]=jColumn;
	      element[1]= -1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  }
	} else {
	  jColumn -= number01; // not 0-1
	  double value2 = solution[jColumn];
	  double lowerValue = lower[jColumn];
	  double upperValue = upper[jColumn];
	  if (entry[j].oneFixed) {
	    double violation = upperValue-value1*(upperValue-lowerValue)-value2;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %g*%d + %d >=%g\n",(upperValue-lowerValue),iColumn,jColumn,upperValue);
	      cut.setLb(upperValue);
	      cut.setUb(COIN_DBL_MAX);
	      column[0]=iColumn;
	      element[0]=upperValue-lowerValue;
	      column[1]=jColumn;
	      element[1]= 1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  } else {
	    double violation = value2-value1*(upperValue-lowerValue)-lowerValue;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %g*%d >= %d -%g\n",(upperValue-lowerValue),iColumn,jColumn,lowerValue);
	      cut.setLb(-lowerValue);
	      cut.setUb(COIN_DBL_MAX);
	      column[0]=iColumn;
	      element[0]=upperValue-lowerValue;
	      column[1]=jColumn;
	      element[1]= -1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  }
	}
      }
      for (int j=toOne[i];j<toZero[i+1];j++) {
	int jColumn=entry[j].sequence;
	if (jColumn<number01) {
	  jColumn=integerVariable[jColumn];
	  assert (jColumn>=0);
	  double value2 = solution[jColumn];
	  if (entry[j].oneFixed) {
	    double violation = value1-value2;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %d <= %d\n",iColumn,jColumn);
	      cut.setLb(-COIN_DBL_MAX);
	      cut.setUb(0.0);
	      column[0]=iColumn;
	      element[0]=1.0;
	      column[1]=jColumn;
	      element[1]= -1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  } else {
	    double violation = value1+value2-1.0;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %d + %d <=1\n",iColumn,jColumn);
	      cut.setLb(-COIN_DBL_MAX);
	      cut.setUb(1.0);
	      column[0]=iColumn;
	      element[0]=1.0;
	      column[1]=jColumn;
	      element[1]= 1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  }
	} else {
	  jColumn -= number01; // not 0-1
	  double value2 = solution[jColumn];
	  double lowerValue = lower[jColumn];
	  double upperValue = upper[jColumn];
	  if (entry[j].oneFixed) {
	    double violation = lowerValue +(upperValue-lowerValue)*value1-value2;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %g*%d <= %d -%g\n",(upperValue-lowerValue),iColumn,jColumn,lowerValue);
	      cut.setLb(-COIN_DBL_MAX);
	      cut.setUb(-lowerValue);
	      column[0]=iColumn;
	      element[0]=upperValue-lowerValue;
	      column[1]=jColumn;
	      element[1]= -1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  } else {
	    double violation = (upperValue-lowerValue)*value1+value2-upperValue;
	    if (violation>requiredViolation_) {
	      //printf("XXX can do %g*%d + %d <=%g\n",(upperValue-lowerValue),iColumn,jColumn,upperValue);
	      cut.setLb(-COIN_DBL_MAX);
	      cut.setUb(upperValue);
	      column[0]=iColumn;
	      element[0]=upperValue-lowerValue;
	      column[1]=jColumn;
	      element[1]= 1.0;
	      cut.setEffectiveness(violation);
	      cut.setRow(2,column,element,false);
	      cs.insert(cut);
	    }
	  }
	}
      }
    }
  }
}
// Add cuts
void 
CglStored::addCut(const OsiCuts & cs)
{
  int numberRowCuts = cs.sizeRowCuts();
  for (int i=0;i<numberRowCuts;i++) {
    cuts_.insert(*cs.rowCutPtr(i));
  }
}
// Add a row cut
void 
CglStored::addCut(const OsiRowCut & cut)
{
  cuts_.insert(cut);
}
// Add a row cut from a packed vector
void 
CglStored::addCut(double lb, double ub, const CoinPackedVector & vector)
{
  OsiRowCut rc;
  rc.setRow(vector);
  rc.mutableRow().setTestForDuplicateIndex(false);
  rc.setLb(lb);
  rc.setUb(ub);   
  cuts_.insert(rc);
}
// Add a row cut from elements
void 
CglStored::addCut(double lb, double ub, int size, const int * colIndices, const double * elements)
{
  OsiRowCut rc;
  rc.setRow(size,colIndices,elements,false);
  rc.setLb(lb);
  rc.setUb(ub);   
  cuts_.insert(rc);
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglStored::CglStored ()
:
  CglCutGenerator(),
  requiredViolation_(1.0e-5),
  probingInfo_(NULL)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglStored::CglStored (const CglStored & source) :
  CglCutGenerator(source),
  requiredViolation_(source.requiredViolation_),
  probingInfo_(NULL),
  cuts_(source.cuts_)
{  
  if (source.probingInfo_)
    probingInfo_ = new CglTreeProbingInfo(*source.probingInfo_);
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglStored::clone() const
{
  return new CglStored(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglStored::~CglStored ()
{
  delete  probingInfo_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglStored &
CglStored::operator=(const CglStored& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    requiredViolation_=rhs.requiredViolation_;
    cuts_=rhs.cuts_;
    delete probingInfo_;
    if (rhs.probingInfo_)
      probingInfo_ = new CglTreeProbingInfo(*rhs.probingInfo_);
    else
      probingInfo_ = NULL;
  }
  return *this;
}
