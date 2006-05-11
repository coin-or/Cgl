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
  rc.setLb(lb);
  rc.setUb(ub);   
  cuts_.insert(rc);
}
// Add a row cut from elements
void 
CglStored::addCut(double lb, double ub, int size, const int * colIndices, const double * elements)
{
  OsiRowCut rc;
  rc.setRow(size,colIndices,elements);
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
requiredViolation_(1.0e-5)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglStored::CglStored (const CglStored & source) :
  CglCutGenerator(source),
  requiredViolation_(source.requiredViolation_),
  cuts_(source.cuts_)
{  
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
  }
  return *this;
}
