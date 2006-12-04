// Name:     CglRedSplitData.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     11/24/06
//---------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.

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

#include "CglRedSplitData.hpp"

/***********************************************************************/
void CglRedSplitData::setSolverPtr(
                                     OsiSolverInterface *givenSolverPtr) {
  solverPtr = givenSolverPtr;
} /* setSolverPtr */

/***********************************************************************/
void CglRedSplitData::setOptimalBasisIsAvailable(
                         const int givenOptBas) {
  optimalBasisIsAvailable = givenOptBas;
} /* setOptimalBasisIsAvailable */

/***********************************************************************/
CglRedSplitData::CglRedSplitData(const OsiSolverInterface *givenSolverPtr,
				 const int givenOptBas) :
  optimalBasisIsAvailable(givenOptBas)
{
  solverPtr = const_cast<OsiSolverInterface *>(givenSolverPtr);
}

/***********************************************************************/
CglRedSplitData::CglRedSplitData(const CglData &source,
				 const OsiSolverInterface *givenSolverPtr,
				 const int givenOptBas) :
  CglData(source),
  optimalBasisIsAvailable(givenOptBas)
{
  solverPtr = const_cast<OsiSolverInterface *>(givenSolverPtr);
}

/***********************************************************************/
CglRedSplitData::CglRedSplitData(const CglRedSplitData &source) :
  CglData(source),
  solverPtr(source.solverPtr),
  optimalBasisIsAvailable(source.optimalBasisIsAvailable)
{}

/***********************************************************************/
CglRedSplitData* CglRedSplitData::clone() const
{
  return new CglRedSplitData(*this);
}

/***********************************************************************/
CglRedSplitData& CglRedSplitData::operator=(const CglRedSplitData &rhs)
{
  if(this != &rhs) {
    CglData::operator=(rhs);

    solverPtr = rhs.solverPtr;
    optimalBasisIsAvailable = rhs.optimalBasisIsAvailable;
  }
  return *this;
}

/***********************************************************************/
CglRedSplitData::~CglRedSplitData()
{}
