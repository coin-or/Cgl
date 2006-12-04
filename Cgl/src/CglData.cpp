// Name:     CglData.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     12/2/06
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

#include "CglData.hpp"

/***********************************************************************/
void CglData::setNrow(const int givenNrow)
{
  nrow = givenNrow;
} /* setNrow */

/***********************************************************************/
void CglData::setNcol(const int givenNcol)
{
  ncol = givenNcol;
} /* setNcol */

/***********************************************************************/
void CglData::setMatrixByColPtr(const CoinPackedMatrix **givenMatrixByColPtr)
{
  matrixByColPtr = givenMatrixByColPtr;
} /* setMatrixByColPtr*/

/***********************************************************************/
void CglData::setMatrixByRowPtr(const CoinPackedMatrix **givenMatrixByRowPtr)
{
  matrixByRowPtr = givenMatrixByRowPtr;
} /* setMatrixByRowPtr*/

/***********************************************************************/
void CglData::setObjPtr(const double **givenObjPtr)
{
  objPtr = givenObjPtr;
} /* setObjPtr*/

/***********************************************************************/
void CglData::setColLowerPtr(const double **givenColLowerPtr)
{
  colLowerPtr = givenColLowerPtr;
} /* setColLowerPtr */

/***********************************************************************/
void CglData::setColUpperPtr(const double **givenColUpperPtr)
{
  colUpperPtr = givenColUpperPtr;
} /* setColUpperPtr */

/***********************************************************************/
void CglData::setRowLowerPtr(const double **givenRowLowerPtr)
{
  rowLowerPtr = givenRowLowerPtr;
} /* setRowLowerPtr */

/***********************************************************************/
void CglData::setRowUpperPtr(const double **givenRowUpperPtr)
{
  rowUpperPtr = givenRowUpperPtr;
} /* setRowUpperPtr */

/***********************************************************************/
void CglData::setRowRhsPtr(const double **givenRowRhsPtr)
{
  rowRhsPtr = givenRowRhsPtr;
} /* setRowUpperPtr */

/***********************************************************************/
void CglData::setRowActivityPtr(const double **givenRowActivityPtr)
{
  rowActivityPtr = givenRowActivityPtr;
} /* setRowActivityPtr */

/***********************************************************************/
void CglData::setColTypePtr(char * const *givenColTypePtr)
{
  colTypePtr = givenColTypePtr;
} /* setColTypePtr */

/***********************************************************************/
void CglData::setSeparateThisPtr(const double **givenSeparateThisPtr)
{
  separateThisPtr = givenSeparateThisPtr;
} /* setSeparateThisPtr */

/***********************************************************************/
void CglData::setDoNotSeparateThisPtr(const 
		  double **givenDoNotSeparateThisPtr)
{
  doNotSeparateThisPtr = givenDoNotSeparateThisPtr;
} /* setDoNotSeparateThisPtr */

/***********************************************************************/
void CglData::setTreeInfoPtr(const CglTreeInfo **givenTreeInfoPtr)
{
  givenTreeInfoPtr = givenTreeInfoPtr;
} /* setTreeInfoPtr */

/***********************************************************************/
CglData::CglData(const int &givenNrow, const int &givenNcol,
		 const CoinPackedMatrix **givenMatrixByColPtr,
		 const CoinPackedMatrix ** givenMatrixByRowPtr,
		 const double **givenObjPtr,
		 const double **givenColLowerPtr, 
		 const double **givenColUpperPtr,   
		 const double **givenRowLowerPtr, 
		 const double **givenRowUpperPtr,
		 const double **givenRowRhsPtr,
		 const double **givenRowActivityPtr,
		 char * const *givenColTypePtr,
		 const double **givenSeparateThisPtr,
		 const CglTreeInfo **givenTreeInfoPtr,
		 const double **givenDoNotSeparateThisPtr) :
  nrow(givenNrow),
  ncol(givenNcol),
  matrixByColPtr(givenMatrixByColPtr),
  matrixByRowPtr(givenMatrixByRowPtr),
  objPtr(givenObjPtr),
  colLowerPtr(givenColLowerPtr),
  colUpperPtr(givenColUpperPtr),
  rowLowerPtr(givenRowLowerPtr),
  rowUpperPtr(givenRowUpperPtr),
  rowRhsPtr(givenRowRhsPtr),
  rowActivityPtr(givenRowActivityPtr),
  colTypePtr(givenColTypePtr),
  separateThisPtr(givenSeparateThisPtr),
  treeInfoPtr(givenTreeInfoPtr),
  doNotSeparateThisPtr(givenDoNotSeparateThisPtr)
{}

/***********************************************************************/
CglData::CglData(const CglData &source) :
  nrow(source.nrow),
  ncol(source.ncol),
  matrixByColPtr(source.matrixByColPtr),
  matrixByRowPtr(source.matrixByRowPtr),
  objPtr(source.objPtr),
  colLowerPtr(source.colLowerPtr),
  colUpperPtr(source.colUpperPtr),
  rowLowerPtr(source.rowLowerPtr),
  rowUpperPtr(source.rowUpperPtr),
  rowRhsPtr(source.rowRhsPtr),
  rowActivityPtr(source.rowActivityPtr),
  colTypePtr(source.colTypePtr),
  separateThisPtr(source.separateThisPtr),
  treeInfoPtr(source.treeInfoPtr),
  doNotSeparateThisPtr(source.doNotSeparateThisPtr)
{}

/***********************************************************************/
CglData* CglData::clone() const
{
  return new CglData(*this);
}

/***********************************************************************/
CglData& CglData::operator=(const CglData &rhs)
{
  if(this != &rhs) {
    nrow = rhs.nrow;
    ncol = rhs.ncol;
    matrixByColPtr = rhs.matrixByColPtr;
    matrixByRowPtr = rhs.matrixByRowPtr;
    objPtr = rhs.objPtr;
    colLowerPtr = rhs.colLowerPtr;
    colUpperPtr = rhs.colUpperPtr;
    rowLowerPtr = rhs.rowLowerPtr;
    rowUpperPtr = rhs.rowUpperPtr;
    rowRhsPtr = rhs.rowRhsPtr;
    rowActivityPtr = rhs.rowActivityPtr;
    colTypePtr = rhs.colTypePtr;
    separateThisPtr = rhs.separateThisPtr; 
    treeInfoPtr = rhs.treeInfoPtr;
    doNotSeparateThisPtr = rhs.doNotSeparateThisPtr; 
  }
  return *this;
}

/***********************************************************************/
CglData::~CglData()
{}
