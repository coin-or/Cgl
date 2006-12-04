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
void CglData::setMatrixByCol(const CoinPackedMatrix *givenMatrixByCol)
{
  matrixByCol = givenMatrixByCol;
} /* setMatrixByCol */

/***********************************************************************/
void CglData::setMatrixByRow(const CoinPackedMatrix *givenMatrixByRow)
{
  matrixByRow = givenMatrixByRow;
} /* setMatrixByRow */

/***********************************************************************/
void CglData::setObj(const double *givenObj)
{
  obj = givenObj;
} /* setObj */

/***********************************************************************/
void CglData::setColLower(const double *givenColLower)
{
  colLower = givenColLower;
} /* setColLower */

/***********************************************************************/
void CglData::setColUpper(const double *givenColUpper)
{
  colUpper = givenColUpper;
} /* setColUpper */

/***********************************************************************/
void CglData::setRowLower(const double *givenRowLower)
{
  rowLower = givenRowLower;
} /* setRowLower */

/***********************************************************************/
void CglData::setRowUpper(const double *givenRowUpper)
{
  rowUpper = givenRowUpper;
} /* setRowUpper */

/***********************************************************************/
void CglData::setRowRhs(const double *givenRowRhs)
{
  rowRhs = givenRowRhs;
} /* setRowUpper */

/***********************************************************************/
void CglData::setRowActivity(const double *givenRowActivity)
{
  rowActivity = givenRowActivity;
} /* setRowActivity */

/***********************************************************************/
void CglData::setColType(const char *givenColType)
{
  colType = givenColType;
} /* setColType */

/***********************************************************************/
void CglData::setSeparateThis(const double *givenSeparateThis)
{
  separateThis = givenSeparateThis;
} /* setSeparateThis */

/***********************************************************************/
void CglData::setDoNotSeparateThis(const 
		  double *givenDoNotSeparateThis)
{
  doNotSeparateThis = givenDoNotSeparateThis;
} /* setDoNotSeparateThis */

/***********************************************************************/
void CglData::setTreeInfo(const CglTreeInfo *givenTreeInfo)
{
  givenTreeInfo = givenTreeInfo;
} /* setTreeInfo */

/***********************************************************************/
CglData::CglData(const int &givenNrow, const int &givenNcol,
		 const CoinPackedMatrix *givenMatrixByCol,
		 const CoinPackedMatrix * givenMatrixByRow,
		 const double *givenObj,
		 const double *givenColLower, 
		 const double *givenColUpper,   
		 const double *givenRowLower, 
		 const double *givenRowUpper,
		 const double *givenRowRhs,
		 const double *givenRowActivity,
		 const char *givenColType,
		 const double *givenSeparateThis,
		 const CglTreeInfo *givenTreeInfo,
		 const double *givenDoNotSeparateThis) :
  nrow(givenNrow),
  ncol(givenNcol),
  matrixByCol(givenMatrixByCol),
  matrixByRow(givenMatrixByRow),
  obj(givenObj),
  colLower(givenColLower),
  colUpper(givenColUpper),
  rowLower(givenRowLower),
  rowUpper(givenRowUpper),
  rowRhs(givenRowRhs),
  rowActivity(givenRowActivity),
  colType(givenColType),
  separateThis(givenSeparateThis),
  treeInfo(givenTreeInfo),
  doNotSeparateThis(givenDoNotSeparateThis)
{}

/***********************************************************************/
CglData::CglData(const CglData &source) :
  nrow(source.nrow),
  ncol(source.ncol),
  matrixByCol(source.matrixByCol),
  matrixByRow(source.matrixByRow),
  obj(source.obj),
  colLower(source.colLower),
  colUpper(source.colUpper),
  rowLower(source.rowLower),
  rowUpper(source.rowUpper),
  rowRhs(source.rowRhs),
  rowActivity(source.rowActivity),
  colType(source.colType),
  separateThis(source.separateThis),
  treeInfo(source.treeInfo),
  doNotSeparateThis(source.doNotSeparateThis)
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
    matrixByCol = rhs.matrixByCol;
    matrixByRow = rhs.matrixByRow;
    obj = rhs.obj;
    colLower = rhs.colLower;
    colUpper = rhs.colUpper;
    rowLower = rhs.rowLower;
    rowUpper = rhs.rowUpper;
    rowRhs = rhs.rowRhs;
    rowActivity = rhs.rowActivity;
    colType = rhs.colType;
    separateThis = rhs.separateThis; 
    treeInfo = rhs.treeInfo;
    doNotSeparateThis = rhs.doNotSeparateThis; 
  }
  return *this;
}

/***********************************************************************/
CglData::~CglData()
{}
