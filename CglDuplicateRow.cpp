// Copyright (C) 2004, International Business Machines
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
//#define PRINT_DEBUG
//#define CGL_DEBUG
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinFinite.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglDuplicateRow.hpp"
//-------------------------------------------------------------------
// Generate duplicate row column cuts
//------------------------------------------------------------------- 
void CglDuplicateRow::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info) const
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n",nPath);
  }
#endif
  int numberColumns = matrix_.getNumCols();
  CoinPackedVector ubs;
  
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  const double * columnLower = si.getColLower();
  const double * columnUpper = si.getColUpper();
  int nFree=0;
  int nOut=0;
  int nFixed=0;
  int i;
  int numberRows=matrix_.getNumRows();
  int * effectiveRhs = CoinCopyOfArray(rhs_,numberRows);
  double * colUpper2 = new double [numberColumns];
  bool infeasible=false;
  for ( i=0;i<numberColumns;i++) {
    colUpper2[i]=columnUpper[i];
    if (columnLower[i]) {
      double value = columnLower[i];
      for (int jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
        int iRow = row[jj];
        nOut += (int) (element[jj]*value);
        effectiveRhs[iRow] -= (int) (element[jj]*value);
      }
    }
  }
  for ( i=0;i<numberColumns;i++) {
    if (columnLower[i]!=columnUpper[i]) {
      bool fixed=false;
      for (int jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
        int iRow = row[jj];
        if (rhs_[iRow]>=0&&element[jj]>effectiveRhs[iRow]) 
          fixed=true;
      }
      if (fixed) {
        nFixed++;
        colUpper2[i]=columnLower[i];
        ubs.insert(i,columnLower[i]);
      } else {
        nFree++;
      }
    }
  }
  // See if anything odd
  char * check = new char[numberColumns];
  memset(check,0,numberColumns);
  int * which2 = new int[numberColumns];
  for (i=0;i<numberRows;i++) {
    if (effectiveRhs[i]>0)
      duplicate_[i]=-1;
    else if (effectiveRhs[i]==0)
      duplicate_[i]=-2;
    else
      duplicate_[i]=-3;
  }
  for (i=0;i<numberRows;i++) {
    // initially just one
    if (effectiveRhs[i]==1&&duplicate_[i]==-1) {
      int nn=0;
      int j,k;
      for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	double value = elementByRow[j];
	if (columnLower[iColumn]!=colUpper2[iColumn]) {
          assert (value==1.0);
          check[iColumn]=1;
          which2[nn++]=iColumn;
        }
      }
      for ( k=i+1;k<numberRows;k++) {
        if (effectiveRhs[k]==1&&duplicate_[i]==-1) {
          int nn2=0;
          int nnsame=0;
          for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
            int iColumn = column[j];
            double value = elementByRow[j];
            if (columnLower[iColumn]!=colUpper2[iColumn]) {
              assert (value==1.0);
              nn2++;
              if (check[iColumn]) 
                nnsame++;
            }
          }
          if (nnsame==nn2) {
            if (nn2<nn) {
              if (logLevel_)
                printf("row %d strict subset of row %d, fix some in row %d\n",
                       k,i,i);
              // treat k as duplicate
              duplicate_[k]=i;
              // zero out check so we can see what is extra
              for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
                int iColumn = column[j];
                check[iColumn]=0; 
              }
              // now redo and fix
              nn=0;
              for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
                int iColumn = column[j];
                if (columnLower[iColumn]!=colUpper2[iColumn]) {
                  if (check[iColumn]) {
                    // fix
                    colUpper2[iColumn]=columnLower[iColumn];
                    nFixed++;
                    ubs.insert(iColumn,columnLower[iColumn]);
                    check[iColumn]=0;
                  } else {
                    check[iColumn]=1;
                    which2[nn++]=iColumn;
                  }
                }
              }
            } else if (nn2==nn) {
              if (logLevel_)
                printf("row %d identical to row %d\n",
                       k,i);
              duplicate_[k]=i;
            } else if (nn2>=nn) {
              abort();
            }
          } else if (nnsame==nn&&nn2>nn) {
            if (logLevel_)
              printf("row %d strict superset of row %d, fix some in row %d\n",
                     k,i,k);
            // treat k as duplicate
            duplicate_[k]=i;
            // set check for k
            for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
              int iColumn = column[j];
              if (columnLower[iColumn]!=colUpper2[iColumn]) 
                check[iColumn]=1; 
            }
            // zero out check so we can see what is extra
            for ( j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
              int iColumn = column[j];
              check[iColumn]=0; 
            }
            //  fix
            for (j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
              int iColumn = column[j];
              if (check[iColumn]) {
                // fix
                colUpper2[iColumn]=columnLower[iColumn];
                nFixed++;
                ubs.insert(iColumn,columnLower[iColumn]);
                check[iColumn]=0;
              }
            }
            // redo
            nn=0;
            for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
              int iColumn = column[j];
              if (columnLower[iColumn]!=colUpper2[iColumn]) {
                check[iColumn]=1;
                which2[nn++]=iColumn;
              }
            }
          }
        }
      }
      for (k=0;k<nn;k++) 
        check[which2[k]]=0;
      
    }
  }
  delete [] check;
  delete [] which2;
  delete [] colUpper2;
  int nRow=0;
  sizeDynamic_=1;
  for (i=0;i<numberRows;i++) {
    if (duplicate_[i]!=-3) {
      if (duplicate_[i]==-1) {
        nRow++;
        int k=effectiveRhs[i];
        while (k) {
          if (sizeDynamic_<1000000000)
            sizeDynamic_ = sizeDynamic_<<1;
          k = k >>1;
        }
      }
    } else {
      duplicate_[i]=-1;
    }
  } 
  delete [] effectiveRhs;
  if (logLevel_)
    printf("%d free (but %d fixed this time), %d out of rhs, DP size %d, %d rows\n",
           nFree,nFixed,nOut,sizeDynamic_,nRow);
  if (nFixed) {
    OsiColCut cc;
    cc.setUbs(ubs);
    cc.setEffectiveness(100.0);
    cs.insert(cc);
  }
  if (infeasible) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(DBL_MAX);
    rc.setUb(0.0);   
    cs.insert(rc);
  }
}
//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglDuplicateRow::CglDuplicateRow ()
:
  CglCutGenerator(),
  rhs_(NULL),
  duplicate_(NULL),
  lower_(NULL),
  maximumRhs_(1),
  sizeDynamic_(INT_MAX),
  logLevel_(0)
{
}
// Useful constructor
CglDuplicateRow::CglDuplicateRow(OsiSolverInterface * solver)
  : CglCutGenerator(),
    rhs_(NULL),
    duplicate_(NULL),
    lower_(NULL),
    maximumRhs_(1),
    sizeDynamic_(INT_MAX),
    logLevel_(0)
{
  refreshSolver(solver);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglDuplicateRow::CglDuplicateRow (  const CglDuplicateRow & rhs)
                                                              :
  CglCutGenerator(rhs),
  matrix_(rhs.matrix_),
  matrixByRow_(rhs.matrixByRow_),
  maximumRhs_(rhs.maximumRhs_),
  sizeDynamic_(rhs.sizeDynamic_),
  logLevel_(rhs.logLevel_)
{  
  int numberRows=matrix_.getNumRows();
  rhs_ = CoinCopyOfArray(rhs.rhs_,numberRows);
  duplicate_ = CoinCopyOfArray(rhs.duplicate_,numberRows);
  lower_ = CoinCopyOfArray(rhs.lower_,numberRows);
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglDuplicateRow::clone() const
{
  return new CglDuplicateRow(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglDuplicateRow::~CglDuplicateRow ()
{
  // free memory
  delete [] rhs_;
  delete [] duplicate_;
  delete [] lower_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglDuplicateRow &
CglDuplicateRow::operator=(
                                         const CglDuplicateRow& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    delete [] rhs_;
    delete [] duplicate_;
    delete [] lower_;
    matrix_=rhs.matrix_;
    matrixByRow_=rhs.matrixByRow_;
    maximumRhs_=rhs.maximumRhs_;
    sizeDynamic_ = rhs.sizeDynamic_;
    logLevel_ = rhs.logLevel_;
    int numberRows=matrix_.getNumRows();
    rhs_ = CoinCopyOfArray(rhs.rhs_,numberRows);
    duplicate_ = CoinCopyOfArray(rhs.duplicate_,numberRows);
    lower_ = CoinCopyOfArray(rhs.lower_,numberRows);
  }
  return *this;
}

/// This can be used to refresh any inforamtion
void 
CglDuplicateRow::refreshSolver(OsiSolverInterface * solver)
{
  delete [] rhs_;
  delete [] duplicate_;
  delete [] lower_;
  matrix_ = *solver->getMatrixByCol();
  matrix_.removeGaps();
  matrixByRow_ = *solver->getMatrixByRow();
  int numberRows=matrix_.getNumRows();
  rhs_ = new int[numberRows];
  duplicate_ = new int[numberRows];
  lower_ = new int[numberRows];
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  int iRow;
  for (iRow=0;iRow<numberRows;iRow++) {
    rhs_[iRow]=-1;
    lower_[iRow]=0;
    duplicate_[iRow]=-1;
    if (rowUpper[iRow]<100) {
      int iRhs= (int) floor(rowUpper[iRow]);
      // check elements
      bool good=true;
      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int iColumn = column[j];
        if (!solver->isInteger(iColumn))
          good=false;
        double value = elementByRow[j];
        if (floor(value)!=value||value<1.0)
          good=false;
      }
      if (good) {
        lower_[iRow] = (int) CoinMax(0.0,ceil(rowLower[iRow]));
        assert (iRhs>=lower_[iRow]);
        rhs_[iRow]=iRhs;
      }
    }
  }
}
