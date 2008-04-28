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
#include "CglStored.hpp"
//-------------------------------------------------------------------
// Generate duplicate row column cuts
//------------------------------------------------------------------- 
void CglDuplicateRow::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info) const
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path\n");
  }
#endif
  // Don't do in tree ?
  if (info.inTree) {
    // but do any stored cuts
    if (storedCuts_)
      storedCuts_->generateCuts(si,cs,info);
    return;
  }
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
  const double * rowLower = si.getRowLower();
  const double * rowUpper = si.getRowUpper();
  int * effectiveRhs = CoinCopyOfArray(rhs_,numberRows);
  int * effectiveLower = CoinCopyOfArray(lower_,numberRows);
  // mark bad rows - also used for domination
  for (i=0;i<numberRows;i++) {
    duplicate_[i]=-1;
    int j;
    for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
      if (elementByRow[j]!=1.0) {
        duplicate_[i]=-3;
        rhs_[i]=-1000000;
        break;
      }
    }
  }
  double * colUpper2 = CoinCopyOfArray(columnUpper,numberColumns);
  if (!info.pass&&(mode_&2)!=0) {
    // First look at duplicate or dominated columns
    double * random = new double[numberRows];
    double * sort = new double[numberColumns+1];
    if (info.randomNumberGenerator) {
      const CoinThreadRandom * randomGenerator = info.randomNumberGenerator;
      for (i=0;i<numberRows;i++) {
	if (rowLower[i]<-1.0e20||rowUpper[i]>1.0e20)
	  random[i]=0.0;
	else
	  random[i] = randomGenerator->randomDouble();
      }
    } else {
      for (i=0;i<numberRows;i++) {
	if (rowLower[i]<-1.0e20||rowUpper[i]>1.0e20)
	  random[i]=0.0;
	else
	  random[i] = CoinDrand48();
      }
    }
    int * which = new int[numberColumns];
    int nPossible=0;
    for ( i=0;i<numberColumns;i++) {
      if (si.isBinary(i)) {
	double value = 0.0;
	for (int jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
	  int iRow = row[jj];
	  value += element[jj]*random[iRow];
	}
	sort[nPossible]=value;
	which[nPossible++]=i;
      }
    }
    sort[nPossible]=COIN_DBL_MAX;
    CoinSort_2(sort,sort+nPossible,which);
    int last=maximumDominated_-1;
    double lastValue=-1.0;
    const double *objective = si.getObjCoefficients() ;
    double direction = si.getObjSense();
    // arrays for checking
    double * elementEqualJ = new double [2*numberRows];
    CoinZeroN(elementEqualJ,numberRows); // for memory checkers
    double * elementGeJ = elementEqualJ + numberRows;
    CoinZeroN(elementGeJ,numberRows);
    int * rowEqualJ = new int[2*numberRows];
    CoinZeroN(rowEqualJ,numberRows); // for memory checkers
    int * rowGeJ = rowEqualJ + numberRows;
    char * mark = new char[numberRows];
    CoinZeroN(mark,numberRows);
    for (i=0;i<nPossible+1;i++) {
      if (sort[i]>lastValue) {
	if (i-last<=maximumDominated_) {
	  // look to see if dominated
	  for (int j=last;j<i;j++) {
	    int jColumn = which[j];
	    // skip if already fixed
	    if (!colUpper2[jColumn])
	      continue;
	    int nGeJ=0;
	    int nEqualJ=0;
	    int jj;
	    for (jj=columnStart[jColumn];jj<columnStart[jColumn]+columnLength[jColumn];jj++) {
	      int iRow = row[jj];
	      if (random[iRow]) {
		elementEqualJ[nEqualJ]=element[jj];
		rowEqualJ[nEqualJ++]=iRow;
	      } else {
		// swap sign so all rows look like G
		elementGeJ[iRow]=(rowUpper[iRow]>1.0e20) ? element[jj] : -element[jj];
		rowGeJ[nGeJ++]=iRow;
	      }
	    }
	    double objValueJ = objective[jColumn]*direction;
	    for (int k=j+1;k<i;k++) {
	      int kColumn = which[k];
	      // skip if already fixed
	      if (!colUpper2[kColumn])
		continue;
	      int nEqualK=0;
	      // -2 no good, -1 J dominates K, 0 unknown or equal, 1 K dominates J
	      int dominate=0;
	      // mark
	      int kk;
	      for (kk=0;kk<nGeJ;kk++)
		mark[rowGeJ[kk]]=1;
	      for (kk=columnStart[kColumn];kk<columnStart[kColumn]+columnLength[kColumn];kk++) {
		int iRow = row[kk];
		if (random[iRow]) {
		  if (iRow!=rowEqualJ[nEqualK]||
		      element[kk]!=elementEqualJ[nEqualK]) {
		    dominate=-2;
		    break;
		  } else {
		    nEqualK++;
		  }
		} else {
		  // swap sign so all rows look like G
		  double valueK = (rowUpper[iRow]>1.0e20) ? element[kk] : -element[kk];
		  double valueJ = elementGeJ[iRow];
		  mark[iRow]=0;
		  if (valueJ==valueK) {
		    // equal
		  } else if (valueJ>valueK) {
		    // J would dominate K
		    if (dominate==1) {
		      // no good
		      dominate=-2;
		      break;
		    } else {
		      dominate=-1;
		    }
		  } else {
		    // K would dominate J
		    if (dominate==-1) {
		      // no good
		      dominate=-2;
		      break;
		    } else {
		      dominate=1;
		    }
		  }
		}
	      }
	      kk=0;
	      if (dominate!=-2) {
		// unmark and check
		for (;kk<nGeJ;kk++) {
		  int iRow = rowGeJ[kk];
		  if (mark[iRow]) {
		    double valueK = 0.0;
		    double valueJ = elementGeJ[iRow];
		    if (valueJ>valueK) {
		      // J would dominate K
		      if (dominate==1) {
			// no good
			dominate=-2;
			break;
		      } else {
			dominate=-1;
		      }
		    } else {
		      // K would dominate J
		      if (dominate==-1) {
			// no good
			dominate=-2;
			break;
		      } else {
			dominate=1;
		      }
		    }
		  }
		  mark[iRow]=0;
		}
	      } 
	      // just unmark rest
	      for (;kk<nGeJ;kk++)
		mark[rowGeJ[kk]]=0;
	      if (nEqualK==nEqualJ&&dominate!=-2) {
		double objValueK = objective[kColumn]*direction;
		if (objValueJ==objValueK) {
		  if (dominate<=0) {
		    // say J dominates
		    assert (colUpper2[kColumn]);
		    dominate=-1;
		  } else {
		    // say K dominates
		    assert (colUpper2[jColumn]);
		    dominate=1;
		  }
		} else if (objValueJ<objValueK&&dominate<=0) {
		  // say J dominates
		  assert (colUpper2[kColumn]);
		  dominate=-1;
		} else if (objValueJ>objValueK&&dominate==1) {
		  // say K dominates
		  assert (colUpper2[jColumn]);
		  dominate=1;
		} else {
		  dominate=0;
		}
		if (dominate) {
		  // see if both can be 1
		  bool canFix=false;
		  for (int jj=0;jj<nEqualJ;jj++) {
		    double value = 2.0*elementEqualJ[jj];
		    int iRow = rowEqualJ[jj];
		    if (duplicate_[iRow]==-1&&rowUpper[iRow]<1.999999) {
		      canFix=true;
		    } else {
		      double minSum=0.0;
		      double maxSum=0.0;
		      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
			int iColumn = column[j];
			if (iColumn!=jColumn&&iColumn!=kColumn) {
			  double elValue = elementByRow[j];
			  double lo = columnLower[iColumn];
			  double up = colUpper2[iColumn];
			  if (elValue>0.0) {
			    minSum += lo*elValue;
			    maxSum += up*elValue;
			  } else {
			    maxSum += lo*elValue;
			    minSum += up*elValue;
			  }
			}
		      }
		      if (minSum+value>rowUpper[iRow]+1.0e-5)
			canFix=true;
		      else if (maxSum+value<rowLower[iRow]-1.0e-5)
			canFix=true;
		    }
		    if (canFix)
		      break;
		  }
		  if (canFix) {
		    int iColumn = (dominate>0) ? jColumn : kColumn;
		    nFixed++;
		    assert (!columnLower[iColumn]);
		    colUpper2[iColumn]=0.0;
		    ubs.insert(iColumn,0.0);
		    if (iColumn==jColumn)
		      break; // no need to carry on on jColumn
		  } else {
		    int iDominated = (dominate>0) ? jColumn : kColumn;
		    int iDominating = (dominate<0) ? jColumn : kColumn;
		    double els[]={1.0,-1.0};
		    int inds[2];
		    inds[0]=iDominating;
		    inds[1]=iDominated;
		    if (!storedCuts_)
		      storedCuts_ = new CglStored();
		    storedCuts_->addCut(0.0,COIN_DBL_MAX,2,inds,els);
		  }
		}
	      }
	    }
	    for (jj=0;jj<nGeJ;jj++) {
	      int iRow = rowGeJ[jj];
	      elementGeJ[iRow]=0.0;
	    }
	  }
	}
	last=i;
	lastValue = sort[i];
      }
    }
    delete [] mark;
    delete [] elementEqualJ;
    delete [] rowEqualJ;
    delete [] random;
    delete [] sort;
    delete [] which;
#ifdef COIN_DEVELOP
    int numberCuts = storedCuts_ ? storedCuts_->sizeRowCuts() : 0;
    if (nFixed||numberCuts) 
      printf("** %d fixed and %d cuts from domination\n",nFixed,numberCuts);
#endif
  }
  bool infeasible=false;
  // if we were just doing columns - mark all as bad
  if ((mode_&1)==0) {
    for (i=0;i<numberRows;i++) {
      duplicate_[i]=-3;
      rhs_[i]=-1000000;
      effectiveLower[i]=-1000000;
    }
  }
  for ( i=0;i<numberColumns;i++) {
    if (columnLower[i]) {
      double value = columnLower[i];
      for (int jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
        int iRow = row[jj];
        nOut += (int) (element[jj]*value);
        effectiveRhs[iRow] -= (int) (element[jj]*value);
        effectiveLower[iRow] -= (int) (element[jj]*value);
      }
    }
  }
  for ( i=0;i<numberColumns;i++) {
    if (columnLower[i]!=colUpper2[i]) {
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
    if (duplicate_[i]==-1) {
      if (effectiveRhs[i]>0) {
	// leave
      } else if (effectiveRhs[i]==0) {
        duplicate_[i]=-2;
      } else {
        duplicate_[i]=-3;
	// leave unless >=1 row
	if (effectiveLower[i]==1&&rhs_[i]<0.0)
	  duplicate_[i]=-4;
      }
    } else {
      effectiveRhs[i]=-1000;
    }
  }
  // Look at <= rows
  for (i=0;i<numberRows;i++) {
    // initially just one
    if (effectiveRhs[i]==1&&duplicate_[i]==-1) {
      int nn=0;
      int j,k;
      for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
          assert (elementByRow[j]==1.0);
#endif
          check[iColumn]=1;
          which2[nn++]=iColumn;
        }
      }
      for ( k=i+1;k<numberRows;k++) {
        if (effectiveRhs[k]==1&&duplicate_[k]==-1) {
          int nn2=0;
          int nnsame=0;
          for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
            int iColumn = column[j];
            if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
              assert (elementByRow[j]==1.0);
#endif
              nn2++;
              if (check[iColumn]) 
                nnsame++;
            }
          }
	  //if (nnsame)
	  //printf("rows %d and %d, %d same - %d %d\n",
	  //   i,k,nnsame,nn,nn2);
          if (nnsame==nn2) {
            if (nn2<nn&&effectiveLower[k]==rhs_[k]&&rhs_[i]==rhs_[k]) {
              if (logLevel_)
                printf("row %d strict subset of row %d, fix some in row %d\n",
                       k,i,i);
              // treat i as duplicate
              duplicate_[i]=k;
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
            } else if (nn2==nn&&effectiveLower[i]==rhs_[i]&&effectiveLower[k]==rhs_[k]) {
              if (logLevel_)
                printf("row %d identical to row %d\n",
                       k,i);
              duplicate_[k]=i;
            } else if (nn2>=nn&&effectiveLower[i]==rhs_[i]&&effectiveLower[k]==rhs_[k]) {
              abort();
            }
          } else if (nnsame==nn&&nn2>nn&&effectiveLower[i]==rhs_[i]&&rhs_[i]<=rhs_[k]) {
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
          } else {
            // may be redundant
            if (nnsame==nn2) {
              // k redundant ?
              if (nn2<nn&&effectiveLower[k]<=0&&rhs_[i]<=rhs_[k]) {
                if (logLevel_)
                  printf("row %d slack subset of row %d, drop row %d\n",
                         k,i,k);
                // treat k as duplicate
                duplicate_[k]=i;
              }
            } else if (nnsame==nn) {
              // i redundant ?
              if (nn2>nn&&effectiveLower[i]<=0&&rhs_[k]<=rhs_[i]) {
                if (logLevel_)
                  printf("row %d slack subset of row %d, drop row %d\n",
                         i,k,i);
                // treat i as duplicate
                duplicate_[i]=k;
              }
            }
          }
        }
      }
      for (k=0;k<nn;k++) 
        check[which2[k]]=0;
      
    }
  }
  // Look at >=1 rows
  for (i=0;i<numberRows;i++) {
    if (duplicate_[i]==-4) {
      int nn=0;
      int j,k;
      for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
          assert (elementByRow[j]==1.0);
#endif
          check[iColumn]=1;
          which2[nn++]=iColumn;
        }
      }
      for ( k=i+1;k<numberRows;k++) {
        if (duplicate_[k]==-4) {
          int nn2=0;
          int nnsame=0;
          for ( j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
            int iColumn = column[j];
            if (columnLower[iColumn]!=colUpper2[iColumn]) {
#ifndef NDEBUG
              assert (elementByRow[j]==1.0);
#endif
              nn2++;
              if (check[iColumn]) 
                nnsame++;
            }
          }
	  // may be redundant
	  if (nnsame==nn||nnsame==nn2) {
	    if (nn2>nn) {
	      // k redundant
	      if (logLevel_) 
		printf("row %d slack superset of row %d, drop row %d\n",
		       k,i,k);
	      // treat k as duplicate
	      duplicate_[k]=i;
	    } else if (nn2<nn) {
	      // i redundant ?
	      if (logLevel_)
		printf("row %d slack superset of row %d, drop row %d\n",
		     i,k,i);
	      // treat i as duplicate
	      duplicate_[i]=k;
	    } else {
	      if (logLevel_) 
		printf("row %d same as row %d, drop row %d\n",
		       k,i,k);
	      // treat k as duplicate
	      duplicate_[k]=i;
	    }
          }
        }
      }
      for (k=0;k<nn;k++) 
        check[which2[k]]=0;
      
    }
  }
  if ((mode_&1)!=0&&true) {
    // look at doubletons
    const double * rowLower = si.getRowLower();
    const double * rowUpper = si.getRowUpper();
    int i;
    int nPossible=0;
    for (i=0;i<numberRows;i++) {
      if (rowLength[i]==2&&(duplicate_[i]<0&&duplicate_[i]!=-2)) {
	bool possible=true;
	int j;
	for (j=rowStart[i];j<rowStart[i]+2;j++) {
	  int iColumn = column[j];
	  if (fabs(elementByRow[j])!=1.0||!si.isInteger(iColumn)) {
	    possible=false;
	    break;
	  }
	}
	if (possible) {
	  int j = rowStart[i];
	  int column0 = column[j];
	  double element0 = elementByRow[j];
	  int column1 = column[j+1];
	  double element1 = elementByRow[j+1];
	  if (element0==1.0&&element1==1.0&&rowLower[i]==1.0&&
	      rowUpper[i]>1.0e30) {
	    if (logLevel_) {
	      printf("Cover row %d %g <= ",i,rowLower[i]);
	      printf("(%d,%g) (%d,%g) ",column0,element0,column1,element1);
	      printf(" <= %g\n",rowUpper[i]);
	    }
	    effectiveRhs[nPossible++]=i;
	  } else {
	    // not at present
	    //printf("NON Cover row %d %g <= ",i,rowLower[i]);
	    //printf("(%d,%g) (%d,%g) ",column0,element0,column1,element1);
	    //printf(" <= %g\n",rowUpper[i]);
	  }
	}
      }
    }
    if (nPossible) {
      int * check2 = new int [numberColumns];
      CoinFillN(check2,numberColumns,-1);
      for (int iPossible=0;iPossible<nPossible;iPossible++) {
#ifndef NDEBUG
	for (i=0;i<numberColumns;i++)
	  assert (check2[i]==-1);
#endif
	i=effectiveRhs[iPossible];
	int j = rowStart[i];
	int column0 = column[j];
	int column1 = column[j+1];
	int k;
	int nMarked=0;
	for (int kPossible=iPossible+1;kPossible<nPossible;kPossible++) {
	  k=effectiveRhs[kPossible];
	  int j = rowStart[k];
	  int columnB0 = column[j];
	  int columnB1 = column[j+1];
	  if (column0==columnB1||column1==columnB1) {
	    columnB1=columnB0;
	    columnB0=column[j+1];
	  }
	  bool good = false;
	  if (column0==columnB0) {
	    if (column1==columnB1) {
	      // probably should have been picked up
	      // safest to ignore
	    } else {
	      good=true;
	    }
	  } else if (column1==columnB0) {
	    if (column0==columnB1) {
	      // probably should have been picked up
	      // safest to ignore
	    } else {
	      good=true;
	    }
	  }
	  if (good) {
	    if (check2[columnB1]<0) {
	      check2[columnB1]=k;
	      which2[nMarked++]=columnB1;
	    } else {
	      // found
#ifndef COIN_DEVELOP
	      if (logLevel_>1)
#endif 
		printf("***Make %d %d %d >=2 and take out rows %d %d %d\n",
		       columnB1,column0,column1,
		       i,k,check2[columnB1]);
	      OsiRowCut rc;
	      rc.setLb(2.0);
	      rc.setUb(DBL_MAX);   
	      int index[3];
	      double element[3]={1.0,1.0,1.0};
	      index[0]=column0;
	      index[1]=column1;
	      index[2]=columnB1;
	      rc.setRow(3,index,element,false);
	      cs.insert(rc);
	      // drop rows
	      duplicate_[i]=-2;
	      duplicate_[k]=-2;
	      duplicate_[check2[columnB1]]=-2;
	    }
	  }
	}
	for (k=0;k<nMarked;k++) {
	  int iColumn = which2[k];
	  check2[iColumn]=-1;
	}
      }
      delete [] check2;
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
  delete [] effectiveLower;
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
  storedCuts_(NULL),
  maximumDominated_(1000),
  maximumRhs_(1),
  sizeDynamic_(COIN_INT_MAX),
  mode_(3),
  logLevel_(0)
{
}
// Useful constructor
CglDuplicateRow::CglDuplicateRow(OsiSolverInterface * solver)
  : CglCutGenerator(),
    rhs_(NULL),
    duplicate_(NULL),
    lower_(NULL),
    storedCuts_(NULL),
    maximumDominated_(1000),
    maximumRhs_(1),
    sizeDynamic_(COIN_INT_MAX),
    mode_(3),
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
  storedCuts_(NULL),
  maximumDominated_(rhs.maximumDominated_),
  maximumRhs_(rhs.maximumRhs_),
  sizeDynamic_(rhs.sizeDynamic_),
  mode_(rhs.mode_),
  logLevel_(rhs.logLevel_)
{  
  int numberRows=matrix_.getNumRows();
  rhs_ = CoinCopyOfArray(rhs.rhs_,numberRows);
  duplicate_ = CoinCopyOfArray(rhs.duplicate_,numberRows);
  lower_ = CoinCopyOfArray(rhs.lower_,numberRows);
  if (rhs.storedCuts_)
    storedCuts_ = new CglStored(*rhs.storedCuts_);
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
  delete storedCuts_;
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
    delete storedCuts_;
    storedCuts_ = NULL;
    matrix_=rhs.matrix_;
    matrixByRow_=rhs.matrixByRow_;
    maximumDominated_ = rhs.maximumDominated_;
    maximumRhs_=rhs.maximumRhs_;
    sizeDynamic_ = rhs.sizeDynamic_;
    mode_ = rhs.mode_;
    logLevel_ = rhs.logLevel_;
    int numberRows=matrix_.getNumRows();
    rhs_ = CoinCopyOfArray(rhs.rhs_,numberRows);
    duplicate_ = CoinCopyOfArray(rhs.duplicate_,numberRows);
    lower_ = CoinCopyOfArray(rhs.lower_,numberRows);
  if (rhs.storedCuts_)
    storedCuts_ = new CglStored(*rhs.storedCuts_);
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
  matrix_.orderMatrix();
  matrixByRow_ = *solver->getMatrixByRow();
  int numberRows=matrix_.getNumRows();
  rhs_ = new int[numberRows];
  duplicate_ = new int[numberRows];
  lower_ = new int[numberRows];
  const double * columnLower = solver->getColLower();
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();
  int iRow;
  int numberGood=0;
  int markBad = -(solver->getNumCols()+1);
  for (iRow=0;iRow<numberRows;iRow++) {
    rhs_[iRow]=markBad;
    lower_[iRow]=markBad;
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
        if (floor(value)!=value||value<1.0) {
          good=false;
        }
      }
      if (good) {
        lower_[iRow] = (int) CoinMax(0.0,ceil(rowLower[iRow]));
        if (iRhs>=lower_[iRow]) {
          rhs_[iRow]=iRhs;
          numberGood++;
        } else {
          // infeasible ?
          lower_[iRow]=markBad;
          rhs_[iRow]=markBad;
        }
      } else {
        lower_[iRow]=markBad;
        rhs_[iRow]=markBad;
      }
    } else if (rowUpper[iRow]>1.0e30&&rowLower[iRow]==1.0) {
      // may be OK to look for dominated in >=1 rows
      // check elements
      bool good=true;
      for (int j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int iColumn = column[j];
        if (!solver->isInteger(iColumn))
          good=false;
        double value = elementByRow[j];
        if (floor(value)!=value||value<1.0) {
          good=false;
        }
	if (columnLower[iColumn]!=0.0)
	  good=false;
      }
      if (good) {
        lower_[iRow] = 1;
      }
    }
  }
}
  /** Fix variables and find duplicate/dominated rows for the model of the 
      solver interface, si.

      This is a very simple minded idea but I (JJF) am using it in a project so thought
      I might as well add it.  It should really be called before first solve and I may
      modify CBC to allow for that.

      This is designed for problems with few rows and many integer variables where the rhs
      are <= or == and all coefficients and rhs are small integers.

      If effective rhs is K then we can fix all variables with coefficients > K to their lower bounds
      (effective rhs just means original with variables with nonzero lower bounds subtracted out).

      If one row is a subset of another and the effective rhs are same we can fix some variables
      and then the two rows are identical.

      This version does deletions and fixings and may return stored cuts for
      dominated columns 
  */
CglStored * 
CglDuplicateRow::outDuplicates( OsiSolverInterface * solver)
{
  
  CglTreeInfo info;
  info.level = 0;
  info.pass = 0;
  int numberRows = solver->getNumRows();
  info.formulation_rows = numberRows;
  info.inTree = false;
  info.strengthenRow= NULL;
  info.pass = 0;
  OsiCuts cs;
  generateCuts(*solver,cs,info);
  // Get rid of duplicate rows
  int * which = new int[numberRows]; 
  int numberDrop=0;
  for (int iRow=0;iRow<numberRows;iRow++) {
    if (duplicate_[iRow]==-2||duplicate_[iRow]>=0) 
      which[numberDrop++]=iRow;
  }
  if (numberDrop) {
    solver->deleteRows(numberDrop,which);
  }
  delete [] which;
  // see if we have any column cuts
  int numberColumnCuts = cs.sizeColCuts() ;
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  for (int k = 0;k<numberColumnCuts;k++) {
    OsiColCut * thisCut = cs.colCutPtr(k) ;
    const CoinPackedVector & lbs = thisCut->lbs() ;
    const CoinPackedVector & ubs = thisCut->ubs() ;
    int j ;
    int n ;
    const int * which ;
    const double * values ;
    n = lbs.getNumElements() ;
    which = lbs.getIndices() ;
    values = lbs.getElements() ;
    for (j = 0;j<n;j++) {
      int iColumn = which[j] ;
      if (values[j]>columnLower[iColumn]) 
        solver->setColLower(iColumn,values[j]) ;
    }
    n = ubs.getNumElements() ;
    which = ubs.getIndices() ;
    values = ubs.getElements() ;
    for (j = 0;j<n;j++) {
      int iColumn = which[j] ;
      if (values[j]<columnUpper[iColumn]) 
        solver->setColUpper(iColumn,values[j]) ;
    }
  }
  return storedCuts_;
}
// Create C++ lines to get to current state
std::string
CglDuplicateRow::generateCpp( FILE * fp) 
{
  CglDuplicateRow other;
  fprintf(fp,"0#include \"CglDuplicateRow.hpp\"\n");
  fprintf(fp,"3  CglDuplicateRow duplicateRow;\n");
  if (logLevel_!=other.logLevel_)
    fprintf(fp,"3  duplicateRow.setLogLevel(%d);\n",logLevel_);
  else
    fprintf(fp,"4  duplicateRow.setLogLevel(%d);\n",logLevel_);
  if (maximumRhs_!=other.maximumRhs_)
    fprintf(fp,"3  duplicateRow.setMaximumRhs(%d);\n",maximumRhs_);
  else
    fprintf(fp,"4  duplicateRow.setMaximumRhs(%d);\n",maximumRhs_);
  if (maximumDominated_!=other.maximumDominated_)
    fprintf(fp,"3  duplicateRow.setMaximumDominated(%d);\n",maximumDominated_);
  else
    fprintf(fp,"4  duplicateRow.setMaximumDominated(%d);\n",maximumDominated_);
  if (mode_!=other.mode_)
    fprintf(fp,"3  duplicateRow.setMode(%d);\n",mode_);
  else
    fprintf(fp,"4  duplicateRow.setMode(%d);\n",mode_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  duplicateRow.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  duplicateRow.setAggressiveness(%d);\n",getAggressiveness());
  return "duplicateRow";
}
