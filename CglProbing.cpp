
// Copyright (C) 2002, International Business Machines
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
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinFinite.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglProbing.hpp"

typedef struct {double infeasibility;int sequence;} double_int_pair;
class double_int_pair_compare {
public:
  bool operator() (double_int_pair x , double_int_pair y) const
  {
    return ( x.infeasibility < y.infeasibility);
  }
}; 

#ifdef CGL_DEBUG
// Checks bounds okay against debugger
static void checkBounds(const OsiRowCutDebugger * debugger,OsiColCut & cut)
{
  if (debugger) {
    // on optimal path
    const double * optimal = debugger->optimalSolution();
    int i;
    int nIndex;
    const double * values;
    const int * index;
    const CoinPackedVector & lbs = cut.lbs();
    values = lbs.getElements();
    nIndex = lbs.getNumElements();
    index = lbs.getIndices();
    for (i=0;i<nIndex;i++) {
      double value=values[i];
      int iColumn = index[i];
      printf("%d optimal %g lower %g\n",iColumn,optimal[iColumn],value);
      assert(value<=optimal[iColumn]+1.0e-5);
    }
    const CoinPackedVector & ubs = cut.ubs();
    values = ubs.getElements();
    nIndex = ubs.getNumElements();
    index = ubs.getIndices();
    for (i=0;i<nIndex;i++) {
      double value=values[i];
      int iColumn = index[i];
      printf("%d optimal %g upper %g\n",iColumn,optimal[iColumn],value);
      assert(value>=optimal[iColumn]-1.0e-5);
    }
  }
}
#endif
// This tightens column bounds (and can declare infeasibility)
// It may also declare rows to be redundant
static int tighten(double *colLower, double * colUpper,
		   const int *column, const double *rowElements, 
		   const int *rowStart, const int * rowLength,
		   double *rowLower, double *rowUpper, 
		   int nRows,int nCols,char * intVar,int maxpass,
		   double tolerance)
{
  int i, j, k, kre;
  int krs;
  int dolrows;
  int iflagu, iflagl;
  int ntotal=0,nchange=1,jpass=0;
  double dmaxup, dmaxdown, dbound;
  int ninfeas=0;
  
  while(nchange) {
    int ilbred = 0; /* bounds reduced */
    int iubred = 0; /* bounds reduced */
    int nrwdrp = 0; /* redundant rows */
    if (jpass==maxpass) break;
    jpass++;
    dolrows = (jpass & 1) == 1;
    
    for (i = 0; i < nRows; ++i) {
      if (rowLower[i]>-1.0e20||rowUpper[i]<1.0e20) {
	iflagu = 0;
	iflagl = 0;
	dmaxup = 0.0;
	dmaxdown = 0.0;
	krs = rowStart[i];
	kre = rowStart[i]+rowLength[i];

	/* ------------------------------------------------------------*/
	/* Compute L(i) and U(i) */
	/* ------------------------------------------------------------*/
	for (k = krs; k < kre; ++k) {
	  double value=rowElements[k];
	  j = column[k];
	  if (value > 0.0) {
	    if (colUpper[j] >= 1e15) {
	      dmaxup = 1e31;
	      ++iflagu;
	    } else {
	      dmaxup += colUpper[j] * value;
	    }
	    if (colLower[j] <= -1e15) {
	      dmaxdown = -1e31;
	      ++iflagl;
	    } else {
	      dmaxdown += colLower[j] * value;
	    }
	  } else if (value<0.0) {
	    if (colUpper[j] >= 1e15) {
	      dmaxdown = -1e31;
	      ++iflagl;
	    } else {
	      dmaxdown += colUpper[j] * value;
	    }
	    if (colLower[j] <= -1e15) {
	      dmaxup = 1e31;
	      ++iflagu;
	    } else {
	      dmaxup += colLower[j] * value;
	    }
	  }
	}
	if (dmaxup <= rowUpper[i] + tolerance && dmaxdown >= rowLower[i] - tolerance) {
	  /*
	   * The sum of the column maxs is at most the row ub, and
	   * the sum of the column mins is at least the row lb;
	   * this row says nothing at all.
	   * I suspect that this corresponds to
	   * an implied column singleton in the paper (viii, on p. 325),
	   * where the singleton in question is the row slack.
	   */
	  ++nrwdrp;
	  rowLower[i]=-DBL_MAX;
	  rowUpper[i]=DBL_MAX;
	} else {
	  if (dmaxup < rowLower[i] -tolerance || dmaxdown > rowUpper[i]+tolerance) {
	    ninfeas++;
	    break;
	  }
	  /*        Finite U(i) */
	  /* -------------------------------------------------------------*/
	  /* below is deliberate mistake (previously was by chance) */
	  /*        never do both */
	  if (iflagu == 0 && rowLower[i] > 0.0 && iflagl == 0 && rowUpper[i] < 1e15)
	    {
	      if (dolrows) {
		iflagu = 1;
	      } else {
		iflagl = 1;
	      }
	    }
	  if (iflagu == 0 && rowLower[i] > -1e15) {
	    for (k = krs; k < kre; ++k) {
	      double value=rowElements[k];
	      j = column[k];
	      if (value > 0.0) {
		if (colUpper[j] < 1e15) {
		  dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
		  if (dbound > colLower[j] + 1.0e-12) {
		    /* we can tighten the lower bound */
		    /* the paper mentions this as a possibility on p. 227 */
		    colLower[j] = dbound;
		    ++ilbred;
		    
		    /* this may have fixed the variable */
		    /* I believe that this roughly corresponds to a
		     * forcing constraint in the paper (p. 226).
		     * If there is a forcing constraint (with respect
		     * to the original, untightened bounds), then in this 
		     * loop we will go through all the columns and fix
		     * each of them to their implied bound, rather than
		     * determining that the row as a whole is forced
		     * and just fixing them without doing computation for
		     * each column (as in the paper).
		     * By doing it this way, we can tighten bounds and
		     * get forcing constraints for free.
		     */
		    if (colUpper[j] - colLower[j] <= tolerance) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      } else {
		if (colLower[j] > -1e15) {
		  dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
		  if (dbound < colUpper[j] - 1.0e-12) {
		    colUpper[j] = dbound;
		    ++iubred;
		    if (colUpper[j] - colLower[j] <= tolerance) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  /* ----------------------------------------------------------------*/
	  /*        Finite L(i) */
	  /* ----------------------------------------------------------------*/
	  if (iflagl == 0 && rowUpper[i] < 1e15) {
	    for (k = krs; k < kre; ++k) {
	      double value=rowElements[k];
	      j = column[k];
	      if (value < 0.0) {
		if (colUpper[j] < 1e15) {
		  dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
		  if (dbound > colLower[j] + 1.0e-12) {
		    colLower[j] = dbound;
		    ++ilbred;
		    if (! (colUpper[j] - colLower[j] > tolerance)) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		} 
	      } else {
		if (colLower[j] > -1e15) {
		  dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
		  if (dbound < colUpper[j] - 1.0e-12) {
		    colUpper[j] = dbound;
		    ++iubred;
		    if (! (colUpper[j] - colLower[j] > tolerance)) {
		      /* --------------------------------------------------*/
		      /*                check if infeasible !!!!! */
		      /* --------------------------------------------------*/
		      if (colUpper[j] - colLower[j] < -100.0*tolerance) {
			ninfeas++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    for (j = 0; j < nCols; ++j) {
      if (intVar[j]) {
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  if (floor(colUpper[j]+1.0e-4)<colUpper[j]) 
	    nchange++;
	  // clean up anyway
	  colUpper[j]=floor(colUpper[j]+1.0e-4);
	  if (ceil(colLower[j]-1.0e-4)>colLower[j]) 
	    nchange++;
	  // clean up anyway
	  colLower[j]=ceil(colLower[j]-1.0e-4);
	  if (colUpper[j]<colLower[j]) {
	    /*printf("infeasible\n");*/
	    ninfeas++;
	  }
	}
      }
    }
    nchange=ilbred+iubred+nrwdrp;
    ntotal += nchange;
    if (ninfeas) break;
  }
  return (ninfeas);
}
// This just sets minima and maxima on rows
static void tighten2(double *colLower, double * colUpper,
		     const int *column, const double *rowElements, 
		     const int *rowStart, const int * rowLength,
		     double *rowLower, double *rowUpper, 
		     double * minR, double * maxR, int * markR,
		     int nRows,int nCols)
{
  int i, j, k, kre;
  int krs;
  int iflagu, iflagl;
  double dmaxup, dmaxdown, value;

  for (i = 0; i < nRows; ++i) {
    if (rowLower[i]>-1.0e20||rowUpper[i]<1.0e20) {
      iflagu = 0;
      iflagl = 0;
      dmaxup = 0.0;
      dmaxdown = 0.0;
      krs = rowStart[i];
      kre = rowStart[i]+rowLength[i];
      
      /* ------------------------------------------------------------*/
      /* Compute L(i) and U(i) */
      /* ------------------------------------------------------------*/
      for (k = krs; k < kre; ++k) {
	value=rowElements[k];
	j = column[k];
	if (value > 0.0) {
	  if (colUpper[j] >= 1e15) {
	    dmaxup = 1e31;
	    ++iflagu;
	  } else {
	    dmaxup += colUpper[j] * value;
	  }
	  if (colLower[j] <= -1e15) {
	    dmaxdown = -1e31;
	    ++iflagl;
	  } else {
	    dmaxdown += colLower[j] * value;
	  }
	} else if (value<0.0) {
	  if (colUpper[j] >= 1e15) {
	    dmaxdown = -1e31;
	    ++iflagl;
	  } else {
	    dmaxdown += colUpper[j] * value;
	  }
	  if (colLower[j] <= -1e15) {
	    dmaxup = 1e31;
	    ++iflagu;
	  } else {
	    dmaxup += colLower[j] * value;
	  }
	}
      }
      if (iflagu)
	maxR[i]=1.0e60;
      else
	maxR[i]=dmaxup;
      if (iflagl) 
	minR[i]=-1.0e60;
      else
	minR[i]=dmaxdown;
      if (minR[i]<-1.0e10&&maxR[i]>1.0e10) {
	markR[i]=-2;
      } else {
	markR[i]=-1;
      }
    } else {
      minR[i]=-1.0e60;
      maxR[i]=1.0e60;
      markR[i]=-2;
    }
  }
}
#ifdef CGL_DEBUG
static int nPath=0;
#endif
//-------------------------------------------------------------------
// Generate disaggregation cuts
//------------------------------------------------------------------- 
void CglProbing::generateCuts(const OsiSolverInterface & si, 
						OsiCuts & cs ) const
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n",nPath);
    nPath++;
    int nCols=si.getNumCols();
    int i;
    const double * solution = si.getColSolution();
    const double * lower = si.getColLower();
    const double * upper = si.getColUpper();
    const double * optimal = debugger->optimalSolution();
    const double * objective = si.getObjCoefficients();
    double objval1=0.0,objval2=0.0;
    for (i=0;i<nCols;i++) {
      printf("%d %g %g %g %g\n",i,lower[i],solution[i],upper[i],optimal[i]);
      objval1 += solution[i]*objective[i];
      objval2 += optimal[i]*objective[i];
      assert(optimal[i]>=lower[i]&&optimal[i]<=upper[i]);
    }
    printf("current obj %g, integer %g\n",objval1,objval2);
  }
#endif
  int nRows=si.getNumRows(); 
  double * rowLower = new double[nRows+1];
  double * rowUpper = new double[nRows+1];

  int nCols=si.getNumCols(); 
  double * colLower = new double[nCols];
  double * colUpper = new double[nCols];

  int ninfeas=gutsOfGenerateCuts(si,cs,rowLower,rowUpper,colLower,colUpper);
  if (ninfeas) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(DBL_MAX);
    rc.setUb(0.0);   
    cs.insert(rc);
#ifdef CGL_DEBUG
    const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
    if (debugger&&debugger->onOptimalPath(si))
      assert(!debugger->invalidCut(rc)); 
#endif
  }
  delete [] rowLower;
  delete [] rowUpper;
  delete [] colLower;
  delete [] colUpper;
}
int CglProbing::generateCutsAndModify(const OsiSolverInterface & si, 
						OsiCuts & cs )
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n",nPath);
    nPath++;
    int nCols=si.getNumCols();
    int i;
    const double * solution = si.getColSolution();
    const double * lower = si.getColLower();
    const double * upper = si.getColUpper();
    const double * optimal = debugger->optimalSolution();
    const double * objective = si.getObjCoefficients();
    double objval1=0.0,objval2=0.0;
    for (i=0;i<nCols;i++) {
      printf("%d %g %g %g %g\n",i,lower[i],solution[i],upper[i],optimal[i]);
      objval1 += solution[i]*objective[i];
      objval2 += optimal[i]*objective[i];
      assert(optimal[i]>=lower[i]&&optimal[i]<=upper[i]);
    }
    printf("current obj %g, integer %g\n",objval1,objval2);
  }
#endif
  int saveMode = mode_;
  if (!mode_)
    mode_=1;
  int nRows=si.getNumRows(); 
  double * rowLower = new double[nRows+1];
  double * rowUpper = new double[nRows+1];

  int nCols=si.getNumCols(); 
  double * colLower = new double[nCols];
  double * colUpper = new double[nCols];

  int ninfeas=gutsOfGenerateCuts(si,cs,rowLower,rowUpper,colLower,colUpper);
  if (ninfeas) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(DBL_MAX);
    rc.setUb(0.0);   
    cs.insert(rc);
#ifdef CGL_DEBUG
    const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
    if (debugger&&debugger->onOptimalPath(si))
      assert(!debugger->invalidCut(rc)); 
#endif
  }

  mode_=saveMode;
  // move bounds so can be used by user
  delete [] rowLower_;
  delete [] rowUpper_;
  delete [] colLower_;
  delete [] colUpper_;
  rowLower_ = rowLower;
  rowUpper_ = rowUpper;
  colLower_	= colLower;
  colUpper_	= colUpper;
  return ninfeas;
}
int CglProbing::gutsOfGenerateCuts(const OsiSolverInterface & si, 
			      OsiCuts & cs ,
			      double * rowLower, double * rowUpper,
			      double * colLower, double * colUpper) const
{
  // Get basic problem information
  int nRows;
  
  CoinPackedMatrix * rowCopy=NULL;

  // get branch and bound cutoff
  double cutoff;
  assert(si.getDblParam(OsiDualObjectiveLimit,cutoff));
  
  int mode=mode_;
  
  int nCols=si.getNumCols(); 

  // get integer variables
  char * intVar = new char[nCols];
  int i;
  int numberIntegers=0;
  for (i=0;i<nCols;i++) {
    if (si.isInteger(i)) {
      intVar[i]=1;
      numberIntegers++;
    } else {
      intVar[i]=0;
    }
  }

  int ninfeas=0;

  // see if using cached copy or not
  if (!numberRows_) {
    // create from current
    nRows=si.getNumRows(); 
    
    // mode==0 is invalid if going from current matrix
    if (mode==0)
      mode=1;
    rowCopy = new CoinPackedMatrix(*si.getMatrixByRow());
    // add in objective if there is a cutoff
    if (cutoff<1.0e30&&usingObjective_) {
      int * columns = new int[nCols];
      double * elements = new double[nCols];
      int n=0;
      const double * objective = si.getObjCoefficients();
      bool maximize = (si.getObjSense()==-1);
      for (i=0;i<nCols;i++) {
	if (objective[i]) {
	  elements[n]= (maximize) ? -objective[i] : objective[i];
	  columns[n++]=i;
	}
      }
      rowCopy->appendRow(n,columns,elements);
      delete [] columns;
      delete [] elements;
      memcpy(rowLower,si.getRowLower(),nRows*sizeof(double));
      memcpy(rowUpper,si.getRowUpper(),nRows*sizeof(double));
      rowLower[nRows]=-DBL_MAX;
      rowUpper[nRows]=cutoff;
      nRows++;
    } else {
      memcpy(rowLower,si.getRowLower(),nRows*sizeof(double));
      memcpy(rowUpper,si.getRowUpper(),nRows*sizeof(double));
    }
    memcpy(colLower,si.getColLower(),nCols*sizeof(double));
    memcpy(colUpper,si.getColUpper(),nCols*sizeof(double));
  } else {
    // use snapshot
    nRows=numberRows_;
    assert(nCols==numberColumns_);
    
    rowCopy = rowCopy_;
    rowLower = new double[nRows];
    rowUpper = new double[nRows];
    memcpy(rowLower,rowLower_,nRows*sizeof(double));
    memcpy(rowUpper,rowUpper_,nRows*sizeof(double));
    rowLower[nRows-1]=-DBL_MAX;
    if (usingObjective_) {
      rowUpper[nRows-1]=cutoff;
    } else {
      rowUpper[nRows-1]=DBL_MAX;
    }
    memcpy(colLower,colLower_,nCols*sizeof(double));
    memcpy(colUpper,colUpper_,nCols*sizeof(double));
    if (mode) {
      // tighten bounds to reflect state of problem
      const double * lower = si.getColLower();
      const double * upper = si.getColUpper();
      for (i=0;i<nCols;i++) {
	if (colLower[i]<lower[i])
	  colLower[i]=lower[i];
	if (colUpper[i]>upper[i])
	  colUpper[i]=upper[i];
      }
    }
  }
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&!debugger->onOptimalPath(si))
    debugger = NULL;
#else
  const OsiRowCutDebugger * debugger = NULL;
#endif
   
  const int * column = rowCopy->getIndices();
  const int * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  
  if (mode) {
    ninfeas= tighten(colLower, colUpper, column, rowElements,
			 rowStart, rowLength, rowLower, rowUpper,
			 nRows, nCols, intVar, 2, primalTolerance_);
    if (!ninfeas) {
      // create column cuts where integer bounds have changed
      {
	const double * lower = si.getColLower();
	const double * upper = si.getColUpper();
	const double * colsol = si.getColSolution();
	int numberChanged=0,ifCut=0;
	CoinPackedVector lbs;
	CoinPackedVector ubs;
	for (i = 0; i < nCols; ++i) {
	  if (intVar[i]) {
	    colUpper[i] = min(upper[i],floor(colUpper[i]+1.0e-4));
	    if (colUpper[i]<upper[i]-1.0e-8) {
	      if (colUpper[i]<colsol[i]-1.0e-8)
		ifCut=1;
	      ubs.insert(i,colUpper[i]);
	      numberChanged++;
	    }
	    colLower[i] = max(lower[i],ceil(colLower[i]-1.0e-4));
	    if (colLower[i]>lower[i]+1.0e-8) {
	      if (colLower[i]>colsol[i]+1.0e-8)
		ifCut=1;
	      lbs.insert(i,colLower[i]);
	      numberChanged++;
	    }
	  }
	}
	if (numberChanged) {
	  OsiColCut cc;
	  cc.setUbs(ubs);
	  cc.setLbs(lbs);
	  if (ifCut) {
	    cc.setEffectiveness(100.0);
	  } else {
	    cc.setEffectiveness(1.0e-5);
	  }
#ifdef CGL_DEBUG
	  checkBounds(debugger,cc);
#endif
	  cs.insert(cc);
	}
      }
      int * markR = new int [nRows];
      double * minR = new double [nRows];
      double * maxR = new double [nRows];
      int * look = new int[nCols];
      int nLook=0;
      // get min max etc for rows
      tighten2(colLower, colUpper, column, rowElements,
	       rowStart, rowLength, rowLower, rowUpper,
	       minR , maxR , markR, nRows, nCols);
      // decide what to look at
      if (mode==1) {
	const double * colsol = si.getColSolution();
	double_int_pair * array = new double_int_pair [nCols];
	for (i=0;i<nCols;i++) {
	  if (intVar[i]&&colUpper[i]-colLower[i]>1.0e-8) {
	    double away = fabs(0.5-(colsol[i]-floor(colsol[i])));
	    if (away<0.49999) {
	      array[nLook].infeasibility=away;
	      array[nLook++].sequence=i;
	    }
	  }
	}
	std::sort(array,array+nLook,double_int_pair_compare());
	nLook=min(nLook,maxProbe_);
	for (i=0;i<nLook;i++) {
	  look[i]=array[i].sequence;
	}
	delete [] array;
      } else {
	for (i=0;i<nCols;i++) {
	  if (intVar[i]&&colUpper[i]-colLower[i]>1.0e-8) {
	    look[nLook++]=i;
	  }
	}
      }
      ninfeas= probe(si, debugger, cs, colLower, colUpper, rowCopy,
		     rowLower, rowUpper,
		     intVar, minR, maxR, markR,
		     look,nLook);
      delete [] look;
      delete [] markR;
      delete [] minR;
      delete [] maxR;
    }
  } else {
    // global cuts from previous calculations
    // could check more thoroughly that integers are correct
    assert(numberIntegers==numberIntegers_);
    // make up list of new variables to look at 
    int * markR = new int [nRows];
    double * minR = new double [nRows];
    double * maxR = new double [nRows];
    int * look = new int[nCols];
    int nLook=0;
    const double * colsol = si.getColSolution();
    for (i=0;i<numberIntegers_;i++) {
      int j=cutVector_[i].sequence;
      if (!cutVector_[i].newValue&&colUpper[j]-colLower[j]>1.0e-8) {
	double away = fabs(0.5-(colsol[j]-floor(colsol[j])));
	if (away<0.4999999) {
	  look[nLook++]=j;
	}
      }
    }
    // get min max etc for rows
    tighten2(colLower, colUpper, column, rowElements,
	     rowStart, rowLength, rowLower, rowUpper,
	     minR , maxR , markR, nRows, nCols);
    OsiCuts csNew;
    ninfeas= probe(si, debugger, csNew, colLower, colUpper, rowCopy,
		   rowLower, rowUpper,
		       intVar, minR, maxR, markR,
		   look,nLook);
    if (!ninfeas) {
      // go through row cuts
      int nCuts = csNew.sizeRowCuts();
      int iCut;
      // need space for backward lookup
      // just for ones being looked at
      int * backward = new int [2*nCols];
      int * onList = backward + nCols;
      for (i=0;i<numberIntegers_;i++) {
	int j=cutVector_[i].sequence;
	backward[j]=i;
	onList[j]=0;
      }
      for (i=0;i<nLook;i++) {
	onList[look[i]]=1;
      }
      // first do counts
      // we know initialized to zero
      for (iCut=0;iCut<nCuts;iCut++) {
	OsiRowCut rcut;
	CoinPackedVector rpv;
	rcut = csNew.rowCut(iCut);
	rpv = rcut.row();
	assert(rpv.getNumElements()==2);
	const int * indices = rpv.getIndices();
	double* elements = rpv.getElements();
	double lb=rcut.lb();
	// find out which integer
        i=backward[indices[1]];
	if (lb==-DBL_MAX) {
	  // UB
	  if (elements[1]<0.0) {
	    cutVector_[i].lastUBWhenAt0++;
	  } else { 
	    cutVector_[i].lastUBWhenAt1++;
	  }
	} else {
	  // LB
	  if (elements[1]>0.0) {
	    cutVector_[i].lastLBWhenAt0++; 
	  } else {
	    cutVector_[i].lastLBWhenAt1++;
	  }
	} 
      }
      // allocate space
      for (i=0;i<numberIntegers_;i++) {
	int j=cutVector_[i].sequence;
	if (onList[j]&&!cutVector_[i].newValue) {
	  disaggregation thisOne=cutVector_[i];
	  int total;
	  // set to start of each type
	  cutVector_[i].lastUBWhenAt0 =0;
	  total = thisOne.lastUBWhenAt0;
	  cutVector_[i].lastUBWhenAt1 = total; 
	  total += thisOne.lastUBWhenAt1;
	  cutVector_[i].lastLBWhenAt0 = total;
	  total += thisOne.lastLBWhenAt0;
	  cutVector_[i].lastLBWhenAt1 = total;
	  total += thisOne.lastLBWhenAt1;
	  cutVector_[i].newValue=new double[total];
	  cutVector_[i].index=new int[total];
	}
      }
      // now put in
      for (iCut=0;iCut<nCuts;iCut++) {
	OsiRowCut rcut;
	CoinPackedVector rpv;
	int iput;
	rcut = csNew.rowCut(iCut);
	rpv = rcut.row();
	assert(rpv.getNumElements()==2);
	const int * indices = rpv.getIndices();
	double* elements = rpv.getElements();
	double lb=rcut.lb();
	// find out which integer
        i=backward[indices[1]];
	if (lb==-DBL_MAX) {
	  // UB
	  if (elements[1]<0.0) {
	    iput=cutVector_[i].lastUBWhenAt0;
	    cutVector_[i].newValue[iput]=elements[1];
	    cutVector_[i].index[iput]=indices[0];
	    cutVector_[i].lastUBWhenAt0++;
	  } else { 
	    iput=cutVector_[i].lastUBWhenAt1;
	    cutVector_[i].newValue[iput]=elements[1];
	    cutVector_[i].index[iput]=indices[0];
	    cutVector_[i].lastUBWhenAt1++;
	  }
	} else {
	  // LB
	  if (elements[1]>0.0) {
	    iput=cutVector_[i].lastLBWhenAt0; 
	    cutVector_[i].newValue[iput]=elements[1];
	    cutVector_[i].index[iput]=indices[0];
	    cutVector_[i].lastLBWhenAt0++; 
	  } else {
	    iput=cutVector_[i].lastLBWhenAt1;
	    cutVector_[i].newValue[iput]=elements[1];
	    cutVector_[i].index[iput]=indices[0];
	    cutVector_[i].lastLBWhenAt1++;
	  }
	} 
      }
      // now see if any disaggregation cuts are violated
      for (i=0;i<numberIntegers_;i++) {
	int j=cutVector_[i].sequence;
	double upperInt=colUpper_[j];
	double lowerInt=colLower_[j];
	double solInt=colsol[j];
	double down = solInt-lowerInt;
	double up = upperInt-solInt;
	double lower, upper, solValue;
	int icol;
	double newSol,value,newValue;
	int index[2];
	double element[2];
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  double away = fabs(0.5-(colsol[j]-floor(colsol[j])));
	  if (away<0.4999999) {
	    disaggregation thisOne=cutVector_[i];
	    int k;
	    OsiRowCut rc;
	    // UB changes when integer goes to lb
	    for (k=0;k<thisOne.lastUBWhenAt0;k++) {
	      icol = thisOne.index[k];
	      value = thisOne.newValue[k]; // new U - old U
	      solValue=colsol[icol];
	      upper=colUpper_[icol];
	      newValue= value + upper;
	      if (solValue > newValue - value * down + primalTolerance_) {
		rc.setLb(-DBL_MAX);
		rc.setUb(newValue + lowerInt*value);
		index[0]=icol;
		element[0]=1.0;
		index[1]=j;
		element[1]= value;
		// effectiveness is how far j moves
		newSol = (rc.ub()-solValue)/element[1];
		if (mode_)
		  assert(newSol>solInt);
		rc.setEffectiveness(newSol-solInt);
		if (rc.effectiveness()>1.0e-3) {
		  rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		  if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		  cs.insert(rc);
		}
	      }
	    }
	    // UB changes when integer goes to ub
	    for (k=thisOne.lastUBWhenAt0;k<thisOne.lastUBWhenAt1;k++) {
	      icol = thisOne.index[k];
	      value = thisOne.newValue[k]; // new U - old U
	      solValue=colsol[icol];
	      upper=colUpper_[icol];
	      newValue= value + upper;
	      if (solValue > newValue - value * up + primalTolerance_) {
		rc.setLb(-DBL_MAX);
		rc.setUb(newValue - upperInt*value);
		index[0]=icol;
		element[0]=1.0;
		index[1]=j;
		element[1]= value;
		// effectiveness is how far j moves
		newSol = (rc.ub()-solValue)/element[1];
		if (mode_)
		  assert(newSol>solInt);
		rc.setEffectiveness(newSol-solInt);
		if (rc.effectiveness()>1.0e-3) {
		  rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		  if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		  cs.insert(rc);
		}
	      }
	    }
	    // LB changes when integer goes to lb
	    for (k=thisOne.lastUBWhenAt1;k<thisOne.lastLBWhenAt0;k++) {
	      icol = thisOne.index[k];
	      value = thisOne.newValue[k]; // new L - old L
	      solValue=colsol[icol];
	      lower=colLower_[icol];
	      newValue= value + lower;
	      if (solValue < newValue - value * down - primalTolerance_) {
		rc.setLb(newValue + lowerInt*value);
		rc.setUb(DBL_MAX);
		index[0]=icol;
		element[0]=1.0;
		index[1]=j;
		element[1]= value;
		// effectiveness is how far j moves
		newSol = (rc.lb()-solValue)/element[1];
		if (mode_)
		  assert(newSol>solInt);
		rc.setEffectiveness(newSol-solInt);
		if (rc.effectiveness()>1.0e-3) {
		  rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		  if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		  cs.insert(rc);
		}
	      }
	    }
	    // LB changes when integer goes to ub
	    for (k=thisOne.lastLBWhenAt0;k<thisOne.lastLBWhenAt1;k++) {
	      icol = thisOne.index[k];
	      value = thisOne.newValue[k]; // new L - old L
	      solValue=colsol[icol];
	      lower=colLower_[icol];
	      newValue= value + lower;
	      if (solValue < newValue - value * up - primalTolerance_) {
		rc.setLb(newValue - upperInt*value);
		rc.setUb(DBL_MAX);
		index[0]=icol;
		element[0]=1.0;
		index[1]=j;
		element[1]= value;
		// effectiveness is how far j moves
		newSol = (rc.lb()-solValue)/element[1];
		if (mode_)
		  assert(newSol>solInt);
		rc.setEffectiveness(newSol-solInt);
		if (rc.effectiveness()>1.0e-3) {
		  rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		  if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		  cs.insert(rc);
		}
	      }
	    }
	  }
	}
      }
      delete [] backward;
    }
    delete [] look;
    delete [] markR;
    delete [] minR;
    delete [] maxR;
  }
  // delete stuff
  if (!numberRows_) {
    delete rowCopy;
  } else {
    delete [] rowLower;
    delete [] rowUpper;
  }
  delete [] intVar;
  return ninfeas;
}
// Does probing and adding cuts
int CglProbing::probe( const OsiSolverInterface & si, 
		       const OsiRowCutDebugger * debugger, 
		       OsiCuts & cs, 
		       double * colLower, double * colUpper, 
		       CoinPackedMatrix *rowCopy,
		       double * rowLower, double * rowUpper,
		       char * intVar, double * minR, double * maxR, 
		       int * markR, 
		       int * look, int nLook) const
{
  int nRows=rowCopy->getNumRows();
  int nCols=rowCopy->getNumCols();
  double * colsol = new double[nCols];
  double * djs = new double[nCols];
  const double * currentColLower = si.getColLower();
  const double * currentColUpper = si.getColUpper();
  double * tempL = new double [nCols];
  double * tempU = new double [nCols];
  int * markC = new int [nCols];
  int * stackC = new int [2*nCols];
  int * stackR = new int [nRows];
  double * saveL = new double [2*nCols];
  double * saveU = new double [2*nCols];
  double * saveMin = new double [nRows];
  double * saveMax = new double [nRows];
  double * element = new double[nCols];
  int * index = new int[nCols];
  CoinPackedMatrix columnCopy = * rowCopy;
  columnCopy.reverseOrdering();
  const int * column = rowCopy->getIndices();
  const int * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  const int * row = columnCopy.getIndices();
  const int * columnStart = columnCopy.getVectorStarts();
  const int * columnLength = columnCopy.getVectorLengths(); 
  const double * columnElements = columnCopy.getElements();
  double movement;
  int i, j, k,kk,jj;
  int jcol,kcol,irow,krow;
  bool anyColumnCuts=false;
  double dbound, value, value2;
  int ninfeas=0;
  int rowCuts;
  double disaggEffectiveness;
  if (mode_) {
    /* clean up djs and solution */
    memcpy(djs,si.getReducedCost(),nCols*sizeof(double));
    memcpy(colsol, si.getColSolution(),nCols*sizeof(double));
    disaggEffectiveness=1.0e-3;
    rowCuts=rowCuts_;
  } else {
    // need to go from a neutral place
    memset(djs,0,nCols*sizeof(double));
    memcpy(colsol, si.getColSolution(),nCols*sizeof(double));
    disaggEffectiveness=-1.0e10;
    rowCuts=1;
  }
  for (i = 0; i < nCols; ++i) {
    /* was if (intVar[i]) */
    if (1) {
      if (colUpper[i]-colLower[i]>1.0e-8) {
	if (colsol[i]<colLower[i]+primalTolerance_) {
	  colsol[i]=colLower[i];
	  djs[i] = max(0.0,djs[i]);
	} else if (colsol[i]>colUpper[i]-primalTolerance_) {
	  colsol[i]=colUpper[i];
	  djs[i] = min(0.0,djs[i]);
	} else {
	  djs[i]=0.0;
	}
	/*if (fabs(djs[i])<1.0e-5) 
	  djs[i]=0.0;*/
      }
    }
  }

  int ipass=0,nfixed=-1;

  double cutoff;
  assert(si.getDblParam(OsiDualObjectiveLimit,cutoff));
  double current = si.getObjValue();
  // make irrelevant if mode is 0
  if (!mode_)
    cutoff=DBL_MAX;
  /* for both way coding */
  int nstackC0=-1;
  int * stackC0 = new int[maxStack_];
  double * lo0 = new double[maxStack_];
  double * up0 = new double[maxStack_];
  int nstackR,nstackC;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3;
    } else {
      markC[i]=0;
    }
  }
  double tolerance = 1.0e1*primalTolerance_;
  while (ipass<maxPass_&&nfixed) {
    int iLook;
    ipass++;
    nfixed=0;
    for (iLook=0;iLook<nLook;iLook++) {
      double solval;
      double down;
      double up;
      int awayFromBound=1;
      j=look[iLook];
      solval=colsol[j];
      down = floor(solval+tolerance);
      up = ceil(solval-tolerance);
      if(colUpper[j]-colLower[j]<1.0e-8) markC[j]=3;
      if (markC[j]||!intVar[j]) continue;
      if (solval>colUpper[j]-tolerance||solval<colLower[j]+tolerance) {
	awayFromBound=0;
	if (solval<colLower[j]+tolerance) {
	  solval = colLower[j]+1.0e-1;
	  down=colLower[j];
	  up=down+1.0;
	} else if (solval>colUpper[j]-tolerance) {
	  solval = colUpper[j]-1.0e-1;
	  up=colUpper[j];
	  down=up-1;
	} 
      }
      if ((solval-down>1.0e-6&&up-solval>1.0e-6)||mode_!=1) {
	int istackC,iway, istackR;
	int way[]={1,2,1};
	int feas[]={1,2,4};
	int feasible=0;
	int notFeasible;
	for (iway=0;iway<3;iway ++) {
	  int fixThis=0;
	  double objVal=current;
	  int goingToTrueBound=0;
	  stackC[0]=j;
	  markC[j]=way[iway];
	  if (way[iway]==1) {
	    movement=down-colUpper[j];
	    assert(movement<-0.99999);
	    if (fabs(down-colLower[j])<1.0e-7) {
	      goingToTrueBound=2;
	      down=colLower[j];
	    }
	  } else {
	    movement=up-colLower[j];
	    assert(movement>0.99999);
	    if (fabs(up-colUpper[j])<1.0e-7) {
	      goingToTrueBound=2;
	      up=colUpper[j];
	    }
	  }
	  if (goingToTrueBound&&(colUpper[j]-colLower[j]>1.5||colLower[j]))
	    goingToTrueBound=1;
	  // switch off disaggregation if not wanted
	  if ((rowCuts&1)==0)
	    goingToTrueBound=0;
#ifdef PRINT_DEBUG
	  if (fabs(movement)>1.01) {
	    printf("big %d %g %g %g\n",j,colLower[j],solval,colUpper[j]);
	  }
#endif
	  if (movement*djs[j]>0.0)
	    objVal += movement*djs[j];
	  nstackC=1;
	  nstackR=0;
	  saveL[0]=colLower[j];
	  saveU[0]=colUpper[j];
	  notFeasible=0;
	  if (movement<0.0) {
	    colUpper[j] += movement;
	    colUpper[j] = floor(colUpper[j]+0.5);
#ifdef PRINT_DEBUG
	    printf("** Trying %d down to 0\n",j);
#endif
	  } else {
	    colLower[j] += movement;
	    colLower[j] = floor(colLower[j]+0.5);
#ifdef PRINT_DEBUG
	    printf("** Trying %d up to 1\n",j);
#endif
	  }
	  if (fabs(colUpper[j]-colLower[j])<1.0e-6)
	    markC[j]=3; // say fixed
	  istackC=0;
	  /* update immediately */
	  for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
	    irow = row[k];
	    value = columnElements[k];
	    if (markR[irow]==-1) {
	      stackR[nstackR]=irow;
	      markR[irow]=nstackR;
	      saveMin[nstackR]=minR[irow];
	      saveMax[nstackR]=maxR[irow];
	      nstackR++;
	    } else if (markR[irow]==-2) {
	      continue;
	    }
	    /* could check immediately if violation */
	    if (movement>0.0) {
	      /* up */
	      if (value>0.0) {
		/* up does not change - down does */
		minR[irow] += value;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      } else {
		/* down does not change - up does */
		maxR[irow] += value;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      }
	    } else {
	      /* down */
	      if (value<0.0) {
		/* up does not change - down does */
		minR[irow] -= value;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      } else {
		/* down does not change - up does */
		maxR[irow] -= value;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      }
	    }
	  }
	  while (istackC<nstackC&&nstackC<maxStack_) {
	    int jway;
	    jcol=stackC[istackC];
	    jway=markC[jcol];
	    // If not first and fixed then skip
	    if (jway==3&&istackC) {
	      istackC++;
	      continue;
	    }
	    for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
	      // break if found not feasible
	      if (notFeasible)
		break;
	      irow = row[k];
	      /*value = columnElements[k];*/
	      if (markR[irow]!=-2) {
		/* see if anything forced */
		for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
		  double moveUp=0.0;
		  double moveDown=0.0;
		  double newUpper=-1.0,newLower=1.0;
		  kcol=column[kk];
		  bool onList = (markC[kcol]!=0);
		  if (markC[kcol]!=3) {
		    value2=rowElements[kk];
		    if (value2 < 0.0) {
		      if (colUpper[kcol] < 1e10 && (markC[kcol]&2)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colUpper[kcol]+
			  (rowUpper[irow]-minR[irow])/value2;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
			    markC[kcol] |= 2;
			    newLower = ceil(dbound-primalTolerance_);
			  } else {
			    markC[kcol]=3;
			    newLower=dbound;
			    if (newLower>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol];
			    }
			  }
			  moveUp = newLower-colLower[kcol];
			}
		      }
		      if (colLower[kcol] > -1e10 && (markC[kcol]&1)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colLower[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markC[kcol] |= 1;
			    newUpper = floor(dbound+primalTolerance_);
			  } else {
			    markC[kcol]=3;
			    newUpper=dbound;
			    if (newUpper<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol];
			    }
			  }
			  moveDown = newUpper-colUpper[kcol];
			}
		      }
		    } else {
		      /* positive element */
		      if (colUpper[kcol] < 1e10 && (markC[kcol]&2)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colUpper[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
			    markC[kcol] |= 2;
			    newLower = ceil(dbound-primalTolerance_);
			  } else {
			    markC[kcol]=3;
			    newLower=dbound;
			    if (newLower>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol];
			    }
			  }
			  moveUp = newLower-colLower[kcol];
			}
		      }
		      if (colLower[kcol] > -1e10 && (markC[kcol]&1)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colLower[kcol] + 
			  (rowUpper[irow]-minR[irow])/value2;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markC[kcol] |= 1;
			    newUpper = floor(dbound+primalTolerance_);
			  } else {
			    markC[kcol]=3;
			    newUpper=dbound;
			    if (newUpper<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol];
			    }
			  }
			  moveDown = newUpper-colUpper[kcol];
			}
		      }
		    }
		    if (moveUp&&nstackC<2*maxStack_) {
		      fixThis++;
#ifdef PRINT_DEBUG
		      printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]);
		      value=0.0;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]);
			} else {
			  value += rowElements[jj]*colLower[ii];
			}
		      }
		      printf(" - fixed %g\n",value);
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
			}
		      }
		      printf("\n");
#endif
		      if (!onList) {
			stackC[nstackC]=kcol;
			saveL[nstackC]=colLower[kcol];
			saveU[nstackC]=colUpper[kcol];
			nstackC++;
			onList=true;
		      }
		      if (newLower>colsol[kcol]) {
			if (djs[kcol]<0.0) {
			  /* should be infeasible */
			  assert (newLower>colUpper[kcol]+primalTolerance_);
			} else {
			  objVal += moveUp*djs[kcol];
			}
		      }
		      if (intVar[kcol])
			newLower = max(colLower[kcol],ceil(newLower-1.0e-4));
		      colLower[kcol]=newLower;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
			markC[kcol]=3; // say fixed
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj];
			value = columnElements[jj];
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow;
			  markR[krow]=nstackR;
			  saveMin[nstackR]=minR[krow];
			  saveMax[nstackR]=maxR[krow];
			  nstackR++;
			} else if (markR[krow]==-2) {
			  continue;
			}
			/* could check immediately if violation */
			/* up */
			if (value>0.0) {
			  /* up does not change - down does */
			  minR[krow] += value*moveUp;
			} else {
			  /* down does not change - up does */
			  maxR[krow] += value*moveUp;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			}
		      }
		    }
		    if (moveDown&&nstackC<2*maxStack_) {
		      fixThis++;
#ifdef PRINT_DEBUG
		      printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]);
		      value=0.0;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]);
			} else {
			  value += rowElements[jj]*colLower[ii];
			}
		      }
		      printf(" - fixed %g\n",value);
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
			}
		      }
		      printf("\n");
#endif
		      if (!onList) {
			stackC[nstackC]=kcol;
			saveL[nstackC]=colLower[kcol];
			saveU[nstackC]=colUpper[kcol];
			nstackC++;
			onList=true;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_);
			} else {
			  objVal += moveDown*djs[kcol];
			}
		      }
		      if (intVar[kcol])
			newUpper = min(colUpper[kcol],floor(newUpper+1.0e-4));
		      colUpper[kcol]=newUpper;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
			markC[kcol]=3; // say fixed
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj];
			value = columnElements[jj];
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow;
			  markR[krow]=nstackR;
			  saveMin[nstackR]=minR[krow];
			  saveMax[nstackR]=maxR[krow];
			  nstackR++;
			} else if (markR[krow]==-2) {
			  continue;
			}
			/* could check immediately if violation */
			/* down */
			if (value<0.0) {
			  /* up does not change - down does */
			  minR[krow] += value*moveDown;
			  if (minR[krow]>rowUpper[krow]+1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			} else {
			  /* down does not change - up does */
			  maxR[krow] += value*moveDown;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			}
		      }
		    }
		    if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
		      notFeasible=1;;
		      k=columnStart[jcol]+columnLength[jcol];
		      istackC=nstackC+1;
#ifdef PRINT_DEBUG
		      printf("** not feasible this way\n");
#endif
		      break;
		    }
		  }
		}
	      }
	    }
	    istackC++;
	  }
	  if (!notFeasible) {
	    if (objVal<=cutoff) {
	      feasible |= feas[iway];
	    } else {
#ifdef PRINT_DEBUG
	      printf("not feasible on dj\n");
#endif
	      notFeasible=1;
	      if (iway==1&&feasible==0) {
		/* not feasible at all */
		ninfeas=1;
		j=nCols-1;
		break;
	      }
	    }
	  } else if (iway==1&&feasible==0) {
	    /* not feasible at all */
	    ninfeas=1;
	    j=nCols-1;
            iLook=nLook;
	    ipass=maxPass_;
	    break;
	  }
	  if (notFeasible)
	    goingToTrueBound=0;
	  if (iway==2||(iway==1&&feasible==2)) {
	    /* keep */
	    iway=3;
	    nfixed++;
            if (mode_) {
	      OsiColCut cc;
	      int nTot=0,nFix=0,nInt=0;
	      bool ifCut=false;
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC];
		if (intVar[icol]) {
		  if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
		    element[nFix]=colUpper[icol];
		    index[nFix++]=icol;
		    nInt++;
		    if (colsol[icol]>colUpper[icol]+primalTolerance_) {
		      ifCut=true;
		      anyColumnCuts=true;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot=nFix;
		cc.setUbs(nFix,index,element);
		nFix=0;
	      }
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC];
		if (intVar[icol]) {
		  if (colLower[icol]>currentColLower[icol]+1.0e-4) {
		    element[nFix]=colLower[icol];
		    index[nFix++]=icol;
		    nInt++;
		    if (colsol[icol]<colLower[icol]-primalTolerance_) {
		      ifCut=true;
		      anyColumnCuts=true;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot+=nFix;
		cc.setLbs(nFix,index,element);
	      }
	      // could tighten continuous as well
	      if (nInt) {
		if (ifCut) {
		  cc.setEffectiveness(100.0);
		} else {
		  cc.setEffectiveness(1.0e-5);
		}
#ifdef CGL_DEBUG
		checkBounds(debugger,cc);
#endif
		cs.insert(cc);
	      }
	    }
	    for (istackC=0;istackC<nstackC;istackC++) {
	      int icol=stackC[istackC];
	      if (colUpper[icol]-colLower[icol]>primalTolerance_) {
		markC[icol]=0;
	      } else {
		markC[icol]=3;
	      }
	    }
	    for (istackR=0;istackR<nstackR;istackR++) {
	      int irow=stackR[istackR];
	      markR[irow]=-1;
	    }
	  } else {
	    /* is it worth seeing if can increase coefficients
	       or maybe better see if it is a cut */
	    if (iway==0) {
	      nstackC0=min(nstackC,maxStack_);
	      double solMove = solval-down;
	      double boundChange;
	      if (notFeasible) {
		nstackC0=0;
	      } else {
		for (istackC=0;istackC<nstackC0;istackC++) {
		  int icol=stackC[istackC];
		  stackC0[istackC]=icol;
		  lo0[istackC]=colLower[icol];
		  up0[istackC]=colUpper[icol];
		}
	      }
	      /* restore all */
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC];
		double oldU=saveU[istackC];
		double oldL=saveL[istackC];
		if(goingToTrueBound==2&&istackC) {
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (jval-down)
		  boundChange = oldU-colUpper[icol];
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(-DBL_MAX);
		    rc.setUb(colUpper[icol]-down*boundChange);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= - boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange;
		    if (mode_) 
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      cs.insert(rc);
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (jval-down)
		  boundChange = oldL-colLower[icol];
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(colLower[icol]-down*boundChange);
		    rc.setUb(DBL_MAX);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]=- boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      cs.insert(rc);
#if 0
		      printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,oldL,oldU,colLower[icol],colsol[icol], j, colsol[j]);
		      printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
			     icol,boundChange,j,rc.ub());
#endif
		    }
		  }
		}
		colUpper[icol]=oldU;
		colLower[icol]=oldL;
		markC[icol]=0;
	      }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR];
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts;
		  double gap = rowUpper[irow]-maxR[irow];
		  double sum=0.0;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]];
		    }
		    if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_) {
		      // can be a cut
		      // subtract gap from upper and integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
		      if (rc.effectiveness()>1.0e-3) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
			cs.insert(rc);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum+gap*colsol[j]<minR[irow]+primalTolerance_) {
		      // can be a cut
		      // add gap to lower and integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap);
		      if (rc.effectiveness()>1.0e-3) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
			cs.insert(rc);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    } else {
	      if (iway==1&&feasible==3) {
		iway=3;
		/* point back to stack */
		for (istackC=nstackC-1;istackC>=0;istackC--) {
		  int icol=stackC[istackC];
		  markC[icol]=istackC+1000;
		}
		if (mode_) {
		  OsiColCut cc;
		  int nTot=0,nFix=0,nInt=0;
		  bool ifCut=false;
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC];
		    int istackC1=markC[icol]-1000;
		    if (istackC1>=0) {
		      if (min(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
			saveL[istackC1]=min(lo0[istackC],colLower[icol]);
			if (intVar[icol]) {
			  element[nFix]=saveL[istackC1];
			  index[nFix++]=icol;
			  nInt++;
			  if (colsol[icol]<saveL[istackC1]-primalTolerance_)
			    ifCut=true;
			}
			nfixed++;
		      }
		    }
		  }
		  if (nFix) {
		    nTot=nFix;
		    cc.setLbs(nFix,index,element);
		    nFix=0;
		  }
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC];
		    int istackC1=markC[icol]-1000;
		    if (istackC1>=0) {
		      if (max(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
			saveU[istackC1]=max(up0[istackC],colUpper[icol]);
			if (intVar[icol]) {
			  element[nFix]=saveU[istackC1];
			  index[nFix++]=icol;
			  nInt++;
			  if (colsol[icol]>saveU[istackC1]+primalTolerance_)
			    ifCut=true;
			}
			nfixed++;
		      }
		    }
		  }
		  if (nFix) {
		    nTot+=nFix;
		    cc.setUbs(nFix,index,element);
		  }
		  // could tighten continuous as well
		  if (nInt) {
		    if (ifCut) {
		      cc.setEffectiveness(100.0);
		    } else {
		      cc.setEffectiveness(1.0e-5);
		    }
#ifdef CGL_DEBUG
		    checkBounds(debugger,cc);
#endif
		    cs.insert(cc);
		  }
		}
	      } else {
		goingToTrueBound=0;
	      }
	      double solMove = up-solval;
	      double boundChange;
	      /* restore all */
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC];
		double oldU=saveU[istackC];
		double oldL=saveL[istackC];
		if(goingToTrueBound==2&&istackC) {
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (up-jval)
		  boundChange = oldU-colUpper[icol];
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(-DBL_MAX);
		    rc.setUb(colUpper[icol]+up*boundChange);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= + boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      cs.insert(rc);
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (up-jval)
		  boundChange = oldL-colLower[icol];
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(colLower[icol]+up*boundChange);
		    rc.setUb(DBL_MAX);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= + boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      cs.insert(rc);
		    }
		  }
		}
		colUpper[icol]=oldU;
		colLower[icol]=oldL;
		markC[icol]=0;
	      }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR];
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts;
		  double gap = rowUpper[irow]-maxR[irow];
		  double sum=0.0;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]];
		    }
		    if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_) {
		      // can be a cut
		      // add gap to integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap);
		      if (rc.effectiveness()>1.0e-3) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
			cs.insert(rc);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum-gap*colsol[j]<rowLower[irow]+primalTolerance_) {
		      // can be a cut
		      // subtract gap from integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]-gap*(colUpper[j]-1));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap);
		      if (rc.effectiveness()>1.0e-3) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
			cs.insert(rc);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    }
	  }
	}
      }
    }
  }
  delete [] stackC0;
  delete [] lo0;
  delete [] up0;
  delete [] tempL;
  delete [] tempU;
  delete [] markC;
  delete [] stackC;
  delete [] stackR;
  delete [] saveL;
  delete [] saveU;
  delete [] saveMin;
  delete [] saveMax;
  delete [] index;
  delete [] element;
  delete [] djs;
  delete [] colsol;
  return (ninfeas);
}
// Create a copy of matrix which is to be used
// this is to speed up process and to give global cuts
// Can give an array with 1 set to select, 0 to ignore
// column bounds are tightened
// If array given then values of 1 will be set to 0 if redundant
void CglProbing::snapshot ( const OsiSolverInterface & si,
		  char * possible)
{
  deleteSnapshot();
  // Get basic problem information
  
  numberColumns_=si.getNumCols(); 
  numberRows_=si.getNumRows();
  colLower_ = new double[numberColumns_];
  colUpper_ = new double[numberColumns_];
  memcpy(colLower_,si.getColLower(),numberColumns_*sizeof(double));
  memcpy(colUpper_,si.getColUpper(),numberColumns_*sizeof(double));
  rowLower_= new double [numberRows_+1];
  rowUpper_= new double [numberRows_+1];
  memcpy(rowLower_,si.getRowLower(),numberRows_*sizeof(double));
  memcpy(rowUpper_,si.getRowUpper(),numberRows_*sizeof(double));

  int i;
  if (possible) {
    for (i=0;i<numberRows_;i++) {
      if (!possible[i]) {
	rowLower_[i]=-DBL_MAX;
	rowUpper_[i]=DBL_MAX;
      }
    }
  }
  

  char * intVar = new char[numberColumns_];
  numberIntegers_=0;
  for (i=0;i<numberColumns_;i++) {
    if (si.isInteger(i)) {
      intVar[i]=1;  
      numberIntegers_++;
    } else {
      intVar[i]=0;
    }
  }
    
  rowCopy_ = new CoinPackedMatrix(*si.getMatrixByRow());

  const int * column = rowCopy_->getIndices();
  const int * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths(); 
  const double * rowElements = rowCopy_->getElements();
  
  int ninfeas= tighten(colLower_, colUpper_, column, rowElements,
			 rowStart, rowLength, rowLower_, rowUpper_,
			 numberRows_, numberColumns_, intVar, 5, primalTolerance_);
  assert (!ninfeas);

  // do integer stuff for mode 0
  cutVector_ = new disaggregation [numberIntegers_];
  memset(cutVector_,0,numberIntegers_*sizeof(disaggregation));
  numberIntegers_=0;
  for (i=0;i<numberColumns_;i++) {
    if (intVar[i]) {
      cutVector_[numberIntegers_++].sequence=i;
    }
  }
  delete [] intVar;

  // now delete rows
  if (possible) {
    for (i=0;i<numberRows_;i++) {
      if (rowLower_[i]<-1.0e30&&rowUpper_[i]>1.0e30) 
	possible[i]=0;
    }
  }
  int * index = new int[numberRows_];
  int nDrop=0,nKeep=0;
  for (i=0;i<numberRows_;i++) {
    if (rowLower_[i]<-1.0e30&&rowUpper_[i]>1.0e30) {
      index[nDrop++]=i;
    } else {
      rowLower_[nKeep]=rowLower_[i];
      rowUpper_[nKeep++]=rowUpper_[i];
    }
  }
  numberRows_=nKeep;
  if (nDrop)
    rowCopy_->deleteRows(nDrop,index);
  delete [] index;
  // add in objective 
  int * columns = new int[numberColumns_];
  double * elements = new double[numberColumns_];
  int n=0;
  const double * objective = si.getObjCoefficients();
  bool maximize = (si.getObjSense()==-1);
  for (i=0;i<numberColumns_;i++) {
    if (objective[i]) {
      elements[n]= (maximize) ? -objective[i] : objective[i];
      columns[n++]=i;
    }
  }
  rowCopy_->appendRow(n,columns,elements);
  delete [] columns;
  delete [] elements;
  numberRows_++;
  
}
// Delete snapshot
void CglProbing::deleteSnapshot()
{
  delete [] rowLower_;
  delete [] rowUpper_;
  delete [] colLower_;
  delete [] colUpper_;
  delete rowCopy_;
  numberRows_=0;
  numberColumns_=0;
  rowCopy_=NULL;
  rowLower_=NULL;
  rowUpper_=NULL;
  colLower_=NULL;
  colUpper_=NULL;
  int i;
  for (i=0;i<numberIntegers_;i++) {
    delete cutVector_[i].newValue;
    delete cutVector_[i].index;
  }
  delete [] cutVector_;
  numberIntegers_=0;
  cutVector_=NULL;
}
// Mode stuff
void CglProbing::setMode(int mode)
{
  if (mode>=0&&mode<3)
    mode_=mode;
}
int CglProbing::getMode() const
{
  return mode_;
}
// Set maximum number of passes per node
void CglProbing::setMaxPass(int value)
{
  if (value>0)
    maxPass_=value;
}
// Get maximum number of passes per node
int CglProbing::getMaxPass() const
{
  return maxPass_;
}
// Set maximum number of unsatisfied variables to look at
void CglProbing::setMaxProbe(int value)
{
  if (value>0)
    maxProbe_=value;
}
// Get maximum number of unsatisfied variables to look at
int CglProbing::getMaxProbe() const
{
  return maxProbe_;
}
// Set maximum number of variables to look at in one probe
void CglProbing::setMaxLook(int value)
{
  if (value>0)
    maxStack_=value;
}
// Get maximum number of variables to look at in one probe
int CglProbing::getMaxLook() const
{
  return maxStack_;
}
// Set whether to use objective
void CglProbing::setUsingObjective(bool yesNo)
{
  usingObjective_=yesNo;
}
// Get whether objective is being used
int CglProbing::getUsingObjective() const
{
  return usingObjective_;
}
// Decide whether to do row cuts
void CglProbing::setRowCuts(int type)
{
  if (type>=0&&type<4)
    rowCuts_=type;
}
// Returns row cuts generation type
int CglProbing::rowCuts() const
{
  return rowCuts_;
}
// Returns tight lower
const double * CglProbing::tightLower() const
{
  return colLower_;
}
// Returns tight upper
const double * CglProbing::tightUpper() const
{
  return colUpper_;
}
// Returns relaxed Row lower
const double * CglProbing::relaxedRowLower() const
{
  return rowLower_;
}
// Returns relaxed Row upper
const double * CglProbing::relaxedRowUpper() const
{
  return rowUpper_;
}


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing ()
:
CglCutGenerator(),
primalTolerance_(1.0e-07),
mode_(1),
rowCuts_(1),
maxPass_(3),
maxProbe_(100),
maxStack_(50),
usingObjective_(false)
{

  numberRows_=0;
  numberColumns_=0;
  rowCopy_=NULL;
  rowLower_=NULL;
  rowUpper_=NULL;
  colLower_=NULL;
  colUpper_=NULL;
  numberIntegers_=0;
  cutVector_=NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing (  const CglProbing & source)
                                                              :
  CglCutGenerator(source),
  primalTolerance_(source.primalTolerance_),
  mode_(source.mode_),
  rowCuts_(source.rowCuts_),
  maxPass_(source.maxPass_),
  maxProbe_(source.maxProbe_),
  maxStack_(source.maxStack_),
  usingObjective_(source.usingObjective_)
{  
  numberRows_=source.numberRows_;
  numberColumns_=source.numberColumns_;
  if (numberRows_>0) {
    rowCopy_= new CoinPackedMatrix(*(source.rowCopy_));
    rowLower_=new double[numberRows_];
    memcpy(rowLower_,source.rowLower_,numberRows_*sizeof(double));
    rowUpper_=new double[numberRows_];
    memcpy(rowUpper_,source.rowUpper_,numberRows_*sizeof(double));
    colLower_=new double[numberColumns_];
    memcpy(colLower_,source.colLower_,numberColumns_*sizeof(double));
    colUpper_=new double[numberColumns_];
    memcpy(colUpper_,source.colUpper_,numberColumns_*sizeof(double));
    int i;
    numberIntegers_=source.numberIntegers_;
    cutVector_=new disaggregation [numberIntegers_];
    memcpy(cutVector_,source.cutVector_,numberIntegers_*sizeof(disaggregation));
    for (i=0;i<numberIntegers_;i++) {
      if (cutVector_[i].newValue) {
	cutVector_[i].newValue = new double [cutVector_[i].lastLBWhenAt1];
	memcpy(cutVector_[i].newValue,source.cutVector_[i].newValue,
	       cutVector_[i].lastLBWhenAt1*sizeof(double));
	cutVector_[i].index = new int [cutVector_[i].lastLBWhenAt1];
	memcpy(cutVector_[i].index,source.cutVector_[i].index,
	       cutVector_[i].lastLBWhenAt1*sizeof(int));
      }
    }
  } else {
    rowCopy_=NULL;
    rowLower_=NULL;
    rowUpper_=NULL;
    colLower_=NULL;
    colUpper_=NULL;
    numberIntegers_=0;
    cutVector_=NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglProbing::~CglProbing ()
{
  // free memory
  delete [] rowLower_;
  delete [] rowUpper_;
  delete [] colLower_;
  delete [] colUpper_;
  delete rowCopy_;
  int i;
  for (i=0;i<numberIntegers_;i++) {
    delete cutVector_[i].newValue;
    delete cutVector_[i].index;
  }
  delete [] cutVector_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglProbing &
CglProbing::operator=(
                                         const CglProbing& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    primalTolerance_=rhs.primalTolerance_;
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    delete [] rowLower_;
    delete [] rowUpper_;
    delete [] colLower_;
    delete [] colUpper_;
    delete rowCopy_;
    mode_=rhs.mode_;
    rowCuts_=rhs.rowCuts_;
    maxPass_=rhs.maxPass_;
    maxProbe_=rhs.maxProbe_;
    maxStack_=rhs.maxStack_;
    usingObjective_=rhs.usingObjective_;
    if (numberRows_>0) {
      rowCopy_= new CoinPackedMatrix(*(rhs.rowCopy_));
      rowLower_=new double[numberRows_];
      memcpy(rowLower_,rhs.rowLower_,numberRows_*sizeof(double));
      rowUpper_=new double[numberRows_];
      memcpy(rowUpper_,rhs.rowUpper_,numberRows_*sizeof(double));
      colLower_=new double[numberColumns_];
      memcpy(colLower_,rhs.colLower_,numberColumns_*sizeof(double));
      colUpper_=new double[numberColumns_];
      memcpy(colUpper_,rhs.colUpper_,numberColumns_*sizeof(double));
      int i;
      numberIntegers_=rhs.numberIntegers_;
      cutVector_=new disaggregation [numberIntegers_];
      memcpy(cutVector_,rhs.cutVector_,numberIntegers_*sizeof(disaggregation));
      for (i=0;i<numberIntegers_;i++) {
	if (cutVector_[i].newValue) {
	  cutVector_[i].newValue = new double [cutVector_[i].lastLBWhenAt1];
	  memcpy(cutVector_[i].newValue,rhs.cutVector_[i].newValue,
		 cutVector_[i].lastLBWhenAt1*sizeof(double));
	  cutVector_[i].index = new int [cutVector_[i].lastLBWhenAt1];
	  memcpy(cutVector_[i].index,rhs.cutVector_[i].index,
		 cutVector_[i].lastLBWhenAt1*sizeof(int));
	}
      }
    } else {
      rowCopy_=NULL;
      rowLower_=NULL;
      rowUpper_=NULL;
      colLower_=NULL;
      colUpper_=NULL;
      numberIntegers_=0;
      cutVector_=NULL;
    }
  }
  return *this;
}

/// This can be used to refresh any inforamtion
void 
CglProbing::refreshSolver(OsiSolverInterface * solver)
{
}
