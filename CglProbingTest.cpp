// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cstdio>
#include <cassert>

#include "OsiPackedMatrix.hpp"
#include "CglProbing.hpp"


#ifdef NDEBUG
#undef NDEBUG
#endif

//--------------------------------------------------------------------------
// test EKKsolution methods.
void
CglProbingUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{
  OsiRelFltEq eq(0.000001);

  // Test default constructor
  {
    CglProbing aGenerator;
  }
  
  // Test copy & assignment
  {
    CglProbing rhs;
    {
      CglProbing bGenerator;
      CglProbing cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  {
    OsiCuts osicuts;
    CglProbing test1;
    OsiSolverInterface  * siP = baseSiP->clone();
    int i;
    int nColCuts;
    int nRowCuts;
    
    std::string fn = mpsDir+"p0033";
    siP->readMps(fn.c_str(),"mps");
    siP->initialSolve();
    // just unsatisfied variables
    test1.generateCuts(*siP,osicuts);
    nColCuts = osicuts.sizeColCuts();
    nRowCuts = osicuts.sizeRowCuts();
    cout<<"There are "<<nRowCuts<<" probing cuts"<<endl;
    {
      cout<<"there are "<<nColCuts<<" probing column cuts"<<endl;
      const double * lo = siP->getColLower();
      const double * up = siP->getColUpper();
      for (i=0; i<nColCuts; i++){
	OsiColCut ccut;
	OsiPackedVector cpv;
	ccut = osicuts.colCut(i);
	cpv = ccut.lbs();
	int n = cpv.getNumElements();
        int j;
	const int * indices = cpv.getIndices();
	double* elements = cpv.getElements();
	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]>lo[icol])
	    cout<<"Can increase lb on "<<icol<<" from "<<lo[icol]<<
	      " to "<<elements[j]<<endl;
	}
	cpv = ccut.ubs();
	n = cpv.getNumElements();
	indices = cpv.getIndices();
	elements = cpv.getElements();
	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]<up[icol])
	    cout<<"Can decrease ub on "<<icol<<" from "<<up[icol]<<
	      " to "<<elements[j]<<endl;
	}
      }
    }
    for (i=0; i<nRowCuts; i++){
      OsiRowCut rcut;
      OsiPackedVector rpv;
      const double * colsol = siP->getColSolution();
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      double lb=rcut.lb();
      double ub=rcut.ub();
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	cout <<endl;
      }
    }
    OsiPackedVector check;
    int index[] = {6,32};
    double el[] = {1,1};
    check.setVector(2,index,el);
    assert (osicuts.sizeRowCuts()==2);
    // sort Elements in increasing order
    OsiPackedVector rpv=osicuts.rowCut(1).row();
    rpv.sortIncrIndex();
    assert (check==rpv);
    assert (osicuts.rowCut(1).lb()==1.0);
    // now all variables
    osicuts=OsiCuts();
    test1.setMode(2);
    test1.generateCuts(*siP,osicuts);
    nColCuts = osicuts.sizeColCuts();
    nRowCuts = osicuts.sizeRowCuts();
    cout<<"There are "<<nRowCuts<<" probing cuts"<<endl;
    {
      cout<<"there are "<<nColCuts<<" probing column cuts"<<endl;
      const double * lo = siP->getColLower();
      const double * up = siP->getColUpper();
      for (i=0; i<nColCuts; i++){
	OsiColCut ccut;
	OsiPackedVector cpv;
	ccut = osicuts.colCut(i);
	cpv = ccut.lbs();
	int n = cpv.getNumElements();
        int j;
	const int * indices = cpv.getIndices();
	double* elements = cpv.getElements();
	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]>lo[icol])
	    cout<<"Can increase lb on "<<icol<<" from "<<lo[icol]<<
	      " to "<<elements[j]<<endl;
	}
	cpv = ccut.ubs();
	n = cpv.getNumElements();
	indices = cpv.getIndices();
	elements = cpv.getElements();
	for (j=0;j<n;j++) {
	  int icol=indices[j];
	  if (elements[j]<up[icol])
	    cout<<"Can decrease ub on "<<icol<<" from "<<up[icol]<<
	      " to "<<elements[j]<<endl;
	}
      }
    }
    for (i=0; i<nRowCuts; i++){
      OsiRowCut rcut;
      OsiPackedVector rpv;
      const double * colsol = siP->getColSolution();
      rcut = osicuts.rowCut(i);
      rpv = rcut.row();
      const int n = rpv.getNumElements();
      const int * indices = rpv.getIndices();
      double* elements = rpv.getElements();
      double sum2=0.0;
      int k=0;
      double lb=rcut.lb();
      double ub=rcut.ub();
      for (k=0; k<n; k++){
	int column=indices[k];
	sum2 += colsol[column]*elements[k];
      }
      if (sum2 >ub + 1.0e-7 ||sum2 < lb - 1.0e-7) {
	cout<<"Cut "<<i<<" lb "<<lb<<" solution "<<sum2<<" ub "<<ub<<endl;
	for (k=0; k<n; k++){
	  int column=indices[k];
	  cout<<"(col="<<column<<",el="<<elements[k]<<",sol="<<
	    colsol[column]<<") ";
	}
	cout <<endl;
      }
    }
    assert (osicuts.sizeRowCuts()==10);
    delete siP;
  }

}

