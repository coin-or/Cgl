// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cstdio>


#include "CglProbing.hpp"
#include "CglProbingDebug.hpp"

namespace {

void printCuts (const OsiSolverInterface *siP, OsiCuts &osicuts)
{
    int nRowCuts = osicuts.sizeRowCuts() ;
    std::cout << "There are " << nRowCuts << " probing row cuts" << std::endl ;
#   ifdef CGL_DEBUG
    for (int i = 0 ; i < nRowCuts ; i++) {
      OsiRowCut rcut ;
      CoinPackedVector rpv ;
      const double *colsol = siP->getColSolution() ;
      rcut = osicuts.rowCut(i) ;
      rpv = rcut.row() ;
      int n = rpv.getNumElements() ;
      const int *indices = rpv.getIndices() ;
      const double *elements = rpv.getElements() ;
      double sum2 = 0.0 ;
      double lb = rcut.lb() ;
      double ub = rcut.ub() ;
      for (int k = 0 ; k < n ; k++) {
	int j = indices[k] ;
	sum2 += colsol[j]*elements[k] ;
      }
      if (sum2 > ub+1.0e-7 || sum2 < lb-1.0e-7) {
	std::cout
	  << "Cut " << i << ": " << lb << " <= " << sum2
	  << " <= " << ub << std::endl ;
	for (int k = 0 ; k < n ; k++) {
	  int j = indices[k] ;
	  std::cout
	    << "(" << elements[k] << ")*("
	    << siP->getColName(j) << "(" << j << ") = "
	    << colsol[j] << ") " ;
	}
	std::cout << std::endl ;
      }
    }
#   endif

    int nColCuts = osicuts.sizeColCuts() ;
    std::cout
      << "There are " << nColCuts << " probing column cuts" << std::endl ;
#   if CGL_DEBUG > 0
    {
      const double *lo = siP->getColLower() ;
      const double *up = siP->getColUpper() ;
      for (int k = 0 ; k < nColCuts ; k++){
	OsiColCut ccut ;
	CoinPackedVector cpv ;
	ccut = osicuts.colCut(k) ;
	cpv = ccut.lbs() ;
	int n = cpv.getNumElements() ;
	const int *indices = cpv.getIndices() ;
	const double *elements = cpv.getElements() ;
	for (int kk = 0 ; kk < n ; kk++) {
	  int j = indices[kk] ;
	  if (elements[kk] > lo[j])
	    std::cout
	      << "Can increase lb on " << j << " from " << lo[j]
	      << " to " << elements[kk] << std::endl ;
	}
	cpv = ccut.ubs() ;
	n = cpv.getNumElements() ;
	indices = cpv.getIndices() ;
	elements = cpv.getElements() ;
	for (int kk = 0 ; kk < n ; kk ++) {
	  int j = indices[kk] ;
	  if (elements[kk] < up[j])
	    std::cout
	      << "Can decrease ub on " << j << " from " << up[j]
	      << " to " << elements[kk] << std::endl ;
	}
      }
    }
#   endif

  return ;
}

}  // end file-local namespace

//--------------------------------------------------------------------------

/*
  The integer return is purely for show, but at least we're headed in the right
  direction.
*/
int CglProbingUnitTest (const OsiSolverInterface *baseSiP,
			const std::string mpsDir)
{
  CoinRelFltEq eq(0.000001) ;

  // Test default constructor
  {
    CglProbing aGenerator ;
  }
  
  // Test copy & assignment
  {
    CglProbing rhs ;
    {
      CglProbing bGenerator ;
      CglProbing cGenerator(bGenerator) ;
      rhs = bGenerator ;
    }
  }

  {
    OsiCuts osicuts ;
    CglProbing test1 ;
    OsiSolverInterface  *siP = baseSiP->clone() ;
    int nColCuts = -1 ;
    int nRowCuts = -1 ;

    siP->setIntParam(OsiNameDiscipline,1) ;

    std::string localDir = "/devel/Coin/Split/Data/Miplib3/" ;
    // std::string probName = "p0033" ;
    // std::string fn = mpsDir+probName ;
    std::string probName = "bell3a" ;
    std::string fn = localDir+probName ;

    siP->readMps(fn.c_str(),"mps") ;
    siP->initialSolve() ;

#   if CGL_DEBUG > 0
    // Activate row cut debugger, if we're doing serious debugging
    siP->activateRowCutDebugger(probName.c_str()) ;
    const OsiRowCutDebugger *debugger = siP->getRowCutDebugger() ;
    if (debugger == 0)
    { std::cout
	<< "ERROR: Could not activate row cut debugger for "
	<< probName << "." << std::endl ;
    }
#   endif

/*
  Test first in the default mode: probing unsatisfied variables only.
*/
    test1.generateCuts(*siP,osicuts) ;
    printCuts(siP,osicuts) ;

/*
  ?? Are we expecting exactly one row cut from p0033 for this test?
  -- lh, 101203 --

  Apparently no, but this is one of the three that are produced. The three
  that are expected are
      -infty <=  x<5> - x<32> <= 0
           1 <=  x<6> + x<32> <= infty
           0 <= -x<5> + x<32> <= infty
*/
    nRowCuts = osicuts.sizeRowCuts() ;
    if (nRowCuts == 1) {
      CoinPackedVector check ;
      int index[] = {6,32} ;
      double el[] = {1,1} ;
      check.setVector(2,index,el) ;
      // sort Elements in increasing order
      CoinPackedVector rpv = osicuts.rowCut(0).row() ;
      assert (rpv.getNumElements() == 2) ;
      rpv.sortIncrIndex() ;
      assert (check == rpv) ;
      assert (osicuts.rowCut(0).lb() == 1.0) ;
    }

/*
  Now try again, probing all variables.
*/
    osicuts = OsiCuts() ;
    test1.setMode(2) ;
    test1.setRowCuts(3) ;
    test1.generateCuts(*siP,osicuts) ;
    printCuts(siP,osicuts) ;
    nColCuts = osicuts.sizeColCuts() ;
    nRowCuts = osicuts.sizeRowCuts() ;

    assert (nRowCuts >= 4) ;
    delete siP ;
  }

  return (0) ;
}

