// $Id: CglZeroHalfTest.cpp 1015 2011-04-29 18:02:51Z stefan $
// Copyright (C) 2010, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "CglZeroHalf.hpp" 
//#include "CglKnapsackCover.hpp" 
#include <stdio.h>

//--------------------------------------------------------------------------
// test the zero half cut generators methods.
void
CglZeroHalfUnitTest(
  const OsiSolverInterface * baseSiP,
  const std::string mpsDir )
{

  // Test default constructor
  {
    CglZeroHalf cg;
  }

  // Test copy & assignment
  {
    CglZeroHalf rhs;
    {
      CglZeroHalf cg;
      CglZeroHalf cgC(cg);
      rhs=cg;
    }
  }




  // Test generate cuts method on p0201
  {
    CglZeroHalf cg;
    
    OsiSolverInterface * siP = baseSiP->clone();
    std::string fn = mpsDir+"lseu";
    siP->readMps(fn.c_str(),"mps");
    siP->initialSolve();
    cg.refreshSolver(siP);
    OsiCuts cuts;
    cg.generateCuts(*siP,cuts);

    // lseu is the optimal solution to lseu
    // Optimal IP solution to lseu    
    int objIndices[13]={0,1,6,13,26,33,38,43,50,52,63,65,85};
    CoinPackedVector lseu(13,objIndices,1.0);

    // test that none of the generated cuts
    // chops off the optimal solution
    int nRowCuts = cuts.sizeRowCuts();
    OsiRowCut rcut;
    CoinPackedVector rpv;
    int i;
    for (i=0; i<nRowCuts; i++){
      rcut = cuts.rowCut(i);
      rpv = rcut.row();
      double lseuSum = (rpv*lseu).sum();
      double rcutub = rcut.ub();
      assert (lseuSum <= rcutub);
    }

    // test that the cuts improve the 
    // lp objective function value
    double lpRelaxBefore=siP->getObjValue();
    OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cuts);
    siP->resolve();
    double lpRelaxAfter=siP->getObjValue(); 
#ifdef CGL_DEBUG
    printf("\n\nOrig LP min=%f\n",lpRelaxBefore);
    printf("Final LP min=%f\n\n",lpRelaxAfter);
#endif
    printf("Zero cuts %d\n",nRowCuts);
    assert( lpRelaxBefore < lpRelaxAfter );
    printf("Good zero %s\n",fn.c_str());

    delete siP;

  }


}
