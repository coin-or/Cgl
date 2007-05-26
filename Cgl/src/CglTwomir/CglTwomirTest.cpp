// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cstdio>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include "CglTwomir.hpp"


void
CglTwomirUnitTest(const OsiSolverInterface *baseSiP,
		  const std::string mpsDir)
{
  // Test default constructor
  {
    CglTwomir aGenerator;
  }
  
  // Test copy & assignment
  {
    CglTwomir rhs;
    {
      CglTwomir bGenerator;
      CglTwomir cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglTwomir getset;
    
    int gtmin = getset.getTmin() + 1;
    int gtmax = getset.getTmax() + 1;
    getset.setMirScale(gtmin, gtmax);
    double gtmin2 = getset.getTmin();
    double gtmax2 = getset.getTmax();
    assert(gtmin == gtmin2);
    assert(gtmax == gtmax2);

    int gamax = 2 * getset.getAmax() + 1;
    getset.setAMax(gamax);
    int gamax2 = getset.getAmax();
    assert(gamax == gamax2);
  }

  // Test generateCuts
  {
    CglTwomir gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"capPlan1";
    std::string fn2 = mpsDir+"capPlan1.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      printf("Can not open file %s;\nSkip test of CglTwomir::generateCuts()\n", fn2.c_str());
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(),"mps");
 
     siP->initialSolve();
     double lpRelax = siP->getObjValue();
 
     OsiCuts cs;
     gct.generateCuts(*siP, cs);
     int nRowCuts = cs.sizeRowCuts();
     std::cout<<"There are "<<nRowCuts<<" Twomir cuts"<<std::endl;
     assert(cs.sizeRowCuts() > 0);
     OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cs);
     
     siP->resolve();
     
     double lpRelaxAfter= siP->getObjValue(); 
     printf("Initial LP value: %f\n", lpRelax);
     printf("LP value with cuts: %f\n", lpRelaxAfter);
     assert( lpRelax < lpRelaxAfter );
    }
    delete siP;
  }

}

