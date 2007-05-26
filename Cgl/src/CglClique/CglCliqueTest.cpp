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
#include "CglClique.hpp"


void
CglCliqueUnitTest(const OsiSolverInterface *baseSiP,
			    const std::string mpsDir)
{
  // Test default constructor
  {
    CglClique aGenerator;
  }
  
  // Test copy & assignment
  {
    CglClique rhs;
    {
      CglClique bGenerator;
      CglClique cGenerator(bGenerator);
      //rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglClique getset;
    // None to test
  }

  // Test generateCuts
  {
    CglClique gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"l152lav";
    std::string fn2 = mpsDir+"l152lav.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      printf("Can not open file %s;\nSkip test of CglClique::generateCuts()\n", fn2.c_str());
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(),"mps");
 
     siP->initialSolve();
     double lpRelax = siP->getObjValue();
 
     OsiCuts cs;
     gct.generateCuts(*siP, cs);
     int nRowCuts = cs.sizeRowCuts();
     std::cout<<"There are "<<nRowCuts<<" Clique cuts"<<std::endl;
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

