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
#include "CglMixedIntegerRounding.hpp"


void
CglMixedIntegerRoundingUnitTest(const OsiSolverInterface *baseSiP,
			    const std::string mpsDir)
{
  // Test default constructor
  {
    CglMixedIntegerRounding aGenerator;
  }
  
  // Test copy & assignment
  {
    CglMixedIntegerRounding rhs;
    {
      CglMixedIntegerRounding bGenerator;
      CglMixedIntegerRounding cGenerator(bGenerator);
      rhs=bGenerator;
    }
  }

  // Test get/set methods
  {
    CglMixedIntegerRounding getset;

    int gagg = 10 * getset.getMAXAGGR_();
    getset.setMAXAGGR_(gagg);
    int gagg2 = getset.getMAXAGGR_();
    assert(gagg == gagg2);

    bool gmult = !getset.getMULTIPLY_();
    getset.setMULTIPLY_(gmult);
    bool gmult2 = getset.getMULTIPLY_();
    assert(gmult == gmult2);

    int gcrit = getset.getCRITERION_();
    gcrit = (gcrit) % 3 + 1;
    getset.setCRITERION_(gcrit);
    int gcrit2 = getset.getCRITERION_();
    assert(gcrit == gcrit2);

    int gpre = getset.getDoPreproc();
    gpre = (gpre + 1) % 3 - 1;
    getset.setDoPreproc(gpre);
    int gpre2 = getset.getDoPreproc();
    assert(gpre == gpre2);
  }

  // Test generateCuts
  {
    CglMixedIntegerRounding gct;
    OsiSolverInterface  *siP = baseSiP->clone();
    std::string fn = mpsDir+"capPlan1";
    std::string fn2 = mpsDir+"capPlan1.mps";
    FILE *in_f = fopen(fn2.c_str(), "r");
    if(in_f == NULL) {
      printf("Can not open file %s;\nSkip test of CglMixedIntegerRounding::generateCuts()\n", fn2.c_str());
    }
    else {
      fclose(in_f);
      siP->readMps(fn.c_str(),"mps");
 
     siP->initialSolve();
     double lpRelax = siP->getObjValue();
 
     OsiCuts cs;
     gct.generateCuts(*siP, cs);
     int nRowCuts = cs.sizeRowCuts();
     std::cout<<"There are "<<nRowCuts<<" MIR cuts"<<std::endl;
     assert(cs.sizeRowCuts() > 0);
     OsiSolverInterface::ApplyCutsReturnCode rc = siP->applyCuts(cs);
     
     siP->resolve();
     
     double lpRelaxAfter= siP->getObjValue(); 
     
#ifdef CGL_DEBUG
     printf("Initial LP value: %f\n", lpRelax);
     printf("LP value with cuts: %f\n", lpRelaxAfter);
#endif
     assert( lpRelax < lpRelaxAfter );
    }
    delete siP;
  }

}

