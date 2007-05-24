// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// Test individual classes or groups of classes

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CglConfig.h"

#include <string>
#include <cstring>
#include <cassert>
#include <iostream>
#include <cstdlib>

#ifdef COIN_HAS_OSL
#include <OsiOslSolverInterface.hpp>
#endif
#ifdef COIN_HAS_CPX
#include <OsiCpxSolverInterface.hpp>
#endif
#ifdef COIN_HAS_XPR
#include <OsiXprSolverInterface.hpp>
#endif
#ifdef COIN_HAS_CLP
#include <OsiClpSolverInterface.hpp>
#endif
#ifdef COIN_HAS_DYLP
#include <OsiDylpSolverInterface.hpp>
#endif

#include "CglSimpleRounding.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglProbing.hpp"
#include "CglGomory.hpp"
#include "CglLandP.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglResidualCapacity.hpp"
//#include "CglFlowCover.hpp"

// Function Prototypes. Function definitions is in this file.
void testingMessage( const char * const msg );

// Command line parameter is directory containing data files.
// If not specified, then "../../Data/Sample/" and
// CglTestData/ are used

int main (int argc, const char *argv[])
{
  // Set directory containing data files.
  std::string mpsDir;
  std::string testDir;
  if (argc >= 2) {
    mpsDir = argv[1];
    testDir = argv[1];
  }
  else {
    mpsDir ="../../Data/Sample/";
    testDir = "CglTestData/";
  }

#ifdef COIN_HAS_OSL
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglGomory with OsiOslSolverInterface\n" );
    CglGomoryUnitTest(&oslSi,mpsDir);
  }  
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglSimpleRounding with OsiOslSolverInterface\n" );
    CglSimpleRoundingUnitTest(&oslSi,mpsDir);
  } 
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglKnapsackCover with OsiOslSolverInterface\n" );
    CglKnapsackCoverUnitTest(&oslSi,mpsDir);
  }  
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglOddHole with OsiOslSolverInterface\n" );
    CglOddHoleUnitTest(&oslSi,mpsDir);
  }  
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglProbing with OsiOslSolverInterface\n" );
    CglProbingUnitTest(&oslSi,mpsDir);
  }  
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglMixedIntegerRounding with OsiOslSolverInterface\n" );
    CglMixedIntegerRoundingUnitTest(&oslSi, testDir);
  }  
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglMixedIntegerRounding2 with OsiOslSolverInterface\n" );
    CglMixedIntegerRounding2UnitTest(&oslSi, testDir);
  }  
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing CglResidualCapacity with OsiOslSolverInterface\n" );
    CglResidualCapacityUnitTest(&oslSi, testDir);
  }  
  {
    OsiOslSolverInterface oslSi;
    //    testingMessage( "Testing CglFlowCover with OsiOslSolverInterface\n" );
    //    CglFlowCoverUnitTest(&oslSi,mpsDir);
  }

#endif
#ifdef COIN_HAS_CPX
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglSimpleRounding with OsiCpxSolverInterface\n" );
    CglSimpleRoundingUnitTest(&cpxSi,mpsDir);
  } 
  {
    OsiCpxSolverInterface cpxSi;
    // Test does not work with Cplex
    //testingMessage( "Testing CglKnapsackCover with OsiCpxSolverInterface\n" );
    //CglKnapsackCoverUnitTest(&cpxSi,mpsDir);
  }  
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglOddHole with OsiCpxSolverInterface\n" );
    CglOddHoleUnitTest(&cpxSi,mpsDir);
  }  
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglProbing with OsiCpxSolverInterface\n" );
    CglProbingUnitTest(&cpxSi,mpsDir);
  }  
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglMixedIntegerRounding with OsiCpxSolverInterface\n" );
    CglMixedIntegerRoundingUnitTest(&cpxSi, testDir);
  }  
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglMixedIntegerRounding2 with OsiCpxSolverInterface\n" );
    CglMixedIntegerRounding2UnitTest(&cpxSi, testDir);
  }  
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglResidualCapacity with OsiCpxSolverInterface\n" );
    CglResidualCapacityUnitTest(&cpxSi, testDir);
  }  
  {
    OsiCpxSolverInterface cpxSi;
    // Test does not work with Cplex
    //    testingMessage( "Testing CglFlowCover with OsiCpxSolverInterface\n" );
    //    CglFlowCoverUnitTest(&cpxSi,mpsDir);
  }

#endif

#ifdef COIN_HAS_XPR
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglSimpleRounding with OsiXprSolverInterface\n" );
    CglSimpleRoundingUnitTest(&xprSi,mpsDir);
  } 
  {
    OsiXprSolverInterface xprSi;
    //testingMessage( "Testing CglKnapsackCover with OsiXprSolverInterface\n" );
    //CglKnapsackCoverUnitTest(&xprSi,mpsDir);
  }  
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglOddHole with OsiXprSolverInterface\n" );
    CglOddHoleUnitTest(&xprSi,mpsDir);
  }  
  {
    OsiXprSolverInterface xprSi;
    // Does not work with Xpress
    //testingMessage( "Testing CglProbing with OsiXprSolverInterface\n" );
    //CglProbingUnitTest(&xprSi,mpsDir);
  }  
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglMixedIntegerRounding with OsiXprSolverInterface\n" );
    CglMixedIntegerRoundingUnitTest(&xprSi, testDir);
  }  
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglMixedIntegerRounding2 with OsiXprSolverInterface\n" );
    CglMixedIntegerRounding2UnitTest(&xprSi, testDir);
  }  
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglResidualCapacity with OsiXprSolverInterface\n" );
    CglResidualCapacityUnitTest(&xprSi, testDir);
  }  
  {
    OsiXprSolverInterface xprSi;
    //    testingMessage( "Testing CglFlowCover with OsiXprSolverInterface\n" );
    //    CglFlowCoverUnitTest(&xprSi,mpsDir);
  }

#endif
#ifdef COIN_HAS_CLP
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglGomory with OsiClpSolverInterface\n" );
    CglGomoryUnitTest(&clpSi,mpsDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglLandp with OsiClpSolverInterface\n" );
    CglLandPUnitTest(&clpSi,mpsDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglSimpleRounding with OsiClpSolverInterface\n" );
    CglSimpleRoundingUnitTest(&clpSi,mpsDir);
  } 
  if (0) {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglKnapsackCover with OsiClpSolverInterface\n" );
    CglKnapsackCoverUnitTest(&clpSi,mpsDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglOddHole with OsiClpSolverInterface\n" );
    CglOddHoleUnitTest(&clpSi,mpsDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglProbing with OsiClpSolverInterface\n" );
    CglProbingUnitTest(&clpSi,mpsDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglMixedIntegerRounding with OsiClpSolverInterface\n" );
    CglMixedIntegerRoundingUnitTest(&clpSi, testDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglMixedIntegerRounding2 with OsiClpSolverInterface\n" );
    CglMixedIntegerRounding2UnitTest(&clpSi, testDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglResidualCapacity with OsiClpSolverInterface\n" );
    CglResidualCapacityUnitTest(&clpSi, testDir);
  }  
  {
    OsiClpSolverInterface clpSi;
    //    testingMessage( "Testing CglFlowCover with OsiClpSolverInterface\n" );
    //    CglFlowCoverUnitTest(&clpSi,mpsDir);
  }

#endif
#ifdef COIN_HAS_DYLP
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglGomory with OsiDylpSolverInterface\n" );
    CglGomoryUnitTest(&dylpSi,mpsDir);
  }  
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglSimpleRounding with OsiDylpSolverInterface\n" );
    CglSimpleRoundingUnitTest(&dylpSi,mpsDir);
  } 
  if (0) {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglKnapsackCover with OsiDylpSolverInterface\n" );
    CglKnapsackCoverUnitTest(&dylpSi,mpsDir);
  }  
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglOddHole with OsiDylpSolverInterface\n" );
    CglOddHoleUnitTest(&dylpSi,mpsDir);
  }  
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglProbing with OsiDylpSolverInterface\n" );
    CglProbingUnitTest(&dylpSi,mpsDir);
  }  
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglMixedIntegerRounding with OsiDylpSolverInterface\n" );
    CglMixedIntegerRoundingUnitTest(&dylpSi, testDir);
  }  
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglMixedIntegerRounding2 with OsiDylpSolverInterface\n" );
    CglMixedIntegerRounding2UnitTest(&dylpSi, testDir);
  }  
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglResidualCapacity with OsiDylpSolverInterface\n" );
    CglResidualCapacityUnitTest(&dylpSi, testDir);
  }  
  {
    OsiDylpSolverInterface dylpSi;
    //    testingMessage( "Testing CglFlowCover with OsiDylpSolverInterface\n" );
    //    CglFlowCoverUnitTest(&dylpSi,mpsDir);
  }

#endif

  testingMessage( "All tests completed successfully\n" );
  return 0;
}

 
// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
  std::cout <<std::endl <<"*****************************************"
            <<std::endl <<msg <<std::endl;
  //std::cerr <<msg;
}

