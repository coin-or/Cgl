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
#ifdef COIN_HAS_GLPK
#include <OsiGlpkSolverInterface.hpp>
#endif
#ifdef COIN_HAS_VOL
#include <OsiVolSolverInterface.hpp>
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
#include "CglRedSplit.hpp"
#include "CglTwomir.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"

// Function Prototypes. Function definitions is in this file.
void testingMessage( const char * const msg );

// Command line parameter is directory containing data files.
// If not specified, then "../../Data/Sample/" and
// "CglTestData/" are used

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
    mpsDir = "../../Data/Sample/";
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
    testingMessage( "Testing CglFlowCover with OsiOslSolverInterface\n" );
    CglFlowCoverUnitTest(&oslSi,testDir);
  }

#endif
#ifdef COIN_HAS_CPX
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglGomory with OsiCpxSolverInterface\n" );
    CglGomoryUnitTest(&cpxSi,mpsDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglSimpleRounding with OsiCpxSolverInterface\n" );
    CglSimpleRoundingUnitTest(&cpxSi,mpsDir);
  }
  if(0) // Test does not work with Cplex
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglKnapsackCover with OsiCpxSolverInterface\n" );
    CglKnapsackCoverUnitTest(&cpxSi,mpsDir);
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
    testingMessage( "Testing CglRedSplit with OsiCpxSolverInterface\n" );
    CglRedSplitUnitTest(&cpxSi, mpsDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglTwomir with OsiCpxSolverInterface\n" );
    CglTwomirUnitTest(&cpxSi, testDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglClique with OsiCpxSolverInterface\n" );
    CglCliqueUnitTest(&cpxSi, testDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing CglFlowCover with OsiCpxSolverInterface\n" );
    CglFlowCoverUnitTest(&cpxSi, testDir);
  }

#endif

#ifdef COIN_HAS_XPR
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglGomory with OsiXprSolverInterface\n" );
    CglGomoryUnitTest(&xprSi,mpsDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglSimpleRounding with OsiXprSolverInterface\n" );
    CglSimpleRoundingUnitTest(&xprSi,mpsDir);
  }
  if(0) 
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglKnapsackCover with OsiXprSolverInterface\n" );
    CglKnapsackCoverUnitTest(&xprSi,mpsDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglOddHole with OsiXprSolverInterface\n" );
    CglOddHoleUnitTest(&xprSi,mpsDir);
  }
  if(0)     // Does not work with Xpress
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglProbing with OsiXprSolverInterface\n" );
    CglProbingUnitTest(&xprSi,mpsDir);
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
    testingMessage( "Testing CglTwomir with OsiXprSolverInterface\n" );
    CglTwomirUnitTest(&xprSi, testDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglClique with OsiXprSolverInterface\n" );
    CglCliqueUnitTest(&xprSi, testDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing CglFlowCover with OsiXprSolverInterface\n" );
    CglFlowCoverUnitTest(&xprSi, testDir);
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
    testingMessage( "Testing CglRedSplit with OsiClpSolverInterface\n" );
    CglRedSplitUnitTest(&clpSi, mpsDir);
  }
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglTwomir with OsiClpSolverInterface\n" );
    CglTwomirUnitTest(&clpSi, testDir);
  }
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglClique with OsiClpSolverInterface\n" );
    CglCliqueUnitTest(&clpSi, testDir);
  }
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing CglFlowCover with OsiClpSolverInterface\n" );
    CglFlowCoverUnitTest(&clpSi, testDir);
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
  if (0)  // needs partial OsiSimplex
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglRedSplit with OsiDylpSolverInterface\n" );
    CglRedSplitUnitTest(&dylpSi, mpsDir);
  }
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglTwomir with OsiDylpSolverInterface\n" );
    CglTwomirUnitTest(&dylpSi, testDir);
  }
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglClique with OsiDylpSolverInterface\n" );
    CglCliqueUnitTest(&dylpSi, testDir);
  }
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing CglFlowCover with OsiDylpSolverInterface\n" );
    CglFlowCoverUnitTest(&dylpSi, testDir);
  }

#endif
#ifdef COIN_HAS_GLPK
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglGomory with OsiGlpkSolverInterface\n" );
    CglGomoryUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglSimpleRounding with OsiGlpkSolverInterface\n" );
    CglSimpleRoundingUnitTest(&glpkSi,mpsDir);
  }
  if (0) {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglKnapsackCover with OsiGlpkSolverInterface\n" );
    CglKnapsackCoverUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglOddHole with OsiGlpkSolverInterface\n" );
    CglOddHoleUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglProbing with OsiGlpkSolverInterface\n" );
    CglProbingUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglMixedIntegerRounding with OsiGlpkSolverInterface\n" );
    CglMixedIntegerRoundingUnitTest(&glpkSi, testDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglMixedIntegerRounding2 with OsiGlpkSolverInterface\n" );
    CglMixedIntegerRounding2UnitTest(&glpkSi, testDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglResidualCapacity with OsiGlpkSolverInterface\n" );
    CglResidualCapacityUnitTest(&glpkSi, testDir);
  }
  if (0)  // needs partial OsiSimplex
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglRedSplit with OsiGlpkSolverInterface\n" );
    CglRedSplitUnitTest(&glpkSi, mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglTwomir with OsiGlpkSolverInterface\n" );
    CglTwomirUnitTest(&glpkSi, testDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglClique with OsiGlpkSolverInterface\n" );
    CglCliqueUnitTest(&glpkSi, testDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing CglFlowCover with OsiGlpkSolverInterface\n" );
    CglFlowCoverUnitTest(&glpkSi, testDir);
  }

#endif

#ifdef COIN_HAS_VOL
  if(0) // p0033: LP not solved to optimality: Finds 2142 versus 2520
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglGomory with OsiVolSolverInterface\n" );
    CglGomoryUnitTest(&volSi,mpsDir);
  }
  if(0) // Not expected number of cuts; might come from different solution?
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglSimpleRounding with OsiVolSolverInterface\n" );
    CglSimpleRoundingUnitTest(&volSi,mpsDir);
  }
  if(0) // tp3: LP not solved to optimality: Finds 97.1842 versus 97.185
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglKnapsackCover with OsiVolSolverInterface\n" );
    CglKnapsackCoverUnitTest(&volSi,mpsDir);
  }
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglOddHole with OsiVolSolverInterface\n" );
    CglOddHoleUnitTest(&volSi,mpsDir);
  }
  if(0) // Not expected number of elements in cut; might come from different solution?
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglProbing with OsiVolSolverInterface\n" );
    CglProbingUnitTest(&volSi,mpsDir);
  }
  if(0) // Throw CoinError since solver can not handle infinite bounds
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglMixedIntegerRounding with OsiVolSolverInterface\n" );
    CglMixedIntegerRoundingUnitTest(&volSi, testDir);
  }
  if(0) // Throw CoinError since solver can not handle infinite bounds
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglMixedIntegerRounding2 with OsiVolSolverInterface\n" );
    CglMixedIntegerRounding2UnitTest(&volSi, testDir);
  }
  if(0) // Throw CoinError since solver can not handle infinite bounds
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglResidualCapacity with OsiVolSolverInterface\n" );
    CglResidualCapacityUnitTest(&volSi, testDir);
  }
  if(0) // Throw CoinError since solver can not handle infinite bounds
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglTwomir with OsiVolSolverInterface\n" );
    CglTwomirUnitTest(&volSi, testDir);
  }
  if(0) // No cuts found
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglClique with OsiVolSolverInterface\n" );
    CglCliqueUnitTest(&volSi, testDir);
  }
  if(0) // Throw CoinError since solver can not handle infinite bounds
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing CglFlowCover with OsiVolSolverInterface\n" );
    CglFlowCoverUnitTest(&volSi, testDir);
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

