// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).


/*
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CglTreeInfo.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "CoinPackedMatrix.hpp"
#include "CglStored.hpp"
*/

#include "CoinPragma.hpp"
#include "OsiRowCut.hpp"
#include "CglTreeInfo.hpp"

// Default constructor 
CglTreeInfo::CglTreeInfo ()
  : level(-1), pass(-1), formulation_rows(-1), options(0), inTree(false),
    strengthenRow(NULL),randomNumberGenerator(NULL) {}

// Copy constructor 
CglTreeInfo::CglTreeInfo (const CglTreeInfo & rhs)
  : level(rhs.level), 
    pass(rhs.pass), 
    formulation_rows(rhs.formulation_rows), 
    options(rhs.options),
    inTree(rhs.inTree),
    strengthenRow(rhs.strengthenRow),
    randomNumberGenerator(rhs.randomNumberGenerator)
{
}
// Clone
CglTreeInfo * 
CglTreeInfo::clone() const
{
  return new CglTreeInfo(*this) ;
}

// Assignment operator 
CglTreeInfo &
CglTreeInfo::operator=(const CglTreeInfo& rhs)
{
  if (this != &rhs) {
    level = rhs.level; 
    pass = rhs.pass; 
    formulation_rows = rhs.formulation_rows; 
    options = rhs.options ;
    inTree = rhs.inTree ;
    strengthenRow = rhs.strengthenRow ;
    randomNumberGenerator = rhs.randomNumberGenerator ;
  }
  return *this ;
}
  
 // Destructor 
CglTreeInfo::~CglTreeInfo ()
{
}
