// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
//#include <stdlib.h>
#include <assert.h>
//#include <float.h>
//#include <iostream>

#include "CglCutGenerator.hpp"
 

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglCutGenerator::CglCutGenerator ()
{
  // nothing to do here
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglCutGenerator::CglCutGenerator (
                  const CglCutGenerator & source)         
{  
  // nothing to do here
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglCutGenerator::~CglCutGenerator ()
{
  // nothing to do here
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglCutGenerator &
CglCutGenerator::operator=(
                   const CglCutGenerator& rhs)
{
  if (this != &rhs) {
    // nothing to do here
  }
  return *this;
}



#ifdef NDEBUG
#undef NDEBUG
#endif

#if 0
//--------------------------------------------------------------------------
// test EKKsolution methods.
//--------------------------------------------------------------------------
void
CglCutGenerator::unitTest()
{
}
#endif
