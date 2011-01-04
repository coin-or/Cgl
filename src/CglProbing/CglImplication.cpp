/*
  Copyright (C) 2002, International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif


#include "CglImplication.hpp"

//-------------------------------------------------------------
void
CglImplication::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
				const CglTreeInfo info) const
{
  if (probingInfo_) {
    //int n1=cs.sizeRowCuts();
    probingInfo_->generateCuts(si,cs,info);
    //int n2=cs.sizeRowCuts();
    //if (n2>n1)
    //printf("added %d cuts\n",n2-n1);
  }
}


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglImplication::CglImplication ()
:
CglCutGenerator(),
probingInfo_(NULL)
{
  // nothing to do here
}


//-------------------------------------------------------------------
// Constructor with info
//-------------------------------------------------------------------
CglImplication::CglImplication (CglTreeProbingInfo * info)
:
CglCutGenerator(),
probingInfo_(info)
{
  // nothing to do here
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglImplication::CglImplication (
                  const CglImplication & source)
:
CglCutGenerator(source),
probingInfo_(source.probingInfo_)
{  
  // Nothing to do here
}


//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglImplication::clone() const
{
  return new CglImplication(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglImplication::~CglImplication ()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglImplication &
CglImplication::operator=(
                   const CglImplication& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    probingInfo_=rhs.probingInfo_;
  }
  return *this;
}
// Create C++ lines to get to current state
std::string
CglImplication::generateCpp( FILE * fp) 
{
  CglImplication other;
  fprintf(fp,"0#include \"CglImplication.hpp\"\n");
  fprintf(fp,"3  CglImplication implication;\n");
  return "implication";
}
