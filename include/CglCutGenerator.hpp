// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglCutGenerator_H
#define CglCutGenerator_H

#include "OsiCuts.hpp"
#include "OsiSolverInterface.hpp"


//-------------------------------------------------------------------
//
// Abstract base class for generating cuts.
//
//-------------------------------------------------------------------
///
/** Cut Generator Base Class

This is an abstract base class for generating cuts.  A specific cut 
generator will inherit from this class.
*/
class CglCutGenerator  {
  
public:
    
  /**@name Generate Cuts */
  //@{
  /** Generate cuts for the model data contained in si.
  The generated cuts are inserted into and returned in the
  collection of cuts cs.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs )const=0; 
  //@}

    
  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglCutGenerator (); 
 
  /// Copy constructor 
  CglCutGenerator ( const CglCutGenerator &);

  /// Assignment operator 
  CglCutGenerator & operator=(const CglCutGenerator& rhs);

  /// Destructor 
  virtual ~CglCutGenerator ();
  //@}
  
  // test this class
  //static void unitTest();
  
// private:
  
 // Presently there is no private member data
};


#endif
