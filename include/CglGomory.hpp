// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglGomory_H
#define CglGomory_H

#include <string>

#include "CglCutGenerator.hpp"

class CoinWarmStartBasis;
/** Gomory Cut Generator Class */
class CglGomory : public CglCutGenerator {
   friend void CglGomoryUnitTest(const OsiSolverInterface * siP,
				  const std::string mpdDir );
 
public:
    
  
  /**@name Generate Cuts */
  //@{
  /** Generate Mixed Integer Gomory cuts for the model of the 
      solver interface, si.

      Insert the generated cuts into OsiCut, cs.

      There is a limit option, which will only generate cuts with
      less than this number of entries.

      We can also only look at 0-1 variables a certain distance
      from integer.
  */
  virtual void generateCuts( const OsiSolverInterface & si, 
			     OsiCuts & cs) const;
  /** Generates cuts given matrix and solution etc,
      returns number of cuts generated */
  int generateCuts( const OsiRowCutDebugger * debugger, 
		    OsiCuts & cs,
		    const CoinPackedMatrix & columnCopy,
		    const double * objective, const double * colsol,
		    const double * colLower, const double * colUpper,
		    const double * rowLower, const double * rowUpper,
		    const char * intVar ,
		    const CoinWarmStartBasis* warm) const;

  //@}

  /**@name Change limit on how many variables in cut (default 50) */
  //@{
  /// Set
  void setLimit(int limit);
  /// Get
  int getLimit() const;
  //@}

  /**@name Change criterion on which variables to look at.  All ones
   more than "away" away from integrality will be investigated 
  (default 0.05) */
  //@{
  /// Set
  void setAway(double value);
  /// Get
  double getAway() const;
  //@}

  /**@name Change criterion on whether a cut is sufficiently good
   i.e. violation (default 0.00001) 
  At present this is not very useful as violation==away
  */
  //@{
  /// Set
  void setViolationTolerance(double value);
  /// Get
  double getViolationTolerance() const;
  //@}



  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglGomory ();
 
  /// Copy constructor 
  CglGomory (
    const CglGomory &);

  /// Assignment operator 
  CglGomory &
    operator=(
    const CglGomory& rhs);
  
  /// Destructor 
  virtual
    ~CglGomory ();
  //@}
      
private:
  
 // Private member methods

  // Private member data

  /**@name Private member data */
  //@{
  double violationTolerance_;
  /// Only investigate if more than this away from integrality
  double away_;
  /// Limit - only generate if fewer than this in cut
  int limit_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglGomory class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void CglGomoryUnitTest(const OsiSolverInterface * siP,
			const std::string mpdDir );
  
#endif
