// Name:     CglRedSplitData.hpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
//           email: fmargot@andrew.cmu.edu
// Date:     12/2/06
//-----------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.

#ifndef CglRedSplitData_H
#define CglRedSplitData_H

#include "CglData.hpp"

  /**@name CglRedSplit Data */
  //@{

  /** Class collecting data used by the Reduced-and-split cut generator.
   */
  //@}

class CglRedSplitData : public CglData {

public:

  /**@name Set/get methods */
  //@{
  /** Set pointer on solver */
  virtual void setSolverPtr(OsiSolverInterface *givenSolverPtr);
  /// Get pointer on solver
  inline OsiSolverInterface* getSolverPtr() const {return solverPtr;};

  /** Set optimalBasisIsAvailable */
  virtual void setOptimalBasisIsAvailable(const int givenOptBas);
  /// Get optimalBasisIsAvailable
  inline int getOptimalBasisIsAvailable() const 
                                     {return optimalBasisIsAvailable;};
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglRedSplitData(const OsiSolverInterface *givenSolverPtr = NULL,
		  const int givenOptBas = 0);

  /// Constructor from CglData.
  CglRedSplitData(const CglData &givenData,
		  const OsiSolverInterface *givenSolverPtr = NULL,
		  const int givenOptBas = 0);

  /// Copy constructor 
  CglRedSplitData(const CglRedSplitData &source);

  /// Clone
  virtual CglRedSplitData* clone() const;

  /// Assignment operator 
  virtual CglRedSplitData& operator=(const CglRedSplitData &rhs);

  /// Destructor 
  virtual ~CglRedSplitData();
  //@}

protected:

  /**@name Parameters */
  //@{
  // Pointer on solver.
  OsiSolverInterface *solverPtr;

  // optimalBasisIsAvailable = 1 if *solverPtr is a pointer on a solver 
  // holding an optimized provle, optimalBasisIsAvailable = 0 otherwise.
  int optimalBasisIsAvailable;
  //@}

};

#endif
