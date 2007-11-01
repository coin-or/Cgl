// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglStored_H
#define CglStored_H

#include <string>

#include "CglCutGenerator.hpp"

class CoinWarmStartBasis;
class CglTreeProbingInfo;
/** Stored Cut Generator Class */
class CglStored : public CglCutGenerator {
 
public:
    
  
  /**@name Generate Cuts */
  //@{
  /** Generate Mixed Integer Stored cuts for the model of the 
      solver interface, si.

      Insert the generated cuts into OsiCut, cs.

      This generator just looks at previously stored cuts
      and inserts any that are violated by enough
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo()) const;
  //@}

  /**@name Change criterion on whether to include cut.
   Violations of more than this will be added to current cut list
  (default 1.0e-5) */
  //@{
  /// Set
  inline void setRequiredViolation(double value)
  { requiredViolation_=value;}
  /// Get
  inline double getRequiredViolation() const
  { return requiredViolation_;}
  /// Takes over ownership of probing info
  inline void setProbingInfo(CglTreeProbingInfo * info)
  { probingInfo_ = info;}
  //@}

  /**@name Cut stuff */
  //@{
  /// Add cuts
  void addCut(const OsiCuts & cs);
  /// Add a row cut
  void addCut(const OsiRowCut & cut);
  /// Add a row cut from a packed vector
  void addCut(double lb, double ub, const CoinPackedVector & vector);
  /// Add a row cut from elements
  void addCut(double lb, double ub, int size, const int * colIndices, const double * elements);
  inline int sizeRowCuts() const
  { return cuts_.sizeRowCuts();}
  const OsiRowCut * rowCutPointer(int index) const
  { return cuts_.rowCutPtr(index);}
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglStored ();
 
  /// Copy constructor 
  CglStored (const CglStored & rhs);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglStored &
    operator=(const CglStored& rhs);
  
  /// Destructor 
  virtual
    ~CglStored ();
  //@}
      
protected:
  
 // Protected member methods

  // Protected member data

  /**@name Protected member data */
  //@{
  /// Only add if more than this requiredViolation
  double requiredViolation_;
  /// Pointer to probing information
  CglTreeProbingInfo * probingInfo_;
  /// Cuts
  mutable OsiCuts cuts_;
  //@}
};
#endif
