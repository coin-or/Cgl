// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglImplication_H
#define CglImplication_H

#include <string>

#include "CglCutGenerator.hpp"


/// This just uses implication info   
class CglImplication : public CglCutGenerator {
 
public:

  /**@name Generate Cuts */
  //@{
  /** Generate cuts from implication table
  Insert generated cuts into the cut set cs.
  */
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info = CglTreeInfo()) const;
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglImplication ();
 
  /// Constructor with info
  CglImplication (CglTreeProbingInfo * info);
 
  /// Copy constructor 
  CglImplication (
    const CglImplication &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglImplication &
    operator=(
    const CglImplication& rhs);
  
  /// Destructor 
  virtual
    ~CglImplication ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  //@}
  /**@name Set implication */
  //@{
  /// Set implication
  inline void setProbingInfo(CglTreeProbingInfo * info)
  { probingInfo_=info;}
  //@}

private:
  /**@name Private member data */
  //@{
  /// Pointer to tree probing info
  CglTreeProbingInfo * probingInfo_;
  //@}
};
#endif
