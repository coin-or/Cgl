// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglStored_H
#define CglStored_H

#include <string>

#include "CglCutGenerator.hpp"

class CoinWarmStartBasis;
// class CglTreeProbingInfo;
/*! \brief A container for storing cuts.

  CglStored is a container for holding previously generated row
  cuts that presents a standard CglCutGenerator interface. Calling
  #generateCuts returns all cuts that meet the specified violation (see
  #setRequiredViolation).

  In addition to a collection of cuts, a CglStored object can hold a small
  amount of miscellaneous information: an objective and solution and upper and
  lower bound arrays.

  More recently, there is an attempt to add a CglTreeProbingInfo object and
  clique-based implication generation. It seems likely that this is not fully
  functional.
*/
class CglStored : public CglCutGenerator {
 
public:
    
  
  /*! \name Cut generation */
  //@{
  /*! \brief Retrieve stored row cuts that satisfy the specified violation.

    Stored cuts that satisfy the specified violation are returned in the
    container \p cs. Violation is evaluated based on the solution held in the
    \p si, not the solution (if any) held by this object.

    There is code that will attempt to generate implication cuts if a
    CglTreeProbingInfo object is present. It seems likely this is not fully
    functional.
    
    The \p info parameter is not used.
  */
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
			    const CglTreeInfo info = CglTreeInfo()) ;
  /// Add a collection of row cuts to the container
  void addCut(const OsiCuts & cs) ;
  /// Add a row cut to the container
  void addCut(const OsiRowCut & cut) ;
  /// Add a row cut to the container
  void addCut(double lb, double ub, const CoinPackedVector &vector) ;
  /// Add a row cut to the container
  void addCut(double lb, double ub, int size,
  	      const int *colIndices, const double *elements) ;
  /// Return the number of row cuts in the container
  inline int sizeRowCuts() const
  { return (cuts_.sizeRowCuts()) ; }
  /// Return a reference to the row cuts in the container
  const OsiRowCut *rowCutPointer(int index) const
  { return (cuts_.rowCutPtr(index)) ; }
  //@}

  /*! \name Cut generation control */
  //@{
  /// Set the violation threshold
  inline void setRequiredViolation(double value)
  { requiredViolation_ = value ; }
  /// Get the violation threshold
  inline double getRequiredViolation() const
  { return (requiredViolation_) ; }
  //@}

  /*! \name Miscellaneous information

    Null is returned if the requested information has not been loaded into
    the object.

    If the column size is <= 0, the solution and bound vectors are not copied
    when the object is copied or cloned.
  */
  //@{
  /// Set column size
  inline void setColumnSize(int nCols)
  { numberColumns_ = nCols ; }
  /// Get column size
  inline int getColumnSize()
  { return (numberColumns_) ; }
  /*! \brief Save miscellaneous information.

    Intended for the objective value, primal solution, and lower and upper
    column bounds. The vectors are copied. You must set the number of columns
    (#setColumnSize or constructor) before calling this method.
  */
  void saveStuff(double bestObjective, const double *bestSolution,
		 const double *lower, const double *upper) ;
  /*! \brief Retrieve the stored objective

    The default value is undefined.
  */
  double bestObjective() const ;
  /// Retrieve the stored solution vector
  inline const double *bestSolution() const
  { return (bestSolution_) ; }
  /// Retrieve the stored lower bounds
  const double *tightLower() const
  { return (bounds_) ; } 
  /// Tight upper bounds
  const double *tightUpper() const
  { return (bounds_+numberColumns_) ; } 
# ifdef CLIQUE_ANALYSIS
  /*! \brief Add a TreeProbingInfo object
  
    The CglStored object takes ownership of CglTreeProbingInfo object.
  */
  inline void setProbingInfo(CglTreeProbingInfo * info)
  { probingInfo_ = info ; }
# endif
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglStored (int numberColumns=0);
 
  /// Copy constructor 
  CglStored (const CglStored & rhs);

  /*! \brief Constructor from file

    A handy way to make an arbitrary collection of cuts available for testing.
    Quietly fails if the file cannot be opened.

    Format? Undocumented, of course. What were you expecting?
  */
  CglStored (const char * fileName);

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
  
  /// Required violation to return a cut
  double requiredViolation_;
# ifdef CLIQUE_ANALYSIS
  /// Pointer to probing information
  CglTreeProbingInfo * probingInfo_;
# endif
  /*! \brief The cut collection
  
    Note that the implementation of CglStored assumes that this collection
    will contain only row cuts.
  */
  OsiCuts cuts_;
  /// Length of primal solution and column bound vectors
  int numberColumns_;
  /// Primal solution
  double * bestSolution_;
  /// Lower (#bounds_) and upper (#bounds_+#numberColumns_) column bounds
  double * bounds_;

};
#endif
