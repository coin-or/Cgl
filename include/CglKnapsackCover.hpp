// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglKnapsackCover_H
#define CglKnapsackCover_H

#include <string>

#include "CglCutGenerator.hpp"

/** Knapsack Cover Cut Generator Class */
class CglKnapsackCover : public CglCutGenerator {
   friend void CglKnapsackCoverUnitTest(const OsiSolverInterface * siP,
					const std::string mpdDir );

public:
  /**@name Generate Cuts */
  //@{
  /** Generate knapsack cover cuts for the model of the solver interface, si. 
      Insert the generated cuts into OsiCut, cs.
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs) const;
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglKnapsackCover ();
 
  /// Copy constructor 
  CglKnapsackCover (
    const CglKnapsackCover &);

  /// Assignment operator 
  CglKnapsackCover &
    operator=(
    const CglKnapsackCover& rhs);
  
  /// Destructor 
  virtual
    ~CglKnapsackCover ();
  //@}


  /**@name Sets and gets
  //@{
  /// Set limit on number in knapsack
  inline void setMaxInKnapsack(int value)
           { if (value>0) maxInKnapsack_ = value;};
  /// get limit on number in knapsack
  inline int getMaxInKnapsack() const
           {return maxInKnapsack_;};
private:
  
 // Private member methods


  /**@name Private methods */
  //@{

  /** deriveAKnapsack 
                 returns 1 if it is able to derive
                 a (canonical) knapsack inequality
                in binary variables of the form ax<=b 
                 from the rowIndex-th  row in the model, 
                returns 0 otherwise.
  */
  int deriveAKnapsack(
    const OsiSolverInterface & si, 
    OsiCuts & cs,
    OsiPackedVector & krow,
    double & b,
    int *  complement,
    double *  xstar,
    int rowIndex,
    const OsiPackedVectorBase & matrixRow) const;

  /** Find a violated minimal cover from 
 a canonical form knapsack inequality by
 solving the -most- violated cover problem
 and postprocess to ensure minimality
  */
  int findExactMostViolatedMinCover(
      int nCols, 
      int row,
      OsiPackedVector & krow,
      double b, 
      double *  xstar, 
      OsiPackedVector & cover,
      OsiPackedVector & remainder) const ;

  /** Find the most violate minimum cover by solving the lp-relaxation of the
      most-violate-min-cover problem 
  */
  int findLPMostViolatedMinCover(
      int nCols,
      int row,
      OsiPackedVector & krow,
      double & b,
      double * xstar, 
      OsiPackedVector & cover,
      OsiPackedVector & remainder) const;  
  
/// find a minimum cover by a simple greedy approach
  int findGreedyCover(
      int row,
      OsiPackedVector & krow,
      double & b,
      double * xstar,
      OsiPackedVector & cover,
      OsiPackedVector & remainder
      ) const;

  /// lift the cover inequality
  int liftCoverCut(
     double & b,
     int nRowElem,
     OsiPackedVector & cover,
     OsiPackedVector & remainder,
     OsiPackedVector & cut ) const;
 
  /// sequence-independent lift and uncomplement and add the resulting cut to the cut set
  int liftAndUncomplementAndAdd(
     double rowub,
     OsiPackedVector & krow,
     double & b,
     int * complement,
     int row,
     OsiPackedVector & cover,
     OsiPackedVector & remainder,
     OsiCuts & cs) const;

  /// sequence-dependent lift, uncomplement and add the resulting cut to the cut set
void seqLiftAndUncomplementAndAdd(
      int nCols,
      double * xstar, 
      int * complement,
      int row,
      int nRowElem,
      double & b,
      OsiPackedVector & cover,      // need not be violated
      OsiPackedVector & remainder,
      OsiCuts & cs) const;

  /// sequence-dependent lift binary variables either up or down, uncomplement and add to the cut set
void liftUpDownAndUncomplementAndAdd(
         int nCols,
         double * xstar, 
         int * complement,
         int row,
         int nRowElem,
         double & b,

         // the following 3 packed vectors partition the krow:
         OsiPackedVector & fracCover, // vars have frac soln values in lp relaxation
                                       // and form cover with the vars atOne
         OsiPackedVector & atOne,     // vars have soln value of 1 in lp relaxation
                                       // and together with fracCover form minimal (?) cover. 
         OsiPackedVector & remainder,
         OsiCuts & cs ) const;

  /// find a cover using a variation of the logic found in OSL (w/o SOS)
  int findPseudoJohnAndEllisCover (
     int row,
     OsiPackedVector & krow,
     double & b,
     double * xstar,                     
     OsiPackedVector & cover,  
     OsiPackedVector & remainder) const;

  /// find a cover using the basic logic found in OSL (w/o SOS)
  int findJohnAndEllisCover (
     int row,
     OsiPackedVector & krow,
     double & b,
     double * xstar,                     
     OsiPackedVector & fracCover,  
     OsiPackedVector & atOnes,  
     OsiPackedVector & remainder) const;


  /** A C-style implementation of the Horowitz-Sahni exact solution 
   procedure for solving knapsack problem. 
   
   (ToDo: implement the more efficient dynamic programming approach)

   (Reference: Martello and Toth, Knapsack Problems, Wiley, 1990, p30.)
  */
  int exactSolveKnapsack(
      int n, 
      double c, 
      double const *pp, 
      double const *ww,
      double & z, 
      int * x) const;

  //@}

  // Private member data

  /**@name Private member data */
  //@{
  /// epsilon
  double epsilon_;  
  /// 1-epsilon
  double onetol_;  
  /// Maximum in knapsack
  int maxInKnapsack_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglKnapsackCover class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void CglKnapsackCoverUnitTest(const OsiSolverInterface * siP,
			      const std::string mpdDir );
  
#endif
