// Last edit: 3/15/07
//
// Name:     CglRedSplit.hpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
//           email: fmargot@andrew.cmu.edu
// Date:     2/6/05
//-----------------------------------------------------------------------------
// Copyright (C) 2005, Francois Margot and others.  All Rights Reserved.

/** Gomory Reduce-and-Split Cut Generator Class; See method generateCuts().
    Based on the paper by K. Anderson, G. Cornuejols, Yanjun Li, 
    "Reduce-and-Split Cuts: Improving the Performance of Mixed Integer 
    Gomory Cuts", Management Science 51 (2005). */

#ifndef CglRedSplit_H
#define CglRedSplit_H

#include "CglCutGenerator.hpp"

class CglRedSplit : public CglCutGenerator {
 
public:
  /**@name Generate Cuts */
  //@{
  /** Generate Reduce-and-Split Mixed Integer Gomory cuts 
      for the model of the solver interface si.

      Insert the generated cuts into OsiCuts cs.

      Warning: This generator currently works only with the Lp solvers Clp or 
      Cplex9.0 or higher. It requires access to the optimal tableau and 
      optimal basis inverse and makes assumptions on the way slack variables 
      are added by the solver. The Osi implementations for Clp and Cplex 
      verify these assumptions.

      When calling the generator, the solver interface si 
      must contain an optimized
      problem and information related to the optimal basis must be available 
      through the OsiSolverInterface methods (si->optimalBasisIsAvailable()
      must return 'true'). It is also essential that the integrality of
      structural variable i can be obtained using si->isInteger(i).

      Reduce-and-Split cuts are variants of Gomory cuts: Starting from
      the current optimal tableau, linear combinations of the rows of 
      the current optimal simplex tableau are used for generating Gomory
      cuts. The choice of the linear combinations is driven by the objective 
      of reducing the coefficients of the non basic continuous variables
      in the resulting row.
      Note that this generator might not be able to generate cuts for some 
      solutions violating integrality constraints. 

      Parameters of the generator are listed below. Modifying the default 
      values for parameters other than the last five might result in 
      invalid cuts.

      - LUB: Value considered large for the absolute value of a lower or upper
             bound on a variable. See method setLUB().
      - MAXDYN: Maximum ratio between largest and smallest non zero 
                coefficients in a cut. See method setMAXDYN().
      - MAXDYN_LUB: Maximum ratio between largest and smallest non zero 
                    coefficients in a cut involving structural variables with
		    lower or upper bound in absolute value larger than LUB.
                    Should logically be larger or equal to MAXDYN.
                    See method setMAXDYN_LUB().
      - EPS: Precision of double computations. See method setEPS().
      - EPS_ELIM: Precision for deciding if a coefficient is zero when 
                  eliminating slack variables. See method setEPS_ELIM().
      - EPS_COEFF: Precision for deciding if a coefficient of a generated cut 
                   is zero. See method setEPS_COEFF().
      - EPS_COEFF_LUB: Precision for deciding if a coefficient of a 
                       generated cut is zero when the corresponding 
                       variable has a lower or upper bound larger than 
                       LUB in absolute value. See method setEPS_COEFF_LUB().
      - EPS_RELAX: Value used to relax slightly the right hand side of each
                   generated cut. See method setEPS_RELAX().

      - MINVIOL: Minimum violation for the current basic solution in 
                 a generated cut. See method setMINVIOL().
      - USE_INTSLACKS: Use integer slacks to generate cuts. (not implemented)
      - USE_CG2: Use alternative formula to generate a mixed integer Gomory
                 cut (see methods generate_cgcut() and generate_cgcut2()).
      - normIsZero: Norm of a vector is considered zero if smaller than
                    this value. See method setNormIsZero().
      - minReduc: Reduction is performed only if the norm of the vector is
                  reduced by this fraction. See method setMinReduc().
      - limit: Generate cuts with at most this number of nonzero entries. 
               See method setLimit().
      - away: Look only at basic integer variables whose current value
              is at least this value from being integer. See method setAway().
      - maxTab: Controls the number of rows selected for the generation. See
                method setMaxTab().
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo());

  /// For compatibility with CglCutGenerator (const method)
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo()) const;

  /// Return true if needs optimal basis to do cuts (will return true)
  virtual bool needsOptimalBasis() const;
  //@}
  
  
  /**@name Parameters */
  //@{
  /** Set limit, the maximum number of non zero coefficients in generated cut;
      Default: 50 */
  void setLimit(int limit);
  /** Get value of limit */
  int getLimit() const;

  /** Set away, the minimum distance from being integer used for selecting 
      rows for cut generation;  all rows whose pivot variable should be 
      integer but is more than away from integrality will be selected; 
      Default: 0.05 */
  void setAway(double value);
  /// Get value of away
  double getAway() const;

  // Compute entries of low_is_lub and up_is_lub.
  void compute_is_lub();

 /** Set the value of LUB, value considered large for the absolute value of
      a lower or upper bound on a variable;
      Default: 1000 */
  void setLUB(double value);
  /** Get the value of LUB */
  double getLUB() const;

  // Set the maximum ratio between largest and smallest non zero 
  // coefficients in a cut. Default: 1e8.
  void setMAXDYN(double value);
  /** Get the value of MAXDYN */
  double getMAXDYN() const;

  // Set the maximum ratio between largest and smallest non zero 
  // coefficient in a cut involving structural variables with
  // lower or upper bound in absolute value larger than LUB.
  // Should logically be larger or equal to MAXDYN. Default: 1e13.
  void setMAXDYN_LUB(double value);
  /** Get the value of MAXDYN_LUB */
  double getMAXDYN_LUB() const;

  /** Set the value of EPS, epsilon for double computations;
      Default: 1e-7 */
  void setEPS(double value);
  /** Get the value of EPS */
  double getEPS() const;

  /** Set the value of EPS_ELIM, epsilon for values of coefficients when 
      eliminating slack variables;
      Default: 1e-12 */
  void setEPS_ELIM(double value);
  /** Get the value of EPS_ELIM */
  double getEPS_ELIM() const;

  /** Set the value of EPS_COEFF, epsilon for values of coefficients;
      Default: 1e-8 */
  void setEPS_COEFF(double value);
  /** Get the value of EPS_COEFF */
  double getEPS_COEFF() const;

  /** Set the value of EPS_COEFF_LUB, epsilon for values of coefficients for 
      variables with absolute value of lower or upper bound larger than LUB;
      Default: 1e-13 */
  void setEPS_COEFF_LUB(double value);
  /** Get the value of EPS_COEFF_LUB */
  double getEPS_COEFF_LUB() const;

  /** Set the value of EPS_RELAX, value used for relaxing the right hand side
      of each generated cut;
      Default: 1e-8 */
  void setEPS_RELAX(double value);
  /** Get the value of EPS_RELAX */
  double getEPS_RELAX() const;

  /** Set the value of MINVIOL, the minimum violation for the current 
      basic solution in a generated cut. Default: 1e-7 */
  void setMINVIOL(double value);
  /** Get the value of MINVIOL */
  double getMINVIOL() const;

  /** Set the value of USE_INTSLACKS. Default: 0 */
  void setUSE_INTSLACKS(int value);
  /** Get the value of USE_INTSLACKS */
  int getUSE_INTSLACKS() const;

  /** Set the value of USE_CG2. Default: 0 */
  void setUSE_CG2(int value);
  /** Get the value of USE_CG2 */
  int getUSE_CG2() const;

  /** Set the value of normIsZero, the threshold for considering a norm to be 
      0; Default: 1e-5 */
  void setNormIsZero(double value);
  /** Get the value of normIsZero */
  double getNormIsZero() const;

  /** Set the value of minReduc, threshold for relative norm improvement for
   performing  a reduction; Default: 0.05 */
  void setMinReduc(double value);
  /// Get the value of minReduc
  double getMinReduc() const;

  /** Set the maximum allowed value for (mTab * mTab * max(mTab, nTab)) where 
      mTab is the number of rows used in the combinations and nTab is the 
      number of continuous non basic variables. The work of the generator is 
      proportional to (mTab * mTab * max(mTab, nTab)). Reducing the value of 
      maxTab makes the generator faster, but weaker. Default: 1e7. */
  void setMaxTab(double value);
  /// Get the value of maxTab
  double getMaxTab() const;

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglRedSplit ();
 
  /// Copy constructor 
  CglRedSplit (const CglRedSplit &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglRedSplit &
    operator=(
    const CglRedSplit& rhs);
  
  /// Destructor 
  virtual
    ~CglRedSplit ();
  /// Create C++ lines to get to current state
  virtual std::string generateCpp( FILE * fp);
  //@}
    
  // Public member methods

/**@name Public member methods */

  //@{
  /** Set given_optsol to the given optimal solution given_sol.
      If given_optsol is set using this method, 
      the code will stop as soon as
      a generated cut is violated by the given solution; exclusively 
      for debugging purposes. */
  void set_given_optsol(const double *given_sol, const int card_sol);

  /** Print some of the data members. */
  void print() const;

  /** Print the current simplex tableau. */
  void printOptTab(OsiSolverInterface *solver) const;
  //@}

private:
  
  // Private member methods

/**@name Private member methods */

  //@{
  /// Compute the fractional part of value, allowing for small error.
  inline double rs_above_integer(double value); 

  /// Perform row r1 of pi := row r1 of pi - step * row r2 of pi.
  void update_pi_mat(int r1, int r2, int step);

  /// Perform row r1 of tab := row r1 of tab - step * row r2 of tab.
  void update_redTab(int r1, int r2, int step);

  /// Find optimal integer step for changing row r1 by adding to it a 
  /// multiple of another row r2.
  void find_step(int r1, int r2, int *step, 
		 double *reduc, double *norm);

  /// Test if an ordered pair of rows yields a reduction. Perform the
  /// reduction if it is acceptable.
  int test_pair(int r1, int r2, double *norm);

  /// Reduce rows of contNonBasicTab.
  void reduce_contNonBasicTab();

  /// Generate a row of the current LP tableau.
  void generate_row(int index_row, double *row);

  /// Generate a mixed integer Chvatal-Gomory cut, when all non basic 
  /// variables are non negative and at their lower bound.
  int generate_cgcut(double *row, double *rhs);

  /// Generate a mixed integer Chvatal-Gomory cut, when all non basic 
  /// variables are non negative and at their lower bound (different formula)
  int generate_cgcut_2(int basic_ind, double *row, double *rhs);

  /// Use multiples of the initial inequalities to cancel out the coefficients
  /// of the slack variables.
  void eliminate_slacks(double *row, 
			const double *elements, 
			const int *start,
			const int *indices,
			const int *rowLength,
			const double *rhs, double *rowrhs);

  /// Change the sign of the coefficients of the continuous non basic
  /// variables at their upper bound.
  void flip(double *row);

  /// Change the sign of the coefficients of the continuous non basic
  /// variables at their upper bound and do the translations restoring
  /// the original bounds. Modify the right hand side
  /// accordingly.
  void unflip(double *row, double *rowrhs, double *slack_val);

  /// Return the scale factor for the row. 
  /// Compute max_coeff: maximum absolute value of the coefficients.
  /// Compute min_coeff: minimum absolute value of the coefficients
  /// larger than EPS_COEFF.
  /// Return -1 if max_coeff < EPS_COEFF or if max_coeff/min_coeff > MAXDYN
  /// or MAXDYN_LUB (depending if the row has a non zero coeff. for a variable
  /// with large lower/upper bound) */.
  double row_scale_factor(double *row);

  /// Generate the packed cut from the row representation.
  int generate_packed_row(const double *xlp, double *row,
			  int *rowind, double *rowelem, 
			  int *card_row, double & rhs);

  /// Check that the generated cuts do not cut a given optimal solution.
  void check_optsol(const OsiSolverInterface *solver, 
		    const int calling_place,
		    const double *xlp, const double *slack_val,
		    const int do_flip);

  /// Check that the generated cuts do not cut a given optimal solution.
  void check_optsol(const OsiSolverInterface *solver, 
		    const int calling_place,
		    const double *xlp, const double *slack_val,
		    const double *ck_row, const double ck_rhs, 
		    const int cut_number, const int do_flip);

  //@}

  
  // Private member data

/**@name Private member data */

  //@{
  /// Number of rows ( = number of slack variables) in the current LP.
  int nrow; 

  /// Number of structural variables in the current LP.
  int ncol;

  /// Lower bounds for structural variables
  const double *colLower;

  /// Upper bounds for structural variables
  const double *colUpper;
  
  /// Lower bounds for constraints
  const double *rowLower;

  /// Upper bounds for constraints
  const double *rowUpper;

  /// Value considered large for the absolute value of lower or upper 
  /// bound on a variable. Default: 1000.
  double LUB;

  // Maximum ratio between largest and smallest non zero 
  // coefficients in a cut. Default: 1e8.
  double MAXDYN;

  // Maximum ratio between largest and smallest non zero 
  // coefficients in a cut involving structural variables with
  // lower or upper bound in absolute value larger than LUB.
  // Should logically be larger or equal to MAXDYN. Default: 1e13.
  double MAXDYN_LUB;

  /// Epsilon for precision. Default: 1e-7.
  double EPS;

  /// Epsilon for value of coefficients when eliminating slack variables. 
  /// Default: 1e-12.
  double EPS_ELIM;

  /// Epsilon for value of coefficients. Default: 1e-8.
  double EPS_COEFF;

  /// Epsilon for value of coefficients for variables with absolute value of
  /// lower or upper bound larger than LUB. Default: 1e-13.
  double EPS_COEFF_LUB;

  /// Epsilon for relaxing the right hand side of each generated constraint. 
  /// Default: 1e-8.
  double EPS_RELAX;

  /** Minimum violation for the current basic solution in a generated cut.
      Default: 1e-7 */
  double MINVIOL;

  /** Use integer slacks to generate cuts if USE_INTSLACKS = 1. Default: 0. */
  int USE_INTSLACKS;

  /** Use second way to generate a mixed integer Gomory cut 
      (see methods generate_cgcut()) and generate_cgcut_2()). Default: 0. */
  int USE_CG2;

  /// Norm of a vector is considered zero if smaller than normIsZero;
  /// Default: 1e-5.
  double normIsZero;

  /// Minimum reduction in percent that must be achieved by a potential 
  /// reduction step in order to be performed; Between 0 and 1, default: 0.05.
  double minReduc;

  /// Number of integer basic structural variables that are fractional in the
  /// current lp solution (at least away_ from being integer).  
  int card_intBasicVar_frac;

  /// Number of integer non basic structural variables in the
  /// current lp solution.  
  int card_intNonBasicVar; 

  /// Number of continuous non basic variables (structural or slack) in the
  /// current lp solution.  
  int card_contNonBasicVar;

  /// Number of non basic variables (structural or slack) at their
  /// upper bound in the current lp solution.
  int card_nonBasicAtUpper; 

  /// Number of non basic variables (structural or slack) at their
  /// lower bound in the current lp solution.
  int card_nonBasicAtLower;

  /// Characteristic vector for integer basic structural variables
  /// with non integer value in the current lp solution.
  int *cv_intBasicVar_frac;  

  /// List of integer structural basic variables 
  /// (in order of pivot in selected rows for cut generation).
  int *intBasicVar_frac;

  /// List of integer structural non basic variables.
  int *intNonBasicVar; 

  /// List of continuous non basic variables (structural or slack). 
  // slacks are considered continuous (no harm if this is not the case).
  int *contNonBasicVar;

  /// List of non basic variables (structural or slack) at their 
  /// upper bound. 
  int *nonBasicAtUpper;

  /// List of non basic variables (structural or slack) at their lower
  /// bound.
  int *nonBasicAtLower;

  /// Number of rows in the reduced tableau (= card_intBasicVar_frac).
  int mTab;

  /// Number of columns in the reduced tableau (= card_contNonBasicVar)
  int nTab;

  /// Tableau of multipliers used to alter the rows used in generation.
  /// Dimensions: mTab by mTab. Initially, pi_mat is the identity matrix.
  int **pi_mat;

  /// Current tableau for continuous non basic variables (structural or slack).
  /// Only rows used for generation.
  /// Dimensions: mTab by nTab.
  double **contNonBasicTab;

  /// Current tableau for integer non basic structural variables.
  /// Only rows used for generation.
  // Dimensions: mTab by card_intNonBasicVar.
  double **intNonBasicTab;

  /// Right hand side of the tableau.
  /// Only rows used for generation.
  double *rhsTab ;

  /// Use row only if pivot variable should be integer but is more 
  /// than away_ from being integer.
  double away_;
  /// Generate cut only if at most limit_ non zero coefficients in cut.
  int limit_;
  
  /// Maximum value for (mTab * mTab * max(mTab, nTab)). See method 
  /// setMaxTab().
  double maxTab_;

  /// Given optimal solution that should not be cut; only for debug. 
  double *given_optsol;

  /// Number of entries in given_optsol.
  int card_given_optsol;

  /// Characteristic vectors of structural integer variables or continuous
  /// variables currently fixed to integer values. 
  int *is_integer;

  /// Characteristic vector of the structural variables whose lower bound 
  /// in absolute value is larger than LUB. 
  int *low_is_lub;

  /// Characteristic vector of the structural variables whose upper bound 
  /// in absolute value is larger than LUB. 
  int *up_is_lub;

  //@}
};
  
#endif
