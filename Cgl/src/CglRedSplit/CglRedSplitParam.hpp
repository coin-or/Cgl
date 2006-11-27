// Name:     CglRedSplitParam.hpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
//           email: fmargot@andrew.cmu.edu
// Date:     11/24/06
//-----------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.

#ifndef CglRedSplitParam_H
#define CglRedSplitParam_H

#include "CglParam.hpp"

/** Class collecting parameters the Reduced-and-split cut generator. */

class CglRedSplitParam : public CglParam {

public:

  /**@name CglRedSplit Parameters */
  /**
      Parameters of the generator are listed below. Modifying the default 
      values for parameters other than the last four might result in 
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
      - EPS_COEFF_LUB: Precision for deciding if a coefficient of a 
                       generated cut is zero when the corresponding 
                       variable has a lower or upper bound larger than 
                       LUB in absolute value. See method setEPS_COEFF_LUB().

      - normIsZero: Norm of a vector is considered zero if smaller than
                    this value. See method setNormIsZero().
      - minReduc: Reduction is performed only if the norm of the vector is
                  reduced by this fraction. See method setMinReduc().
      - away: Look only at basic integer variables whose current value
              is at least this value from being integer. See method setAway().
      - maxTab: Controls the number of rows selected for the generation. See
                method setMaxTab().
  */

  /**@name Set/get methods */
  //@{
  /** Set away, the minimum distance from being integer used for selecting 
      rows for cut generation;  all rows whose pivot variable should be 
      integer but is more than away from integrality will be selected; 
      Default: 0.05 */
  virtual void setAway(const double value);
  /// Get value of away
  inline double getAway() const {return away_;};

 /** Set the value of LUB, value considered large for the absolute value of
      a lower or upper bound on a variable;
      Default: 1000 */
  virtual void setLUB(const double value);
  /** Get the value of LUB */
  inline double getLUB() const {return LUB;};

  // Set the maximum ratio between largest and smallest non zero 
  // coefficients in a cut. Default: 1e8.
  virtual void setMAXDYN(double value);
  /** Get the value of MAXDYN */
  inline double getMAXDYN() const {return MAXDYN_LUB;};

  // Set the maximum ratio between largest and smallest non zero 
  // coefficient in a cut involving structural variables with
  // lower or upper bound in absolute value larger than LUB.
  // Should logically be larger or equal to MAXDYN. Default: 1e13.
  virtual void setMAXDYN_LUB(double value);
  /** Get the value of MAXDYN_LUB */
  inline double getMAXDYN_LUB() const {return MAXDYN_LUB;};

  /** Set the value of EPS_COEFF_LUB, epsilon for values of coefficients for 
      variables with absolute value of lower or upper bound larger than LUB;
      Default: 1e-13 */
  virtual void setEPS_COEFF_LUB(const double value);
  /** Get the value of EPS_COEFF_LUB */
  inline double getEPS_COEFF_LUB() const {return EPS_COEFF_LUB;};

  /** Set the value of normIsZero, the threshold for considering a norm to be 
      0; Default: 1e-5 */
  virtual void setNormIsZero(const double value);
  /** Get the value of normIsZero */
  inline double getNormIsZero() const {return normIsZero;};

  /** Set the value of minReduc, threshold for relative norm improvement for
   performing  a reduction; Default: 0.05 */
  virtual void setMinReduc(const double value);
  /// Get the value of minReduc
  inline double getMinReduc() const {return minReduc;};

  /** Set the maximum allowed value for (mTab * mTab * max(mTab, nTab)) where 
      mTab is the number of rows used in the combinations and nTab is the 
      number of continuous non basic variables. The work of the generator is 
      proportional to (mTab * mTab * max(mTab, nTab)). Reducing the value of 
      maxTab makes the generator faster, but weaker. Default: 1e7. */
  virtual void setMaxTab(const double value);
  /// Get the value of maxTab
  inline double getMaxTab() const {return maxTab_;};
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglRedSplitParam(const double inf = DBL_MAX, 
		   const double eps = 1e-7,
		   const double eps_coeff = 1e-8,
		   const double eps_relax_abs = 1e-8,
		   const double eps_relax_rel = 0.0,
		   const int max_support = 50,
		   const double lub = 1000.0,
		   const double max_dyn = 1e8,
		   const double max_dyn_lub = 1e13,
		   const double eps_coeff_lub = 1e-13,
		   const double norm_zero = 1e-5,
		   const double min_reduc = 0.05,
		   const double away = 0.05,
		   const double max_tab = 1e7);

   /// Constructor from CglParam
  CglRedSplitParam(const CglParam &source,
		   const double lub = 1000.0,
		   const double max_dyn = 1e8,
		   const double max_dyn_lub = 1e13,
		   const double eps_coeff_lub = 1e-13,
		   const double norm_zero = 1e-5,
		   const double min_reduc = 0.05,
		   const double away = 0.05,
		   const double max_tab = 1e7);

  /// Copy constructor 
  CglRedSplitParam(const CglRedSplitParam &source);

  /// Clone
  virtual CglRedSplitParam* clone() const;

  /// Assignment operator 
  virtual CglRedSplitParam& operator=(const CglRedSplitParam &rhs);

  /// Destructor 
  virtual ~CglRedSplitParam();
  //@}

protected:

  /**@name Parameters */
  //@{

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

  /// Epsilon for value of coefficients for variables with absolute value of
  /// lower or upper bound larger than LUB. Default: 1e-13.
  double EPS_COEFF_LUB;

  /// Norm of a vector is considered zero if smaller than normIsZero;
  /// Default: 1e-5.
  double normIsZero;

  /// Minimum reduction in percent that must be achieved by a potential 
  /// reduction step in order to be performed; Between 0 and 1, default: 0.05.
  double minReduc;

  /// Use row only if pivot variable should be integer but is more 
  /// than away_ from being integer.
  double away_;
  
  /// Maximum value for (mTab * mTab * max(mTab, nTab)). See method 
  /// setMaxTab().
  double maxTab_;

  //@}
};

#endif
