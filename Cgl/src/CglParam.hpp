// Name:     CglParam.hpp
// Author:   Francois Margot
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
//           email: fmargot@andrew.cmu.edu
// Date:     11/24/06
//-----------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.

#ifndef CglParam_H
#define CglParam_H

/** Class collecting parameters for all cut generators. Each generator
    may have a derived class to add parameters. Each generator might
    also set different default values for the parameters in CglParam.  */

class CglParam {

public:

  /**@name Public Set/get methods */
  //@{

  /** Set INFINIT */
  virtual void setINFINIT(const double inf);
  /** Get value of INFINIT */
  inline double getINFINIT() const {return INFINIT;};

  /** Set EPS */
  virtual void setEPS(const double eps);
  /** Get value of EPS */
  inline double getEPS() const {return EPS;};

  /** Set EPS_COEFF */
  virtual void setEPS_COEFF(const double eps_c);
  /** Get value of EPS_COEFF */
  inline double getEPS_COEFF() const {return EPS_COEFF;};

  /** Set MAX_SUPPORT */
  virtual void setMAX_SUPPORT(const int max_s);
  /** Get value of MAX_SUPPORT */
  inline int getMAX_SUPPORT() const {return MAX_SUPPORT;};
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglParam(const double inf = DBL_MAX, const double eps = 0, 
	   const double eps_rr = 0, const int max_s = INT_MAX);
 
  /// Copy constructor 
  CglParam(const CglParam&);

  /// Clone
  virtual CglParam* clone() const;

  /// Assignment operator 
  virtual CglParam& operator=(const CglParam &rhs);

  /// Destructor 
  virtual ~CglParam();
  //@}

protected:

  // Protected member data

  /**@name Protected member data */

  //@{
  // Value for infinity. Default: DBL_MAX. 
  double INFINIT;

  // EPSILON for double comparisons. Default: 0.
  double EPS;

  // Returned cuts do not have coefficients with absolute value smaller 
  // than EPS_COEFF. Default: 0.
  double EPS_COEFF;

  /** Maximum number of non zero coefficients in a generated cut;
      Default: INT_MAX */ 
  int MAX_SUPPORT; 
  //@}

};

#endif
