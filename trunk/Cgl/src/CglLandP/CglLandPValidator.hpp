// Copyright (C) 2005, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     11/22/05
//---------------------------------------------------------------------------

#ifndef CglValidator_H
#define CglValidator_H
#include "OsiSolverInterface.hpp"
#include "CglParam.hpp"
#include <vector>

/** constants describing rejection codes*/
//[5] = {"Accepted", "violation too small", "small coefficient too small", "big dynamic","too dense"}



/** Class to validate or reject a cut */
class CglValidator
{
public:
    /** Reasons for rejecting a cut */
  enum RejectionsReasons{
    NoneAccepted=0 /**Cut was accepted*/,
    SmallViolation /** Violation of the cut is too small */,
    SmallCoefficient /** There is a small coefficient we can not get rid off.*/,
    BigDynamic /** Dynamic of coefficinet is too important. */,
    DenseCut/**cut is too dense */,
    EmptyCut/**After cleaning cut has become empty*/,
    DummyEnd/** dummy*/
  };

  /** Constructor with default values */
  CglValidator(double maxFillIn = 1.,
	       double maxRatio = 1e8, 
	       double minViolation = 0,
          bool scale = false);

  /** Clean an OsiCut */
  int cleanCut(OsiRowCut & aCut, const double * solCut,const OsiSolverInterface &si, const CglParam & par) const;
    /** Clean an OsiCut by another method */
  int cleanCut2(OsiRowCut & aCut, const double * solCut, const OsiSolverInterface &si, const CglParam & par) const;
  /** Call the cut cleaner */
  int operator()(OsiRowCut & aCut, const double * solCut,const OsiSolverInterface &si, const CglParam & par) const
  {return cleanCut2(aCut, solCut, si, par);}
  /** @name set functions */
  /** @{ */
  void setMaxFillIn(double value) { maxFillIn_ = value;}
  void setMaxRatio(double value) { maxRatio_ = value;}
  void setMinViolation(double value) {minViolation_ = value;}
  /** @} */
  /** @name get functions */
  /** @{ */
  double getMaxFillIn() {return maxFillIn_;}
  double getMaxRatio() { return maxRatio_;}
  double getMinViolation() {return minViolation_;}
  /** @} */

  const std::string& failureString(RejectionsReasons code) const {return rejections_[(int) code];} 
  const std::string& failureString(int code) const {return rejections_[ code];} 
  int numRejected(RejectionsReasons code)const{ return numRejected_[(int) code];}
  int numRejected(int code)const{ return numRejected_[ code];}
private:
  static void fillRejectionReasons();
  /** max percentage of given formulation fillIn should be accepted for cut fillin.*/
  double maxFillIn_;
  /** max ratio between smallest and biggest coefficient */
  double maxRatio_;
  /** minimum violation for accepting a cut */
  double minViolation_;
  /** Do we do scaling? */
  bool scale_;
  /** Strings explaining reason for rejections */
  static std::vector<std::string> rejections_;
  /** Number of cut rejected for each of the reasons.*/
  mutable std::vector<int> numRejected_;
};
#endif
