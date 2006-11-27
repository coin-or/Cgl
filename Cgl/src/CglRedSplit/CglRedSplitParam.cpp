// Name:     CglRedSplitParam.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     11/24/06
//---------------------------------------------------------------------------
// Copyright (C) 2006, Francois Margot and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CglRedSplitParam.hpp"

/***********************************************************************/
void CglRedSplitParam::setAway(const double value)
{
  if (value > 0.0 && value <= 0.5)
    away_ = value;
}

/***********************************************************************/
void CglRedSplitParam::setMaxTab(const double value)
{
  if (value > 10) {
    maxTab_ = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setMaxTab(): value: %f ignored\n", 
	   value);
  }
}

/***********************************************************************/
void CglRedSplitParam::setLUB(const double value)
{
  if (value > 0.0) {
    LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setLUB(): value: %f ignored\n", value);
  }
} /* setLUB */

/***********************************************************************/
void CglRedSplitParam::setMAXDYN(double value)
{
    if (value > 1.0) {
    MAXDYN = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setMAXDYN(): value: %f ignored\n", 
	   value);
  }
} /* setMAXDYN */

/***********************************************************************/
void CglRedSplitParam::setMAXDYN_LUB(double value)
{
  if (value > 1.0) {
    MAXDYN_LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setMAXDYN_LUB(): value: %f ignored\n", 
	   value);
  }
} /* setMAXDYN_LUB */

/***********************************************************************/
void CglRedSplitParam::setEPS_COEFF_LUB(const double value)
{
  if (value > 0.0 && value <= 0.1) {
    EPS_COEFF_LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setEPS_COEFF_LUB(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_COEFF_LUB */

/***********************************************************************/
void CglRedSplitParam::setNormIsZero(const double value)
{
  if (value > 0.0 && value <= 1) {
    normIsZero = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::setNormIsZero(): value: %f ignored\n",
	   value);
  }
} /* setNormIsZero */

/***********************************************************************/
void CglRedSplitParam::setMinReduc(const double value)
{
  if (value > 0.0 && value <= 1) {
    minReduc = value;
  }
  else {
    printf("### WARNING: CglRedSplitParam::MinReduc(): value: %f ignored\n",
	   value);
  }
} /* setMinReduc */

/***********************************************************************/
CglRedSplitParam::CglRedSplitParam(const double inf, 
				   const double eps,
				   const double eps_coeff,
				   const double eps_relax_abs,
				   const double eps_relax_rel,
				   const int max_support,
				   const double lub,
				   const double max_dyn,
				   const double max_dyn_lub,
				   const double eps_coeff_lub,
				   const double norm_zero,
				   const double min_reduc,
				   const double away,
				   const double max_tab) :
  CglParam(inf, eps, eps_coeff, eps_relax_abs, eps_relax_rel, max_support),
  LUB(lub),
  MAXDYN(1e8),
  MAXDYN_LUB(1e13),
  EPS_COEFF_LUB(eps_coeff_lub),
  normIsZero(norm_zero),
  minReduc(min_reduc),
  away_(away),
  maxTab_(max_tab)
{}

/***********************************************************************/
CglRedSplitParam::CglRedSplitParam(const CglParam &source,
				   const double lub,
				   const double max_dyn,
				   const double max_dyn_lub,
				   const double eps_coeff_lub,
				   const double norm_zero,
				   const double min_reduc,
				   const double away,
				   const double max_tab) :

  CglParam(source), 
	   //  CglParam(source.getINFINIT(), source.getEPS(), source.getEPS_COEFF(), 
	   //source.getEPS_RELAX_ABS(), source.getEPS_RELAX_REL(), 
	   //source.getMAX_SUPPORT()),
  LUB(100),
  MAXDYN(1e8),
  MAXDYN_LUB(1e13),
  EPS_COEFF_LUB(1e-13),
  normIsZero(1e-05),
  minReduc(0.05),
  away_(0.05),
  maxTab_(1e7)
{}

/***********************************************************************/
CglRedSplitParam::CglRedSplitParam(const CglRedSplitParam &source) :
  CglParam(source.INFINIT, source.EPS, source.EPS_COEFF, 
	   source.EPS_RELAX_ABS, source.EPS_RELAX_REL, 
	   source.MAX_SUPPORT),
  LUB(source.LUB),
  MAXDYN(source.MAXDYN),
  MAXDYN_LUB(source.MAXDYN_LUB),
  EPS_COEFF_LUB(source.EPS_COEFF_LUB),
  normIsZero(source.normIsZero),
  minReduc(source.minReduc),
  away_(source.away_),
  maxTab_(source.maxTab_)
{}

/***********************************************************************/
CglRedSplitParam* CglRedSplitParam::clone() const
{
  return new CglRedSplitParam(*this);
}

/***********************************************************************/
CglRedSplitParam& CglRedSplitParam::operator=(const CglRedSplitParam &rhs)
{
  if(this != &rhs) {
    INFINIT = rhs.INFINIT;
    EPS = rhs.EPS; 
    EPS_COEFF = rhs.EPS_COEFF; 
    EPS_RELAX_ABS = rhs.EPS_RELAX_ABS;
    EPS_RELAX_REL = rhs.EPS_RELAX_REL;
    MAX_SUPPORT = rhs.MAX_SUPPORT;

    LUB = rhs.LUB;
    MAXDYN = rhs.MAXDYN;
    MAXDYN_LUB = rhs.MAXDYN_LUB;
    EPS_COEFF_LUB = rhs.EPS_COEFF_LUB;
    normIsZero = rhs.normIsZero;
    minReduc = rhs.minReduc;
    away_ = rhs.away_;
    maxTab_ = rhs.maxTab_;
  }
  return *this;
}

/***********************************************************************/
CglRedSplitParam::~CglRedSplitParam()
{}
