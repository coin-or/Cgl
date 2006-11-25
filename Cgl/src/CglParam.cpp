// Name:     CglParam.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     24/11/06
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

#include "CglParam.hpp"

/***********************************************************************/
void CglParam::setINFINIT(const double inf)
{
  if(inf > 0)
    INFINIT = inf;
} /* setINFINIT */

/***********************************************************************/
void CglParam::setEPS(const double eps)
{
  if(eps >= 0)
    EPS = eps;
} /* setEPS */

/***********************************************************************/
void CglParam::setEPS_COEFF(const double eps_c)
{
  if(eps_c >= 0)
    EPS_COEFF = eps_c;
} /* setEPS_COEFF */

/***********************************************************************/
void CglParam::setEPS_RELAX_ABS(const double eps_ra)
{
  if(eps_ra >= 0)
    EPS_RELAX_ABS = eps_ra;
} /* setEPS_RELAX_ABS */

/***********************************************************************/
void CglParam::setEPS_RELAX_REL(const double eps_rr)
{
  if(eps_rr >= 0)
    EPS_RELAX_REL = eps_rr;
} /* setEPS_RELAX_REL */

/***********************************************************************/
void CglParam::setMAX_SUPPORT(const int max_s)
{
  if(max_s > 0)
    MAX_SUPPORT = max_s;
} /* setMAX_SUPPORT */

/***********************************************************************/
CglParam::CglParam(const double inf, const double eps, const double eps_c, 
		   const double eps_ra, const double eps_rr, const int max_s) :
  INFINIT(inf),
  EPS(eps),
  EPS_COEFF(eps_c),
  EPS_RELAX_ABS(eps_ra),
  EPS_RELAX_REL(eps_rr),
  MAX_SUPPORT(max_s)
{}

/***********************************************************************/
CglParam::CglParam(const CglParam &source) :
  INFINIT(source.INFINIT),
  EPS(source.EPS),
  EPS_COEFF(source.EPS_COEFF),
  EPS_RELAX_ABS(source.EPS_RELAX_ABS),
  EPS_RELAX_REL(source.EPS_RELAX_REL),
  MAX_SUPPORT(source.MAX_SUPPORT)
{}

/***********************************************************************/
CglParam* CglParam::clone() const
{
  return new CglParam(*this);
}

/***********************************************************************/
CglParam& CglParam::operator=(const CglParam &rhs)
{
  if(this != &rhs) {
    INFINIT = rhs.INFINIT;
    EPS = rhs.EPS;
    EPS_COEFF = rhs.EPS_COEFF;
    EPS_RELAX_ABS = rhs.EPS_RELAX_ABS;
    EPS_RELAX_REL = rhs.EPS_RELAX_REL;
    MAX_SUPPORT = rhs.MAX_SUPPORT;
  }
  return *this;
}

/***********************************************************************/
CglParam::~CglParam()
{}
