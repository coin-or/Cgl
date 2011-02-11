/*
  Copyright (C) 2010, Lou Hafer, International Business Machines Corporation,
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

#ifndef CglProbingDebug_H
#define CglProbingDebug_H

/*! \file CglProbingDebug.hpp
    \brief Definitions for debugging functions.
*/

/*
  CGL_DEBUG controls debugging effort. If defined to zero, you'll get asserts
  and nothing else. Values greater than zero enable debugging output, the
  methods declared below, and increasingly expensive checks.
*/
#define CGL_DEBUG 0

#if defined(CGL_DEBUG)
# ifdef NDEBUG
#   undef NDEBUG
# endif
# include <assert.h>
#endif


#if CGL_DEBUG > 0

#include "OsiSolverInterface.hpp"
#include "OsiRowCutDebugger.hpp"

namespace CglProbingDebug {

/*! \brief Check column bounds against optimal solution

  Check the column bounds passed as parameters against the optimal solution
  known to the debugger.
*/
void checkBounds(const OsiSolverInterface &si, const OsiColCut &cut) ;

/*! \brief Print a summary of the primal variables (bounds, value, solution)

  For each primal variable, print a line of the form
  \verbatim
  name(index) <lower bound> <= <value> <= <upper bound>; opt = <optimal value>
  \endverbatim
  The optimal value is available only if a row cut debugger is active.
*/
void dumpSoln (const OsiSolverInterface &si) ;

/*! \brief Print detailed information about a row

  Starts with the line
  \verbatim
  rhsLow <= lhsLow <= name (i) <= lhsHi <= rhsHi
  \end{verbatim}
  where lhsLow and lhsHi are the lower and upper bounds on the lhs calculated
  using column lower and upper bounds.
  If \p unfixed is true, coefficients of unfixed variables are listed with the
  format
  \verbatim
  [ name(index) coeff, val ]
  \end{verbatim}
  If \p fixed is true, coefficients of fixed variables are listed with the
  same format. If \p si is null, variable names and the current solution value
  are not printed.
*/
void dump_row (int i, double rhsLow, double rhsHi, double lhsLow, double lhsHi,
	       const OsiSolverInterface *si,
	       bool unfixed, bool fixed, int indent,
	       int rowLen, const int *indices, const double *coeffs,
	       double tol, const double *vlbs, const double *vubs) ;

}  // end namespace CglProbingDebug

#endif  // CGL_DEBUG > 0

#endif
