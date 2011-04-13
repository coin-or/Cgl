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

/*! \brief Check for an active row cut debugger

  If \p paranoia = 0, the call is a noop and returns 0. If \p paranoia > 0
  the return value is
    0: no debugger available
    1: debugger active, not on optimal path
    2: debugger active, on optimal path
  If verbosity >= 1, you'll get a message for case 0. If verbosity >= 3, you'll
  get a message for cases 1 and 2.
*/
int checkForRowCutDebugger(const OsiSolverInterface &si,
			   int paranoia, int verbosity) ;

/*! \brief Check column bounds against optimal solution

  Check the column bounds passed as parameters against the optimal solution
  known to the debugger. \p Verbosity > 0 allows messages when there's a
  violation. \p Verbosity = 4 prints a summary, 5 prints each check.
*/
int checkBounds(const OsiSolverInterface &si, const OsiColCut &cut,
		int verbosity) ;

/*! \brief Print a summary of the primal variables (bounds, value, solution)

  For each primal variable, print a line of the form
  \verbatim
  x(index) lb <= x <= ub ; opt = x*
  \endverbatim
  The optimal value is available only if a row cut debugger is active. When
  the optimal solution is available, the method checks that
  \verbatim
  lb-tol < x* < ub+tol
  \endverbatim
*/
void dumpSoln (const OsiSolverInterface &si, double tol) ;

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
