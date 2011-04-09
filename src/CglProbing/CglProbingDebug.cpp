/*
  Copyright (C) 2010, Lou Hafer, International Business Machines Corporation,
  and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

/*
  This file contains debugging functions. None of these should be used unless
  CGL_DEBUG is defined to be > 0.
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <sstream>
#include <iomanip>

#include "CglProbing.hpp"
#include "CglProbingDebug.hpp"

#include "CoinPackedMatrix.hpp"

#if CGL_DEBUG > 0

namespace CglProbingDebug {

int nPath = 0 ;

/*
  Check for the presence of an active row cut debugger.
*/
int checkForRowCutDebugger (const OsiSolverInterface &si,
			    int paranoia, int verbosity)
{
  static bool printedWarning = false ;
  const OsiRowCutDebugger *debugger = 0 ;
  int retval = 0 ;
  if (paranoia > 0) {
    debugger = si.getRowCutDebugger() ;
    if (debugger) {
      if (verbosity >= 3)
	std::cout << "      On optimal path " << nPath << std::endl ;
      nPath++ ;
      retval = 2 ;
    } else {
      debugger = si.getRowCutDebuggerAlways() ;
      if (debugger) {
        if (verbosity >= 3)
	  std::cout << "      Not on optimal path." << std::endl ;
	retval = 1 ;
      } else {
        if (verbosity > 0 && !printedWarning) {
	  std::cout
	    << "      WARNING: row cut debugger not active!" << std::endl ;
	  printedWarning = true ;
	}
      }
    }
  }
  return (retval) ;
}

/*
  This function checks that the column bounds passed in cut do not cut off
  the optimal solution known to the debugger. It returns the number of
  violations.
*/
int checkBounds (const OsiSolverInterface &si, const OsiColCut &cut,
		 int verbosity)
{
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  int infeasCnt = 0 ;
  if (debugger) {
    if (verbosity >= 4)
      std::cout
	<< "  checkBounds: checking col cuts against optimal solution."
	<< std::endl ;
    const double *xopt = debugger->optimalSolution() ;
    const int *indices ;

    const CoinPackedVector &lbVec = cut.lbs() ;
    int numlbs = lbVec.getNumElements() ;
    const double *lbs = lbVec.getElements() ;
    indices = lbVec.getIndices() ;
    for (int k = 0 ; k < numlbs ; k++) {
      int j = indices[k] ;
      const double &lj = lbs[k] ;
      const double &xoptj = xopt[j] ;
      const double err = lj-xoptj ;
      if (err > 1.0e-5) infeasCnt++ ;
      if (verbosity >= 5 || (verbosity > 0 && err > 1.0e-5)) {
	std::cout << "    x<" << j << "> l = " << lj << ", x* = " << xoptj ;
	if (err > 1.0e-5) std::cout << "; infeas = " << err ;
	std::cout << "." << std::endl ;
      }
    }

    const CoinPackedVector & ubVec = cut.ubs() ;
    int numubs = ubVec.getNumElements() ;
    const double *ubs = ubVec.getElements() ;
    indices = ubVec.getIndices() ;
    for (int k = 0 ; k < numubs ; k++) {
      int j = indices[k] ;
      const double &uj = ubs[k] ;
      const double &xoptj = xopt[j] ;
      const double err = xoptj-uj ;
      if (err > 1.0e-5) infeasCnt++ ;
      if (verbosity >= 5 || (verbosity > 0 && err > 1.0e-5)) {
	std::cout << "    x<" << j << "> u = " << uj << ", x* = " << xoptj ;
	if (err > 1.0e-5) std::cout << "; infeas = " << err ;
	std::cout << "." << std::endl ;
      }
    }
  } else {
    if (verbosity >= 4)
      std::cout
	<< "  checkBounds: not on optimal path or no optimal solution."
	<< std::endl ;
  }
  if (infeasCnt && verbosity > 0)
    std::cout
      << "  checkBounds: " << infeasCnt
      << " bounds cut off the optimal solution." << std::endl ;
  return (infeasCnt) ;
}

/*
  Debugging utility to print row information. Far too many parameters, but
  doing it right requires too much code to replicate at multiple locations.
*/

void dump_row (int i, double rhsLow, double rhsHi, double lhsLow, double lhsHi,
	       const OsiSolverInterface *si,
	       bool unfixed, bool fixed, int indent,
	       int rowLen, const int *indices, const double *coeffs,
	       double tol, const double *vlbs, const double *vubs)

{ int ndx, j ;
  double lbj, ubj, aj ;

  int pfxLen = indent+8 ;	// indent+length("unfixed:")
  int lineLen,varLen ;

/*
  rhsLow <= lhsLow <= name (i) <= lhsHi <= rhsHi
*/
  std::cout
    << std::setw(indent) << " "
    << rhsLow << " <= " << lhsLow << " <= " ;
  if (si)
  { std::cout << si->getRowName(i) << " " ; }
  std::cout << "(" << i << ") <= " << lhsHi << " <= " << rhsHi ;
/*
  If rowLen is 0 but one of fixed or unfixed is true, warn the user (who
  perhaps wasn't expecting a zero-length row. But if we have 0, false, false,
  take it as a sign nothing more was wanted.
*/
  if (rowLen == 0)
  { if (fixed || unfixed) {
      std::cout
	<< std::setw(indent) << " "
	<< "Information requested for row " << i
	<< " but length is " << rowLen << "." << std::endl ;
    }
    return ;
  }
/*
  Scan the row and calculate the contribution of fixed variables.
*/
  double sumFixed = 0 ;
  int nFixed = 0 ;
  for (ndx = 0 ; ndx < rowLen ; ndx++)
  { j = indices[ndx] ;
    if (fabs(vubs[j]-vlbs[j]) < tol)
    { sumFixed += coeffs[ndx]*vubs[j] ;
      nFixed++ ; } }
  if (nFixed > 0)
  { std::cout
      << ";  " << rowLen << " vars, "
      << nFixed << " fixed, sum = " << sumFixed ; }
/*
  Details of unfixed variables, if requested.
*/
  if (unfixed)
  { std::cout << std::endl << std::setw(indent) << " " << "unfixed:" ;
    lineLen = pfxLen ;
    const double *soln = 0 ;
    if (si)
    { soln = si->getColSolution() ; }
    for (ndx = 0 ; ndx < rowLen ; ndx++)
    { j = indices[ndx] ;
      lbj = vlbs[j] ;
      ubj = vubs[j] ;
      aj = coeffs[ndx] ;

      if (fabs(ubj-lbj) >= tol)
      { std::ostringstream oneVar ;
	oneVar << " [" ;
	if (si)
	{ oneVar << si->getColName(j) << " " ; }
	oneVar << " (" << j << "), " << aj ;
	if (si)
	{ oneVar << ", " << soln[j] ; }
	oneVar << "]" ;
	varLen = static_cast<int>(oneVar.str().length()) ;
	if (lineLen+varLen > 80)
	{ std::cout << std::endl << std::setw(pfxLen) << " " ;
	  lineLen = pfxLen ; }
	std::cout << oneVar.str() ;
	lineLen += varLen ; } } }
/*
  Details of fixed variables, if requested.
*/
  if (fixed)
  { std::cout << std::endl << std::setw(indent) << " " << "  fixed:" ;
    lineLen = pfxLen ;
    for (ndx = 0 ; ndx < rowLen ; ndx++)
    { j = indices[ndx] ;
      lbj = vlbs[j] ;
      ubj = vubs[j] ;
      aj = coeffs[ndx] ;

      if (fabs(ubj-lbj) < tol)
      { std::ostringstream oneVar ;
	oneVar << " [" ;
	if (si)
	{ oneVar << si->getColName(j) << " " ; }
	oneVar << " (" << j << "), " << aj << ", " << lbj << "]" ;
	varLen = static_cast<int>(oneVar.str().length()) ;
	if (lineLen+varLen > 80)
	{ std::cout << std::endl << std::setw(pfxLen) << " " ;
	  lineLen = pfxLen ; }
	std::cout << oneVar.str() ;
	lineLen += varLen ; } } }

  std::cout << std::endl ;

  return ; }

/*
  Debugging utility to dump the current solution and the optimal solution (if
  available).
*/
void dumpSoln (const OsiSolverInterface &si)
{
  const OsiRowCutDebugger *debugger = si.getRowCutDebuggerAlways() ;
  const double *xopt = 0 ;

  if (debugger) xopt = debugger->optimalSolution() ;

  int n = si.getNumCols() ; 
  const double *x = si.getColSolution() ;
  const double *lbs = si.getColLower() ;
  const double *ubs = si.getColUpper() ;
  const double *c = si.getObjCoefficients() ;
  double z = 0.0 ;
  double zopt = 0.0 ;

  std::cout << "  " << "Solution:" << std::endl ;
  for (int j = 0 ; j < n ; j++) {
    std::cout
      << "      " << "x<" << j << ">: " 
      << lbs[j] << " <= " << x[j] << " <= " << ubs[j] ;
    if (xopt) {
      std::cout << "; x* = " << xopt[j] ;
    }
    z += x[j]*c[j] ;
    if (xopt) {
      zopt += xopt[j]*c[j] ;
      if (lbs[j]-xopt[j] > 1.0e-8)
        std::cout << "; lb violation " << lbs[j]-xopt[j] ;
      if (xopt[j]-ubs[j] > 1.0e-8)
        std::cout << "; ub violation " << xopt[j]-ubs[j] ;
    }
    std::cout << "." << std::endl ;
  }

  std::cout << "    z = " << z ;
  if (xopt) {
    std::cout << ", z* = " << zopt ;
  }
  std::cout << "." << std::endl ;
  
  return ;
}

}  // end namespace CglProbingDebug

#endif  // CGL_DEBUG > 0
