/*
  Copyright (C) 2010, Lou Hafer, International Business Machines Corporation,
  and others.
  All Rights Reserved.
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

#if CGL_DEBUG > 0

/*
  This function checks that the column bounds passed in cut do not cut off
  the optimal solution known to the debugger.
*/
void checkBounds (const OsiSolverInterface &si, const OsiColCut &cut)
{
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  std::cout
    << "  checkBounds: checking bounds against optimal solution."
    << std::endl ;
  if (debugger) {
    const double *optimal = debugger->optimalSolution() ;
    int i ;
    int nIndex ;
    const double *values ;
    const int *index ;
    const CoinPackedVector &lbs = cut.lbs() ;
    values = lbs.getElements() ;
    nIndex = lbs.getNumElements() ;
    index = lbs.getIndices() ;
    for (i = 0 ; i < nIndex ; i++) {
      int iColumn = index[i] ;
      std::cout
        << "    " << si.getColName(iColumn) << " (" << iColumn
	<< ") lb = " << values[i] << ", opt = " << optimal[iColumn]
	<< "." << std::endl ;
      assert(values[i] <= optimal[iColumn]+1.0e-5) ;
    }
    const CoinPackedVector & ubs = cut.ubs() ;
    values = ubs.getElements() ;
    nIndex = ubs.getNumElements() ;
    index = ubs.getIndices() ;
    for (i = 0 ; i < nIndex ; i++) {
      int iColumn = index[i];
      std::cout
        << "    " << si.getColName(iColumn) << " (" << iColumn
	<< ") ub = " << values[i] << ", opt = " << optimal[iColumn]
	<< "." << std::endl ;
      assert(values[i] >= optimal[iColumn]-1.0e-5) ;
    }
  }
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
  Debugging utility to dump current solution, and optimal solution if
  available via the RowCutDebugger.
*/
void dumpSoln (const OsiSolverInterface &si)

{ const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
  const double *optimal = 0 ;

  if (debugger) optimal = debugger->optimalSolution() ;

  int i ;
  const double *solution = si.getColSolution() ;
  const double *lower = si.getColLower() ;
  const double *upper = si.getColUpper() ;
  const double *objective = si.getObjCoefficients() ;
  double objval1 = 0.0 ;
  double objval2 = 0.0 ;
  int nCols = si.getNumCols() ; 

  for (i = 0 ; i < nCols ; i++)
  {
    std::cout
      << "    " << si.getColName(i) << " (" << i << ") "
      << lower[i] << " <= " << solution[i] << " <= " << upper[i] ;
    if (optimal) {
      std::cout << "; opt = " << optimal[i] ;
    }
    std::cout << "." << std::endl ;
    objval1 += solution[i]*objective[i] ;
/*
  TODO: In the original code, generateCuts had this test, but in
	generateCutsAndModify, the test was weakened by a hardwired tolerance
	of 1e-5. The two tests should be the same (which is now the case). If
	we need to reintroduce a tolerance, it should not be hardwired.
	-- lh, 101203 --
*/
    if (optimal)
    { objval2 += optimal[i]*objective[i] ;
      assert(optimal[i] >= lower[i] && optimal[i] <= upper[i]) ; } }

  std::cout << "  current obj " << objval1 ;
  if (optimal)
  { std::cout << ", integer " << objval2 ; }
  std::cout << "." << std::endl ;
  
  return ; }

#endif  // CGL_DEBUG > 0
