// Last edit: 11/4/06
//
// Name:     CglRedSplit.cpp
// Author:   Francois Margot                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     2/6/05
//---------------------------------------------------------------------------
// Copyright (C) 2005, Francois Margot and others.  All Rights Reserved.

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

//#define RS_TRACE
//#define RS_TRACEALL
//#define RS_TRACETAB

#include "OsiSolverInterface.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglRedSplit.hpp"
#include "CoinFinite.hpp"

//-------------------------------------------------------------------
// Generate Gomory Reduce-and-Split cuts
//------------------------------------------------------------------- 

/***************************************************************************/
// Returns (value - floor) but allowing for small errors
inline double CglRedSplit::rs_above_integer(double value) 
{
  double value2=floor(value);
  double value3=floor(value+0.5);
  if (fabs(value3-value)< EPS * (fabs(value3)+1.0))
    return 0.0;
  return value-value2;
} /* rs_above_integer */

/**********************************************************/
void rs_allocmatINT(int ***v, const int m, const int n)
{
  int i;

  *v = (int **) calloc (m, sizeof(int *));
  if (*v == NULL) {
    printf("###ERROR: INTEGER matrix allocation failed\n");
    exit(1);
  }

  for(i=0; i<m; i++) {
    (*v)[i] = (int *) calloc (n, sizeof(int));
    if ((*v)[i] == NULL) {
      printf("###ERROR: INTEGER matrix allocation failed\n");
      exit(1);
    }
  }
} /* rs_allocmatINT */

/**********************************************************/
void rs_deallocmatINT(int ***v, const int m, const int n)
{
  int i;

  for(i=0; i<m; i++) {
    free((void *) (*v)[i]);
  }
  free((void *) (*v));
} /* rs_deallocmatINT */

/**********************************************************/
void rs_allocmatDBL(double ***v, const int m, const int n)
{
  int i;

  *v = (double **) calloc (m, sizeof(double *));
  if (*v == NULL) {
    printf("###ERROR: DOUBLE matrix allocation failed\n");
    exit(1);
  }

  for(i=0; i<m; i++) {
    (*v)[i] = (double *) calloc (n, sizeof(double));
    if ((*v)[i] == NULL) {
      printf("###ERROR: DOUBLE matrix allocation failed\n");
      exit(1);
    }
  }
} /* rs_allocmatDBL */

/**********************************************************/
void rs_deallocmatDBL(double ***v, const int m, const int n)
{
  int i;

  for(i=0; i<m; i++) {
    free((void *) (*v)[i]);
  }
  free((void *) (*v));
} /* rs_deallocmatDBL */

/**********************************************************/
void rs_printvecINT(char const *vecstr, const int *x, const int n)
{
  int num, fromto, upto, j, i;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10*j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      printf(" %4d", x[i]);
    printf("\n");
  }
  printf("\n");
} /* rs_printvecINT */

/**********************************************************/
void rs_printvecDBL(char const *vecstr, const double *x, const int n)
{
  int num, fromto, upto, j, i;

  num = (n/10) + 1;
  printf("%s :\n", vecstr);
  for(j=0; j<num; j++){
    fromto = 10*j;
    upto = 10 * (j+1);
    if(n <= upto) upto = n;
    for(i=fromto; i<upto; i++)
      printf(" %7.3f", x[i]);
    printf("\n");
  }
  printf("\n");
} /* rs_printvecDBL */

/**********************************************************/
void rs_printmatINT(char const *vecstr, const int **x, 
		    const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4d", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* rs_printmatINT */

/**********************************************************/
void rs_printmatINT(char const *vecstr, int **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %4d", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* rs_printmatINT */

/**********************************************************/
void rs_printmatDBL(char const *vecstr, double **x, const int m, const int n)
{
  int i, j;

  printf("%s :\n", vecstr);

  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf(" %7.3f", x[i][j]);
    }
      printf("\n");
  }
  printf("\n");
} /* rs_printmatDBL */

/***************************************************************************/
double rs_dotProd(const double *u, const double *v, const int dim) {

  double result = 0;
  for(int i=0; i<dim; i++) {
    result += u[i] * v[i];
  }
  return(result);
} /* rs_dotProd */

/***************************************************************************/
double rs_dotProd(const int *u, const double *v, const int dim) {

  double result = 0;
  for(int i=0; i<dim; i++) {
    result += u[i] * v[i];
  }
  return(result);
} /* rs_dotProd */

/***************************************************************************/
double rs_genalea (int *x0)
{
  int m = 2147483647;
  int a = 16807 ;
  int b = 127773 ;
  int c = 2836 ;
  int x1, k;

  k = (int) ((*x0)/b) ;
  x1 = a*(*x0 - k*b) - k*c ;
  if(x1 < 0) x1 = x1 + m;
  *x0 = x1;

  return((double)x1/(double)m);
} /* rs_genalea */

/***********************************************************************/
void CglRedSplit::update_pi_mat(int r1, int r2, int step) {

  for(int j=0; j<mTab; j++) {
    pi_mat[r1][j] = pi_mat[r1][j] - step * pi_mat[r2][j];
  }
} /* update_pi_mat */

/***********************************************************************/
void CglRedSplit::update_redTab(int r1, int r2, int step) {

  for(int j=0; j<nTab; j++) {
    contNonBasicTab[r1][j] = 
                  contNonBasicTab[r1][j] - step * contNonBasicTab[r2][j];
  }
} /* update_redTab */

/***********************************************************************/
void CglRedSplit::find_step(int r1, int r2, int *step, 
			    double *reduc, double *norm) {

   double btb_val = rs_dotProd(contNonBasicTab[r1], contNonBasicTab[r2], nTab);
   double opt_step = btb_val/norm[r2];

   int f_step= (int) floor(opt_step);
   int c_step = f_step + 1;

   double val_f = norm[r1] + f_step * f_step * norm[r2] - 2 * btb_val * f_step;
   double val_c = norm[r1] + c_step * c_step * norm[r2] - 2 * btb_val * c_step;

   if( val_f <= val_c ) {
     (*step) = f_step;
     (*reduc) = norm[r1] - val_f;
   }
   else {
     (*step) = c_step;
     (*reduc) = norm[r1] - val_c;
   }
} /* find_step */

/***************************************************************************/
int CglRedSplit::test_pair(int r1, int r2, double *norm) {

   int step;
   double reduc;
   
   find_step(r1, r2, &step, &reduc, norm); 

   if(reduc/norm[r1] >= minReduc ) {
     update_pi_mat(r1, r2, step);
     update_redTab(r1, r2, step);
     norm[r1] = rs_dotProd(contNonBasicTab[r1], contNonBasicTab[r1], nTab);

#ifdef RS_TRACEALL
     printf("Use %d and %d for reduction (step: %d)\n", r1, r2, step);
#endif

     return(1);
   }

   return(0);
} /* test_pair */

/***************************************************************************/
void CglRedSplit::reduce_contNonBasicTab() {

  int i, j;

  double *norm = new double[mTab];
  for(i=0; i<mTab; i++) {
    norm[i] = rs_dotProd(contNonBasicTab[i], contNonBasicTab[i], nTab);
  }

  double sum_norms = 0;
  for(i=0; i<mTab; i++) {
    sum_norms += norm[i];
  }

#ifdef RS_TRACE
  printf("CglRedSplit::reduce_contNonBasicTab():Initial sum of  norms: %f\n", 
	 sum_norms);
#endif

  int iter = 0, done = 0;
  int *changed = new int[mTab]; // changed[i]: last iter where row i updated
  int **checked; // checked[i][j]: last iter where pair (i, j) checked
  
  rs_allocmatINT(&checked, mTab, mTab);
  for(i=0; i<mTab; i++) {
    changed[i] = 0;
    for(j=0; j<mTab; j++) {
      checked[i][j] = -1;
    }
    checked[i][i] = 0;
  }

  while(!done) {
    done = 1;

#ifdef RS_TRACEALL
    rs_printmatDBL("contNonBasicTab", contNonBasicTab, mTab, nTab);
    rs_printmatINT("checked", checked, mTab, mTab);
    rs_printvecINT("changed", changed, mTab);
    rs_printvecDBL("norm", norm, mTab);
    rs_printmatINT("pi_mat", pi_mat, mTab, mTab);
#endif

    for(i=0; i<mTab; i++) {
      if(norm[i] > normIsZero) {
	for(j=i+1; j<mTab; j++) {
	  if(norm[j] > normIsZero) {
	    if((checked[i][j] < changed[i]) || (checked[i][j] < changed[j])) {
	      if(test_pair(i, j, norm)) {
		changed[i] = iter+1;
		done = 0;
	      }
	      checked[i][j] = iter;

	      if((checked[j][i] < changed[i]) || 
		 (checked[j][i] < changed[j])) {
		if(test_pair(j, i, norm)) {
		  changed[j] = iter+1;
		  done = 0;
		}
		checked[j][i] = iter;
	      }
	    }
	  }
	}
      }
    }
    iter++;
  }

#ifdef RS_TRACEALL
  rs_printmatDBL("contNonBasicTab", contNonBasicTab, mTab, nTab);
  rs_printmatINT("checked", checked, mTab, mTab);
  rs_printvecINT("changed", changed, mTab);
  rs_printvecDBL("norm", norm, mTab);
  rs_printmatINT("pi_mat", pi_mat, mTab, mTab);
#endif

  sum_norms = 0;
  for(i=0; i<mTab; i++) {
    sum_norms += norm[i];
  }

#ifdef RS_TRACE
  printf("CglRedSplit::reduce_contNonBasicTab():Final sum of norms: %f\n", sum_norms);
#endif

  delete[] norm;
  delete[] changed;
  rs_deallocmatINT(&checked, mTab, mTab);

} /* reduce_contNonBasicTab */

/************************************************************************/
void CglRedSplit::generate_row(int index_row, double *row) {

  int i;
  for( i=0; i<ncol+nrow; i++) {
    row[i] = 0;
  }
  for( i=0; i<card_intBasicVar_frac; i++) {
    row[intBasicVar_frac[i]] += pi_mat[index_row][i];
  }
  for( i=0; i<card_intNonBasicVar; i++) {
    int locind = intNonBasicVar[i];
    row[locind] = 0;
    for(int j =0; j<mTab; j++) {
      row[locind] += pi_mat[index_row][j] * intNonBasicTab[j][i];
    }
  }
  for( i=0; i<card_contNonBasicVar; i++) {
    row[contNonBasicVar[i]] = contNonBasicTab[index_row][i];
  }
} /* generate_row */

/************************************************************************/
int CglRedSplit::generate_cgcut(double *row, double *rhs) {
  
  double f0 = rs_above_integer(*rhs);
  double f0compl = 1 - f0;

#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::generate_cgcut(): starting row", 
		 row, ncol+nrow);
  printf("rhs: %f\n", *rhs);
#endif

  // See Wolsey "Integer Programming" (1998), p. 130, second line of proof of 
  // Proposition 8.8

  if((f0 < away_) || (f0compl < away_)) {
    return(0);
  }

  int i;
  for( i=0; i<card_intNonBasicVar; i++) {
    int locind = intNonBasicVar[i];
    double f = rs_above_integer(row[locind]);
    row[locind] -= f;
    if(f > f0) {
      row[locind] += (f-f0)/f0compl;
    }
  }

  for( i=0; i<card_contNonBasicVar; i++) {
    if(row[contNonBasicVar[i]] < 0) {
      row[contNonBasicVar[i]] /= f0compl;
    }
    else {
      row[contNonBasicVar[i]] = 0;
    }
  }
  (*rhs) -= f0;
  
#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::generate_cgcut: row", row, ncol+nrow);
  printf("rhs: %f\n", *rhs);
#endif

  return(1);
} /* generate_cgcut */

/************************************************************************/
void CglRedSplit::eliminate_slacks(double *row, 
				   const double *elements, 
				   const int *rowStart,
				   const int *indices,
				   const int *rowLength,
				   const double *rhs, double *rowrhs) {


  for(int i=0; i<nrow; i++) {
    if(fabs(row[ncol+i]) > EPS) {
      int upto = rowStart[i] + rowLength[i];
      for(int j=rowStart[i]; j<upto; j++) {
	row[indices[j]] -= row[ncol+i] * elements[j];
      }
      *rowrhs -= row[ncol+i] * rhs[i];
    }
  }

#ifdef RS_TRACEALL
  rs_printvecDBL("CglRedSplit::eliminate_slacks: row", row, ncol+nrow);
  printf("rhs: %f\n", *rowrhs);
#endif

} /* eliminate_slacks */

/************************************************************************/
void CglRedSplit::flip(double *row) {
  
  for(int i=0; i<card_nonBasicAtUpper; i++) {
    row[nonBasicAtUpper[i]] = -row[nonBasicAtUpper[i]];
  }
} /* flip */

/************************************************************************/
void CglRedSplit::unflip(double *row, double *rowrhs,
			 const double *colLower, const double *colUpper,
			 double *slack_val) {
  
  int i;
  for( i=0; i<card_nonBasicAtLower; i++) {
    int locind = nonBasicAtLower[i];
    if(locind < ncol) {
      *rowrhs += row[locind] * colLower[locind];
    }
    else {
      *rowrhs += row[locind] * slack_val[locind-ncol];
    }
  }
  for( i=0; i<card_nonBasicAtUpper; i++) {
    int locind = nonBasicAtUpper[i];
    row[locind] = -row[locind];
    if(locind < ncol) {
      *rowrhs += row[locind] * colUpper[locind];
    }
    else {
      *rowrhs += row[locind] * slack_val[locind-ncol];
    }
  }

#ifdef RS_TRACEALL
  rs_printvecDBL("After unflip: row", row, ncol+nrow);
  printf("rhs: %f\n", *rowrhs);
#endif

} /* unflip */

/************************************************************************/
int CglRedSplit::generate_packed_row( const OsiSolverInterface *solver,
                                             double *row,
                                             int *rowind, double *rowelem, 
                                             int *card_row, double & rhs) {
  *card_row = 0;
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();
  for(int i=0; i<ncol; i++) {
    double value = row[i];
    if(fabs(value) > EPS_COEFF) {
      rowind[*card_row] = i;
      rowelem[*card_row] = value;
      (*card_row)++;
      if(*card_row > limit_) {
#ifdef RS_TRACE
	printf("Cut discarded since too many non zero coefficients\n");
#endif	
	return(0);
      }
    } else {
      // small - adjust rhs if rhs reasonable
      if (value>0.0&&colLower[i] > -LUB) {
        rhs -= value*colLower[i];
      } else if (value<0.0&&colUpper[i] < LUB) {
        rhs -= value*colUpper[i];
      } else if (fabs(value) > EPS_COEFF_LUB) {
        // take anyway
        rowind[*card_row] = i;
        rowelem[*card_row] = value;
        (*card_row)++;
        if(*card_row > limit_) {
#ifdef RS_TRACE
          printf("Cut discarded since too many non zero coefficients\n");
#endif	
          return(0);
        }
      }
    }
  }
  return(1);
} /* generate_packed_row */

/************************************************************************/
void CglRedSplit::set_given_optsol(const double *given_sol, int card_sol) {
  if(given_optsol != NULL) {
    delete[] given_optsol;
    given_optsol = NULL;
  }
  given_optsol =  new double[card_sol];
  card_given_optsol = card_sol;
  for(int i=0; i<card_sol; i++) {
    given_optsol[i] = given_sol[i];
  }
} /* set_given_optsol */

/************************************************************************/
void CglRedSplit::check_optsol(const OsiSolverInterface *solver, 
			       const int calling_place,
			       const double *xlp, const double *slack_val,
			       const int do_flip) {

  if(card_given_optsol != ncol) {
    printf("### ERROR: CglRedSplit(): card_given_optsol: %d  ncol: %d\n", 
	   card_given_optsol, ncol);
    exit(1);
  }

  int i;
  const CoinPackedMatrix *ck_byRow = solver->getMatrixByRow();
  const double *rhs = solver->getRightHandSide();

  double *ck_slack = new double[nrow];

#ifdef RS_TRACEALL
  print();
#endif

  ck_byRow->timesMinor(given_optsol, ck_slack);
  int irow;
  for( irow=0; irow<nrow; irow++) {
    ck_slack[irow] = rhs[irow] - ck_slack[irow];  
                                        // slack values for optimal solution
  }
  
  double *ck_row = new double[ncol+nrow];
  
  for( irow=0; irow<mTab; irow++) {
    for(i=0; i<ncol+nrow; i++) {
      ck_row[i] = 0;
    }
    for(i=0; i<card_intBasicVar_frac; i++) {
      ck_row[intBasicVar_frac[i]] = pi_mat[irow][i];
    }
    for(i=0; i<card_intNonBasicVar; i++) {
      ck_row[intNonBasicVar[i]] = 0;
      for(int j=0; j<mTab; j++) {
	ck_row[intNonBasicVar[i]] += pi_mat[irow][j] * intNonBasicTab[j][i];
      }
    }
    for(i=0; i<card_contNonBasicVar; i++) {
      ck_row[contNonBasicVar[i]] = contNonBasicTab[irow][i];
    }

    if(do_flip) {
      for(i=0; i<card_nonBasicAtUpper; i++) {
	int locind = nonBasicAtUpper[i];
	ck_row[locind] = -ck_row[locind];
      }
    }

#ifdef RS_TRACEALL
    rs_printvecDBL("ck_row", ck_row, ncol);
    rs_printvecDBL("given_optsol", given_optsol, ncol);
    rs_printvecDBL("ck_row(slacks)", &(ck_row[ncol]), nrow);
    rs_printvecDBL("ck_slack", ck_slack, nrow);
#endif

    double ck_lhs = rs_dotProd(ck_row, given_optsol, ncol);
    ck_lhs += rs_dotProd(&(ck_row[ncol]), ck_slack, nrow);
    
    double ck_rhs = rs_dotProd(ck_row, xlp, ncol);
    ck_rhs += rs_dotProd(&(ck_row[ncol]), slack_val, nrow);
    
    if((ck_lhs < ck_rhs - EPS) || (ck_lhs > ck_rhs + EPS)) {
      printf("### ERROR: CglRedSplit::check_optsol(): Cut %d cuts given_optsol\n", 
	     irow);
      rs_printvecDBL("ck_row", ck_row, ncol+nrow);
      printf("lhs: %f  rhs: %f    calling_place: %d\n", 
	     ck_lhs, ck_rhs, calling_place);
      exit(1);
    }
  }
  delete[] ck_slack;
  delete[] ck_row;
} /* check_optsol */

/************************************************************************/
void CglRedSplit::check_optsol(const OsiSolverInterface *solver, 
			       const int calling_place,
			       const double *ck_row, const double ck_rhs,
			       const int cut_number, const int do_flip) {

  if(card_given_optsol != ncol) {
    printf("### ERROR: CglRedSplit(): card_given_optsol: %d  ncol: %d\n", 
	   card_given_optsol, ncol);
    exit(1);
  }

  const CoinPackedMatrix *ck_byRow = solver->getMatrixByRow();
  const double *rhs = solver->getRightHandSide();
  double *cpy_row = new double[ncol+nrow];
  double *ck_slack = new double[nrow];

#ifdef RS_TRACEALL
  print();
#endif

  for(int i=0; i<ncol+nrow; i++) {
    cpy_row[i] = ck_row[i];
  }

  ck_byRow->timesMinor(given_optsol, ck_slack);
  for(int irow=0; irow<nrow; irow++) {
    ck_slack[irow] = rhs[irow] - ck_slack[irow];  
                                       // slack values for optimal solution
  }
  
  if(do_flip) {
    for(int i=0; i<card_nonBasicAtUpper; i++) {
      int locind = nonBasicAtUpper[i];
      cpy_row[locind] = -cpy_row[locind];
    }
  }  printf("Use const method\n");


  double ck_lhs = rs_dotProd(cpy_row, given_optsol, ncol);
  ck_lhs += rs_dotProd(&(cpy_row[ncol]), ck_slack, nrow);
    
  if(ck_lhs > ck_rhs + EPS) {
    printf("### ERROR: CglRedSplit::check_optsol(): Cut %d cuts given_optsol\n", 
	   cut_number);
    rs_printvecDBL("cpy_row", cpy_row, ncol+nrow);
    printf("lhs: %f  rhs: %f    calling_place: %d\n", 
	   ck_lhs, ck_rhs, calling_place);
    exit(1);
  }
  delete[] cpy_row;
  delete[] ck_slack;
} /* check_optsol */

/************************************************************************/
void CglRedSplit::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			       const CglTreeInfo info) const
{
  if(nrow) {
    // should never get here
    printf("CglRedSplit::generateCuts(): bad status i.e. nrow > 0\n");
    exit(1);
  } else {

    // kludge to be able to modify the CglRedSplit object if it is const

    CglRedSplit temp(*this);
    temp.generateCuts(si, cs, info);
  }

} /* generateCuts */

/************************************************************************/
void CglRedSplit::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			       const CglTreeInfo info)
{
  if (nrow) {
    // should never get here
    printf("CglRedSplit::generateCuts(): bad status i.e. nrow > 0\n");
    exit(1);
  }
  OsiSolverInterface *solver = const_cast<OsiSolverInterface *>(&si);
  if(!solver->optimalBasisIsAvailable()) {
    printf("### WARNING: CglRedSplit::generateCuts(): no optimal basis available.\n");
    return;
  }
  solver->enableFactorization();
  
  int i;

  // Get basic problem information
  ncol = solver->getNumCols(); 
  nrow = solver->getNumRows(); 

  const double *xlp = solver->getColSolution();
  const double *row_act = solver->getRowActivity();
  const double *rhs = solver->getRightHandSide();
  const double *colLower = solver->getColLower();
  const double *colUpper = solver->getColUpper();

  int *cstat = new int[ncol];
  int *rstat = new int[nrow];
  solver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
                                          // 2: upper 3: lower

  int *basis_index = new int[nrow]; // basis_index[i] = 
                                    //        index of pivot var in row i
                                    //        (slack if number >= ncol) 

#ifdef RS_TRACETAB
  rs_printvecINT("cstat", cstat, ncol);
  rs_printvecINT("rstat", rstat, nrow);
#endif

  solver->getBasics(basis_index);

  cv_intBasicVar_frac = new int[ncol];  
  intBasicVar_frac = new int[ncol];                                 
  intNonBasicVar = new int[ncol];       
  contNonBasicVar = new int[ncol+nrow]; 
  nonBasicAtUpper = new int[ncol+nrow]; 
  nonBasicAtLower = new int[ncol+nrow]; 
  double dist_int;

  for(i=0; i<ncol; i++) {
    cv_intBasicVar_frac[i] = 0;

    switch(cstat[i]) {
    case 1: // basic variable
      
      dist_int = rs_above_integer(xlp[i]);
      if(solver->isInteger(i) && 
	  (dist_int > away_) && (dist_int < 1 - away_)) {
	cv_intBasicVar_frac[i] = 1;
	card_intBasicVar_frac++;

	// intBasicVar_frac computed below, 
	// as order must be according to selected rows

      }
      break;

    case 2: // Non basic at upper bound: must be flipped and shifted
            // so that it becomes non negative with lower bound 0 
            // It is assumed that bounds for integer variables have been
            // tightend so that non basic integer structural variables 
            // have integer values

      nonBasicAtUpper[card_nonBasicAtUpper] = i;
      card_nonBasicAtUpper++;

      if(solver->isInteger(i)) {
	intNonBasicVar[card_intNonBasicVar] = i;
	card_intNonBasicVar++;
      }
      else {
	contNonBasicVar[card_contNonBasicVar] = i;
	card_contNonBasicVar++;
      }
      break;

    case 3 : // non basic at lower bound: must be shifted so that it becomes
             // non negative with lower bound 0
            // It is assumed that bounds for integer variables have been
            // tightend so that they are integer

      nonBasicAtLower[card_nonBasicAtLower] = i;
      card_nonBasicAtLower++;

      if(solver->isInteger(i)) {
	intNonBasicVar[card_intNonBasicVar] = i;
	card_intNonBasicVar++;
      }
      else {
	contNonBasicVar[card_contNonBasicVar] = i;
	card_contNonBasicVar++;
      }
      break;

    default: // free variable ? Don't know how to handle 
      printf("### ERROR: CglRedSplit::generateCuts(): cstat[%d]: %d\n",
	     i, cstat[i]);
      exit(1);
      break;
    } 
  }

  for(i=0; i<nrow; i++) {
    switch(rstat[i]) {
    case 1: // basic slack
      break;

    case 2: // non basic slack at upper; flipped and shifted
      nonBasicAtUpper[card_nonBasicAtUpper] = ncol+i;
      card_nonBasicAtUpper++;

      contNonBasicVar[card_contNonBasicVar] = ncol+i;
      card_contNonBasicVar++;
      break;

    case 3: // non basic slack at lower; shifted
      nonBasicAtLower[card_nonBasicAtLower] = ncol+i;
      card_nonBasicAtLower++;

      contNonBasicVar[card_contNonBasicVar] = ncol+i;
      card_contNonBasicVar++;
      break;

    default: 
      printf("### ERROR: CglRedSlpit::generateCuts(): rstat[%d]: %d\n",
	     i, rstat[i]);
      exit(1);
      break;
    }
  }

  if((card_contNonBasicVar == 0) || (card_intBasicVar_frac == 0)) {
    delete[] cstat;
    delete[] rstat;
    delete[] basis_index;

    delete[] cv_intBasicVar_frac;  
    delete[] intBasicVar_frac;
    delete[] intNonBasicVar;
    delete[] contNonBasicVar;
    delete[] nonBasicAtUpper;
    delete[] nonBasicAtLower;
    solver->disableFactorization();

    return; // no cuts can be generated
  }

  /* Loop is mTab * mTab * max(mTab, nTab) so may be very expensive. 
     Reduce mTab if the above value is larger than maxTab_ */

  int new_mTab = card_intBasicVar_frac;
  double nc = (double) card_contNonBasicVar;
  double nc3 = nc * nc * nc;

  if(nc3 > maxTab_) {
    new_mTab = (int) (sqrt(maxTab_/card_contNonBasicVar));
  }
  else {
    new_mTab = (int) (CoinCbrt(maxTab_));
  }

  if(new_mTab == 0) {
    delete[] cstat;
    delete[] rstat;
    delete[] basis_index;
    
    delete[] cv_intBasicVar_frac;  
    delete[] intBasicVar_frac;
    delete[] intNonBasicVar;
    delete[] contNonBasicVar;
    delete[] nonBasicAtUpper;
    delete[] nonBasicAtLower;
    solver->disableFactorization();

    return; // no cuts can be generated
  }

  int start = 0;  // first row for selecting intBasicVar_frac, if too many.
                  // Removing rows used for generation whose pivot is some 
                  // var in intBasicVar_frac is valid since the corresponding 
                  // column in the optimal tableau is a column of the 
                  // identity matrix 

  if(new_mTab < card_intBasicVar_frac) {
    // Remove some of the rows used for generation.
    // Poor randomness; could do better if needed

    int seed = card_intBasicVar_frac;
    start = (int) (nrow * rs_genalea(&seed));
    
#ifdef RS_TRACE
    printf("CglRedSlpit::generateCuts(): mTab: %d  new_mTab: %d\n", 
	   card_intBasicVar_frac, new_mTab);
#endif

    card_intBasicVar_frac = new_mTab;
  }
  
  double *slack_val = new double[nrow];

  for(i=0; i<nrow; i++) {
    slack_val[i] = rhs[i] - row_act[i];
  }
  double *z = new double[ncol];  // workspace to get row of the tableau
  double *slack = new double[nrow];  // workspace to get row of the tableau

#ifdef RS_TRACETAB
  printOptTab(solver);
#endif

  mTab = card_intBasicVar_frac;
  nTab = card_contNonBasicVar;

  rhsTab = new double[mTab];
  int card_rowTab = 0;

  rs_allocmatDBL(&contNonBasicTab, mTab, nTab);
  rs_allocmatDBL(&intNonBasicTab, mTab, card_intNonBasicVar);

  card_intBasicVar_frac = 0; // recompute in pivot order

  for(i=0; i<nrow; i++) {

    if(mTab <= card_intBasicVar_frac) { 
      break;
    }

    int ind_row = start + i;
    if(ind_row > nrow) {
      ind_row -= nrow;
    }

    if(basis_index[ind_row] >= ncol) {
      continue;
    } 

    if(cv_intBasicVar_frac[basis_index[ind_row]] == 1) { 
                                                  // row used in generation
      intBasicVar_frac[card_intBasicVar_frac] = basis_index[ind_row];
      card_intBasicVar_frac++;
      rhsTab[card_rowTab] = xlp[basis_index[ind_row]];
      solver->getBInvARow(ind_row, z, slack);
      int ii;
      for( ii=0; ii<card_contNonBasicVar; ii++) {
	int locind = contNonBasicVar[ii];
	if(locind < ncol) {
	  contNonBasicTab[card_rowTab][ii] = z[locind];
	}
	else {
	  contNonBasicTab[card_rowTab][ii] = slack[locind - ncol];
	}
      }

      for( ii=0; ii<card_intNonBasicVar; ii++) {
	int locind = intNonBasicVar[ii];
	if(locind < ncol) {
	  intNonBasicTab[card_rowTab][ii] = z[locind];
	}
	else {
	  printf("### ERROR: CglRedSplit::generateCuts(): integer slack unexpected\n");
	  exit(1);
	}
      }

      card_rowTab++;
    }
  }

  rs_allocmatINT(&pi_mat, mTab, mTab);
  for(i=0; i<mTab; i++) {
    for(int ii=0; ii<mTab; ii++) {
      pi_mat[i][ii] = 0;
    }
    pi_mat[i][i] = 1;
  }

#ifdef RS_TRACE
  printf("intBasicVar_frac:\n");
  for(i=0; i<card_intBasicVar_frac; i++) {
    printf("%d ", intBasicVar_frac[i]);
  }
  printf("\n");
  printf("intNonBasicVar:\n");
  for(i=0; i<card_intNonBasicVar; i++) {
    printf("%d ", intNonBasicVar[i]);
  }
  printf("\n");
  printf("contNonBasicVar:\n");
  for(i=0; i<card_contNonBasicVar; i++) {
    printf("%d ", contNonBasicVar[i]);
  }
  printf("\n");
#endif

  if(given_optsol) {
    check_optsol(solver, 1, xlp, slack_val, 0);
  }

  reduce_contNonBasicTab();

  if(given_optsol) {
    check_optsol(solver, 2, xlp, slack_val, 0);
  }

  int card_row;
  double *row = new double[ncol+nrow];
  int *rowind = new int[ncol];
  double *rowelem = new double[ncol];

  const CoinPackedMatrix *byRow = solver->getMatrixByRow();
  const double *elements = byRow->getElements();
  const int *rowStart = byRow->getVectorStarts();
  const int *indices = byRow->getIndices();
  const int *rowLength = byRow->getVectorLengths(); 

  for(i=0; i<mTab; i++) {
    generate_row(i, row);
    flip(row);

    // RHS of equalities after flipping/translating non basic variables
    // is given by the current LP solution (re-ordered according to 
    // basis_index), i.e. what is now in rhsTab. 
    // RHS of row i of (pi_mat * contNonBasicTab) is then simply 
    // pi_mat[i] * rhsTab

    double rowrhs = rs_dotProd(pi_mat[i], rhsTab, mTab); 

    if(generate_cgcut(row, &rowrhs)) {
      unflip(row, &rowrhs, colLower, colUpper, slack_val);

      if(given_optsol) {
	check_optsol(solver, 3, row, rowrhs, i, 0);
      }

      eliminate_slacks(row, elements, rowStart, indices, 
		       rowLength, rhs, &rowrhs);

      if(given_optsol) {
	check_optsol(solver, 4, row, rowrhs, i, 0);
      }

      if(generate_packed_row(solver,row, rowind, rowelem, &card_row,rowrhs)) {
      	OsiRowCut rc;
	rc.setRow(card_row, rowind, rowelem);
	rc.setLb(-DBL_MAX);
	rc.setUb(rowrhs + EPS_RELAX);   // relax the constraint slightly
	cs.insertIfNotDuplicate(rc, CoinAbsFltEq(EPS_COEFF));
      }
    }
  }

  delete[] cstat;
  delete[] rstat;
  delete[] basis_index;
  delete[] slack;
  delete[] z;
  delete[] slack_val;
  delete[] row;
  delete[] rowind;
  delete[] rowelem;

  delete[] cv_intBasicVar_frac;  
  delete[] intBasicVar_frac;
  delete[] intNonBasicVar;
  delete[] contNonBasicVar;
  delete[] nonBasicAtUpper;
  delete[] nonBasicAtLower;
  rs_deallocmatDBL(&contNonBasicTab, mTab, nTab);
  rs_deallocmatDBL(&intNonBasicTab, mTab, card_intNonBasicVar);
  rs_deallocmatINT(&pi_mat, mTab, mTab);
  delete[] rhsTab;
  solver->disableFactorization();

  //return(cs.sizeRowCuts());
} /* generateCuts */

/***********************************************************************/
void CglRedSplit::setLimit(int limit)
{
  if (limit>0)
    limit_=limit;
} /* setLimit */

/***********************************************************************/
int CglRedSplit::getLimit() const
{
  return limit_;
} /* getLimit */

/***********************************************************************/
void CglRedSplit::setAway(double value)
{
  if (value>0.0 && value<=0.5)
    away_ = value;
}

/***********************************************************************/
double CglRedSplit::getAway() const
{
  return away_;
}

/***********************************************************************/
void CglRedSplit::setMaxTab(double value)
{
  if (value > 10) {
    maxTab_ = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setMaxTab(): value: %f ignored\n", 
	   value);
  }
}

/***********************************************************************/
double CglRedSplit::getMaxTab() const
{
  return maxTab_;
}

/***********************************************************************/
void CglRedSplit::setLUB(double value)
{
  if (value>0.0) {
    LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setLUB(): value: %f ignored\n", value);
  }
} /* setLUB */

/***********************************************************************/
double CglRedSplit::getLUB() const
{
  return LUB;
} /* getLUB */

/***********************************************************************/
void CglRedSplit::setEPS(double value)
{
  if (value>0.0 && value<=0.1) {
    EPS = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS(): value: %f ignored\n", value);
  }
} /* setEPS */

/***********************************************************************/
double CglRedSplit::getEPS() const
{
  return EPS;
} /* getEPS */

/***********************************************************************/
void CglRedSplit::setEPS_COEFF(double value)
{
  if (value > 0.0 && value <= 0.1) {
    EPS_COEFF = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS_COEFF(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_COEFF */

/***********************************************************************/
double CglRedSplit::getEPS_COEFF() const
{
  return EPS_COEFF;
} /* getEPS_COEFF */

/***********************************************************************/
void CglRedSplit::setEPS_COEFF_LUB(double value)
{
  if (value > 0.0 && value <= 0.1) {
    EPS_COEFF_LUB = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS_COEFF_LUB(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_COEFF_LUB */

/***********************************************************************/
double CglRedSplit::getEPS_COEFF_LUB() const
{
  return EPS_COEFF_LUB;
} /* getEPS_COEFF_LUB */

/***********************************************************************/
void CglRedSplit::setEPS_RELAX(double value)
{
  if (value > 0.0 && value <= 0.1) {
    EPS_RELAX = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setEPS_RELAX(): value: %f ignored\n", 
	   value);
  }
} /* setEPS_RELAX */

/***********************************************************************/
double CglRedSplit::getEPS_RELAX() const
{
  return EPS_RELAX;
} /* getEPS_RELAX */

/***********************************************************************/
void CglRedSplit::setNormIsZero(double value)
{
  if (value>0.0 && value<=1) {
    normIsZero = value;
  }
  else {
    printf("### WARNING: CglRedSplit::setNormIsZero(): value: %f ignored\n",
	   value);
  }
} /* setNormIsZero */

/***********************************************************************/
double CglRedSplit::getNormIsZero() const
{
  return normIsZero;
} /* getNormIsZero */

/***********************************************************************/
void CglRedSplit::setMinReduc(double value)
{
  if (value>0.0 && value<=1) {
    minReduc = value;
  }
  else {
    printf("### WARNING: CglRedSplit::MinReduc(): value: %f ignored\n",
	   value);
  }
} /* setMinReduc */

/***********************************************************************/
double CglRedSplit::getMinReduc() const
{
  return minReduc;
} /* getMinReduc */

/***********************************************************************/
void CglRedSplit::print() const
{
  rs_printvecINT("intBasicVar_frac", intBasicVar_frac, card_intBasicVar_frac);
  rs_printmatINT("pi_mat", pi_mat, card_intBasicVar_frac, 
		 card_intBasicVar_frac);
  rs_printvecINT("intNonBasicVar", intNonBasicVar, card_intNonBasicVar);
  rs_printmatDBL("intNonBasicTab", intNonBasicTab, card_intBasicVar_frac, 
		 card_intNonBasicVar);
  rs_printvecINT("contNonBasicVar", contNonBasicVar, card_contNonBasicVar);
  rs_printmatDBL("contNonBasicTab", contNonBasicTab, card_intBasicVar_frac, 
		 card_contNonBasicVar);
  rs_printvecINT("nonBasicAtLower", nonBasicAtLower, card_nonBasicAtLower);
  rs_printvecINT("nonBasicAtUpper", nonBasicAtUpper, card_nonBasicAtUpper);

} /* print */

/***********************************************************************/
void CglRedSplit::printOptTab(OsiSolverInterface *solver) const
{
  int i;
  const double *row_act = solver->getRowActivity();
  const double *rhs = solver->getRightHandSide();
  int *cstat = new int[ncol];
  int *rstat = new int[nrow];

  solver->getBasisStatus(cstat, rstat);   // 0: free  1: basic  
                                          // 2: upper 3: lower

  int *basis_index = new int[nrow]; // basis_index[i] = 
                                    //        index of pivot var in row i
                                    //        (slack if number >= ncol) 
  solver->getBasics(basis_index);

  double *z = new double[ncol];  // workspace to get row of the tableau
  double *slack = new double[nrow];  // workspace to get row of the tableau
  double *slack_val = new double[nrow];

  for(i=0; i<nrow; i++) {
    slack_val[i] = rhs[i] - row_act[i];
  }

  const double *rc = solver->getReducedCost();
  const double *dual = solver->getRowPrice();
  const double *solution = solver->getColSolution();

  rs_printvecINT("cstat", cstat, ncol);
  rs_printvecINT("rstat", rstat, nrow);
  rs_printvecINT("basis_index", basis_index, nrow);

  rs_printvecDBL("solution", solution, ncol);
  rs_printvecDBL("slack_val", slack_val, nrow);
  rs_printvecDBL("reduced_costs", rc, ncol);
  rs_printvecDBL("dual solution", dual, nrow);

  printf("Optimal Tableau:\n");

  for(i=0; i<nrow; i++) {
    solver->getBInvARow(i, z, slack);
    int ii;
    for( ii=0; ii<ncol; ii++) {
      printf("%5.2f ", z[ii]);
    }
    printf(" | ");
    for( ii=0; ii<nrow; ii++) {
      printf("%5.2f ", slack[ii]);
    }
    printf(" | ");
    if(basis_index[i] < ncol) {
      printf("%5.2f ", solution[basis_index[i]]);
    }
    else {
      printf("%5.2f ", slack_val[basis_index[i]-ncol]);
    }
    printf("\n");
  }
  int ii;
  for( ii=0; ii<7*(ncol+nrow+1); ii++) {
    printf("-");
  }
  printf("\n");

  for( ii=0; ii<ncol; ii++) {
    printf("%5.2f ", rc[ii]);    
  }
  printf(" | ");
  for( ii=0; ii<nrow; ii++) {
    printf("%5.2f ", -dual[ii]);
  }
  printf(" | ");
  printf("%5.2f\n", -solver->getObjValue());

  delete[] cstat;
  delete[] rstat;
  delete[] basis_index;
  delete[] slack;
  delete[] z;
  delete[] slack_val;
} /* printOptTab */

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglRedSplit::CglRedSplit () :
CglCutGenerator(),
nrow(0),
ncol(0),
LUB(1000.0),
EPS(1e-07),
EPS_COEFF(1e-08),
EPS_COEFF_LUB(1e-13),
EPS_RELAX(1e-8),
normIsZero(1e-05),
minReduc(0.05),
card_intBasicVar_frac(0),
card_intNonBasicVar(0),
card_contNonBasicVar(0),
card_nonBasicAtUpper(0),
card_nonBasicAtLower(0),
cv_intBasicVar_frac(0),
intBasicVar_frac(0),
intNonBasicVar(0), 
contNonBasicVar(0),
nonBasicAtUpper(0),
nonBasicAtLower(0),
mTab(0),
nTab(0),
pi_mat(0),
contNonBasicTab(0),
intNonBasicTab(0),
rhsTab(0),
away_(0.05),
limit_(50),
maxTab_(1e7),
given_optsol(0),
card_given_optsol(0)
{}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglRedSplit::CglRedSplit (const CglRedSplit & source) :
  CglCutGenerator(source),
  nrow(0),
  ncol(0),
  LUB(source.LUB),
  EPS(source.EPS),
  EPS_COEFF(source.EPS_COEFF),
  EPS_COEFF_LUB(source.EPS_COEFF_LUB),
  EPS_RELAX(source.EPS_RELAX),
  normIsZero(source.normIsZero),
  minReduc(source.minReduc),
  card_intBasicVar_frac(0),
  card_intNonBasicVar(0),
  card_contNonBasicVar(0),
  card_nonBasicAtUpper(0),
  card_nonBasicAtLower(0),
  cv_intBasicVar_frac(NULL),
  intBasicVar_frac(NULL),
  intNonBasicVar(NULL), 
  contNonBasicVar(NULL),
  nonBasicAtUpper(NULL),
  nonBasicAtLower(NULL),
  mTab(0),
  nTab(0),
  pi_mat(NULL),
  contNonBasicTab(NULL),
  intNonBasicTab(NULL),
  rhsTab(NULL),
  away_(source.away_),
  limit_(source.limit_),
  maxTab_(source.maxTab_),
  given_optsol(NULL),
  card_given_optsol(source.card_given_optsol)
{  
  if (source.nrow) {
    printf("### ERROR: CglRedSplit::CglRedSplit(): copy not implemented\n");
    exit(1);
  } else {
    // we are in good shape
    if (card_given_optsol) {
      assert (source.given_optsol);
      given_optsol = CoinCopyOfArray(source.given_optsol,card_given_optsol);
    }
  }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglRedSplit::clone() const
{
  return new CglRedSplit(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglRedSplit::~CglRedSplit ()
{}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglRedSplit &
CglRedSplit::operator=(const CglRedSplit& rhs)
{  

  printf("### ERROR: CglRedSplit::CglRedSplit(): operator '=' not implemented\n");
  exit(1);

  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    away_=rhs.away_;
    limit_=rhs.limit_;
  }
  return *this;
}
// Returns true if needs optimal basis to do cuts
bool 
CglRedSplit::needsOptimalBasis() const
{
  return true;
}
// Create C++ lines to get to current state
std::string
CglRedSplit::generateCpp( FILE * fp) 
{
  CglRedSplit other;
  fprintf(fp,"0#include \"CglRedSplit.hpp\"\n");
  fprintf(fp,"3  CglRedSplit redSplit;\n");
  if (getLimit()!=other.getLimit())
    fprintf(fp,"3  redSplit.setLimit(%d);\n",getLimit());
  else
    fprintf(fp,"4  redSplit.setLimit(%d);\n",getLimit());
  if (getAway()!=other.getAway())
    fprintf(fp,"3  redSplit.setAway(%g);\n",getAway());
  else
    fprintf(fp,"4  redSplit.setAway(%g);\n",getAway());
  if (getLUB()!=other.getLUB())
    fprintf(fp,"3  redSplit.setLUB(%g);\n",getLUB());
  else
    fprintf(fp,"4  redSplit.setLUB(%g);\n",getLUB());
  if (getEPS()!=other.getEPS())
    fprintf(fp,"3  redSplit.setEPS(%g);\n",getEPS());
  else
    fprintf(fp,"4  redSplit.setEPS(%g);\n",getEPS());
  if (getEPS_COEFF()!=other.getEPS_COEFF())
    fprintf(fp,"3  redSplit.setEPS_COEFF(%g);\n",getEPS_COEFF());
  else
    fprintf(fp,"4  redSplit.setEPS_COEFF(%g);\n",getEPS_COEFF());
  if (getEPS_COEFF_LUB()!=other.getEPS_COEFF_LUB())
    fprintf(fp,"3  redSplit.setEPS_COEFF_LUB(%g);\n",getEPS_COEFF_LUB());
  else
    fprintf(fp,"4  redSplit.setEPS_COEFF_LUB(%g);\n",getEPS_COEFF_LUB());
  if (getEPS_RELAX()!=other.getEPS_RELAX())
    fprintf(fp,"3  redSplit.setEPS_RELAX(%g);\n",getEPS_RELAX());
  else
    fprintf(fp,"4  redSplit.setEPS_RELAX(%g);\n",getEPS_RELAX());
  if (getNormIsZero()!=other.getNormIsZero())
    fprintf(fp,"3  redSplit.setNormIsZero(%g);\n",getNormIsZero());
  else
    fprintf(fp,"4  redSplit.setNormIsZero(%g);\n",getNormIsZero());
  if (getMinReduc()!=other.getMinReduc())
    fprintf(fp,"3  redSplit.setMinReduc(%g);\n",getMinReduc());
  else
    fprintf(fp,"4  redSplit.setMinReduc(%g);\n",getMinReduc());
  if (getMaxTab()!=other.getMaxTab())
    fprintf(fp,"3  redSplit.setMaxTab(%g);\n",getMaxTab());
  else
    fprintf(fp,"4  redSplit.setMaxTab(%g);\n",getMaxTab());
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  redSplit.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  redSplit.setAggressiveness(%d);\n",getAggressiveness());
  return "redSplit";
}
