// 
// Copyright (C) 2001, International Business Machines
// Corporation and others.  All Rights Reserved.
//----------------------------------------------------- 
// Simple example usage of doing "cut-and-branch" with the
// lift-and-project cut generator
//
// This sample program iteratively tightens a 
// given formulation by adding LAP cuts, then calls 
// branch-and-bound to solve the tightened
// formulation.
//
// usage:
//   cgl2 mpsFileName objectiveSense
// where:
//   mpsFileName: Name of an mps file (without the
//                file extension)
//   objectiveSense: min for minimization, 
//                   max for maximization. 
// example:
//   cgl2 ../../Data/Sample/p0033 min              
//----------------------------------------------------- 
#include <cassert>
#include <iostream>

using std::cout;
using std::endl;

#include "OsiCuts.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CglLiftAndProject.hpp"

using std::cerr;
using std::cout;
using std::endl;

int convert(OsiSolverInterface * stdSi, OsiSolverInterface * origSi);

int main(int argc, const char *argv[])
{
  // Make sure a file name and objective sense
  // were specified
  if ( argc<3 ) {
    cerr <<"Incorrect number of command line parameters." <<endl;
    cerr <<"  usage:" <<endl;
    cerr <<"    "<<argv[0] <<" mpsFileName objectiveSense" <<endl;
    cerr <<"  where:" <<endl;
    cerr <<"    mpsFileName: Name of an mps file" <<endl;
    cerr <<"                 without \".mps\" file extension" <<endl;
    cerr <<"    objectiveSense: min for minimization," <<endl; 
    cerr <<"                    max for maximization." <<endl; 
    return 1;
  }  

  // Make sure valid objective sense was specified
  if( strcmp(argv[2],"min") && strcmp(argv[2],"max" ) ){
    cerr <<"Unrecognized objective sense specifed" <<endl;
    cerr <<"  specified value: \"" <<argv[2] <<"\"" <<endl;

    cerr <<"  valid values: \"min\" for minimization," <<endl; 
    cerr <<"                \"max\" for maximization." <<endl; 
    return 1;
  }  

  // Instantiate a specific solver interface
  OsiClpSolverInterface si;

  // Read file describing problem
  si.readMps(argv[1],"mps"); 

  // Create an (empty) clone of the solver interface
  // and populate it with a canonical version 
  // of the original probelm. The canonical form is:
  //    min cx
  //    s.t.
  //    Ax>=b, x>=0, x_i={0,1} for i in B.
  //

  OsiSolverInterface * stdSi=si.clone(false);
  int retCode = convert(stdSi, &si);
  if (retCode > 0){
    cout << "Unable to convert problem into canonical form" << endl;
    return 0;
  }
  
  // Solve continuous problem
  stdSi->initialSolve();

  // Save the orig lp relaxation value for 
  // comparisons later
  double origLpObj = stdSi->getObjValue();

  // Track the total number of cuts applied
  int totalNumberApplied = 0;

  // Instantiate cut generators
  CglLiftAndProject lapCg;
  lapCg.setBeta(1);

  //---------------------------------------------------
  // Keep applying cuts until 
  //   1. no more cuts are generated
  // or
  //   2. the objective function value doesn't change
  //---------------------------------------------------
  bool equalObj;
  CoinRelFltEq eq(0.0001);
  OsiSolverInterface::ApplyCutsReturnCode acRc;
  double obj;

  do {
    // Get current solution value
    obj = stdSi->getObjValue();

    // Generate and apply cuts
    OsiCuts cuts;
    lapCg.generateCuts(*stdSi,cuts);
    acRc = stdSi->applyCuts(cuts,0.0);

    // Print applyCuts return code
    cout << endl << endl;
    cout <<cuts.sizeCuts() <<" cuts were generated" <<endl;
    cout <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<endl;
    cout <<"  " <<acRc.getNumInconsistentWrtIntegerModel() 
         <<" were inconsistent for this problem" <<endl;
    cout <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<endl;
    cout <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<endl;
    cout <<"  " <<acRc.getNumApplied() <<" were applied" <<endl;
    cout << endl << endl;

    // Increment the counter of total cuts applied
    totalNumberApplied += acRc.getNumApplied();

    // If no cuts were applied, then done
    if ( acRc.getNumApplied()==0 ) break;

    // Resolve
    stdSi->resolve();

    cout << endl;
    cout <<"After applying cuts objective value changed from "
         <<obj <<" to " <<stdSi->getObjValue() <<endl << endl;
    
    // -----------------------------------------------
    // Set Boolean flag to true if new objective is 
    // almost equal to prior value.
    //
    // The test is:
    // abs(oldObj-newObj) <= 0.0001*(max(abs(oldObj),abs(newObj))+1.);
    // see CoinRelFloatEqual.h 
    // -----------------------------------------------
    equalObj = eq( stdSi->getObjValue(), obj );
  } while( !equalObj );
  
  // Print total number of cuts applied, 
  // and total improvement in the lp objective value
  cout << endl << endl;
  cout << "----------------------------------------------------------" 
       << endl;
  cout << "Cut generation phase completed:" << endl;
  cout << "   " << totalNumberApplied << " cuts were applied in total," 
       << endl;
  cout << "   changing the lp objective value from " << origLpObj 
       << " to " << stdSi->getObjValue() << endl;
  cout << "----------------------------------------------------------" 
       << endl;
  cout << endl << endl;

  // If you want to solve problem, change "#if 0" to "#if 1"
#if 0
  // Solve MIP Problem
  stdSi->branchAndBound();

  // Print the solution
  cout <<"The objective function value: " <<stdSi->objvalue() <<endl;
  const double * soln = stdSi->colsol();
  int i;
  for ( i=0; i<stdSi->numcols(); i ++ ) {
    cout <<" x[" <<i <<"] = " <<soln[i] <<endl;
  }
#endif

  return 0;
}

//#if 0

//=============================================================

int convert(OsiSolverInterface * stdSi, OsiSolverInterface * origSi)
{
  // Convert routine
  // Converts the origSi's lp into canonical form:
  //    min cx
  //    s.t.
  //    Ax>=b, x>=0, x_i={0,1} for i in B.
  // and stores the convert problem in stdSi.
  // The routine returns a nonzero the conversion fails.
  
  // Osi get/set methods return FALSE 
  // if method not implemented for solver.
  // IMPROVEME: test return codes of get/sets. 
  
  CoinRelFltEq eq;

  int origNumRows = origSi->getNumRows();
  int origNumCols = origSi->getNumCols();
  
  // Check for binary variables
  //  if any general integer variables, reject (in future: bin expansion)
  //  if no bin variables, then reject (print reason).
  int numBin=0;
  int numGenInt=0;
  int numCont=0;
  int colIndex=0;

  // stdNumCols will be equal to
  // origNumCols + number of unbounded cols
  // which is less than 2*origNumCols
  const int stdMaxNumCols = 2*origNumCols; 
  
  // Note: Empty clones don't contain binary information.
  // The loadProblem method doesn't pass binary information.
  // Binary information must be tracked and loaded 
  // separately. 
  int * binIndices = new int[origNumCols];
  int * contIndices = new int[origNumCols]; 
  while (numGenInt==0 && colIndex<origNumCols){
    if(origSi->isContinuous(colIndex)){
      contIndices[numCont]=colIndex;
      numCont++;
    }
    else if(origSi->isIntegerNonBinary(colIndex))numGenInt++;
    else if(origSi->isBinary(colIndex)){
      binIndices[numBin]=colIndex;
      numBin++;
    }
    colIndex++; 
  }
  if (numGenInt>0){
    delete [] binIndices;
    delete [] contIndices;
    return 1; 
  }
  if (numBin==0) {  
    delete [] binIndices;
    delete [] contIndices;
    return 1; 
  }
  
  // stdNumRows will be equal to
  // origNumRows +num'E' rows + num'R'rows 
  //   + 2*(origNumCols) - (num cols with exactly one finite bound)
  // which is less than
  //   2*origNumRows + 2*origNumCols
  const int stdMaxNumRows = 2*origNumRows + 2*origNumCols;
  const double * origRhs = origSi->getRightHandSide();
  double * stdRhs = new double[stdMaxNumRows]; 
  memset(stdRhs,0,(stdMaxNumRows)*sizeof(double));
  memcpy(stdRhs,origRhs,origNumRows*sizeof(double));

  const double * origObjCoefs = origSi->getObjCoefficients();
  double * stdObjCoefs = new double[stdMaxNumCols]; 
  memset(stdObjCoefs,0,(stdMaxNumCols)*sizeof(double));
  memcpy(stdObjCoefs,origObjCoefs,origNumCols*sizeof(double)); 

  const double * origColUppers = origSi->getColUpper();
  double * stdColUppers = new double[stdMaxNumCols]; 
  memset(stdColUppers,0,(stdMaxNumCols)*sizeof(double));
  memcpy(stdColUppers,origColUppers,origNumCols*sizeof(double)); 
  
  // Transform objective to a minimization problem.
  // Set the objective sense after loading the problem
  // into the solver interface.
  if ( origSi->getObjSense()!=1){
    // if max, negate objective coeffients
    // assuming that the code will exit gracefully if numCols<=0;
    if ( origObjCoefs!=NULL ){
      int i;
      for (i=0; i<origNumCols; i++) stdObjCoefs[i]= -origObjCoefs[i]; 
    }
  }

  // Here's the story plot from here down.
  // If this looks ugly - it's because it is :-)
  // 
  // Get a col order copy of the constraint matrix.
  // Convert the columns to std form. 
  // Build up the "partial"-ly converted matrix
  // Create a packed matrix from the partially converted info
  // Transpose the packed matrix into row order
  // Get a row order copy of the matrix
  // Convert the rows to std form
  // Create a packed matrix
  // Load the packed matrix into the solver interface

  // Get a column Ordered packed Matrix and pointers
  const CoinPackedMatrix * origAByCol =  origSi->getMatrixByCol();
  const double * origColElements =  origAByCol->getElements();
  const int * origColIndices =  origAByCol->getIndices();
  const int * origColStarts = origAByCol->getVectorStarts();
  const int * origColLengths = origAByCol->getVectorLengths();
  
  // parElements, parIndices, parStarts
  // IMPROVEME: don't need so much space for just column transformation
  int stdMaxNumElems = stdMaxNumCols*stdMaxNumRows;
  double * parElements = new double[stdMaxNumElems];
  int * parIndices = new int[stdMaxNumElems];
  int * parStarts = new int[stdMaxNumCols+1];
  int * parLengths = new int[stdMaxNumCols];
  int parNumRows=origNumRows;
  int parNumCols=origNumCols;
  // If no gaps, this would be numElements, 
  // but it's really the end position of the parElements's vector
  int parNumElements = origColStarts[origNumCols];

  // orig matrix may have gaps so copy it all
  int fullLength = origAByCol->getVectorLast(origNumCols-1);
  memset(parElements,0,stdMaxNumElems*sizeof(double));
  memcpy(parElements,origColElements,fullLength*sizeof(double));
  memset(parIndices,0,stdMaxNumElems*sizeof(int));
  memcpy(parIndices,origColIndices,fullLength*sizeof(int));
  memset(parStarts,0,(stdMaxNumCols+1)*sizeof(int));
  memcpy(parStarts,origColStarts,(origNumCols+1)*sizeof(int));
  memset(parLengths,0,(stdMaxNumCols)*sizeof(int));
  memcpy(parLengths,origColLengths,(origNumCols)*sizeof(int));

  //
  // Handle Columns
  //    Just make transformations on cols
  //    Don't actually add constraints arising 
  //    from explicit bounds until "Handle Rows" section
  //    Number of rows remains constant.
  //    New cols will be added in transforming free 
  //    variables into standard form.
  //

  const double inf = origSi->getInfinity();
  
  const double * origColLowers = origSi->getColLower();

  //  for i not Binary:
  //     lb <= x_j <= ub
  //     lb = -inf, neg, 0, pos (if inf then prob doesn't make sense)
  //     ub = neg, 0, pos, inf (if -inf then prob doesn't make sense)
  int i;
  for (i=0; i<numCont; i++){
    int origContIndex = contIndices[i];
    double origColLb = origColLowers[origContIndex];
    double origColUb = origColUppers[origContIndex];
    
    // if (lb >= inf) then there's a problem(?)so exitGracefully;
    if (origColLb >= inf) {
      delete [] binIndices;
      delete [] contIndices;
      delete [] stdObjCoefs;
      delete [] stdRhs;
      delete [] parElements;
      delete [] parIndices;
      delete [] parStarts;
      return 1;
    }

    // IMPROVEME:Could test if |ub-lb|<epsilon, then fix var.
    
    // if (lb==0) then 
    //     ultimately add constraint x_j >= 0.0
    //     if (ub < inf) 
    //       ultimately add constraint -x_j >= -ub
    //     no transformations needed here
    
    // else if (lb > -inf) then
    //     // replace x_j >= lb with y_j (= x_j - lb) >= 0
    //     // (changes optimal objValue, but not optimal col solution)
    //     for every row i, such that a_ij<>0 
    // ***    rhs[row_i]-=a_ij*lb
    //     ultimately add constraint x_j >= 0
    //     if(ub < inf) then
    //        ub=ub-lb
    //        ultimately add constraint -x_j >= -ub
    if (origColLb > -inf){
      int k;
      int origColEnd = origColStarts[origContIndex]+origColLengths[origContIndex];
      for(k=origColStarts[origContIndex]; k<origColEnd; k++){
	int rowI = origColIndices[k];
	double elemIJ = origColElements[k]; 
	stdRhs[rowI] -= elemIJ * origColLb;
      }
      if (origColUb < inf){
	stdColUppers[origContIndex] -= origColLb;
      }
    }

    // else if (lb = -inf) then
    else if (origColLb <= -inf){

    //     if (ub = 0) then
    //      // essentially, replace x_j with y_j=-x_j
    //         for every row i, such that a_ij<>0
    //  ***        a_ij = -a_ij
    //         ultimately add constraint x_j >= 0
    //         ub = inf
    //         negate objCoef_j
      if(eq(origColUb,0)){
	int k;
	int origColEnd = origColStarts[origContIndex]+origColLengths[origContIndex];
	for(k=origColStarts[origContIndex]; k<origColEnd; k++){
	  parElements[k]= -parElements[k];
	  stdColUppers[origContIndex]= inf;
	  stdObjCoefs[origContIndex]=-stdObjCoefs[origContIndex];
	}
      }

    //     else if (ub<inf) then
    //         // essentially replace x_j with y_j=-x_j+ub (x_j=ub-y_j)
    //         for every row i, such that a_ij<>0
    // ***         a_ij = -a_ij
    //             rhs[row_i]-=a_ij*ub
    //         ultimately add constraint, x_j >= 0
    //         ub = inf
    //         negate objCoef
      else if(origColUb < inf){
	int k;
	int origColEnd = origColStarts[origContIndex]+origColLengths[origContIndex];
	for(k=origColStarts[origContIndex]; k<origColEnd; k++){
	  parElements[k]= -parElements[k];
          stdRhs[origColIndices[k]] -= parElements[k]*origColUb;
	  stdColUppers[origContIndex]= inf;
	  stdObjCoefs[origContIndex]=-stdObjCoefs[origContIndex];
	}
      }

    //     else if (ub=inf) then
    //          // replace x_i free by x_i=x_i+ - x_i-, each nonneg
    //          copy the column assoc with x_i, 
    //             Colindex=numCols, numCols=numCols+1
    //          negate the sign of every entry
    //   ***    add the new col. (Add length to parLengths)
    //          ultimately add contraints x_j>=0, and x_numCols>=0
    //          ub_j = inf and ub_numCols = inf
    //          objCoef x_i- = - objCoef of x_i
    //    Note:at this point contIndices is no longer in sync with parIndices 
      else if (origColUb >= inf) {
	int k;
	int origColEnd = origColStarts[origContIndex]+origColLengths[origContIndex];
	for(k=origColStarts[origContIndex]; k<origColEnd; k++){
	  parElements[parNumElements] = - origColElements[k];
	  parIndices[parNumElements] = origColIndices[k];
	  parNumElements++;
	}
	parLengths[parNumCols]= origColLengths[origContIndex];
	stdColUppers[parNumCols]=inf;
	stdObjCoefs[parNumCols]= -origObjCoefs[origContIndex];
	parNumCols++;
	parStarts[parNumCols] = parNumElements;
      }
    }
  }

  // Create packed matrix from partially transormed information

  // There's a bug. When last paramater of the constructor is 0, the 
  // packedMatrix still has gaps.
  // When you feed in the lengths explicitly, it still gives you gaps.
  // It'd be cleaner not to have gaps (I don't need them)..
  // but this is it for now.

  // Laci suggested to try writing a nonconst version of getElements
  // declaring an CoinPackedMatrix, then
  // Using the "copyOf" method, setting extraMajor=0, but extraGap=%
  // (or is that vice-versa? Anyway gap within rows=0, roomToAddRows=%)
  // And working within the packed matrix.
  // Append rows by new-ing once, and reusing the memory.
  // For now...we use what's currently in -- and tested -- in Osi.
  
  // This gives a warning about taking the address of temp memory..
  // IMPROVEME: Need to make parAByRow solid and change -> to .
  CoinPackedMatrix * parAByRow =
     new CoinPackedMatrix(true, parNumRows, parNumCols, parNumElements,
			  parElements, parIndices, parStarts,parLengths);
  parAByRow->reverseOrdering();
  delete [] parElements;
  delete [] parIndices;
  delete [] parStarts;
  delete [] parLengths;

  // Get a Row Ordered vectors
  const double * parRowElements =  parAByRow->getElements();
  const int * parRowIndices =  parAByRow->getIndices();
  const int * parRowStarts = parAByRow->getVectorStarts();
  const int * parRowLengths = parAByRow->getVectorLengths();

  // stdElements, stdIndices, stdStarts
  // IMPROVEME: should have better estimate now for maxnumelems...
  double * stdElements = new double[stdMaxNumElems];
  int * stdIndices = new int[stdMaxNumElems];
  int * stdStarts = new int[stdMaxNumRows+1];
  int * stdLengths = new int[stdMaxNumRows];
  int stdNumRows=parNumRows;
  int stdNumCols=parNumCols;

  // If no gaps, then this would be numElements, but
  // otherwise, it's the end position 
  int stdNumElements = parRowStarts[parNumRows]; 

  // there may be gaps (until the bug is fixed), so copy full length
  fullLength = parAByRow->getVectorLast(parNumRows-1);
  memset(stdElements,0,stdMaxNumElems*sizeof(double));
  memcpy(stdElements,parRowElements,fullLength*sizeof(double));
  memset(stdIndices,0,stdMaxNumElems*sizeof(int));
  memcpy(stdIndices,parRowIndices,fullLength*sizeof(int));
  memset(stdStarts,0,(stdMaxNumRows+1)*sizeof(int));
  memcpy(stdStarts,parRowStarts,(parNumRows+1)*sizeof(int));
  memset(stdLengths,0,(stdMaxNumRows)*sizeof(int));
  memcpy(stdLengths,parRowLengths,(parNumRows)*sizeof(int));

  //
  // Handle Rows - number of columns remains constant
  //
  
  // Change all row senses to ">="
  //  if <= then change sign of row coef and rhs, and sense
  //  if =  then create 2 rows one >= and one neg >= neg
  //  if "range", handle similiarily to "=" constraints
  const char * origRowSense = origSi->getRowSense();
  const double * origRowRange = origSi->getRowRange(); // unaffected by transformation
  int origRowIndex;
  for (origRowIndex=0; origRowIndex<origNumRows; origRowIndex++){
    // if (origRowSense[origRowIndex]=='G')
      // row in std format -- do nothing

    if (origRowSense[origRowIndex]=='L'){
      // negate row coeff and rhs; change sense
      stdRhs[origRowIndex] = -stdRhs[origRowIndex];
      int k;
      int parRowEnd = parRowStarts[origRowIndex]+ parRowLengths[origRowIndex];
      for(k=parRowStarts[origRowIndex]; k<parRowEnd; k++){
	  stdElements[k]= - stdElements[k];
      }
    }

    else if (origRowSense[origRowIndex]=='R'){
      // row <= ub, range = ub-lb (lb<=row<=ub, lb=ub-range)
      // change sense to "G", replace rhs with lb
      // create a duplicate row
      // negate duplicate row, change sense to "G", new rhs = -ub 
      // add duplicate row to problem, increment stdNumRows
      // (i.e. row >=ub-range, -row>=-ub)
      // Note: stdRhs <> origRhs at this point
      stdRhs[stdNumRows]= -stdRhs[origRowIndex];
      stdRhs[origRowIndex] = stdRhs[origRowIndex]-origRowRange[origRowIndex];	
      // stdStarts[stdNumRows]=stdStarts[stdNumRows-1]+stdLengths[stdNumRows-1];
      int k;
      int parRowEnd = parRowStarts[origRowIndex]+parRowLengths[origRowIndex];
      for(k=parRowStarts[origRowIndex]; k<parRowEnd; k++){
	stdElements[stdNumElements] = -parRowElements[k];
	stdIndices[stdNumElements] = parRowIndices[k];
	stdNumElements++;
      }
      stdLengths[stdNumRows]=parRowLengths[origRowIndex];
      stdNumRows++;
      stdStarts[stdNumRows] = stdStarts[stdNumRows-1]+stdLengths[stdNumRows-1];
    }

    else if (origRowSense[origRowIndex]=='E'){
      // change row sense, 
      // duplicate row, negate sense of duplicate and duplicate rhs, 
      // change sense, add new row to the problem
      // row >= rhs  -row >= -rhs
      stdRhs[stdNumRows] = - stdRhs[origRowIndex];
      // stdStarts[stdNumRows]=stdStarts[stdNumRows-1]+stdLengths[stdNumRows-1];
      int k;
      int parRowEnd = parRowStarts[origRowIndex]+parRowLengths[origRowIndex];
      for(k=parRowStarts[origRowIndex]; k<parRowEnd; k++){
	stdElements[stdNumElements] = -parRowElements[k];
	stdIndices[stdNumElements] = parRowIndices[k];
	stdNumElements++;
      }
      stdLengths[stdNumRows]=parRowLengths[origRowIndex];
      stdNumRows++;
      stdStarts[stdNumRows] = stdStarts[stdNumRows-1]+stdLengths[stdNumRows-1];
    } 
    // N: else if (origRowSense[origRowIndex]=='N') do nothing; ignore
  } 


  //
  // Make Bounds Explicit - number of cols remains constant
  //

  // IMPROVEME: check for lb=ub?  
  // At this point all var lbs 
  // have been transformed to 0.
 
  // For every variable, 
  //    add constriaint x_i>=0
  //    increment stdNumRows
  //    if ub < inf, 
  //       add -x_i >= -ub
  //       increment stdNumRows
  int j;
  for (j=0; j<stdNumCols; j++){
    stdRhs[stdNumRows]=0;
    stdLengths[stdNumRows]=1;
    stdElements[stdNumElements]=1;
    stdIndices[stdNumElements]=j;
    stdNumElements++;
    stdNumRows++;
    stdStarts[stdNumRows]= stdNumElements;
    if(stdColUppers[j]<inf){
      stdRhs[stdNumRows]= -stdColUppers[j];
      stdLengths[stdNumRows]=1;
      stdElements[stdNumElements]=-1;
      stdIndices[stdNumElements]=j;
      stdNumElements++;
      stdNumRows++;
      stdStarts[stdNumRows]=stdNumElements;
    }
  }
  stdStarts[stdNumRows]=stdNumElements;

  // senses have been transformed to "G"
  char * stdRowSense = new char[stdNumRows];
  memset(stdRowSense,'G',(stdNumRows)*sizeof(char));

  // all column lower bounds 
  // have been transformed to  zero
  double * stdColLowers = new double[stdNumCols];
  memset(stdColLowers,0,(stdNumCols)*sizeof(double));

  // all row upper bounds have been 
  // transformed to inf
  double * stdRowUppers = new double[stdNumRows];
  CoinFillN(stdRowUppers,stdNumRows,inf);

  // load the standardized problem into stdSi    
  CoinPackedMatrix * stdAByRow =
     new CoinPackedMatrix(false, stdNumCols, stdNumRows, stdNumElements,
			  stdElements,stdIndices, stdStarts,stdLengths);

  stdSi->loadProblem (*stdAByRow, stdColLowers, stdColUppers, stdObjCoefs,
		      stdRhs, stdRowUppers);

  stdSi->setObjSense(1.0);

  stdSi->setInteger(binIndices, numBin);


  return 0;  
  
}

