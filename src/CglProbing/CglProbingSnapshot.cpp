/*
  Copyright (C) 2002 -- 2010
  International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id#
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinPackedMatrix.hpp"
#include "CglProbing.hpp"

/*
  Create row- and column-major copies of the matrix passed in via the solver
  interface for use by the CglProbing object in further probing work.

  The parameter possible allows for selective use of constraints. Set entry i
  to 1 to select constraint i, 0 to ignore it.
  
  Column bounds are tightened (call to tighten()) and row bounds are
  reevaluated before the column-major copy is made.  If the possible array is
  present, values of 1 will be changed to 0 if the corresponding row is found
  to be redundant. The point here is that it's painful to excise rows from the
  column-major copy.

  The purpose is to speed up probing (we don't have to regenerate the row-
  and column-major copies at each call) and also to make it easier to do
  global cuts (because we can use the possible array to specify that only
  original constraints are to be considered) or impose other restrictions.

  TODO: Now that I've gone through the setup code in gutsOfGenerateCuts and
  	pulled it out into setupRowCopy and groomModel, I can really see the
	commonality. It should be possible to merge these. Also notice that
	we call tighten here, which is something I've been thinking is a good
	idea.  -- lh, 110212 --
*/
int CglProbing::snapshot ( const OsiSolverInterface & si,
			   char * possible, bool withObjective)
{
  deleteSnapshot() ;

/*
  Get the basic problem information and create modifiable copies of the row and
  column bounds arrays.
*/
  numberColumns_ = si.getNumCols(); 
  numberRows_ = si.getNumRows() ;
  colLower_ = new double[numberColumns_] ;
  colUpper_ = new double[numberColumns_] ;
  CoinMemcpyN(si.getColLower(),numberColumns_,colLower_) ;
  CoinMemcpyN(si.getColUpper(),numberColumns_,colUpper_) ;
  rowLower_ = new double [numberRows_+1] ;
  rowUpper_ = new double [numberRows_+1] ;
  CoinMemcpyN(si.getRowLower(),numberRows_,rowLower_) ;
  CoinMemcpyN(si.getRowUpper(),numberRows_,rowUpper_) ;
/*
  Disable rows specified in possible.
*/
  if (possible != 0) {
    for (int i = 0 ; i < numberRows_ ; i++) {
      if (!possible[i]) {
	rowLower_[i] = -DBL_MAX ;
	rowUpper_[i] = DBL_MAX ;
      }
    }
  }
/*
  Get a modifiable copy of the column type array and count integer variables.

  TODO: See below. As far as I can tell, the only reason we need a modifiable
	copy is because tighten() takes a non-const parameter. The modified
	vector is discarded.  -- lh, 100917 --
*/
  const char *intVarOriginal = si.getColType(true) ;
  char *intVar = CoinCopyOfArray(intVarOriginal,numberColumns_) ;
  numberIntegers_ = 0 ;
  number01Integers_ = 0 ;
  for (int i = 0 ; i < numberColumns_ ; i++) {
    if (intVar[i]) {
      numberIntegers_++ ;
      if (intVar[i] == 1) {
        number01Integers_++ ;
      }
    }
  }
/*
  Grab the row-major copy and massage it so that all negative coefficients in a
  constraint precede all positive coefficients. Remember where the positive
  coefficients start in each row.

  TODO: rowStartPos isn't deleted until after the call to tighten(), for no
	apparent reason. At a guess, John might have intended to add this as
	a parameter to tighten() and either never got to it, or concluded it
	wasn't necessary.  -- lh, 100917 --
*/
  rowCopy_ = new CoinPackedMatrix(*si.getMatrixByRow()) ;
  int *column = rowCopy_->getMutableIndices() ;
  const CoinBigIndex *rowStart = rowCopy_->getVectorStarts() ;
  const int *rowLength = rowCopy_->getVectorLengths(); 
  double *rowElements = rowCopy_->getMutableElements() ;
  int *column2 = new int[numberColumns_] ;
  double *elements2 = new double[numberColumns_] ;
  CoinBigIndex *rowStartPos = new CoinBigIndex [numberRows_] ;
  for (int i = 0 ; i < numberRows_ ; i++) {
    CoinBigIndex start = rowStart[i] ;
    CoinBigIndex end = start+rowLength[i] ;
    int nOther = 0 ;
    for (CoinBigIndex j = start ; j < end ; j++) { 
      int iColumn = column[j] ;
      double value = rowElements[j] ;
      if (value < 0.0) {
	rowElements[start] = value ;
	column[start++] = iColumn ;
      } else {
	elements2[nOther] = value ;
	column2[nOther++] = iColumn ;
      }
    }
    rowStartPos[i] = start ;
    for (int k = 0 ; k < nOther ; k++) {
      rowElements[start] = elements2[k] ;
      column[start++] = column2[k] ;
    }
  }
  delete [] column2 ;
  delete [] elements2 ;
/*
  Tighten column and row bounds, if possible. This can result in infeasibility.

  TODO: If ninfeas > 1 (one or more variables infeasible), shouldn't we
	just clean up and bail out now?  -- lh, 100917 --
  
  TODO: It seems likely that tighten() could reduce general integers to the
	point where they are limited to [0,1]. Does it rewrite intVar? Seems
	likely. Check.  -- lh, 100917 --

  TODO: Seems like it'd make sense to have a version of tighten() that takes
	a CglProbing object as a parameter. A bit further: why not package the
	snapshot as its own class? That'd simplify parameter passing and we
	could just have a mutable snapshot.
*/
  int returnCode = 0 ;
  int ninfeas = tighten(colLower_,colUpper_,column,rowElements,
  			rowStart,NULL,rowLength,rowLower_,rowUpper_,
			numberRows_,numberColumns_,intVar,5,primalTolerance_) ;
  delete [] rowStartPos ;
  if (ninfeas) {
    returnCode = 1 ;
  }
/*
  jjf: do integer stuff for mode 0

  TODO: We throw our mutable intVar away after this loop. If it's really
	different from the nonmutable intVarOriginal, seems like we're in
	trouble. Not the least because if it has more binary vars, we've
	overrun the end of cutVector_. -- lh, 100917 --

  Ok, now I think I see the point here. We need to have this in place in the
  snapshot for the benefit of the mode 0 code that populates lookedAt_. It
  will scan this vector and add variables for which there is no probing info
  (cutVector_[].index == 0) to lookedAt_. So here we're arranging for all
  binary variables to be probed.  -- lh, 110212 --
*/
  cutVector_ = new disaggregation [number01Integers_] ;
  memset(cutVector_,0,number01Integers_*sizeof(disaggregation)) ;
  number01Integers_ = 0 ;
  for (int i = 0 ; i < numberColumns_ ; i++) {
    if (intVar[i] == 1)
      cutVector_[number01Integers_++].sequence = i ;
  }
  delete [] intVar ;
/*
  Check for redundant rows and update possible. Physically delete the
  redundant rows from the row-major copy.

  TODO: Yet another magic number for an `infinite' bound. That's three so far.
	-- lh, 100917 --

  TODO: Seems like we could avoid running two loops when possible is present.
	Rewrite the if with one loop that does possible and index, etc. and
	an else loop that just does index, etc.  -- lh, 100917 --
*/
  if (possible != 0) {
    for (int i = 0 ; i < numberRows_ ; i++) {
      if (rowLower_[i] < -1.0e30 && rowUpper_[i] > 1.0e30) 
	possible[i] = 0 ;
    }
  }
  int *index = new int[numberRows_] ;
  int nDrop = 0 ;
  int nKeep = 0 ;
  for (int i = 0 ; i < numberRows_ ; i++) {
    if (rowLower_[i] < -1.0e30 && rowUpper_[i] > 1.0e30) {
      index[nDrop++] = i ;
    } else {
      rowLower_[nKeep] = rowLower_[i] ;
      rowUpper_[nKeep++] = rowUpper_[i] ;
    }
  }
  numberRows_ = nKeep ;
  if (nDrop)
    rowCopy_->deleteRows(nDrop,index) ;
  delete [] index ;
/*
  Add the objective to the matrix, if requested.

  TODO: What's the point? No row bounds are set; the row bound arrays aren't
	enlarged to accommodate it. We'll have to deal with this as a special
	case everywhere. Far worse, gutsOfGenerateCut *does* assume that the
	row bound arrays have an extra entry for the objective bounds. Check
	this for consistency!

	And the answer, looking further in gutsOfGenerateCut, is that we
	don't actually use this cached row copy. No, we create a new copy
	of the cached copy, and new copies of the cached row bounds, and
	work on them. The row bounds are resized appropriately there.
	
	-- lh, 100917 --

  TODO: Most places, negating the objective is phrased as sense*coefficient.
	Why is it an if here?  -- lh, 100917 --

  TODO: Maybe we could just keep one of the previous column-length arrays
	hanging around for use here. column2 and elements2 would do fine.
*/
  if (withObjective) {
    int *columns = new int[numberColumns_] ;
    double *elements = new double[numberColumns_] ;
    int n = 0 ;
    const double *objective = si.getObjCoefficients() ;
    bool maximize = (si.getObjSense() == -1) ;
    for (int i = 0 ; i < numberColumns_ ; i++) {
      if (objective[i]) {
        elements[n] = (maximize) ? -objective[i] : objective[i] ;
        columns[n++] = i ;
      }
    }
    rowCopy_->appendRow(n,columns,elements) ;
    delete [] columns ;
    delete [] elements ;
    numberRows_++ ;
  }
/*
  Create the column-major copy (no extra space). Allow that the matrix may be
  empty by this point.

  Make sure that we have the full column dimension. It can happen that we've
  dropped all the constraints that reference a variable. You'd think we could
  resize the row copy and the column copy would then be correct by
  construction, but the CoinPackedMatrix copy methods are too smart for our
  own good here; empty columns don't register.
*/
  if (rowCopy_->getNumElements()) {
    columnCopy_ = new CoinPackedMatrix(*rowCopy_,0,0,true) ;
  } else {
    columnCopy_ = new CoinPackedMatrix() ;
  }
  columnCopy_->setDimensions(numberRows_,numberColumns_) ;
  rowCopy_->setDimensions(numberRows_,numberColumns_) ;

  return returnCode ;
}


// Delete snapshot
void CglProbing::deleteSnapshot()
{
  delete [] rowLower_ ;
  delete [] rowUpper_ ;
  delete [] colLower_ ;
  delete [] colUpper_ ;
  delete rowCopy_ ;
  delete columnCopy_ ;
  rowCopy_ = NULL ;
  columnCopy_ = NULL ;
  rowLower_ = NULL ;
  rowUpper_ = NULL ;
  colLower_ = NULL ;
  colUpper_ = NULL ;
  for (int i = 0 ; i < number01Integers_ ; i++) {
    delete [] cutVector_[i].index ;
  }
  delete [] cutVector_ ;
  numberIntegers_ = 0 ;
  number01Integers_ = 0 ;
  cutVector_ = NULL ;
}

