// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
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

#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#ifdef FULL_DEBUG
#include "OsiSimplexSolution.hpp"
#endif
#include "CoinWarmStartBasis.hpp"
#include "CglGomory.hpp"
#include "CoinFinite.hpp"

//-------------------------------------------------------------------
// Generate disaggregation cuts
//------------------------------------------------------------------- 
void CglGomory::generateCuts(const OsiSolverInterface & si, 
						OsiCuts & cs ) const
{
  // Get basic problem information
  int numberColumns=si.getNumCols(); 

  // get integer variables and basis
  char * intVar = new char[numberColumns];
  int i;
  CoinWarmStart * warmstart = si.getWarmStart();
  const CoinWarmStartBasis* warm =
    dynamic_cast<const CoinWarmStartBasis*>(warmstart);
  const double * colUpper = si.getColUpper();
  const double * colLower = si.getColLower();
  for (i=0;i<numberColumns;i++) {
    if (si.isInteger(i)) {
      if (colUpper[i]>colLower[i]+0.5) {
	if (fabs(colUpper[i]-1.0)<1.0e-12&&fabs(colLower[i])<1.0e-12) {
	  intVar[i]=1; //0-1
	} else  if (colLower[i]>=0.0) {
	  intVar[i] = 2; // other
	} else {
	  // negative bounds - I am not sure works
	  intVar[i] = 0;
	}
      }
    } else {
      intVar[i]=0;
    }
  }
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&!debugger->onOptimalPath(si))
    debugger = NULL;
#else
  const OsiRowCutDebugger * debugger = NULL;
#endif

  generateCuts(debugger, cs, *si.getMatrixByCol(),
	   si.getObjCoefficients(), si.getColSolution(),
	   si.getColLower(), si.getColUpper(), 
	   si.getRowLower(), si.getRowUpper(),
	   intVar,warm);

  delete warmstart;
  delete [] intVar;
}

// Returns value - floor but allowing for small errors
inline double above_integer(double value) 
{
  double value2=floor(value);
  if (value-value2>0.9999999) {
    value2 = floor(value+0.5);
    assert (value-value2<=0.00000011);
  }
  return value-value2;
}
//-------------------------------------------------------------------
// Returns the greatest common denominator of two 
// positive integers, a and b, found using Euclid's algorithm 
//-------------------------------------------------------------------
int gcd(int a, int b) 
{
  int remainder = -1;
#if CGL_DEBUG>1
  printf("gcd of %d and %d\n",a,b);
  int nLoop=0;
#endif
  // make sure a<=b (will always remain so)
  if(a > b) {
    // Swap a and b
    int temp = a;
    a = b;
    b = temp;
  }
  // if zero then gcd is nonzero (zero may occur in rhs of packed)
  if (!a) {
    if (b) {
      return b;
    } else {
      printf("**** gcd given two zeros!!\n");
      abort();
    }
  }
  while (remainder) {

#if CGL_DEBUG>1
    nLoop++;
    if (nLoop>50) {
      abort();
      return -1;
    }
#endif
    remainder = b % a;
    b = a;
    a = remainder;
  }
#if CGL_DEBUG>1
  printf("=> %d\n",b);
#endif
  return b;
}

//-------------------------------------------------------------------
// Returns the nearest rational with denominator < maxDenominator
//-------------------------------------------------------------------
typedef struct {
  int numerator;
  int denominator;
} Rational;
inline Rational nearestRational(double value, int maxDenominator) 
{
  Rational tryThis;
  Rational tryA;
  Rational tryB;
  double integerPart;

#if CGL_DEBUG>1
  printf("Rational of %g is ",value);
#endif
  int nLoop=0;

  tryA.numerator=0;
  tryA.denominator=1;
  tryB.numerator=1;
  tryB.denominator=0;

  if (fabs(value)<1.0e-10)
    return tryA;
  integerPart = floor(value);
  value -= integerPart;
  tryThis.numerator = tryB.numerator* (int ) integerPart + tryA.numerator;
  tryThis.denominator = tryB.denominator* (int ) integerPart + tryA.denominator;
  tryA = tryB;
  tryB = tryThis;

  while (value>1.0e-10 && tryB.denominator <=maxDenominator) {
    nLoop++;
    if (nLoop>50) {
      Rational bad;
      bad.numerator=-1;
      bad.denominator=-1;
#if CGL_DEBUG>1
      printf(" *** bad rational\n");
#endif
      return bad;
    }
    value = 1.0/value;
    integerPart = floor(value+1.0e-10);
    value -= integerPart;
    tryThis.numerator = tryB.numerator* (int ) integerPart + tryA.numerator;
    tryThis.denominator = tryB.denominator* (int ) integerPart + tryA.denominator;
    tryA = tryB;
    tryB = tryThis;
  }
  if (tryB.denominator <= maxDenominator) {
#if CGL_DEBUG>1
    printf("%d/%d\n",tryB.numerator,tryB.denominator);
#endif
    return tryB;
  } else {
#if CGL_DEBUG>1
    printf("%d/%d\n",tryA.numerator,tryA.denominator);
#endif
    return tryA;
  }
}

// Does actual work - returns number of cuts
int
CglGomory::generateCuts( const OsiRowCutDebugger * debugger, 
		     OsiCuts & cs,
		     const CoinPackedMatrix & columnCopy,
		     const double * objective, const double * colsol,
		     const double * colLower, const double * colUpper,
		     const double * rowLower, const double * rowUpper,
			 const char * intVar,
		    const CoinWarmStartBasis* warm) const
{
  int numberRows=columnCopy.getNumRows();
  int numberColumns=columnCopy.getNumCols(); 
  
  // check factorization is okay
  CoinFactorization factorization;
  // We can either set increasing rows so ...IsBasic gives pivot row
  // or we can just increment iBasic one by one
  // for now let ...iBasic give pivot row
  factorization.increasingRows(2);
  int status=-100;
  // probably could use pivotVariables from OsiSimplexModel
  int * rowIsBasic = new int[numberRows];
  int * columnIsBasic = new int[numberColumns];
  int i;
  int numberBasic=0;
  for (i=0;i<numberRows;i++) {
    if (warm->getArtifStatus(i) == CoinWarmStartBasis::basic) {
      rowIsBasic[i]=1;
      numberBasic++;
    } else {
      rowIsBasic[i]=-1;
    }
  }
  for (i=0;i<numberColumns;i++) {
    if (warm->getStructStatus(i) == CoinWarmStartBasis::basic) {
      columnIsBasic[i]=1;
      numberBasic++;
    } else {
      columnIsBasic[i]=-1;
    }
  }
  //returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
  while (status<-98) {
    status=factorization.factorize(columnCopy,
				   rowIsBasic, columnIsBasic);
    if (status==-99) factorization.areaFactor(factorization.areaFactor() * 2.0);
  } 
  if (status) {
    std::cout<<"Bad factorization of basis - status "<<status<<std::endl;
    return -1;
  }
  
  double bounds[2]={-DBL_MAX,0.0};
  int iColumn,iRow;

  // get row copy of matrix
  CoinPackedMatrix rowCopy =  columnCopy;
  rowCopy.reverseOrdering();
  const int * column = rowCopy.getIndices();
  const int * rowStart = rowCopy.getVectorStarts();
  const int * rowLength = rowCopy.getVectorLengths(); 
  const double * rowElements = rowCopy.getElements();
  const int * row = columnCopy.getIndices();
  const int * columnStart = columnCopy.getVectorStarts();
  const int * columnLength = columnCopy.getVectorLengths(); 
  const double * columnElements = columnCopy.getElements();

  // we need to do book-keeping for variables at ub
  double tolerance = 1.0e-7;
  bool * swap= new bool [numberColumns];
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (columnIsBasic[iColumn]<0&&
	colUpper[iColumn]-colsol[iColumn]<=tolerance) {
      swap[iColumn]=true;
    } else {
      swap[iColumn]=false;
    }
  }

  // get row activities (could use solver but lets do here )
  double * rowActivity = new double [numberRows];
  CoinFillN(rowActivity,numberRows,0.0);
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double value = colsol[iColumn];
    int k;
    for (k=columnStart[iColumn];k<columnStart[iColumn]+columnLength[iColumn];k++) {
      iRow = row[k];
      rowActivity[iRow] += columnElements[k]*value;
    }
  }
#ifdef FULL_DEBUG
  // check solution matches
  {
    OsiSimplexSolution test;
    test.borrowModel(columnCopy, colLower, colUpper,
		     objective, 
		     rowLower, rowUpper);
    test.factorize(factorization, *warm);
    test.getSolution(rowActivity, colsol);
    assert (test.largestSolutionError()<1.0e-5);
    assert (test.primalFeasible());
  }
#endif
  /* we need to mark rows:
     0) equality
     1) slack at lb (activity at ub)
     2) slack at ub (activity at lb)
     and 4 bit is set if slack must be integer

  */
  int * rowType = new int[numberRows];
  for (iRow=0;iRow<numberRows;iRow++) {
    if (rowIsBasic[iRow]<0&&rowUpper[iRow]>rowLower[iRow]+1.0e-7) {
      int type=0;
      double rhs=0.0;
      if (rowActivity[iRow]>rowUpper[iRow]-1.0e-7) {
	type=1;
	rhs=rowUpper[iRow];
      } else if (rowActivity[iRow]<rowLower[iRow]+1.0e-7) {
	type=2;
	rhs=rowLower[iRow];
      } else {
	// wrong - but probably large rhs
#ifdef CGL_DEBUG
	assert (min(rowUpper[iRow]-rowActivity[iRow],
		    rowActivity[iRow]-rowUpper[iRow])<1.0e-7);
#else
	continue;
#endif
      }
      if (above_integer(rhs)<1.0e-10) {
	// could be integer slack
	bool allInteger=true;
	int k;
	for (k=rowStart[iRow];
	     k<rowStart[iRow]+rowLength[iRow];k++) {
	  int iColumn=column[k];
	  if (!intVar[iColumn]||above_integer(rowElements[k])>1.0e-10) {
	    // not integer slacks
	    allInteger=false;
	    break;
	  }
	}
	if (allInteger) {
	  type |= 4;
	}
      }
      rowType[iRow]=type;
    } else {
      // row is equality or basic
      rowType[iRow]=0;
    }
  }

  // two vectors for updating (one is work)
  CoinIndexedVector work;
  CoinIndexedVector array;
  // make sure large enough
  work.reserve(numberRows);
  array.reserve(numberRows);
  int * arrayRows = array.getIndices();
  double * arrayElements = array.denseVector();
  int numberAdded=0;
  // we also need somewhere to accumulate cut
  CoinIndexedVector cutVector;
  cutVector.reserve(numberColumns+1);
  int * cutIndex = cutVector.getIndices();
  double * cutElement = cutVector.denseVector(); 
  // and for packed form (as not necessarily in order)
  double * packed = new double[numberColumns+1];

  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    double reducedValue=above_integer(colsol[iColumn]);;
    int iBasic=columnIsBasic[iColumn];
    double ratio=reducedValue/(1.0-reducedValue);
    if (iBasic>=0) {
      int j;
#ifdef CGL_DEBUG
      {
	// put column into array
	array.setVector(columnLength[iColumn],row+columnStart[iColumn],
			columnElements+columnStart[iColumn]);
	// get column in tableau
	factorization.updateColumn ( &work, &array );
	int nn=0;
	int numberInArray=array.getNumElements();
	for (j=0;j<numberInArray;j++) {
	  int indexValue=arrayRows[j];
	  double value=arrayElements[indexValue];
	  if (fabs(value)>1.0e-6) {
	    assert (fabs(value-1.0)<1.0e-7);
	    assert (indexValue==iBasic);
	    nn++;
	  }
	}
	assert (nn==1);
      }
#endif
      if(intVar[iColumn]&&reducedValue<1.0-away_&&reducedValue>away_) {
#ifdef CGL_DEBUG
	cutVector.checkClear();
#endif
	// get row of tableau
	double one =1.0;
	array.setVector(1,&iBasic,&one);
	int numberNonInteger=0;
	// get pi
	factorization.updateColumnTranspose ( &work, &array );
	int numberInArray=array.getNumElements();
	// check pivot on iColumn
	{
	  double value=0.0;
	  int k;
	  // add in row of tableau
	  for (k=columnStart[iColumn];
	       k<columnStart[iColumn]+columnLength[iColumn];k++) {
	    iRow = row[k];
	    value += columnElements[k]*arrayElements[iRow];
	  }
	  // should be 1
#ifdef CGL_DEBUG
	  assert (fabs(value-1.0) < 1.0e-7);
#endif
	}
	//reducedValue=colsol[iColumn];
	// coding from pg 130 of Wolsey 
	// adjustment to rhs
	double rhs=0.0;
	int number=0;
	// columns
	for (j=0;j<numberColumns;j++) {
	  double value=0.0;
	  if (columnIsBasic[j]<0&&
	      colUpper[j]>colLower[j]+1.0e-8) {
	    int k;
	    // add in row of tableau
	    for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
	      iRow = row[k];
	      value += columnElements[k]*arrayElements[iRow];
	    }
	    if (fabs(value)<1.0e-12) {
	      // small value
	      continue;
	    } else {
#if CGL_DEBUG>1
	      printf("for basic %d, column %d has alpha %g, colsol %g\n",
		     iColumn,j,value,colsol[j]);
#endif
	      // deal with bounds
	      if (swap[j]) {
		//reducedValue -= value*colUpper[j];
		// negate
		value = - value;
	      } else {
		//reducedValue -= value*colLower[j];
	      }
	      double coefficient;
	      if (intVar[j]) {
		// integer
		coefficient = above_integer(value);
		if (coefficient > reducedValue) {
		  coefficient = ratio * (1.0-coefficient);
		} 
	      } else {
		// continuous
		numberNonInteger++;
		if (value > 0.0) {
		  coefficient = value;
		} else {
		  //??? sign wrong in book
		  coefficient = -ratio*value;
		}
	      }
	      if (swap[j]) {
		// negate
		coefficient = - coefficient;
		rhs += colUpper[j]*coefficient;
	      } else {
		rhs += colLower[j]*coefficient;
	      }
	      if (fabs(coefficient)>= COIN_INDEXED_TINY_ELEMENT) {
		 cutElement[j] = coefficient;
	      } else {
		 cutElement[j] = 1.0e-100;
	      }
	      //	      cutElement[j]=coefficient;
	      cutIndex[number++]=j;
	    }
	  } else {
	    // basic
	    continue;
	  }
	}
	cutVector.setNumElements(number);
	//check will be cut
	//reducedValue=above_integer(reducedValue);
	rhs += reducedValue;
	double violation = reducedValue;
	if (violation>violationTolerance_) {
#ifdef CGL_DEBUG
	  std::cout<<"cut has violation of "<<violation
		   <<" value "<<colsol[iColumn]<<std::endl;
#endif
	  // now do slacks part
	  for (j=0;j<numberInArray;j++) {
	    iRow=arrayRows[j];
	    double value = arrayElements[iRow];
	    int type=rowType[iRow];
	    if (type&&fabs(value)>=1.0e-12) {
	      if ((type&1)==0) {
		// negate to get correct coefficient
		value = - value;
	      }
	      double coefficient;
	      if ((type&4)!=0) {
		// integer
		coefficient = above_integer(value);
		if (coefficient > reducedValue) {
		  coefficient = ratio * (1.0-coefficient);
		} 
	      } else {
		numberNonInteger++;
		// continuous
		if (value > 0.0) {
		  coefficient = value;
		} else {
		  coefficient = -ratio*value;
		}
	      }
	      if ((type&1)!=0) {
		// slack at ub - treat as +1.0
		rhs -= coefficient*rowUpper[iRow];
	      } else {
		// negate yet again
		coefficient = - coefficient;
		rhs -= coefficient*rowLower[iRow];
	      }
	      int k;
	      for (k=rowStart[iRow];
		   k<rowStart[iRow]+rowLength[iRow];k++) {
		int jColumn=column[k];
		double value=rowElements[k];
		cutVector.add(jColumn,-coefficient*value);
	      }
	    }
	  }
	  //check again and pack down
	  // also change signs
	  // also zero cutElement
	  double sum=0.0;
	  rhs = - rhs;
	  int n = cutVector.getNumElements();
	  number=0;
	  for (j=0;j<n;j++) {
	    int jColumn =cutIndex[j];
	    double value=-cutElement[jColumn];
	    cutElement[jColumn]=0.0;
	    if (fabs(value)>1.0e-12) {
	      sum+=value*colsol[jColumn];
	      packed[number]=value;
	      cutIndex[number++]=jColumn;
	    }
	  }
	  // say zeroed out
	  cutVector.setNumElements(0);
	  if (sum >rhs+0.9*violationTolerance_&&
	      fabs((sum-rhs)-violation)<1.0e-6) {
#ifdef CGL_DEBUG
#if CGL_DEBUG==1
	    if (number<10) {
#endif
	      for (j=0;j<number;j++) {
		std::cout<<" ["<<cutIndex[j]<<","<<packed[j]<<"]";
	      }
	      std::cout<<" <= "<<rhs<<std::endl;
#if CGL_DEBUG==1
	    }
#endif
#endif
	    if (!numberNonInteger&&number) {
#ifdef CGL_DEBUG
	      assert (sizeof(Rational)==sizeof(double));
#endif
	      Rational * cleaned = (Rational *) cutElement;
	      int * xInt = (int *) cutElement;
	      // cut should have an integer slack so try and simplify
	      // add in rhs and put in cutElements (remember to zero later)
	      cutIndex[number]=numberColumns+1;
	      packed[number]=rhs;
	      int lcm = 1;
	      
	      for (j=0;j<number+1;j++) {
		double value=above_integer(fabs(packed[j]));
		cleaned[j]=nearestRational(value,100000);
		if (cleaned[j].denominator<0) {
		  // bad rational
		  lcm=-1;
		  break;
		}
		int thisGcd = gcd(lcm,cleaned[j].denominator);
		// may need to check for overflow
		lcm /= thisGcd;
		lcm *= cleaned[j].denominator;
	      }
	      if (lcm>0) {
		double multiplier = lcm;
		int nOverflow = 0; 
		for (j=0; j<number+1;j++) {
		  double value = fabs(packed[j]);
		  double dxInt = value*multiplier;
		  xInt[j]= (int) (dxInt+0.5); 
#if CGL_DEBUG>1
		  printf("%g => %g   \n",value,dxInt);
#endif
		  if (dxInt>1.0e9||fabs(dxInt-xInt[j])> 1.0e-8) {
		    nOverflow++;
		    break;
		  }
		}

		if (nOverflow){
#ifdef CGL_DEBUG
		  printf("Gomory Scaling: Warning: Overflow detected \n");
#endif
		  numberNonInteger=-1;
		} else {
		  
		  // find greatest common divisor of the elements
		  int thisGcd = gcd(xInt[0],xInt[1]);
		  for (j=2;j<number+1;j++) {
		    thisGcd = gcd(thisGcd,xInt[j]);
		  }
		  
#if CGL_DEBUG>1
		  printf("The gcd of xInt is %i\n",thisGcd);    
#endif
		  
		  // construct new cut by dividing through by gcd and 
		  double minMultiplier=1.0e100;
		  double maxMultiplier=0.0;
		  for (j=0; j<number+1;j++) {
		    double old=packed[j];
		    if (old>0.0) {
		      packed[j]=xInt[j]/thisGcd;
		    } else {
		      packed[j]=-xInt[j]/thisGcd;
		    }
#if CGL_DEBUG>1
		    printf("%g => %g   \n",old,packed[j]);
#endif
		    if (packed[j]) {
		      if (fabs(packed[j])>maxMultiplier*fabs(old))
			maxMultiplier = packed[j]/old;
		      if (fabs(packed[j])<minMultiplier*fabs(old))
			minMultiplier = packed[j]/old;
		    }
		  }
		  rhs = packed[number];
		  double ratio=maxMultiplier/minMultiplier;
#ifdef CGL_DEBUG
		  printf("min, max multipliers - %g, %g\n",
			 minMultiplier,maxMultiplier);
#endif
		  assert(ratio>0.9999&&ratio<1.0001);
		}
	      }
	      // erase cutElement
	      CoinFillN(cutElement,number+1,0.0);
	    } else {
	      // relax rhs a tiny bit
	      rhs += 1.0e-8;
	      // relax if lots of elements for mixed gomory
	      if (number>=20) {
		rhs  += 1.0e-7*((double) (number/20));
	      }
	    }
	    if (number<limit_||!numberNonInteger) {
	      bounds[1]=rhs;
	      {
		OsiRowCut rc;
		rc.setRow(number,cutIndex,packed);
		rc.setLb(bounds[0]);
		rc.setUb(bounds[1]);   
		cs.insert(rc);
#ifdef CGL_DEBUG
	      if (!number)
		std::cout<<"******* Empty cut - infeasible"<<std::endl;
		if (debugger) 
		  assert(!debugger->invalidCut(rc)); 
#endif
	      }
	      numberAdded++;
	    } else {
#ifdef CGL_DEBUG
	      std::cout<<"cut has "<<number<<" entries - skipped"
		     <<std::endl;
	      if (!number)
		std::cout<<"******* Empty cut - infeasible"<<std::endl;
#endif
	    }
	  } else {
	    // why dropped?
#ifdef CGL_DEBUG
	    std::cout<<"first violation "<<violation<<" now "
		     <<sum-rhs<<" why?, rhs= "<<rhs<<std::endl;
	    
	    for (j=0;j<number;j++) {
	      int jColumn =cutIndex[j];
	      double value=packed[j];
	      std::cout<<"("<<jColumn<<","<<value<<","<<colsol[jColumn]
		       <<") ";
	    }
	    std::cout<<std::endl;
	    abort();
#endif
	  }
	} else {
	  // clean anyway
	  cutVector.clear();
	}
      }
    } else {
      // not basic
#if CGL_DEBUG>1
      {
	// put column into array
	array.setVector(columnLength[iColumn],row+columnStart[iColumn],
			columnElements+columnStart[iColumn]);
	// get column in tableau
	factorization.updateColumn ( &work, &array );
	int numberInArray=array.getNumElements();
	printf("non-basic %d\n",iColumn);
	for (j=0;j<numberInArray;j++) {
	  int indexValue=arrayRows[j];
	  double value=arrayElements[indexValue];
	  if (fabs(value)>1.0e-6) {
	    printf("%d %g\n",indexValue,value);
	  }
	}
      }
#endif
    }
  }

  delete [] rowActivity;
  delete [] swap;
  delete [] rowType;
  delete [] packed;
  delete [] rowIsBasic;
  delete [] columnIsBasic;
  return numberAdded;
}
// Limit stuff
void CglGomory::setLimit(int limit)
{
  if (limit>0)
    limit_=limit;
}
int CglGomory::getLimit() const
{
  return limit_;
}

// Away stuff
void CglGomory::setAway(double value)
{
  if (value>0.0&&value<=0.5)
    away_=value;
}
double CglGomory::getAway() const
{
  return away_;
}

// ViolationTolerance stuff
void CglGomory::setViolationTolerance(double value)
{
  if (value>0.0&&value<=0.5)
    violationTolerance_=value;
}
double CglGomory::getViolationTolerance() const
{
  return violationTolerance_;
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglGomory::CglGomory ()
:
CglCutGenerator(),
violationTolerance_(0.00001),
away_(0.05),
limit_(50)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglGomory::CglGomory (const CglGomory & source) :
  CglCutGenerator(source),
  violationTolerance_(source.violationTolerance_),
  away_(source.away_),
  limit_(source.limit_)
{  
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglGomory::~CglGomory ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglGomory &
CglGomory::operator=(const CglGomory& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    violationTolerance_=rhs.violationTolerance_;
    away_=rhs.away_;
    limit_=rhs.limit_;
  }
  return *this;
}
