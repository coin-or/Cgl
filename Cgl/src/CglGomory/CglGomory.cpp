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
//#define CGL_DEBUG 1
//#ifdef NDEBUG
//#undef NDEBUG
//#endif
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinIndexedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinFactorization.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglGomory.hpp"
#include "CoinFinite.hpp"
//-------------------------------------------------------------------
// Generate Gomory cuts
//------------------------------------------------------------------- 
void CglGomory::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			     const CglTreeInfo info) const
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
      } else {
	intVar[i] = 0;
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

  generateCuts(debugger, cs, *si.getMatrixByCol(), *si.getMatrixByRow(),
	   si.getObjCoefficients(), si.getColSolution(),
	   si.getColLower(), si.getColUpper(), 
	   si.getRowLower(), si.getRowUpper(),
	   intVar,warm,info);

  delete warmstart;
  delete [] intVar;
}

// Returns value - floor but allowing for small errors
inline double above_integer(double value) 
{
  double value2=floor(value);
  double value3=floor(value+0.5);
  if (fabs(value3-value)<1.0e-7*(fabs(value3)+1.0))
    return 0.0;
  return value-value2;
}
//-------------------------------------------------------------------
// Returns the greatest common denominator of two 
// positive integers, a and b, found using Euclid's algorithm 
//-------------------------------------------------------------------
static int gcd(int a, int b) 
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
                         const CoinPackedMatrix & rowCopy,
                         const double * objective, const double * colsol,
                         const double * colLower, const double * colUpper,
                         const double * rowLower, const double * rowUpper,
			 const char * intVar,
                         const CoinWarmStartBasis* warm,
                         const CglTreeInfo info) const
{
  // get limit on length of cut
  int limit = info.inTree ? limit_ : CoinMax(limit_,limitAtRoot_);
  int numberRows=columnCopy.getNumRows();
  int numberColumns=columnCopy.getNumCols(); 
  // Allow bigger length on initial matrix (if special setting)
  if (limit==512&&!info.inTree&&!info.pass)
    limit=1024;
  
  // Start of code to create a factorization from warm start (A) ====
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
#ifdef COIN_DEVELOP
    std::cout<<"Bad factorization of basis - status "<<status<<std::endl;
#endif
    return -1;
  }
  // End of creation of factorization (A) ====
  
  double relaxation = factorization.conditionNumber();
#ifdef COIN_DEVELOP
  if (relaxation>1.0e12)
    printf("condition %g\n",relaxation);
#endif
  relaxation *= 1.0e-18;
  double bounds[2]={-DBL_MAX,0.0};
  int iColumn,iRow;

  const int * column = rowCopy.getIndices();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  const int * rowLength = rowCopy.getVectorLengths(); 
  const double * rowElements = rowCopy.getElements();
  const int * row = columnCopy.getIndices();
  const CoinBigIndex * columnStart = columnCopy.getVectorStarts();
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
      if (rowActivity[iRow]>=rowUpper[iRow]-1.0e-7) {
	type=1;
	rhs=rowUpper[iRow];
      } else if (rowActivity[iRow]<=rowLower[iRow]+1.0e-7) {
	type=2;
	rhs=rowLower[iRow];
      } else {
	// wrong - but probably large rhs
	rowType[iRow]=type;
#ifdef CGL_DEBUG
	assert (min(rowUpper[iRow]-rowActivity[iRow],
		    rowActivity[iRow]-rowUpper[iRow])<1.0e-7);
	//abort();
        continue;
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

  // Start of code to create work arrays for factorization (B) ====
  // two vectors for updating (one is work)
  CoinIndexedVector work;
  CoinIndexedVector array;
  // make sure large enough
  work.reserve(numberRows);
  array.reserve(numberRows);
  int * arrayRows = array.getIndices();
  double * arrayElements = array.denseVector();
  // End of code to create work arrays (B) ====

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
    // This returns pivot row for columns or -1 if not basic (C) ====
    int iBasic=columnIsBasic[iColumn];
    double ratio=reducedValue/(1.0-reducedValue);
    if (iBasic>=0) {
  // Debug code below computes tableau column of basic ====
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
	  if (fabs(value)>1.0e-5) {
	    assert (fabs(value-1.0)<1.0e-7);
	    assert (indexValue==iBasic);
	    nn++;
	  }
	}
	assert (nn==1);
	array.clear();
	work.checkClear();
      }
#endif
      array.clear();
      if(intVar[iColumn]&&reducedValue<1.0-away_&&reducedValue>away_) {
#ifdef CGL_DEBUG
	cutVector.checkClear();
#endif
	// get row of tableau
	double one =1.0;
	array.setVector(1,&iBasic,&one);
	int numberNonInteger=0;
	//Code below computes tableau row ====
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
	    // value is entry in tableau row end (C) ====
	    if (fabs(value)<1.0e-16) {
	      // small value
	      continue;
	    } else {
#if CGL_DEBUG>1
	      if (iColumn==52) printf("for basic %d, column %d has alpha %g, colsol %g\n",
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
#if CGL_DEBUG>1
	      if (iColumn==52) printf("%d value %g reduced %g int %d rhs %g swap %d\n",
		     j,value,reducedValue,intVar[j],rhs,swap[j]);
#endif
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
		 cutIndex[number++]=j;
		 // If too many - break from loop
		 if (number>limit) 
		   break;
	      }
	    }
	  } else {
	    // basic
	    continue;
	  }
	}
	cutVector.setNumElements(number);
	// If too many - just clear vector and skip
	if (number>limit) {
	  cutVector.clear();
	  continue;
	}
	//check will be cut
	//reducedValue=above_integer(reducedValue);
	rhs += reducedValue;
	double violation = reducedValue;
#ifdef CGL_DEBUG
	std::cout<<"cut has violation of "<<violation
		 <<" value "<<colsol[iColumn]<<std::endl;
#endif
	// now do slacks part
	for (j=0;j<numberInArray;j++) {
	  iRow=arrayRows[j];
	  double value = arrayElements[iRow];
	  int type=rowType[iRow];
	  if (type&&fabs(value)>=1.0e-16) {
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
	      // negate yet again ?
	      coefficient = - coefficient;
	      rhs -= coefficient*rowLower[iRow];
	    }
	    int k;
	    for (k=rowStart[iRow];
		 k<rowStart[iRow]+rowLength[iRow];k++) {
	      int jColumn=column[k];
	      double value=rowElements[k];
	      cutVector.quickAdd(jColumn,-coefficient*value);
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
	  if (fabs(value)>1.0e-8) {
	    sum+=value*colsol[jColumn];
	    packed[number]=value;
	    cutIndex[number++]=jColumn;
          } else {
            // small - adjust rhs if rhs reasonable
            if (value>0.0&&colLower[jColumn]>-1000.0) {
              rhs -= value*colLower[jColumn];
            } else if (value<0.0&&colUpper[jColumn]<1000.0) {
              rhs -= value*colUpper[jColumn];
            } else if (fabs(value)>1.0e-13) {
              // take anyway
              sum+=value*colsol[jColumn];
              packed[number]=value;
              cutIndex[number++]=jColumn;
            } 
          }
	}
	// say zeroed out
	cutVector.setNumElements(0);
	if (sum >rhs+0.9*away_&&
	    fabs((sum-rhs)-violation)<1.0e-6) {
	  //#ifdef CGL_DEBUG
#ifdef CGL_DEBUG
#if CGL_DEBUG<=1
	  if (number<=-10) {
#endif
	    for (j=0;j<number;j++) {
	      std::cout<<" ["<<cutIndex[j]<<","<<packed[j]<<"]";
	    }
	    std::cout<<" <= "<<rhs<<std::endl;
#if CGL_DEBUG<=1
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
	    int numberNonSmall=0;
	    int lcm = 1;
	    
	    for (j=0;j<number+1;j++) {
	      double value=above_integer(fabs(packed[j]));
	      if (fabs(value)<1.0e-4) {
		// too small
		continue;
	      } else {
		numberNonSmall++;
	      }
	      
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
	    if (lcm>0&&numberNonSmall) {
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
		j=0;
		while (!xInt[j])
		  j++; // skip zeros
		int thisGcd = gcd(xInt[j],xInt[j+1]);
		j++;
		for (;j<number+1;j++) {
		  if (xInt[j])
		    thisGcd = gcd(thisGcd,xInt[j]);
		}
#if 0
		// Check nothing too illegal - FIX this
		for (j=0;j<number+1;j++) {
		  double old = lcm*packed[j];
		  int newOne;
		  if (old>0.0)
		    newOne=xInt[j]/thisGcd;
		  else
		    newOne=-xInt[j]/thisGcd;
		  if (fabs(((double) newOne)-old)>
		      1.0e-10*(fabs(newOne)+1.0)) {
		    // say no good - first see if happens
		    printf("Fix this test 456 - just skip\n");
		    abort();
		  }
		} 
#endif		  
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
#ifdef CGL_DEBUG
		printf("min, max multipliers - %g, %g\n",
		       minMultiplier,maxMultiplier);
#endif
		assert(maxMultiplier/minMultiplier>0.9999&&maxMultiplier/minMultiplier<1.0001);
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
	  // Take off tiny elements
	  // for first pass reject
#define TINY_ELEMENT 1.0e-12
	  {
	    int i,number2=number;
	    number=0;
	    double largest=0.0;
	    double smallest=1.0e30;
	    for (i=0;i<number2;i++) {
	      double value=fabs(packed[i]);
	      if (value<TINY_ELEMENT) {
		int iColumn = cutIndex[i];
		if (colUpper[iColumn]-colLower[iColumn]<10.0) {
		  // weaken cut
		  if (packed[i]>0.0) 
		    rhs -= value*colLower[iColumn];
		  else
		    rhs += value*colUpper[iColumn];
		} else {
		  // throw away
		  number=limit+1;
		  numberNonInteger=1;
		  break;
		}
	      } else {
		int iColumn = cutIndex[i];
		if (colUpper[iColumn]!=colLower[iColumn]) {
		  largest=max(largest,value);
		  smallest=min(smallest,value);
		  cutIndex[number]=cutIndex[i];
		  packed[number++]=packed[i];
		} else {
		  // fixed so subtract out
		  rhs -= packed[i]*colLower[iColumn];
		}
	      }
	    }
	    if (largest>1.0e9*smallest) {
	      number=limit+1; //reject
	      numberNonInteger=1;
	    }
	  }
	  if (number<limit||!numberNonInteger) {
	    bounds[1]=rhs;
	    if (number>50&&numberNonInteger)
	      bounds[1] = bounds[1]+1.0e-6+1.0e-8*fabs(rhs); // weaken
	    if (number>5&&numberNonInteger&&relaxation>1.0e-20) {
	      //printf("relaxing rhs by %g\n",relaxation*fabs(rhs));
	      bounds[1] = bounds[1]+relaxation*fabs(rhs); // weaken
	    }
	    {
	      OsiRowCut rc;
	      rc.setRow(number,cutIndex,packed,false);
	      rc.setLb(bounds[0]);
	      rc.setUb(bounds[1]);   
	      cs.insert(rc);
              //rc.print();
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
	  //abort();
#endif
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
	for (int j=0;j<numberInArray;j++) {
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
// Limit stuff at root
void CglGomory::setLimitAtRoot(int limit)
{
  if (limit>0)
    limitAtRoot_=limit;
}
int CglGomory::getLimitAtRoot() const
{
  return limitAtRoot_;
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

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglGomory::CglGomory ()
:
CglCutGenerator(),
away_(0.05),
limit_(50),
limitAtRoot_(0)
{

}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglGomory::CglGomory (const CglGomory & source) :
  CglCutGenerator(source),
  away_(source.away_),
  limit_(source.limit_),
  limitAtRoot_(source.limitAtRoot_)
{  
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglGomory::clone() const
{
  return new CglGomory(*this);
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
    away_=rhs.away_;
    limit_=rhs.limit_;
    limitAtRoot_=rhs.limitAtRoot_;
  }
  return *this;
}
// Returns true if needs optimal basis to do cuts
bool 
CglGomory::needsOptimalBasis() const
{
  return true;
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
                         const CoinWarmStartBasis* warm,
                         const CglTreeInfo info) const
{
  CoinPackedMatrix rowCopy;
  rowCopy.reverseOrderedCopyOf(columnCopy);
  return generateCuts( debugger, cs, columnCopy, rowCopy,
		       objective, colsol, colLower, colUpper,
		       rowLower, rowUpper, intVar, warm, info);
}
// Create C++ lines to get to current state
std::string
CglGomory::generateCpp( FILE * fp) 
{
  CglGomory other;
  fprintf(fp,"0#include \"CglGomory.hpp\"\n");
  fprintf(fp,"3  CglGomory gomory;\n");
  if (limit_!=other.limit_)
    fprintf(fp,"3  gomory.setLimit(%d);\n",limit_);
  else
    fprintf(fp,"4  gomory.setLimit(%d);\n",limit_);
  if (limitAtRoot_!=other.limitAtRoot_)
    fprintf(fp,"3  gomory.setLimitAtRoot(%d);\n",limitAtRoot_);
  else
    fprintf(fp,"4  gomory.setLimitAtRoot(%d);\n",limitAtRoot_);
  if (away_!=other.away_)
    fprintf(fp,"3  gomory.setAway(%g);\n",away_);
  else
    fprintf(fp,"4  gomory.setAway(%g);\n",away_);
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  gomory.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  gomory.setAggressiveness(%d);\n",getAggressiveness());
  return "gomory";
}
