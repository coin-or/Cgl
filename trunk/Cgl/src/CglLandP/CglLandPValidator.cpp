// Copyright (C) 2005, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     11/22/05
//---------------------------------------------------------------------------
#include "CglLandPValidator.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiRowCut.hpp"

std::vector<std::string> CglValidator::rejections_;

/** Clean an OsiCut 
\return 1 if min violation is too small
\return 2 if small coefficient can not be removed
\return 3 if dynamic is too big
\return 4 if too many non zero element*/
int 
CglValidator::cleanCut(OsiRowCut & aCut, const double * solCut, const OsiSolverInterface &si, const CglParam& par) const
{
  /** Compute fill-in in si */
  int numcols = si.getNumCols();

  int maxNnz = (int) (maxFillIn_ * (double) numcols);
  
  double rhs = aCut.lb();
  assert (aCut.ub()> 1e50);
  
  CoinPackedVector *vec = const_cast<CoinPackedVector *>(&aCut.row());
  int * indices = vec->getIndices();
  double * elems = vec->getElements();
  int n = vec->getNumElements();

  /** First compute violation if it is too small exit */
  double violation = aCut.violated(solCut);
  if(violation < minViolation_)
    return 1;
  
  /** Now relax get dynamic and remove tiny elements */
  const double * colLower = si.getColLower();
  const double * colUpper = si.getColUpper();
  int offset = 0;
  rhs -= 1e-8;
  double smallest = 1e100;
  double biggest = 0;
  for(int i = 0 ; i < n ; i++)
  {
    double val = fabs(elems[i]);
    if(val <= par.getEPS())//try to remove coef
    {
      if(val>0 && val<1e-40) {offset++; continue; throw;}    
      if(val==0) {offset++; continue;}    

      int & iCol = indices[i];
      if(elems[i]>0. && colUpper[iCol] < 10000.)
	{
	  offset++;
	  rhs -= elems[i] * colUpper[iCol];
	  elems[i]=0;
	}
      else if (elems[i]<0. && colLower[iCol] > -10000.)
	{
	  offset++;
	  rhs -= elems[i] * colLower[iCol];
	  elems[i]=0.;
	}
      else
	{
#ifdef DEBUG
	  std::cout<<"Small coefficient : "<<elems[i]<<" bounds : ["<<colLower[iCol]<<", "<<colUpper[iCol]<<std::endl;
#endif
	  numRejected_[SmallCoefficient]++;
	  return SmallCoefficient;
	}
    }
  
    else //Not a small coefficient keep it
    {
      smallest = min(val,smallest);
      biggest = max (val,biggest);
      if(biggest > maxRatio_ * smallest)
	{
#ifdef DEBUG
	  std::cout<<"Whaooo "<<biggest/smallest<<std::endl;
#endif
	  numRejected_[BigDynamic]++;
	  return BigDynamic;
	}
      if(offset)//if offset is zero current values are ok otherwise translate
      {
        int i2 = i - offset;
        indices[i2] = indices[i];
        elems[i2] = elems[i];
      }
    }
  }
   if((n - offset) > maxNnz) 
     {
       numRejected_[DenseCut] ++;
       return DenseCut;
     }
  if(offset == n)
    {
      numRejected_[EmptyCut]++;
      return EmptyCut;
    }

  if(offset)
    vec->truncate(n - offset);
  
  indices = vec->getIndices();
  elems = vec->getElements();
  n = vec->getNumElements();
  
  aCut.setLb(rhs);
  violation = aCut.violated(solCut);
  if(violation < minViolation_)
    {
      numRejected_[SmallViolation]++;
      return SmallViolation;
    }

  return NoneAccepted;
}

/**Clean cut 2, different algorithm. First check the dynamic of the cut if < maxRatio scale to a biggest coef of 1
   otherwise scale it so that biggest coeff is 1 and try removing tinys ( < 1/maxRatio) either succeed or fail */
int 
CglValidator::cleanCut2(OsiRowCut & aCut, const double * solCut, const OsiSolverInterface &si, const CglParam &par) const
{
  /** Compute fill-in in si */
  int numcols = si.getNumCols();
  // int numrows = si.getNumRows();
  
  
  int maxNnz = (int) ( maxFillIn_ * (double) numcols);
  
  double rhs = aCut.lb();
  assert (aCut.ub()> 1e50);
  
  CoinPackedVector *vec = const_cast<CoinPackedVector *>(&aCut.row());
//  vec->sortIncrIndex();

  int * indices = vec->getIndices();
  double * elems = vec->getElements();
  int n = vec->getNumElements();
  if(n==0)
    {
      numRejected_[EmptyCut]++;
      return EmptyCut;
    }
  /** First compute violation if it is too small exit */
  double violation = aCut.violated(solCut);
  if(violation < minViolation_)
    return 1;
  
  /** Now relax get dynamic and remove tiny elements */
  const double * colLower = si.getColLower();
  const double * colUpper = si.getColUpper();
  int offset = 0;
  rhs -= 1e-10;
  double smallest = fabs(rhs);
  double biggest = smallest;
  double veryTiny = 1e-15;
  for(int i = 0 ; i < n ; i++)
  {
    double val = fabs(elems[i]);
    if (val > veryTiny)//tiny should be very very small
    {
      smallest = min(val,smallest);
      biggest = max (val,biggest);
    }
  }

  if(biggest > 1e9)
	{
#ifdef DEBUG
	  std::cout<<"Whaooo "<<biggest/smallest<<std::endl;
#endif
	  numRejected_[BigDynamic]++;
	  return BigDynamic;
	}

  //rescale the cut so that biggest is 1e1.
  double toBeBiggest = 1.;
  rhs *= (toBeBiggest / biggest);
  toBeBiggest /= biggest;
  for(int i = 0 ; i < n ; i++)
  {
    elems[i] *= toBeBiggest;
  }
  if(biggest > maxRatio_ * smallest)//we have to remove some small coefficients
  {
    double myTiny = biggest * toBeBiggest / maxRatio_;
    veryTiny *= toBeBiggest ;
    for(int i = 0 ; i < n ; i++)
    {
      double val = fabs(elems[i]);
      if(val < myTiny)
	    {
 	      if(val< veryTiny) {offset++; continue;}
 	      int & iCol = indices[i];
 	      if(elems[i]>0. && colUpper[iCol] < 1000.)
        {
          offset++;
          rhs -= elems[i] * colUpper[iCol];
          elems[i]=0;
        }
      else if (elems[i]<0. && colLower[iCol] > -1000.)
        {
          offset++;
          rhs -= elems[i] * colLower[iCol];
          elems[i]=0.;
        }
 	      else
        {
          numRejected_[SmallCoefficient]++;
          return SmallCoefficient;
        }
 	    }
      else //Not a small coefficient keep it
 	    {
 	      if(offset)//if offset is zero current values are ok
		{
		  int i2 = i - offset;
		  indices[i2] = indices[i];
		  elems[i2] = elems[i];
		}
 	    }
    }
  }
   if((n - offset) > maxNnz) 
   {
     numRejected_[DenseCut] ++;
     return DenseCut;
   }


  if(offset)
    vec->truncate(n - offset);

  if(vec->getNumElements() == 0 )
    {
      numRejected_[EmptyCut]++;
      return EmptyCut;
    }

  /** recheck violation */
  aCut.setLb(rhs);
  violation = aCut.violated(solCut);
  if(violation < minViolation_)
  {
    numRejected_[SmallViolation]++;
    return SmallViolation;
  }
  assert(fabs(rhs)<1e09);
  return NoneAccepted;
}


/** Constructor with default values */
CglValidator::CglValidator(double maxFillIn,
                           double maxRatio, 
                           double minViolation,
                           bool scale):
  maxFillIn_(maxFillIn), 
  maxRatio_(maxRatio), 
  minViolation_(minViolation),
  scale_(scale),
  numRejected_(DummyEnd,0)
{
  fillRejectionReasons();
}


void 
CglValidator::fillRejectionReasons()
{
  if (rejections_.size() == 0){
    rejections_.resize(DummyEnd) ;
    rejections_[NoneAccepted] = "Cut was accepted";
    rejections_[SmallViolation] = "Violation of the cut is too small ";
    rejections_[SmallCoefficient] = "There is a small coefficient we can not get rid off.";
    rejections_[BigDynamic] = "Dynamic of coefficinet is too important. ";
    rejections_[DenseCut] = "Cut is too dense.";
    rejections_[EmptyCut] = "Cleaned cut is empty";
  }
}

