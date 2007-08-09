
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
//#define PRINT_DEBUG
//#define CGL_DEBUG 1
//#undef NDEBUG
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinFinite.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CglProbing.hpp"
//#define PROBING_EXTRA_STUFF true
#define PROBING_EXTRA_STUFF false
typedef struct {double infeasibility;int sequence;} double_int_pair;
class double_int_pair_compare {
public:
  bool operator() (double_int_pair x , double_int_pair y) const
  {
    return ( x.infeasibility < y.infeasibility);
  }
};
#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__CYGWIN32__)|| (__GNUC__ <3)
// Seems to bug in cygwin
#define GNU_OLDWAY
#endif
#ifndef GNU_OLDWAY
// Use a hash set to find duplicates
#include <ext/hash_map>
#include <ext/hash_set>
using namespace __gnu_cxx;
class row_cut_compare {
public:
  bool operator() (const OsiRowCut2 & x, const OsiRowCut2 & y) const
  {
    int xN =x.row().getNumElements();
    int yN =y.row().getNumElements();
    bool identical=false;
    if (xN==yN) {
      double xLb = x.lb();
      double xUb = x.ub();
      double yLb = y.lb();
      double yUb = y.ub();
      if (fabs(xLb-yLb)<1.0e-8&&fabs(xUb-yUb)<1.0e-8) {
        const int * xIndices = x.row().getIndices();
        const double * xElements = x.row().getElements();
        const int * yIndices = y.row().getIndices();
        const double * yElements = y.row().getElements();
        int j;
        for( j=0;j<xN;j++) {
          if (xIndices[j]!=yIndices[j])
            break;
          if (fabs(xElements[j]-yElements[j])>1.0e-12)
            break;
        }
        identical =  (j==xN);
      }
    }
    return identical;
  }
};
static double multiplier[] = {1.23456789e2,-9.87654321};
class row_cut_hash {
public:
  size_t operator () (const OsiRowCut2 & x) const
  {
    int xN =x.row().getNumElements();
    double xLb = x.lb();
    double xUb = x.ub();
    const int * xIndices = x.row().getIndices();
    const double * xElements = x.row().getElements();
    size_t hashValue;
    double value=1.0;
    if (xLb>-1.0e10)
      value += xLb*multiplier[0];
    if (xUb<1.0e10)
      value += xUb*multiplier[1];
    for( int j=0;j<xN;j++) {
      int xColumn = xIndices[j];
      double xValue = xElements[j];
      int k=(j&1);
      value += (j+1)*multiplier[k]*(xColumn+1)*xValue;
    }
    // should be compile time but too lazy for now
    if (sizeof(value)>sizeof(hashValue)) {
      assert (sizeof(value)==2*sizeof(hashValue));
      unsigned int * xx = (unsigned int *) (&value);
      xx[0] += xx[1];
      hashValue = (size_t) xx[0];
    } else {
      assert (sizeof(value)==sizeof(hashValue));
      hashValue = (size_t) value;
    }
    return hashValue;
  }
};
class row_cut {
public:

  row_cut(int nRows )
  {
    numberCuts_=0;
    int maxRowCuts = 2*nRows + 200;
    size_=maxRowCuts;
  };
  ~row_cut()
  {
    //hash_set<OsiRowCut2  , row_cut_hash, row_cut_compare>::const_iterator i;
    // already done by addCuts
    //for (i=rowCut_.begin();i!=rowCut_.end();i++)
    //delete *i;
    rowCut_.clear();
  }
  int numberCuts() const
  { return numberCuts_;};
  inline bool outOfSpace() const
  { return size_==numberCuts_;};
  hash_set<OsiRowCut2  , row_cut_hash, row_cut_compare> rowCut_;
  int size_;
  int numberCuts_;
  // Return 0 if added, 1 if not, -1 if not added because of space
  int addCutIfNotDuplicate(OsiRowCut & cut,int whichRow=-1)
  {
    if (numberCuts_<size_) {
      numberCuts_++;
      double newLb = cut.lb();
      double newUb = cut.ub();
      CoinPackedVector vector = cut.row();
      int numberElements =vector.getNumElements();
      int * newIndices = vector.getIndices();
      double * newElements = vector.getElements();
      CoinSort_2(newIndices,newIndices+numberElements,newElements);
      OsiRowCut2 newCut(whichRow);
      newCut.setLb(newLb);
      newCut.setUb(newUb);
      newCut.setRow(vector);
      hash_set<OsiRowCut2 , row_cut_hash, 
        row_cut_compare>::iterator find;
      find =  rowCut_.find(newCut);
      if (find == rowCut_.end()) {
        rowCut_.insert(newCut);
        return 0;
      } else {
        return 1;
      }
    } else {
      return -1;
    }
  }
  void addCuts(OsiCuts & cs, OsiRowCut ** whichRow)
  {
    if (!numberCuts_)
      return;
    hash_set<OsiRowCut2  , row_cut_hash, row_cut_compare>::const_iterator i;
    int numberCuts=cs.sizeRowCuts();
    int nRows = (size_-200)/2;
    if (numberCuts_<nRows) {
      for (i=rowCut_.begin();i!=rowCut_.end();i++) {
        cs.insert (*i);
        if (whichRow) {
          //cs.rowCutPtr(numberCuts)->print();
          int iRow= i->whichRow();
          if (iRow>=0&&iRow<nRows&&!whichRow[iRow]) {
            whichRow[iRow]=cs.rowCutPtr(numberCuts);
          }
        }
        numberCuts++;
      }
    } else {
      // just best
      double * effectiveness = new double[numberCuts_];
      int iCut=0;
      for (i=rowCut_.begin();i!=rowCut_.end();i++) {
        effectiveness[iCut++]=-i->effectiveness();
      }
      std::sort(effectiveness,effectiveness+iCut);
      double threshold = -1.0e20;
      if (iCut>nRows)
        threshold = effectiveness[nRows];
      for (i=rowCut_.begin();i!=rowCut_.end();i++) {
        if (i->effectiveness()>threshold) {
          cs.insert (*i);
          if (whichRow) {
            int iRow= i->whichRow();
            if (iRow>=0&&!whichRow[iRow])
              whichRow[iRow]=cs.rowCutPtr(numberCuts);;
          }
          numberCuts++;
        }
      }
      delete [] effectiveness;
    }
    rowCut_.clear();
    numberCuts_=0;
  }
};
#else
class row_cut {
public:

  row_cut(int nRows )
  {
    numberCuts_=0;
    int maxRowCuts = 2*nRows + 200;
    size_=maxRowCuts;
    rowCut_ = new  OsiRowCut2 * [size_];
    numberCuts_=0;
  };
  ~row_cut()
  {
    for (int i=0;i<numberCuts_;i++)
      delete rowCut_[i];
    delete [] rowCut_;
  }
  OsiRowCut2 * cut(int i) const
  { return rowCut_[i];};
  int numberCuts() const
  { return numberCuts_;};
  inline bool outOfSpace() const
  { return size_==numberCuts_;};
  OsiRowCut2 ** rowCut_;
  int size_;
  int numberCuts_;
  // Return 0 if added, 1 if not, -1 if not added because of space
  int addCutIfNotDuplicate(OsiRowCut & cut,int whichRow=-1)
  {
    if (numberCuts_<size_) {
      double newLb = cut.lb();
      double newUb = cut.ub();
      CoinPackedVector vector = cut.row();
      int numberElements =vector.getNumElements();
      int * newIndices = vector.getIndices();
      double * newElements = vector.getElements();
      CoinSort_2(newIndices,newIndices+numberElements,newElements);
      bool notDuplicate=true;
      for ( int i =0; i<numberCuts_;i++) {
        const OsiRowCut2 * cutPtr = rowCut_[i];
        if (cutPtr->row().getNumElements()!=numberElements)
          continue;
        if (fabs(cutPtr->lb()-newLb)>1.0e-12)
          continue;
        if (fabs(cutPtr->ub()-newUb)>1.0e-12)
          continue;
        const CoinPackedVector * thisVector = &(cutPtr->row());
        const int * indices = thisVector->getIndices();
        const double * elements = thisVector->getElements();
        int j;
        for(j=0;j<numberElements;j++) {
          if (indices[j]!=newIndices[j])
            break;
          if (fabs(elements[j]-newElements[j])>1.0e-12)
            break;
        }
        if (j==numberElements) {
          notDuplicate=false;
          break;
        }
      }
      if (notDuplicate) {
        OsiRowCut2 * newCutPtr = new OsiRowCut2(whichRow);
        newCutPtr->setLb(newLb);
        newCutPtr->setUb(newUb);
        newCutPtr->setRow(vector);
        rowCut_[numberCuts_++]=newCutPtr;
        return 0;
      } else {
        return 1;
      }
    } else {
      return -1;
    }
  }
  void addCuts(OsiCuts & cs, OsiRowCut ** whichRow)
  {
    int numberCuts=cs.sizeRowCuts();
    int nRows = (size_-200)/2;
    int i ;
    if (numberCuts_<nRows) {
      for (i=0;i<numberCuts_;i++) {
        cs.insert(*rowCut_[i]);
        if (whichRow) {
          int iRow= rowCut_[i]->whichRow();
          if (iRow>=0&&!whichRow[iRow])
            whichRow[iRow]=cs.rowCutPtr(numberCuts);;
        }
        numberCuts++;
      }
    } else {
      // just best
      double * effectiveness = new double[numberCuts_];
      int iCut=0;
      for (i=0;i<numberCuts_;i++) {
        effectiveness[iCut++]=-rowCut_[i]->effectiveness();
      }
      std::sort(effectiveness,effectiveness+numberCuts_);
      double threshold = effectiveness[nRows];
      for ( i=0;i<numberCuts_;i++) {
        if (rowCut_[i]->effectiveness()>threshold) {
          cs.insert(*rowCut_[i]);
          if (whichRow) {
            int iRow= rowCut_[i]->whichRow();
            if (iRow>=0&&!whichRow[iRow])
              whichRow[iRow]=cs.rowCutPtr(numberCuts);;
          }
          numberCuts++;
        }
      }
      delete[] effectiveness ;
    }
    for (i = 0 ; i < numberCuts_ ; i++)
    { delete rowCut_[i] ;
      rowCut_[i] = 0 ; }
    numberCuts_=0;
  }
};
#endif
// Adds in cut to list 
#ifdef CGL_DEBUG
// Checks bounds okay against debugger
static void checkBounds(const OsiRowCutDebugger * debugger,OsiColCut & cut)
{
  if (debugger) {
    // on optimal path
    const double * optimal = debugger->optimalSolution();
    int i;
    int nIndex;
    const double * values;
    const int * index;
    const CoinPackedVector & lbs = cut.lbs();
    values = lbs.getElements();
    nIndex = lbs.getNumElements();
    index = lbs.getIndices();
    for (i=0;i<nIndex;i++) {
      double value=values[i];
      int iColumn = index[i];
      printf("%d optimal %g lower %g\n",iColumn,optimal[iColumn],value);
      assert(value<=optimal[iColumn]+1.0e-5);
    }
    const CoinPackedVector & ubs = cut.ubs();
    values = ubs.getElements();
    nIndex = ubs.getNumElements();
    index = ubs.getIndices();
    for (i=0;i<nIndex;i++) {
      double value=values[i];
      int iColumn = index[i];
      printf("%d optimal %g upper %g\n",iColumn,optimal[iColumn],value);
      assert(value>=optimal[iColumn]-1.0e-5);
    }
  }
}
#endif
// This tightens column bounds (and can declare infeasibility)
// It may also declare rows to be redundant
int 
CglProbing::tighten(double *colLower, double * colUpper,
                    const int *column, const double *rowElements, 
                    const CoinBigIndex *rowStart, const int * rowLength,
                    double *rowLower, double *rowUpper, 
                    int nRows,int nCols,char * intVar,int maxpass,
                    double tolerance) const
{
  int i, j, k, kre;
  int krs;
  int dolrows;
  int iflagu, iflagl;
  int ntotal=0,nchange=1,jpass=0;
  double dmaxup, dmaxdown, dbound;
  int ninfeas=0;
  // For clique stuff
  double * cliqueMin=NULL;
  double * cliqueMax=NULL;
  // And second best ones
  double * cliqueMin2 = NULL;
  double * cliqueMax2 = NULL;
  if (cliqueRowStart_&&numberRows_&&cliqueRowStart_[numberRows_]) {
    cliqueMin = new double[nCols];
    cliqueMax = new double[nCols];
    cliqueMin2 = new double[nCols];
    cliqueMax2 = new double[nCols];
  }
  
  while(nchange) {
    int ilbred = 0; /* bounds reduced */
    int iubred = 0; /* bounds reduced */
    int nrwdrp = 0; /* redundant rows */
    if (jpass==maxpass) break;
    jpass++;
    dolrows = (jpass & 1) == 1;
    
    for (i = 0; i < nRows; ++i) {
      bool cliqueChanges=false;
      if (rowLower[i]>-1.0e20||rowUpper[i]<1.0e20) {
	iflagu = 0;
	iflagl = 0;
	dmaxup = 0.0;
	dmaxdown = 0.0;
	krs = rowStart[i];
	kre = rowStart[i]+rowLength[i];

	/* ------------------------------------------------------------*/
	/* Compute L(i) and U(i) */
	/* ------------------------------------------------------------*/
        if (!cliqueMin||i>=numberRows_||cliqueRowStart_[i]==cliqueRowStart_[i+1]) {
          // without cliques
          for (k = krs; k < kre; ++k) {
            double value=rowElements[k];
            j = column[k];
            if (value > 0.0) {
              if (colUpper[j] >= 1e15) {
                dmaxup = 1e31;
                ++iflagu;
              } else {
                dmaxup += colUpper[j] * value;
              }
              if (colLower[j] <= -1e15) {
                dmaxdown = -1e31;
                ++iflagl;
              } else {
                dmaxdown += colLower[j] * value;
              }
            } else if (value<0.0) {
              if (colUpper[j] >= 1e15) {
                dmaxdown = -1e31;
                ++iflagl;
              } else {
                dmaxdown += colUpper[j] * value;
              }
              if (colLower[j] <= -1e15) {
                dmaxup = 1e31;
                ++iflagu;
              } else {
                dmaxup += colLower[j] * value;
              }
            }
          }
        } else {
          // with cliques
          int nClique=0;
          int bias = cliqueRowStart_[i]-krs;
          double dmaxup2=0.0;
          double dmaxdown2=0.0;
          double sumZeroFixes=0.0;
          for (k = krs; k < kre; ++k) {
            double value=rowElements[k];
            j = column[k];
            int iClique = (int) (cliqueRow_[k+bias].sequence);
            bool oneFixes = (cliqueRow_[k+bias].oneFixes!=0);
            if (iClique>=numberColumns_||colUpper[j]==colLower[j]) {
              if (value > 0.0) {
                if (colUpper[j] >= 1e15) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colLower[j] * value;
                }
              } else if (value<0.0) {
                if (colUpper[j] >= 1e15) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colLower[j] * value;
                }
              }
            } else {
              // clique may help
              if (iClique>=nClique) {
                //zero out
                for (int j=nClique;j<=iClique;j++) {
                  cliqueMin[j]=0.0;
                  cliqueMax[j]=0.0;
                  cliqueMin2[j]=0.0;
                  cliqueMax2[j]=0.0;
                }
                nClique=iClique+1;
              }
              //  Update best and second best
              if (oneFixes) {
                if (value > 0.0) {
                  dmaxup2 += value;
                  cliqueMax2[iClique] = cliqueMax[iClique];
                  cliqueMax[iClique] = CoinMax(cliqueMax[iClique],value);
                } else if (value<0.0) {
                  dmaxdown2 +=  value;
                  cliqueMin2[iClique] = cliqueMin[iClique];
                  cliqueMin[iClique] = CoinMin(cliqueMin[iClique],value);
                }
              } else {
                sumZeroFixes += value;
                if (value > 0.0) {
                  dmaxup2 += value;
                  cliqueMin2[iClique] = cliqueMin[iClique];
                  cliqueMin[iClique] = CoinMin(cliqueMin[iClique],-value);
                } else if (value<0.0) {
                  dmaxdown2 +=  value;
                  cliqueMax2[iClique] = cliqueMax[iClique];
                  cliqueMax[iClique] = CoinMax(cliqueMax[iClique],-value);
                }
              }
            }
          }
          double dmaxup3 = dmaxup + sumZeroFixes;
          double dmaxdown3 = dmaxdown + sumZeroFixes;
          for (int iClique=0;iClique<nClique;iClique++) {
            dmaxup3 += cliqueMax[iClique];
            dmaxdown3 += cliqueMin[iClique];
          }
          dmaxup += dmaxup2;
          dmaxdown += dmaxdown2;
          assert (dmaxup3<=dmaxup+1.0e-8);
          assert (dmaxdown3>=dmaxdown-1.0e-8);
          if (dmaxup3<dmaxup-1.0e-8||dmaxdown3>dmaxdown+1.0e-8) {
            cliqueChanges=true;
            //printf("normal min/max %g , %g clique %g , %g\n",
            //     dmaxdown,dmaxup,dmaxdown3,dmaxup3);
            dmaxdown=dmaxdown3;
            dmaxup=dmaxup3;
          }
        }
	if (dmaxup <= rowUpper[i] + tolerance && dmaxdown >= rowLower[i] - tolerance) {
	  /*
	   * The sum of the column maxs is at most the row ub, and
	   * the sum of the column mins is at least the row lb;
	   * this row says nothing at all.
	   * I suspect that this corresponds to
	   * an implied column singleton in the paper (viii, on p. 325),
	   * where the singleton in question is the row slack.
	   */
	  ++nrwdrp;
	  rowLower[i]=-DBL_MAX;
	  rowUpper[i]=DBL_MAX;
	} else {
	  if (dmaxup < rowLower[i] -tolerance || dmaxdown > rowUpper[i]+tolerance) {
	    ninfeas++;
            assert (!cliqueChanges);
	    break;
	  }
	  /*        Finite U(i) */
	  /* -------------------------------------------------------------*/
	  /* below is deliberate mistake (previously was by chance) */
	  /*        never do both */
	  if (iflagu == 0 && rowLower[i] > 0.0 && iflagl == 0 && rowUpper[i] < 1e15) {
            if (dolrows) {
              iflagu = 1;
            } else {
              iflagl = 1;
            }
          }
          if (!cliqueChanges) {
            // without cliques
            if (iflagu == 0 && rowLower[i] > -1e15) {
              for (k = krs; k < kre; ++k) {
                double value=rowElements[k];
                j = column[k];
                if (value > 0.0) {
                  if (colUpper[j] < 1e15) {
                    dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
                    if (dbound > colLower[j] + 1.0e-8) {
                      /* we can tighten the lower bound */
                      /* the paper mentions this as a possibility on p. 227 */
                      colLower[j] = dbound;
                      ++ilbred;
                      
                      /* this may have fixed the variable */
                      /* I believe that this roughly corresponds to a
                       * forcing constraint in the paper (p. 226).
                       * If there is a forcing constraint (with respect
                       * to the original, untightened bounds), then in this 
                       * loop we will go through all the columns and fix
                       * each of them to their implied bound, rather than
                       * determining that the row as a whole is forced
                       * and just fixing them without doing computation for
                       * each column (as in the paper).
                       * By doing it this way, we can tighten bounds and
                       * get forcing constraints for free.
                       */
                      if (colUpper[j] - colLower[j] <= tolerance) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  }
                } else {
                  if (colLower[j] > -1e15) {
                    dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
                    if (dbound < colUpper[j] - 1.0e-8) {
                      colUpper[j] = dbound;
                      ++iubred;
                      if (colUpper[j] - colLower[j] <= tolerance) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  }
                }
              }
            }
            
            /* ----------------------------------------------------------------*/
            /*        Finite L(i) */
            /* ----------------------------------------------------------------*/
            if (iflagl == 0 && rowUpper[i] < 1e15) {
              for (k = krs; k < kre; ++k) {
                double value=rowElements[k];
                j = column[k];
                if (value < 0.0) {
                  if (colUpper[j] < 1e15) {
                    dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
                    if (dbound > colLower[j] + 1.0e-8) {
                      colLower[j] = dbound;
                      ++ilbred;
                      if (! (colUpper[j] - colLower[j] > tolerance)) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  } 
                } else {
                  if (colLower[j] > -1e15) {
                    dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
                    if (dbound < colUpper[j] - 1.0e-8) {
                      colUpper[j] = dbound;
                      ++iubred;
                      if (! (colUpper[j] - colLower[j] > tolerance)) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            // with cliques
            int bias = cliqueRowStart_[i]-krs;
            if (iflagu == 0 && rowLower[i] > -1e15) {
              for (k = krs; k < kre; ++k) {
                double value=rowElements[k];
                j = column[k];
                int iClique = (int) (cliqueRow_[k+bias].sequence);
                //bool oneFixes = (cliqueRow_[k+bias].oneFixes!=0);
                if (iClique>=numberColumns_) {
                  if (value > 0.0) {
                    if (colUpper[j] < 1e15) {
                      dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
                      if (dbound > colLower[j] + 1.0e-8) {
                        /* we can tighten the lower bound */
                        /* the paper mentions this as a possibility on p. 227 */
                        colLower[j] = dbound;
                        ++ilbred;
                        
                        /* this may have fixed the variable */
                        /* I believe that this roughly corresponds to a
                         * forcing constraint in the paper (p. 226).
                         * If there is a forcing constraint (with respect
                         * to the original, untightened bounds), then in this 
                         * loop we will go through all the columns and fix
                         * each of them to their implied bound, rather than
                         * determining that the row as a whole is forced
                         * and just fixing them without doing computation for
                         * each column (as in the paper).
                         * By doing it this way, we can tighten bounds and
                         * get forcing constraints for free.
                         */
                        if (colUpper[j] - colLower[j] <= tolerance) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowUpper[i]>1.0e20) {
                        // can we modify coefficient
                        if (dmaxdown+value>rowLower[i]+1.0e-8) {
                          assert (dmaxdown<rowLower[i]+1.0e-8);
                          double change = dmaxdown+value - rowLower[i];
                          double newValue = value - change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g\n",
                                 value,newValue);
                          // dmaxup -= change; 
                        }
#endif
                      }
                    }
                  } else {
                    if (colLower[j] > -1e15) {
                      dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
                      if (dbound < colUpper[j] - 1.0e-8) {
                        colUpper[j] = dbound;
                        ++iubred;
                        if (colUpper[j] - colLower[j] <= tolerance) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowUpper[i]>1.0e20) {
                        // can we modify coefficient
                        if (dmaxdown-value>rowLower[i]+1.0e-8) {
                          assert (dmaxdown<rowLower[i]+1.0e-8);
                          double change = dmaxdown-value-rowLower[i];
                          double newValue = value+change;
                          double newLower = rowLower[i]+change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g and lorow from %g to %g\n",
                                 value,newValue,rowLower[i],newLower);
                         // dmaxdown += change                          
                        }
#endif
                      }
                    }
                  }
                } else if (colUpper[j]>colLower[j]) {
                  // in clique
                  // adjustment
                  double dmaxup2=dmaxup;
                  assert (cliqueMax[iClique]>=0);
                  assert (cliqueMax2[iClique]>=0);
                  /* get max up if at other bound
                     May not go down at all but will not go up */
                  if (fabs(value)==fabs(cliqueMax[iClique]))
                    dmaxup2 -= cliqueMax[iClique]-cliqueMax2[iClique];
                  if (dmaxup2<rowLower[i]-1.0e-8) {
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if ( dmaxup<rowLower[i]-1.0e-8) {
                      ninfeas++;
                    } else {
                      if (value > 0.0) {
                        colLower[j] = 1.0;
                        ++ilbred;
                      } else {
                        colUpper[j] = 0.0;
                        ++iubred;
                      }
                    }
                  }
                }
              }
            }
            
            /* ----------------------------------------------------------------*/
            /*        Finite L(i) */
            /* ----------------------------------------------------------------*/
            if (iflagl == 0 && rowUpper[i] < 1e15) {
              for (k = krs; k < kre; ++k) {
                double value=rowElements[k];
                j = column[k];
                int iClique = (int) (cliqueRow_[k+bias].sequence);
                //bool oneFixes = (cliqueRow_[k+bias].oneFixes!=0);
                if (iClique>=numberColumns_) {
                  if (value < 0.0) {
                    if (colUpper[j] < 1e15) {
                      dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
                      if (dbound > colLower[j] + 1.0e-8) {
                        colLower[j] = dbound;
                        ++ilbred;
                        if (! (colUpper[j] - colLower[j] > tolerance)) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowLower[i]<-1.0e20) {
                        // can we modify coefficient
                        if (dmaxup+value<rowUpper[i]-1.0e-8) {
                          assert (dmaxup>rowUpper[i]-1.0e-8);
                          double change = dmaxup+value - rowUpper[i];
                          double newValue = value - change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g b\n",
                                 value,newValue);
                          // dmaxdown -= change; 
                        }
#endif
                      }
                    } 
                  } else {
                    if (colLower[j] > -1e15) {
                      dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
                      if (dbound < colUpper[j] - 1.0e-8) {
                        colUpper[j] = dbound;
                        ++iubred;
                        if (! (colUpper[j] - colLower[j] > tolerance)) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0*tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowLower[i]<-1.0e20) {
                        // can we modify coefficient
                        if (dmaxup-value<rowUpper[i]-1.0e-8) {
                          assert (dmaxup>rowUpper[i]-1.0e-8);
                          double change = dmaxup-value-rowUpper[i];
                          double newValue = value+change;
                          double newUpper = rowUpper[i]+change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g and uprow from %g to %g b\n",
                                 value,newValue,rowLower[i],newUpper);
                         // dmaxup += change                          
                        }
#endif
                      }
                    }
                  }
                } else if (colUpper[j]>colLower[j]) {
                  // in clique
                  // adjustment
                  double dmaxdown2=dmaxdown;
                  assert (cliqueMin[iClique]<=0);
                  assert (cliqueMin2[iClique]<=0);
                  /* get max down if this is at other bound
                     May not go up at all but will not go down */
                  if (fabs(value)==fabs(cliqueMin[iClique]))
                    dmaxdown2 -= cliqueMin[iClique]-cliqueMin2[iClique];
                  if (dmaxdown2>rowUpper[i]+1.0e-8) {
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if ( dmaxdown>rowUpper[i]+1.0e-8) {
                      ninfeas++;
                    } else {
                      if (value < 0.0) {
                        colLower[j] = 1.0;
                        ++ilbred;
                      } else {
                        colUpper[j] = 0.0;
                        ++iubred;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    for (j = 0; j < nCols; ++j) {
      if (intVar[j]) {
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  if (floor(colUpper[j]+1.0e-4)<colUpper[j]) 
	    nchange++;
	  // clean up anyway
	  colUpper[j]=floor(colUpper[j]+1.0e-4);
	  if (ceil(colLower[j]-1.0e-4)>colLower[j]) 
	    nchange++;
	  // clean up anyway
	  colLower[j]=ceil(colLower[j]-1.0e-4);
	  if (colUpper[j]<colLower[j]) {
	    /*printf("infeasible\n");*/
	    ninfeas++;
	  }
	}
      }
    }
    nchange=ilbred+iubred+nrwdrp;
    ntotal += nchange;
    if (ninfeas) break;
  }
  delete [] cliqueMin;
  delete [] cliqueMax;
  delete [] cliqueMin2;
  delete [] cliqueMax2;
  return (ninfeas);
}
// This just sets minima and maxima on rows
void 
CglProbing::tighten2(double *colLower, double * colUpper,
		     const int *column, const double *rowElements, 
		     const CoinBigIndex *rowStart, const int * rowLength,
		     double *rowLower, double *rowUpper, 
		     double * minR, double * maxR, int * markR,
		     int nRows,int nCols) const
{
  int i, j, k, kre;
  int krs;
  int iflagu, iflagl;
  double dmaxup, dmaxdown, value;

  for (i = 0; i < nRows; ++i) {
    if (rowLower[i]>-1.0e20||rowUpper[i]<1.0e20) {
      iflagu = 0;
      iflagl = 0;
      dmaxup = 0.0;
      dmaxdown = 0.0;
      krs = rowStart[i];
      kre = rowStart[i]+rowLength[i];
      
      /* ------------------------------------------------------------*/
      /* Compute L(i) and U(i) */
      /* ------------------------------------------------------------*/
      for (k = krs; k < kre; ++k) {
	value=rowElements[k];
	j = column[k];
	if (value > 0.0) {
	  if (colUpper[j] >= 1e15) {
	    dmaxup = 1e31;
	    ++iflagu;
	  } else {
	    dmaxup += colUpper[j] * value;
	  }
	  if (colLower[j] <= -1e15) {
	    dmaxdown = -1e31;
	    ++iflagl;
	  } else {
	    dmaxdown += colLower[j] * value;
	  }
	} else if (value<0.0) {
	  if (colUpper[j] >= 1e15) {
	    dmaxdown = -1e31;
	    ++iflagl;
	  } else {
	    dmaxdown += colUpper[j] * value;
	  }
	  if (colLower[j] <= -1e15) {
	    dmaxup = 1e31;
	    ++iflagu;
	  } else {
	    dmaxup += colLower[j] * value;
	  }
	}
      }
      if (iflagu)
	maxR[i]=1.0e60;
      else
	maxR[i]=dmaxup;
      if (iflagl) 
	minR[i]=-1.0e60;
      else
	minR[i]=dmaxdown;
      if (minR[i]<-1.0e10&&maxR[i]>1.0e10) {
	markR[i]=-2;
      } else {
	markR[i]=-1;
      }
    } else {
      minR[i]=-1.0e60;
      maxR[i]=1.0e60;
      markR[i]=-2;
    }
  }
}
#ifdef CGL_DEBUG
static int nPath=0;
#endif
//-------------------------------------------------------------------
// Generate disaggregation cuts
//------------------------------------------------------------------- 
void CglProbing::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			      const CglTreeInfo info) const
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n",nPath);
    nPath++;
    int nCols=si.getNumCols();
    int i;
    const double * solution = si.getColSolution();
    const double * lower = si.getColLower();
    const double * upper = si.getColUpper();
    const double * optimal = debugger->optimalSolution();
    const double * objective = si.getObjCoefficients();
    double objval1=0.0,objval2=0.0;
    for (i=0;i<nCols;i++) {
#if CGL_DEBUG>1
      printf("%d %g %g %g %g\n",i,lower[i],solution[i],upper[i],optimal[i]);
#endif
      objval1 += solution[i]*objective[i];
      objval2 += optimal[i]*objective[i];
      assert(optimal[i]>=lower[i]&&optimal[i]<=upper[i]);
    }
    printf("current obj %g, integer %g\n",objval1,objval2);
  }
#endif
  int saveRowCuts=rowCuts_;
  if (rowCuts_<0) {
    if (info.inTree)
      rowCuts_=4;
    else
      rowCuts_=-rowCuts_;
  }
  int nRows=si.getNumRows(); 
  double * rowLower = new double[nRows+1];
  double * rowUpper = new double[nRows+1];

  int nCols=si.getNumCols();
  // Set size if not set
  if (!rowCopy_) {
    numberRows_=nRows;
    numberColumns_=nCols;
  }
  double * colLower = new double[nCols];
  double * colUpper = new double[nCols];

  int ninfeas=gutsOfGenerateCuts(si,cs,rowLower,rowUpper,colLower,colUpper,info);
  if (ninfeas) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(DBL_MAX);
    rc.setUb(0.0);   
    cs.insert(rc);
#ifdef CGL_DEBUG
    const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
    if (debugger&&debugger->onOptimalPath(si))
      assert(!debugger->invalidCut(rc)); 
#endif
  }
  delete [] rowLower;
  delete [] rowUpper;
  delete [] colLower;
  delete [] colUpper;
  rowCuts_=saveRowCuts;
}
int CglProbing::generateCutsAndModify(const OsiSolverInterface & si, 
				      OsiCuts & cs,
				      CglTreeInfo * info2) 
{
  CglTreeInfo info = *info2;
  return generateCutsAndModify(si,cs,info);
}
int CglProbing::generateCutsAndModify(const OsiSolverInterface & si, 
				      OsiCuts & cs,
				      const CglTreeInfo info) 
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n",nPath);
    nPath++;
    int nCols=si.getNumCols();
    int i;
    const double * solution = si.getColSolution();
    const double * lower = si.getColLower();
    const double * upper = si.getColUpper();
    const double * optimal = debugger->optimalSolution();
    const double * objective = si.getObjCoefficients();
    double objval1=0.0,objval2=0.0;
    for (i=0;i<nCols;i++) {
#if CGL_DEBUG>1
      printf("%d %g %g %g %g\n",i,lower[i],solution[i],upper[i],optimal[i]);
#endif
      objval1 += solution[i]*objective[i];
      objval2 += optimal[i]*objective[i];
      assert(optimal[i]>=lower[i]-1.0e-5&&optimal[i]<=upper[i]+1.0e-5);
    }
    printf("current obj %g, integer %g\n",objval1,objval2);
  }
#endif
  int saveRowCuts=rowCuts_;
  if (rowCuts_<0) {
    if (info.inTree)
      rowCuts_=4;
    else
      rowCuts_=-rowCuts_;
  }
  int saveMode = mode_;
  bool rowCliques=false;
  if (!mode_) {
    if (info.pass!=4||info.inTree) {
      mode_=1;
    } else {
      saveMode=1; // make sure do just once
      rowCliques=true;
    }
  }
  int nRows=si.getNumRows(); 
  double * rowLower = new double[nRows+1];
  double * rowUpper = new double[nRows+1];

  int nCols=si.getNumCols(); 
  double * colLower = new double[nCols];
  double * colUpper = new double[nCols];

  int ninfeas=gutsOfGenerateCuts(si,cs,rowLower,rowUpper,colLower,colUpper,info);
  if (ninfeas) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(DBL_MAX);
    rc.setUb(0.0);   
    cs.insert(rc);
#ifdef CGL_DEBUG
    const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
    if (debugger&&debugger->onOptimalPath(si))
      assert(!debugger->invalidCut(rc)); 
#endif
  }
  rowCuts_=saveRowCuts;
  mode_=saveMode;
  // move bounds so can be used by user
  if (mode_==3) {
    delete [] rowLower_;
    delete [] rowUpper_;
    rowLower_ = rowLower;
    rowUpper_ = rowUpper;
  } else {
    delete [] rowLower;
    delete [] rowUpper;
  }
  delete [] colLower_;
  delete [] colUpper_;
  colLower_	= colLower;
  colUpper_	= colUpper;
  // Setup information 
  if (rowCliques&&numberRows_&&numberColumns_)
    setupRowCliqueInformation(si);
  return ninfeas;
}
bool analyze(const OsiSolverInterface * solverX, char * intVar,
             double * lower, double * upper)
{
  OsiSolverInterface * solver = solverX->clone();
  const double *objective = solver->getObjCoefficients() ;
  int numberColumns = solver->getNumCols() ;
  int numberRows = solver->getNumRows();
  double direction = solver->getObjSense();
  int iRow,iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix  matrixByCol(*solver->getMatrixByCol());
  const double * element = matrixByCol.getElements();
  const int * row = matrixByCol.getIndices();
  const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
  const int * columnLength = matrixByCol.getVectorLengths();

  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();

  char * ignore = new char [numberRows];
  int * which = new int[numberRows];
  double * changeRhs = new double[numberRows];
  memset(changeRhs,0,numberRows*sizeof(double));
  memset(ignore,0,numberRows);
  int numberChanged=0;
  bool finished=false;
  while (!finished) {
    int saveNumberChanged = numberChanged;
    for (iRow=0;iRow<numberRows;iRow++) {
      int numberContinuous=0;
      double value1=0.0,value2=0.0;
      bool allIntegerCoeff=true;
      double sumFixed=0.0;
      int jColumn1=-1,jColumn2=-1;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (upper[jColumn] > lower[jColumn]+1.0e-8) {
          if (!intVar[jColumn]) {
            if (numberContinuous==0) {
              jColumn1=jColumn;
              value1=value;
            } else {
              jColumn2=jColumn;
              value2=value;
            }
            numberContinuous++;
          } else {
            if (fabs(value-floor(value+0.5))>1.0e-12)
              allIntegerCoeff=false;
          }
        } else {
          sumFixed += lower[jColumn]*value;
        }
      }
      double low = rowLower[iRow];
      if (low>-1.0e20) {
        low -= sumFixed;
        if (fabs(low-floor(low+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      double up = rowUpper[iRow];
      if (up<1.0e20) {
        up -= sumFixed;
        if (fabs(up-floor(up+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      if (!allIntegerCoeff)
        continue; // can't do
      if (numberContinuous==1) {
        // see if really integer
        // This does not allow for complicated cases
        if (low==up) {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // integer
              numberChanged++;
              intVar[jColumn1]=77;
            }
          }
        } else {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // This constraint will not stop it being integer
              ignore[iRow]=1;
            }
          }
        }
      } else if (numberContinuous==2) {
        if (low==up) {
          /* need general theory - for now just look at 2 cases -
             1 - +- 1 one in column and just costs i.e. matching objective
             2 - +- 1 two in column but feeds into G/L row which will try and minimize
             (take out 2 for now - until fixed)
          */
          if (fabs(value1)==1.0&&value1*value2==-1.0&&!lower[jColumn1]
              &&!lower[jColumn2]&&columnLength[jColumn1]==1&&columnLength[jColumn2]==1) {
            int n=0;
            int i;
            double objChange=direction*(objective[jColumn1]+objective[jColumn2]);
            double bound = CoinMin(upper[jColumn1],upper[jColumn2]);
            bound = CoinMin(bound,1.0e20);
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                which[n++]=jRow;
                changeRhs[jRow]=value;
              }
            }
            for ( i=columnStart[jColumn2];i<columnStart[jColumn2]+columnLength[jColumn2];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                if (!changeRhs[jRow]) {
                  which[n++]=jRow;
                  changeRhs[jRow]=value;
                } else {
                  changeRhs[jRow]+=value;
                }
              }
            }
            if (objChange>=0.0) {
              // see if all rows OK
              bool good=true;
              for (i=0;i<n;i++) {
                int jRow = which[i];
                double value = changeRhs[jRow];
                if (value) {
                  value *= bound;
                  if (rowLength[jRow]==1) {
                    if (value>0.0) {
                      double rhs = rowLower[jRow];
                      if (rhs>0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    } else {
                      double rhs = rowUpper[jRow];
                      if (rhs<0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    }
                  } else if (rowLength[jRow]==2) {
                    if (value>0.0) {
                      if (rowLower[jRow]>-1.0e20)
                        good=false;
                    } else {
                      if (rowUpper[jRow]<1.0e20)
                        good=false;
                    }
                  } else {
                    good=false;
                  }
                }
              }
              if (good) {
                // both can be integer
                numberChanged++;
                intVar[jColumn1]=77;
                numberChanged++;
                intVar[jColumn2]=77;
              }
            }
            // clear
            for (i=0;i<n;i++) {
              changeRhs[which[i]]=0.0;
            }
          }
        }
      }
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8&&!intVar[iColumn]) {
        double value;
        value = upper[iColumn];
        if (value<1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        value = lower[iColumn];
        if (value>-1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        bool integer=true;
        for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          if (!ignore[iRow]) {
            integer=false;
            break;
          }
        }
        if (integer) {
          // integer
          numberChanged++;
          intVar[iColumn]=77;
        }
      }
    }
    finished = numberChanged==saveNumberChanged;
  }
  bool feasible=true;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (intVar[iColumn]==77) {
      if (upper[iColumn]>1.0e20) {
        upper[iColumn] = 1.0e20;
      } else {
        upper[iColumn] = floor(upper[iColumn]+1.0e-5);
      }
      if (lower[iColumn]<-1.0e20) {
        lower[iColumn] = -1.0e20;
      } else {
        lower[iColumn] = ceil(lower[iColumn]-1.0e-5);
        if (lower[iColumn]>upper[iColumn])
          feasible=false;
      }
      if (lower[iColumn]==0.0&&upper[iColumn]==1.0)
        intVar[iColumn]=1;
      else if (lower[iColumn]==upper[iColumn])
        intVar[iColumn]=0;
      else
        intVar[iColumn]=2;
    }
  }
  delete [] which;
  delete [] changeRhs;
  delete [] ignore;
  //if (numberChanged)
  //printf("%d variables could be made integer\n",numberChanged);
  delete solver;
  return feasible;
}
int CglProbing::gutsOfGenerateCuts(const OsiSolverInterface & si, 
                                   OsiCuts & cs ,
                                   double * rowLower, double * rowUpper,
                                   double * colLower, double * colUpper,
                                   const CglTreeInfo info) const
{
  // Get basic problem information
  int nRows;
  
  CoinPackedMatrix * rowCopy=NULL;

  // get branch and bound cutoff
  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff);
  if (!cutoff_available) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  cutoff *= si.getObjSense();
  if (fabs(cutoff)>1.0e30)
    assert (cutoff>1.0e30);
  int mode=mode_;
  
  int nCols=si.getNumCols(); 

  // get integer variables
  char * intVar = new char[nCols];
  int i;
  int numberIntegers=0;
  memcpy(colLower,si.getColLower(),nCols*sizeof(double));
  memcpy(colUpper,si.getColUpper(),nCols*sizeof(double));
  const double * colsol =si.getColSolution();
  // and put reasonable bounds on integer variables
  for (i=0;i<nCols;i++) {
    if (si.isInteger(i)) {
      intVar[i]=2;
      numberIntegers++;
      if (si.isBinary(i)) {
        intVar[i]=1;  
      } else {
	// make sure reasonable bounds
	if (colsol[i]<1.0e10&&colUpper[i]>1.0e12) 
	  colUpper[i] = 1.23456789e11;
	if (colsol[i]>-1.0e10&&colLower[i]<-1.0e12) 
	  colLower[i] = -1.23456789e11;
      }
    } else {
      intVar[i]=0;
    }
  }
  bool feasible=true;
  if (!info.inTree) {
    // make more integer
    feasible = analyze(&si,intVar,colLower,colUpper);
  }
  if (feasible&&PROBING_EXTRA_STUFF) {
    // tighten bounds on djs
    // should be in CbcCutGenerator and check if basic 
    const double * djs =si.getReducedCost();
    const double * colsol =si.getColSolution();
    double direction = si.getObjSense();
    double cutoff;
    bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff);
    if (!cutoff_available) { // cut off isn't set or isn't valid
      cutoff = si.getInfinity();
    }
    cutoff *= direction;
    if (fabs(cutoff)>1.0e30)
      assert (cutoff>1.0e30);
    double current = si.getObjValue();
    current *= direction;
    double gap=CoinMax(cutoff-current,1.0e-1);
    for (int i = 0; i < nCols; ++i) {
      double djValue = djs[i]*direction;
      if (colUpper[i]-colLower[i]>1.0e-8) {
        if (colsol[i]<colLower[i]+primalTolerance_) {
          if (djValue>gap) {
            if (si.isInteger(i)) {
              printf("why int %d not fixed at lb\n",i);
              colUpper[i]= colLower[i];
            } else {
              double newUpper = colLower[i] + gap/djValue;
              if (newUpper<colUpper[i]) {
                //printf("%d ub from %g to %g\n",i,colUpper[i],newUpper);
                colUpper[i]= CoinMax(newUpper,colLower[i]+1.0e-5);
              }
            }
          }
        } else if (colsol[i]>colUpper[i]-primalTolerance_) {
          if (-djValue>gap) {
            if (si.isInteger(i)) {
              printf("why int %d not fixed at ub\n",i);
              colLower[i]= colUpper[i];
            } else {
              double newLower = colUpper[i] + gap/djValue;
              if (newLower>colLower[i]) {
                //printf("%d lb from %g to %g\n",i,colLower[i],newLower);
                colLower[i]= CoinMin(newLower,colUpper[i]-1.0e-5);
              }
            }
          }
        }
      }
    }
  }
  int ninfeas=0;
  // Set up maxes
  int maxProbe = info.inTree ? maxProbe_ : maxProbeRoot_;
  int maxElements = info.inTree ? maxElements_ : maxElementsRoot_;

  // Get objective offset
  double offset;
  si.getDblParam(OsiObjOffset,offset);
  // see if using cached copy or not
  if (!rowCopy_) {
    // create from current
    nRows=si.getNumRows(); 
    
    // mode==0 is invalid if going from current matrix
    if (mode==0)
      mode=1;
    rowCopy = new CoinPackedMatrix(*si.getMatrixByRow());
    // add in objective if there is a cutoff
    if (cutoff<1.0e30&&usingObjective_) {
      int * columns = new int[nCols];
      double * elements = new double[nCols];
      int n=0;
      const double * objective = si.getObjCoefficients();
      bool maximize = (si.getObjSense()==-1);
      for (i=0;i<nCols;i++) {
	if (objective[i]) {
	  elements[n]= (maximize) ? -objective[i] : objective[i];
	  columns[n++]=i;
	}
      }
      rowCopy->appendRow(n,columns,elements);
      delete [] columns;
      delete [] elements;
      memcpy(rowLower,si.getRowLower(),nRows*sizeof(double));
      memcpy(rowUpper,si.getRowUpper(),nRows*sizeof(double));
      rowLower[nRows]=-DBL_MAX;
      rowUpper[nRows]=cutoff+offset;
      nRows++;
    } else {
      memcpy(rowLower,si.getRowLower(),nRows*sizeof(double));
      memcpy(rowUpper,si.getRowUpper(),nRows*sizeof(double));
    }
  } else {
    // use snapshot
    nRows=numberRows_;
    assert(nCols==numberColumns_);
    
    rowCopy = rowCopy_;
    assert (rowCopy_->getNumRows()==numberRows_);
    rowLower = new double[nRows];
    rowUpper = new double[nRows];
    memcpy(rowLower,rowLower_,nRows*sizeof(double));
    memcpy(rowUpper,rowUpper_,nRows*sizeof(double));
    if (usingObjective_) {
      rowLower[nRows-1]=-DBL_MAX;
      rowUpper[nRows-1]=cutoff+offset;
    }
  }
  int nRowsSafe=CoinMin(nRows,si.getNumRows());
#ifdef CGL_DEBUG
  const OsiRowCutDebugger * debugger = si.getRowCutDebugger();
  if (debugger&&!debugger->onOptimalPath(si))
    debugger = NULL;
#else
  const OsiRowCutDebugger * debugger = NULL;
#endif
   
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info.inTree ? nRowsSafe/3 : nRowsSafe;
  row_cut rowCut(nRowsFake);
  int * markR = new int [nRows];
  double * minR = new double [nRows];
  double * maxR = new double [nRows];
  if (mode) {
    ninfeas= tighten(colLower, colUpper, column, rowElements,
			 rowStart, rowLength, rowLower, rowUpper,
			 nRows, nCols, intVar, 2, primalTolerance_);
    if (!feasible)
      ninfeas=1;
    if (!ninfeas) {
      // create column cuts where integer bounds have changed
      {
	const double * lower = si.getColLower();
	const double * upper = si.getColUpper();
	const double * colsol = si.getColSolution();
	int numberChanged=0,ifCut=0;
	CoinPackedVector lbs;
	CoinPackedVector ubs;
	for (i = 0; i < nCols; ++i) {
	  if (intVar[i]) {
	    colUpper[i] = CoinMin(upper[i],floor(colUpper[i]+1.0e-4));
	    if (colUpper[i]<upper[i]-1.0e-8) {
	      if (colUpper[i]<colsol[i]-1.0e-8)
		ifCut=1;
	      ubs.insert(i,colUpper[i]);
	      numberChanged++;
	    }
	    colLower[i] = CoinMax(lower[i],ceil(colLower[i]-1.0e-4));
	    if (colLower[i]>lower[i]+1.0e-8) {
	      if (colLower[i]>colsol[i]+1.0e-8)
		ifCut=1;
	      lbs.insert(i,colLower[i]);
	      numberChanged++;
	    }
	  }
	}
	if (numberChanged) {
	  OsiColCut cc;
	  cc.setUbs(ubs);
	  cc.setLbs(lbs);
	  if (ifCut) {
	    cc.setEffectiveness(100.0);
	  } else {
	    cc.setEffectiveness(1.0e-5);
	  }
#ifdef CGL_DEBUG
	  checkBounds(debugger,cc);
#endif
	  cs.insert(cc);
	}
      }
      if (maxProbe>0) {
        int * look = new int[nCols];
        int nLook=0;
        // get min max etc for rows
        tighten2(colLower, colUpper, column, rowElements,
                 rowStart, rowLength, rowLower, rowUpper,
                 minR , maxR , markR, nRows, nCols);
        // decide what to look at
        if (mode==1) {
          const double * colsol = si.getColSolution();
          double_int_pair * array = new double_int_pair [nCols];
#	  ifdef ZEROFAULT
	  std::memset(array,0,sizeof(double_int_pair)*nCols) ;
#	  endif
          for (i=0;i<nCols;i++) {
            if (intVar[i]&&colUpper[i]-colLower[i]>1.0e-8) {
              double away = fabs(0.5-(colsol[i]-floor(colsol[i])));
              if (away<0.49999||!info.inTree) {
                array[nLook].infeasibility=away;
                array[nLook++].sequence=i;
              }
            }
          }
          std::sort(array,array+nLook,double_int_pair_compare());
          nLook=CoinMin(nLook,maxProbe);
          for (i=0;i<nLook;i++) {
            look[i]=array[i].sequence;
          }
          delete [] array;
        } else {
          for (i=0;i<nCols;i++) {
            if (intVar[i]&&colUpper[i]-colLower[i]>1.0e-8) {
              look[nLook++]=i;
            }
          }
        }
        // Only look at short rows
        for (i=0;i<nRows;i++) {
          if (rowLength[i]>maxElements)
            markR[i]=-2;
        }
        if (!numberCliques_) {
          ninfeas= probe(si, debugger, cs, colLower, colUpper, rowCopy,
                         rowLower, rowUpper,
                         intVar, minR, maxR, markR,
                         look,nLook,info);
        } else {
          ninfeas= probeCliques(si, debugger, cs, colLower, colUpper, rowCopy,
                                rowLower, rowUpper,
                                intVar, minR, maxR, markR,
                                look,nLook,info);
        }
        delete [] look;
      }
    }
  } else if (maxProbe>0) {
    // global cuts from previous calculations
    // could check more thoroughly that integers are correct
    assert(numberIntegers==numberIntegers_);
    // make up list of new variables to look at 
    int * look = new int[nCols];
    int nLook=0;
    const double * colsol = si.getColSolution();
    double_int_pair * array = new double_int_pair [nCols];
#   ifdef ZEROFAULT
    std::memset(array,0,sizeof(double_int_pair)*nCols) ;
#   endif
    for (i=0;i<number01Integers_;i++) {
      int j=cutVector_[i].sequence;
      if (!cutVector_[i].index&&colUpper[j]-colLower[j]>1.0e-8) {
	double away = fabs(0.5-(colsol[j]-floor(colsol[j])));
        array[nLook].infeasibility=away;
        array[nLook++].sequence=i;
      }
    }
    std::sort(array,array+nLook,double_int_pair_compare());
    nLook=CoinMin(nLook,maxProbe);
    for (i=0;i<nLook;i++) {
      look[i]=array[i].sequence;
    }
    delete [] array;
    // get min max etc for rows
    tighten2(colLower, colUpper, column, rowElements,
	     rowStart, rowLength, rowLower, rowUpper,
	     minR , maxR , markR, nRows, nCols);
    OsiCuts csNew;
    // don't do cuts at all if 0 (i.e. we are just checking bounds)
    if (rowCuts_) {
      // Only look at short rows
      for (i=0;i<nRows;i++) {
        if (rowLength[i]>maxElements)
          markR[i]=-2;
      }
      ninfeas= probeCliques(si, debugger, csNew, colLower, colUpper, rowCopy,
		     rowLower, rowUpper,
		     intVar, minR, maxR, markR,
		     look,nLook,info);
    }
    if (!ninfeas) {
      // go through row cuts
      int nCuts = csNew.sizeRowCuts();
      int iCut;
      // need space for backward lookup
      // just for ones being looked at
      int * backward = new int [2*nCols];
      int * onList = backward + nCols;
      for (i=0;i<nCols;i++) {
        backward[i]=-1;
        onList[i]=0;
      }
      for (i=0;i<number01Integers_;i++) {
	int j=cutVector_[i].sequence;
	backward[j]=i;
	onList[j]=1;
      }
      // first do counts
      // we know initialized to zero
      for (iCut=0;iCut<nCuts;iCut++) {
	OsiRowCut rcut;
	CoinPackedVector rpv;
	rcut = csNew.rowCut(iCut);
	rpv = rcut.row();
	assert(rpv.getNumElements()==2);
	const int * indices = rpv.getIndices();
	double* elements = rpv.getElements();
	double lb=rcut.lb();
	// find out which integer
        int which=0;
        i=backward[indices[0]];
        if (i<0||!onList[indices[0]]) {
          which=1;
          i=backward[indices[1]];
          // Just possible variable was general integer but now 0-1
          if (!onList[indices[which]])
            continue;
        }
        int other = indices[1-which];
	if (lb==-DBL_MAX) {
          if (!rcut.ub()) {
            // UB
            if (elements[which]<0.0) {
              //assert (elements[1-which]>0.0);
              // delta to 0 => x to 0.0
              cutVector_[i].length++;
            } else {
              if (elements[1-which]<0.0&&fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other];
                    cutVector_[i].length++;
                    cutVector_[j].length++;
                  } else {
                    abort();
                  }
                }
              }
            }
          } else {
            if (onList[other]) {
              if (elements[0]==1.0&&elements[1]==1.0&&rcut.ub()==1.0) {
                // can do something ?
                int j=backward[other];
                cutVector_[i].length++;
                cutVector_[j].length++;
              } else {
                abort();
              }
            }
          }
	} else {
          assert(rcut.ub()==DBL_MAX);
          if (!lb) {
            // LB
            if (elements[which]>0.0) {
              //assert (elements[1-which]<0.0);
              // delta to 0 => x to 0.0
              // flip so same as UB
              cutVector_[i].length++; 
            } else {
              if (elements[1-which]<0.0&&fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other];
                    cutVector_[i].length++;
                    cutVector_[j].length++;
                  } else {
                    abort();
                  }
                }
              }
            }
          }
	} 
      }
      // allocate space
      for (i=0;i<number01Integers_;i++) {
	int j=cutVector_[i].sequence;
	if (onList[j]&&!cutVector_[i].index) {
	  disaggregation thisOne=cutVector_[i];
	  cutVector_[i].index=new disaggregationAction [thisOne.length];
          cutVector_[i].length=0;
	}
      }
      // now put in
      for (iCut=0;iCut<nCuts;iCut++) {
	OsiRowCut rcut;
	CoinPackedVector rpv;
	int iput;
	rcut = csNew.rowCut(iCut);
	rpv = rcut.row();
	assert(rpv.getNumElements()==2);
	const int * indices = rpv.getIndices();
	double* elements = rpv.getElements();
	double lb=rcut.lb();
	// find out which integer
	// find out which integer
        int which=0;
        i=backward[indices[0]];
        if (i<0||!onList[indices[0]]) {
          which=1;
          i=backward[indices[1]];
          // Just possible variable was general integer but now 0-1
          if (!onList[indices[which]])
            continue;
        }
        int other = indices[1-which];
        int j = other ? backward[other] : -1;
	if (lb==-DBL_MAX) {
          if (!rcut.ub()) {
            // UB
            if (elements[which]<0.0) {
              iput=cutVector_[i].length;
              if (j>=0)
                cutVector_[i].index[iput].affected=j;
              else
                cutVector_[i].index[iput].affected=other;
              cutVector_[i].index[iput].whenAtUB=0;
              cutVector_[i].index[iput].affectedToUB=0;
              cutVector_[i].index[iput].zeroOne = onList[other] ? 1 : 0;
              cutVector_[i].length++;
            } else { 
              if (elements[1-which]<0.0&&fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                // delta to 1 => x to upper bound
                iput=cutVector_[i].length;
                if (j>=0)
                  cutVector_[i].index[iput].affected=j;
                else
                  cutVector_[i].index[iput].affected=other;
                cutVector_[i].index[iput].whenAtUB=1;
                cutVector_[i].index[iput].affectedToUB=1;
                cutVector_[i].index[iput].zeroOne = onList[other] ? 1 : 0;
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other];
                    assert (j>=0);
                    // flip so value0 1.0
                    if (value1==1.0) {
                      j=i;
                      i=backward[other];
                      value1=value0;
                      value0=1.0;
                    }
                    assert (value0==1.0);
                    assert (value1==-1.0);
                    iput=cutVector_[i].length;
                    cutVector_[i].index[iput].affected=j;
                    cutVector_[i].index[iput].whenAtUB=1;
                    cutVector_[i].index[iput].affectedToUB=1;
                    cutVector_[i].index[iput].zeroOne = 1;
                    cutVector_[i].length++;
                    iput=cutVector_[j].length;
                    cutVector_[j].index[iput].affected=i;
                    cutVector_[j].index[iput].whenAtUB=0;
                    cutVector_[j].index[iput].affectedToUB=0;
                    cutVector_[j].index[iput].zeroOne = 1;
                    cutVector_[j].length++;
                  }
                }
              }
            }
          } else {
            if (onList[other]) {
              if (elements[0]==1.0&&elements[1]==1.0&&rcut.ub()==1.0) {
                // can do something ?
                int j=backward[other];
                assert (j>=0);
                iput=cutVector_[i].length;
                cutVector_[i].index[iput].affected=j;
                cutVector_[i].index[iput].whenAtUB=1;
                cutVector_[i].index[iput].affectedToUB=0;
                cutVector_[i].index[iput].zeroOne = 1;
                cutVector_[i].length++;
                iput=cutVector_[j].length;
                cutVector_[j].index[iput].affected=i;
                cutVector_[j].index[iput].whenAtUB=1;
                cutVector_[j].index[iput].affectedToUB=0;
                cutVector_[j].index[iput].zeroOne = 1;
                cutVector_[j].length++;
              } else {
                abort();
              }
            }
          }
	} else {
          assert(rcut.ub()==DBL_MAX);
          if (!lb) {
            // LB
            if (elements[which]>0.0) {
              iput=cutVector_[i].length;
              if (j>=0)
                cutVector_[i].index[iput].affected=j;
              else
                cutVector_[i].index[iput].affected=other;
              cutVector_[i].index[iput].whenAtUB=0;
              cutVector_[i].index[iput].affectedToUB=0;
              cutVector_[i].index[iput].zeroOne = onList[other] ? 1 : 0;
              cutVector_[i].length++;
            } else { 
              if (elements[1-which]<0.0&&fabs(elements[which]/elements[1-which]-
                                              colUpper[other])<1.0e-5) {
                iput=cutVector_[i].length;
                if (j>=0)
                  cutVector_[i].index[iput].affected=j;
                else
                  cutVector_[i].index[iput].affected=other;
                cutVector_[i].index[iput].whenAtUB=1;
                cutVector_[i].index[iput].affectedToUB=1;
                cutVector_[i].index[iput].zeroOne = onList[other] ? 1 : 0;
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0*value1==-1.0) {
                    // can do something ?
                    int j=backward[other];
                    assert (j>=0);
                    // flip so value0 -1.0
                    if (value1==-1.0) {
                      j=i;
                      i=backward[other];
                      value1=value0;
                      value0=-1.0;
                    }
                    assert (value0==-1.0);
                    assert (value1==1.0);
                    iput=cutVector_[i].length;
                    cutVector_[i].index[iput].affected=j;
                    cutVector_[i].index[iput].whenAtUB=1;
                    cutVector_[i].index[iput].affectedToUB=1;
                    cutVector_[i].index[iput].zeroOne = 1;
                    cutVector_[i].length++;
                    iput=cutVector_[j].length;
                    cutVector_[j].index[iput].affected=i;
                    cutVector_[j].index[iput].whenAtUB=0;
                    cutVector_[j].index[iput].affectedToUB=0;
                    cutVector_[j].index[iput].zeroOne = 1;
                    cutVector_[j].length++;
                  }
                }
              }
            }
          }
	} 
      }
      delete [] backward;
      // Now sort and get rid of duplicates
      // could also see if any are cliques
      int longest=0;
      for (i=0;i<number01Integers_;i++) 
        longest = CoinMax(longest, cutVector_[i].length);
      unsigned int * sortit = new unsigned int[longest];
      for (i=0;i<number01Integers_;i++) {
        disaggregation & thisOne=cutVector_[i];
        int k;
        int number = thisOne.length;
        for (k=0;k<number;k++) {
          unsigned int affected = thisOne.index[k].affected;
          unsigned int zeroOne = thisOne.index[k].zeroOne;
          unsigned int whenAtUB = thisOne.index[k].whenAtUB;
          unsigned int affectedToUB = thisOne.index[k].affectedToUB;
          sortit[k]=(affected<<3)|(zeroOne<<2)|(whenAtUB<<1)|affectedToUB;
        }
        std::sort(sortit,sortit+number);
        unsigned int affectedLast = 0xffffffff;
        unsigned int zeroOneLast = 0;
        unsigned int whenAtUBLast = 0;
        unsigned int affectedToUBLast = 0; 
        int put=0;
        for (k=0;k<number;k++) {
          unsigned int affected = sortit[k]>>3;
          unsigned int zeroOne = (sortit[k]&4)>>2;
          unsigned int whenAtUB = (sortit[k]&2)>>1;
          unsigned int affectedToUB = sortit[k]&1;
          disaggregationAction action;
          action.affected=affected;
          action.zeroOne=zeroOne;
          action.whenAtUB=whenAtUB;
          action.affectedToUB=affectedToUB;
          if (affected!=affectedLast||zeroOne!=zeroOneLast) {
            // new variable
            thisOne.index[put++]=action;
          } else if (whenAtUB!=whenAtUBLast||affectedToUB!=affectedToUBLast) {
            // new action - what can we discover
            thisOne.index[put++]=action;
            int j=cutVector_[i].sequence;
            int k=affected;
            if (zeroOne) {
              k=cutVector_[k].sequence;
              if (logLevel_>1)
                printf("For %d %d 0-1 pair",j,k) ;
            } else {
              if (logLevel_>1)
                printf("For %d %d pair",j,k) ;
            }
            if (logLevel_>1)
              printf(" old whenAtUB, affectedToUB %d %d, new whenAtUB, affectedToUB %d %d\n",
                     whenAtUBLast, affectedToUBLast,whenAtUB, affectedToUB);
          }
          affectedLast=affected;
          zeroOneLast=zeroOne;
          whenAtUBLast=whenAtUB;
          affectedToUBLast=affectedToUB;
        }
        if (put<number) {
          //printf("%d reduced from %d to %d\n",i,number,put);
          thisOne.length=put;
        }
      }
      // And look at all where two 0-1 variables involved
      for (i=0;i<number01Integers_;i++) {
        disaggregation & thisOne=cutVector_[i];
        int k;
        int number = thisOne.length;
        for (k=0;k<number;k++) {
          unsigned int affected = thisOne.index[k].affected;
          unsigned int zeroOne = thisOne.index[k].zeroOne;
          if (zeroOne&&(int)affected>i) {
            unsigned int whenAtUB = thisOne.index[k].whenAtUB;
            unsigned int affectedToUB = thisOne.index[k].affectedToUB;
            disaggregation otherOne=cutVector_[affected];
            int numberOther = otherOne.length;
            // Could do binary search if a lot
            int lastAction=-1;
            for (int j=0;j<numberOther;j++) {
              if ((int) otherOne.index[j].affected==i) {
                unsigned int whenAtUBOther = otherOne.index[j].whenAtUB;
                unsigned int affectedToUBOther = otherOne.index[j].affectedToUB;
                /* action -
                   0 -> x + y <=1 (1,1 impossible)
                   1 -> x - y <=0 (1,0 impossible)
                   2 -> -x + y <=0 (0,1 impossible)
                   3 -> -x -y <= -1 (0,0 impossible)
                  10 -> x == y
                  11 -> x + y == 1
                  20 -> x == 0
                  21 -> x == 1
                  22 -> y == 0
                  23 -> y == 1
                */
                int action=-1;
                if (whenAtUB) {
                  if (affectedToUB) {
                    // x -> 1 => y -> 1
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=10; // x,y must be same
                      } else {
                        // y -> 1 => x -> 0 
                        action=20; // If x is 1 then contradiction
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=23; // if y is 0 then contradiction
                      } else {
                        // y -> 0 => x -> 0
                        action=1; // x,y 1,0 impossible 
                      }
                    }
                  } else {
                    // x -> 1 => y -> 0
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=22; // If y is 1 then contradiction
                      } else {
                        // y -> 1 => x -> 0 
                        action=0; 
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=11; // x,y with same values impossible
                      } else {
                        // y -> 0 => x -> 0
                        action=20; // If x is 1 then contradiction
                      }
                    }
                  }
                } else {
                  if (affectedToUB) {
                    // x -> 0 => y -> 1
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=21; // If x is 0 then contradiction
                      } else {
                        // y -> 1 => x -> 0 
                        action=11; // x,y must be different
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=3; // one of x,y must be 1
                      } else {
                        // y -> 0 => x -> 0
                        action=23; // if y is 0 then contradiction
                      }
                    }
                  } else {
                    // x -> 0 => y -> 0
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action=2; // x,y 0,1 impossible
                      } else {
                        // y -> 1 => x -> 0 
                        action=22; // If y is 1 then contradiction
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1 
                        action=21; // if x is 0 then contradiction
                      } else {
                        // y -> 0 => x -> 0
                        action=10; // x,y must be same
                      }
                    }
                  }
                }
                assert (action>=0);
                if (action<4) {
                  // clique - see if there
                  if (oneFixStart_) {
                    switch (action) {
                    case 0:
                      break;
                    case 1:
                      break;
                    case 2:
                      break;
                    case 3:
                      break;
                    }
                    // If not can we add or strengthen
                  }
                  // check last action
                  if (lastAction>=0) {
                    if (logLevel_>1)
                      printf("XX lastAction %d, this %d\n",lastAction,action);
                  }
                } else if (action<12) {
                  if (logLevel_>1)
                    printf("XX Could eliminate one of %d %d 0-1 variables %c\n",i,affected,
                           (lastAction>=0) ? '*' : ' ');
                  if (info.strengthenRow) {
                    OsiRowCut rc;
                    int index[2];
                    double element[2];
                    index[0]=cutVector_[i].sequence;
                    element[0]=1.0;
                    index[1]=cutVector_[affected].sequence;
                    if (action==10) {
                      // 10 -> x == y
                      rc.setLb(0.0);
                      rc.setUb(0.0);   
                      element[1]= -1.0;
                    } else {
                      // 11 -> x + y == 1
                      rc.setLb(1.0);
                      rc.setUb(1.0);   
                      element[1]= 1.0;
                    }
                    rc.setRow(2,index,element);
                    cs.insert(rc);
                  }
                } else {
                  if (action<22)
                    if (logLevel_>1)
                      printf("XX Could fix a 0-1 variable %d\n",i);
                  else
                    if (logLevel_>1)
                      printf("XX Could fix a 0-1 variable %d\n",affected);
                }
                //printf("%d when %d forces %d to %d , %d when %d forces %d to %d\n",
                //     i,whenAtUB,affected,affectedToUB, 
                //     affected, whenAtUBOther,i, affectedToUBOther);
              }
            }
          }
        }
      }
      delete [] sortit;
    }
    delete [] look;
    if (cutVector_) {
      // now see if any disaggregation cuts are violated
      for (i=0;i<number01Integers_;i++) {
	int j=cutVector_[i].sequence;
	double solInt=colsol[j];
	double  upper, solValue;
	int icol;
	int index[2];
	double element[2];
	if (colUpper[j]-colLower[j]>1.0e-8) {
	  double away = fabs(0.5-(solInt-floor(solInt)));
	  if (away<0.4999999) {
	    disaggregation thisOne=cutVector_[i];
	    int k;
	    OsiRowCut rc;
	    for (k=0;k<thisOne.length;k++) {
	      icol = thisOne.index[k].affected;
              if (thisOne.index[k].zeroOne)
                icol = cutVector_[icol].sequence;
	      solValue=colsol[icol];
	      upper=colUpper_[icol];
              double infeasibility=0.0;
              if (!thisOne.index[k].whenAtUB) {
                if (!thisOne.index[k].affectedToUB) {
                  // delta -> 0 => x to lb (at present just 0)
                  infeasibility = solValue - upper * solInt;
                  if (infeasibility > 1.0e-3) {
                    rc.setLb(-DBL_MAX);
                    rc.setUb(0.0);
                    index[0]=icol;
                    element[0]=1.0;
                    index[1]=j;
                    element[1]= -upper;
                  } else {
                    infeasibility=0.0;
                  }
                } else {
                  // delta -> 0 => x to ub
                  abort();
                }
              } else {
                if (thisOne.index[k].affectedToUB) {
                  // delta -> 1 => x to ub (?)
                  icol = thisOne.index[k].affected;
                  if (thisOne.index[k].zeroOne)
                    icol = cutVector_[icol].sequence;
                  solValue=colsol[icol];
                  upper=colUpper_[icol];
                  if (!colLower[icol]) {
                    infeasibility = upper * solInt - solValue;
                    if (infeasibility > 1.0e-3) {
                      rc.setLb(-DBL_MAX);
                      rc.setUb(0.0);
                      index[0]=icol;
                      element[0]=-1.0;
                      index[1]=j;
                      element[1]= upper;
                    } else {
                      infeasibility=0.0;
                    }
                  } else {
                    assert (upper==colLower[icol]);
                    infeasibility=0.0;
                  }
                } else {
                  // delta + delta2 <= 1
                  assert (thisOne.index[k].zeroOne);
                  // delta -> 1 => delta2 -> 0
                  icol = thisOne.index[k].affected;
                  icol = cutVector_[icol].sequence;
                  // only do if icol > j
                  if (icol >j && colUpper[icol] ) {
                    solValue=colsol[icol];
                    if (!colLower[icol]) {
                      infeasibility = solInt + solValue - 1.0;
                      if (infeasibility > 1.0e-3) {
                        rc.setLb(-DBL_MAX);
                        rc.setUb(1.0);
                        index[0]=icol;
                        element[0]=1.0;
                        index[1]=j;
                        element[1]= 1.0;
                      } else {
                        infeasibility=0.0;
                      }
                    } else {
                      assert (upper==colLower[icol]);
                      infeasibility=0.0;
                    }
                  }
                }
              }
              if (infeasibility) {
                rc.setEffectiveness(infeasibility);
                rc.setRow(2,index,element);
                if (logLevel_>1)
                  printf("%g <= %g * x%d + %g * x%d <= %g\n",
                         rc.lb(),element[0],index[0],element[1],index[1],rc.ub());
#ifdef CGL_DEBUG
                if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                rowCut.addCutIfNotDuplicate(rc);
              }
	    }
	  }
	}
      }
    }
  }
  delete [] markR;
  delete [] minR;
  delete [] maxR;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info.strengthenRow);
  }
  // delete stuff
  if (!rowCopy_) {
    delete rowCopy;
  } else {
    delete [] rowLower;
    delete [] rowUpper;
  }
  delete [] intVar;
  // and put back unreasonable bounds on integer variables
  const double * trueLower = si.getColLower();
  const double * trueUpper = si.getColUpper();
  for (i=0;i<nCols;i++) {
    if (si.isInteger(i)&&!si.isBinary(i)) {
      if (colUpper[i] == 1.23456789e11) 
	colUpper[i] = trueUpper[i];
      if (colLower[i] == -1.23456789e11) 
	colLower[i] = trueLower[i];
    }
  }
  return ninfeas;
}
// Does probing and adding cuts
int CglProbing::probe( const OsiSolverInterface & si, 
		       const OsiRowCutDebugger * debugger, 
		       OsiCuts & cs, 
		       double * colLower, double * colUpper, 
		       CoinPackedMatrix *rowCopy,
		       double * rowLower, double * rowUpper,
		       char * intVar, double * minR, double * maxR, 
		       int * markR, 
		       int * look, int nLook,
                       const CglTreeInfo info) const
{
  int nRows=rowCopy->getNumRows();
  int nRowsSafe=CoinMin(nRows,si.getNumRows());
  int nCols=rowCopy->getNumCols();
  double * colsol = new double[nCols];
  double * djs = new double[nCols];
  const double * objective = si.getObjCoefficients();
  const double * currentColLower = si.getColLower();
  const double * currentColUpper = si.getColUpper();
  double * tempL = new double [nCols];
  double * tempU = new double [nCols];
  int * markC = new int [nCols];
  int * stackC = new int [2*nCols];
  int * stackR = new int [nRows];
  double * saveL = new double [2*nCols];
  double * saveU = new double [2*nCols];
  double * saveMin = new double [nRows];
  double * saveMax = new double [nRows];
  double * element = new double[nCols];
  int * index = new int[nCols];
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info.inTree ? nRowsSafe/3 : nRowsSafe;
  row_cut rowCut(nRowsFake);
  CoinPackedMatrix * columnCopy=NULL;
  // Set up maxes
  int maxStack = info.inTree ? maxStack_ : maxStackRoot_;
  int maxPass = info.inTree ? maxPass_ : maxPassRoot_;
  if (!rowCopy_) {
    columnCopy = new CoinPackedMatrix(*rowCopy);
    columnCopy->reverseOrdering();
  } else {
    columnCopy = columnCopy_;
  }
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements();
  double movement;
  int i, j, k,kk,jj;
  int jcol,kcol,irow,krow;
  bool anyColumnCuts=false;
  double dbound, value, value2;
  int ninfeas=0;
  int rowCuts;
  double disaggEffectiveness;
  /* clean up djs and solution */
  memcpy(djs,si.getReducedCost(),nCols*sizeof(double));
  memcpy(colsol, si.getColSolution(),nCols*sizeof(double));
  disaggEffectiveness=1.0e-3;
  rowCuts=rowCuts_;
  
  double direction = si.getObjSense();
  for (i = 0; i < nCols; ++i) {
    double djValue = djs[i]*direction;
    if (colUpper[i]-colLower[i]>1.0e-8) {
      if (colsol[i]<colLower[i]+primalTolerance_) {
        colsol[i]=colLower[i];
        djs[i] = max(0.0,djValue);
      } else if (colsol[i]>colUpper[i]-primalTolerance_) {
        colsol[i]=colUpper[i];
        djs[i] = min(0.0,djValue);
      } else {
        djs[i]=0.0;
      }
    }
  }

  int ipass=0,nfixed=-1;

  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff);
  if (!cutoff_available) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  cutoff *= direction;
  if (fabs(cutoff)>1.0e30)
    assert (cutoff>1.0e30);
  double current = si.getObjValue();
  current *= direction;
  /* for both way coding */
  int nstackC0=-1;
  int * stackC0 = new int[maxStack];
  double * lo0 = new double[maxStack];
  double * up0 = new double[maxStack];
  int nstackR,nstackC;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3;
    } else {
      markC[i]=0;
    }
  }
  double tolerance = 1.0e1*primalTolerance_;
  // If we are going to replace coefficient then we don't need to be effective
  double needEffectiveness = info.strengthenRow ? -1.0e10 : 1.0e-3;
  if (PROBING_EXTRA_STUFF) {
    int nCut=0;
    for (int iRow=0;iRow<nRows;iRow++) {
      int numberInt=0;
      int whichInt=-1;
      int numberNeg=0;
      double sumFixed=0.0;
      double intValue=0.0;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int jColumn = column[j];
        double value = rowElements[j];
        if (colUpper[jColumn] > colLower[jColumn]+1.0e-8) {
          if (intVar[jColumn]) {
            numberInt++;
            whichInt=jColumn;
            intValue=value;
          } else if (value<0) {
            numberNeg++;
          }
        } else {
          sumFixed += colLower[jColumn]*value;
        }
      }
      if (numberInt==1&&numberNeg==0&&intValue<0.0&&!rowUpper[iRow]&&rowLower[iRow]<-1.0e30&&!sumFixed) {
        double intSol = colsol[whichInt];
        for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
          int jColumn = column[j];
          //double value = rowElements[j];
          if (colUpper[jColumn] > colLower[jColumn]+1.0e-8) {
            if (!intVar[jColumn]) {
              if (colLower[jColumn]||colUpper[jColumn]>1.0)
                continue;;
              double upper = colUpper[jColumn];
              if (colsol[jColumn]>intSol*upper+1.0e-4) {
                nCut++;
                OsiRowCut rc;
                rc.setLb(-COIN_DBL_MAX);
                rc.setUb(0.0);
                rc.setEffectiveness(1.0e-5);
                int index[2];
                double element[2];
                index[0]=jColumn;
                index[1]=whichInt;
                element[0]=1.0;
                element[1]=-upper;
                rc.setRow(2,index,element);
                cs.insert(rc);
              }
            }
          }
        }
      }
    }
    if (nCut)
      printf("%d possible cuts\n",nCut);
  }
  while (ipass<maxPass&&nfixed) {
    int iLook;
    ipass++;
    nfixed=0;
    for (iLook=0;iLook<nLook;iLook++) {
      double solval;
      double down;
      double up;
      if (rowCut.outOfSpace())
        break;
      int awayFromBound=1;
      j=look[iLook];
      solval=colsol[j];
      down = floor(solval+tolerance);
      up = ceil(solval-tolerance);
      if(colUpper[j]-colLower[j]<1.0e-8) markC[j]=3;
      if (markC[j]||!intVar[j]) continue;
      double saveSolval = solval;
      if (solval>=colUpper[j]-tolerance||solval<=colLower[j]+tolerance||up==down) {
	awayFromBound=0;
	if (solval<=colLower[j]+2.0*tolerance) {
	  solval = colLower[j]+1.0e-1;
	  down=colLower[j];
	  up=down+1.0;
	} else if (solval>=colUpper[j]-2.0*tolerance) {
	  solval = colUpper[j]-1.0e-1;
	  up=colUpper[j];
	  down=up-1;
	} else {
          // odd
          up=down+1.0;
          solval = down+1.0e-1;
        }
      }
      assert (up<=colUpper[j]);
      assert (down>=colLower[j]);
      assert (up>down);
      int istackC,iway, istackR;
      int way[]={1,2,1};
      int feas[]={1,2,4};
      int feasible=0;
      int notFeasible;
      for (iway=0;iway<3;iway ++) {
        int fixThis=0;
        double objVal=current;
        int goingToTrueBound=0;
        stackC[0]=j;
        markC[j]=way[iway];
        double solMovement;
        if (way[iway]==1) {
          movement=down-colUpper[j];
          solMovement = down-colsol[j];
          assert(movement<-0.99999);
          if (fabs(down-colLower[j])<1.0e-7) {
            goingToTrueBound=2;
            down=colLower[j];
          }
        } else {
          movement=up-colLower[j];
          solMovement = up-colsol[j];
          assert(movement>0.99999);
          if (fabs(up-colUpper[j])<1.0e-7) {
            goingToTrueBound=2;
            up=colUpper[j];
          }
        }
        if (goingToTrueBound&&(colUpper[j]-colLower[j]>1.5||colLower[j]))
          goingToTrueBound=1;
        // switch off disaggregation if not wanted
        if ((rowCuts&1)==0)
          goingToTrueBound=0;
#ifdef PRINT_DEBUG
        if (fabs(movement)>1.01) {
          printf("big %d %g %g %g\n",j,colLower[j],solval,colUpper[j]);
        }
#endif
        if (solMovement*djs[j]>0.0)
          objVal += solMovement*djs[j];
        nstackC=1;
        nstackR=0;
        saveL[0]=colLower[j];
        saveU[0]=colUpper[j];
        assert (saveU[0]>saveL[0]);
        notFeasible=0;
        if (movement<0.0) {
          colUpper[j] += movement;
          colUpper[j] = floor(colUpper[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d down to 0\n",j);
#endif
        } else {
          colLower[j] += movement;
          colLower[j] = floor(colLower[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d up to 1\n",j);
#endif
        }
        if (fabs(colUpper[j]-colLower[j])<1.0e-6)
          markC[j]=3; // say fixed
        istackC=0;
        /* update immediately */
        for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
          irow = row[k];
          value = columnElements[k];
          if (markR[irow]==-1) {
            stackR[nstackR]=irow;
            markR[irow]=nstackR;
            saveMin[nstackR]=minR[irow];
            saveMax[nstackR]=maxR[irow];
            nstackR++;
          } else if (markR[irow]==-2) {
            continue;
          }
          /* could check immediately if violation */
          if (movement>0.0) {
            /* up */
            if (value>0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] += value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] += value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          } else {
            /* down */
            if (value<0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] -= value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] -= value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          }
        }
        while (istackC<nstackC&&nstackC<maxStack) {
          int jway;
          jcol=stackC[istackC];
          jway=markC[jcol];
          // If not first and fixed then skip
          if (jway==3&&istackC) {
            //istackC++;
            //continue;
            //printf("fixed %d on stack\n",jcol);
          }
          for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
            // break if found not feasible
            if (notFeasible)
              break;
            irow = row[k];
            /*value = columnElements[k];*/
            if (markR[irow]!=-2) {
              /* see if anything forced */
              for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
                double moveUp=0.0;
                double moveDown=0.0;
                double newUpper=-1.0,newLower=1.0;
                kcol=column[kk];
                bool onList = (markC[kcol]!=0);
                if (markC[kcol]!=3) {
                  value2=rowElements[kk];
                  int markIt=markC[kcol];
                  if (value2 < 0.0) {
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colUpper[kcol]+
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            //markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            //markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  } else {
                    /* positive element */
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colUpper[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            //markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            //markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  }
                  if (nstackC<2*maxStack) 
                    markC[kcol] = markIt;
                  if (moveUp&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      assert (nstackC<nCols);
                      nstackC++;
                      onList=true;
                    }
                    if (newLower>colsol[kcol]) {
                      if (djs[kcol]<0.0) {
                        /* should be infeasible */
                        assert (newLower>colUpper[kcol]+primalTolerance_);
                      } else {
                        objVal += moveUp*djs[kcol];
                      }
                    }
                    if (intVar[kcol])
                      newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4));
                    colLower[kcol]=newLower;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
                      markC[kcol]=3; // say fixed
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* up */
                      if (value>0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveUp;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveUp;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (moveDown&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      assert (nstackC<nCols);
                      nstackC++;
                      onList=true;
                    }
                    if (newUpper<colsol[kcol]) {
                      if (djs[kcol]>0.0) {
                        /* should be infeasible */
                        assert (colLower[kcol]>newUpper+primalTolerance_);
                      } else {
                        objVal += moveDown*djs[kcol];
                      }
                    }
                    if (intVar[kcol])
                      newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4));
                    colUpper[kcol]=newUpper;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
                      markC[kcol]=3; // say fixed
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* down */
                      if (value<0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveDown;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveDown;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
                    notFeasible=1;;
                    k=columnStart[jcol]+columnLength[jcol];
                    istackC=nstackC+1;
#ifdef PRINT_DEBUG
                    printf("** not feasible this way\n");
#endif
                    break;
                  }
                }
              }
            }
          }
          istackC++;
        }
        if (!notFeasible) {
          if (objVal<=cutoff) {
            feasible |= feas[iway];
          } else {
#ifdef PRINT_DEBUG
            printf("not feasible on dj\n");
#endif
            notFeasible=1;
            if (iway==1&&feasible==0) {
              /* not feasible at all */
              ninfeas=1;
              j=nCols-1;
              break;
            }
          }
        } else if (iway==1&&feasible==0) {
          /* not feasible at all */
          ninfeas=1;
          j=nCols-1;
          iLook=nLook;
          ipass=maxPass;
          break;
        }
        if (notFeasible)
          goingToTrueBound=0;
        if (iway==2||(iway==1&&feasible==2)) {
          /* keep */
          iway=3;
          nfixed++;
          OsiColCut cc;
          int nTot=0,nFix=0,nInt=0;
          bool ifCut=false;
          for (istackC=0;istackC<nstackC;istackC++) {
            int icol=stackC[istackC];
            if (intVar[icol]) {
              if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
                element[nFix]=colUpper[icol];
                index[nFix++]=icol;
                nInt++;
                if (colsol[icol]>colUpper[icol]+primalTolerance_) {
                  ifCut=true;
                  anyColumnCuts=true;
                }
              }
            }
          }
          if (nFix) {
            nTot=nFix;
            cc.setUbs(nFix,index,element);
            nFix=0;
          }
          for (istackC=0;istackC<nstackC;istackC++) {
            int icol=stackC[istackC];
            if (intVar[icol]) {
              if (colLower[icol]>currentColLower[icol]+1.0e-4) {
                element[nFix]=colLower[icol];
                index[nFix++]=icol;
                nInt++;
                if (colsol[icol]<colLower[icol]-primalTolerance_) {
                  ifCut=true;
                  anyColumnCuts=true;
                }
              }
            }
          }
          if (nFix) {
            nTot+=nFix;
            cc.setLbs(nFix,index,element);
          }
          // could tighten continuous as well
          if (nInt) {
            if (ifCut) {
              cc.setEffectiveness(100.0);
            } else {
              cc.setEffectiveness(1.0e-5);
            }
#ifdef CGL_DEBUG
            checkBounds(debugger,cc);
#endif
            cs.insert(cc);
          }
          for (istackC=0;istackC<nstackC;istackC++) {
            int icol=stackC[istackC];
            if (colUpper[icol]-colLower[icol]>primalTolerance_) {
              markC[icol]=0;
            } else {
              markC[icol]=3;
            }
          }
          for (istackR=0;istackR<nstackR;istackR++) {
            int irow=stackR[istackR];
            markR[irow]=-1;
          }
        } else {
          /* is it worth seeing if can increase coefficients
             or maybe better see if it is a cut */
          if (iway==0) {
            nstackC0=CoinMin(nstackC,maxStack);
            double solMove = saveSolval-down;
            double boundChange;
            if (notFeasible) {
              nstackC0=0;
            } else {
              for (istackC=0;istackC<nstackC0;istackC++) {
                int icol=stackC[istackC];
                stackC0[istackC]=icol;
                lo0[istackC]=colLower[icol];
                up0[istackC]=colUpper[icol];
              }
            }
            /* restore all */
            assert (iway==0);
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              if(goingToTrueBound==2&&istackC) {
                // upper disaggregation cut would be
                // xval < upper + (old_upper-upper) (jval-down)
                boundChange = oldU-colUpper[icol];
                if (boundChange>0.0&&oldU<1.0e10&&
                    (colsol[icol]>colUpper[icol]
                     + boundChange*solMove+primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  rc.setLb(-DBL_MAX);
                  rc.setUb(colUpper[icol]-down*boundChange);
                  index[0]=icol;
                  element[0]=1.0;
                  index[1]=j;
                  element[1]= - boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol]-colUpper[icol])/
                    boundChange;
                  assert(newSol>solMove);
                  rc.setEffectiveness(newSol-solMove);
                  if (rc.effectiveness()>disaggEffectiveness) {
                    rc.setRow(2,index,element);
#ifdef CGL_DEBUG
                    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                    rowCut.addCutIfNotDuplicate(rc);
                  }
                }
                // lower disaggregation cut would be
                // xval > lower + (old_lower-lower) (jval-down)
                boundChange = oldL-colLower[icol];
                if (boundChange<0.0&&oldL>-1.0e10&&
                    (colsol[icol]<colLower[icol]
                     + boundChange*solMove-primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  rc.setLb(colLower[icol]-down*boundChange);
                  rc.setUb(DBL_MAX);
                  index[0]=icol;
                  element[0]=1.0;
                  index[1]=j;
                  element[1]=- boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol]-colLower[icol])/
                    boundChange;
                  assert(newSol>solMove);
                  rc.setEffectiveness(newSol-solMove);
                  if (rc.effectiveness()>disaggEffectiveness) {
                    rc.setRow(2,index,element);
#ifdef CGL_DEBUG
                    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                    rowCut.addCutIfNotDuplicate(rc);
#if 0
                    printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,oldL,oldU,colLower[icol],colsol[icol], j, colsol[j]);
                    printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
                           icol,boundChange,j,rc.ub());
#endif
                  }
                }
              }
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              markC[icol]=0;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0&&goingToTrueBound) {
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  // also see if singletons can go to good objective
                  bool moveSingletons=true;
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    int iColumn = column[kk];
                    double value = rowElements[kk];
                    sum += value*colsol[iColumn];
                    if (moveSingletons&&j!=iColumn) {
                      if (colUpper[iColumn]>colLower[iColumn]) {
                        if (columnLength[iColumn]!=1) {
                          moveSingletons=false;
                        }
                      }
                    }
                  }
                  if (moveSingletons) {
                    // can fix any with good costs
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      int iColumn = column[kk];
                      if (j!=iColumn) {
                        if (colUpper[iColumn]>colLower[iColumn]) {
                          if (columnLength[iColumn]==1) {
                            double value = rowElements[kk];
                            if (direction*objective[iColumn]*value<0.0&&!markC[iColumn]) {
                              // Fix
                              if (nstackC0+1<maxStack) {
                                stackC0[nstackC0]=iColumn;
                                if (value>0.0) {
                                  lo0[nstackC0]=colUpper[iColumn];
                                  up0[nstackC0]=colUpper[iColumn];
                                } else {
                                  lo0[nstackC0]=colLower[iColumn];
                                  up0[nstackC0]=colLower[iColumn];
                                }
                                nstackC0++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_||(info.strengthenRow&&rowLower[irow]<-1.0e20)) {
                    // can be a cut
                    // subtract gap from upper and integer coefficient
                    // saveU and saveL spare
                    int * index = (int *)saveL;
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]-gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=-gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(-DBL_MAX);
                    rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0));
                    // effectiveness is how far j moves
                    rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      // If strengthenRow point to row
                      //if(info.strengthenRow)
                      //printf("a point to row %d\n",irow);
                      rowCut.addCutIfNotDuplicate(rc,rowLower[irow]<-1.0e20 ? irow :-1);
                    }
                  }
                }
                gap = minR[irow]-rowLower[irow];
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  sum =0.0;
                  // also see if singletons can go to good objective
                  bool moveSingletons=true;
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    int iColumn = column[kk];
                    double value = rowElements[kk];
                    sum += value*colsol[iColumn];
                    if (moveSingletons&&j!=iColumn) {
                      if (colUpper[iColumn]>colLower[iColumn]) {
                        if (columnLength[iColumn]!=1) {
                          moveSingletons=false;
                        }
                      }
                    }
                  }
                  if (moveSingletons) {
                    // can fix any with good costs
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      int iColumn = column[kk];
                      if (j!=iColumn) {
                        if (colUpper[iColumn]>colLower[iColumn]) {
                          if (columnLength[iColumn]==1) {
                            double value = rowElements[kk];
                            if (direction*objective[iColumn]*value>0.0&&!markC[iColumn]) {
                              // Fix
                              if (nstackC0+1<maxStack) {
                                stackC0[nstackC0]=iColumn;
                                if (value<0.0) {
                                  lo0[nstackC0]=colUpper[iColumn];
                                  up0[nstackC0]=colUpper[iColumn];
                                } else {
                                  lo0[nstackC0]=colLower[iColumn];
                                  up0[nstackC0]=colLower[iColumn];
                                }
                                nstackC0++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  if (sum+gap*colsol[j]<minR[irow]-primalTolerance_||(info.strengthenRow&&rowUpper[irow]>1.0e20)) {
                    // can be a cut
                    // add gap to lower and integer coefficient
                    // saveU and saveL spare
                    int * index = (int *)saveL;
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]+gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0));
                    rc.setUb(DBL_MAX);
                    // effectiveness is how far j moves
                    rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      //if(info.strengthenRow)
                      //printf("b point to row %d\n",irow);
                      rowCut.addCutIfNotDuplicate(rc,rowUpper[irow]>1.0e20 ? irow : -1);
                    }
                  }
                }
              }
              minR[irow]=saveMin[istackR];
              maxR[irow]=saveMax[istackR];
              markR[irow]=-1;
            }
          } else {
            if (iway==1&&feasible==3) {
              iway=3;
              // look for singletons that can move (just at root)
              if ((rowCuts&2)!=0&&goingToTrueBound&&info.strengthenRow) {
                for (istackR=0;istackR<nstackR;istackR++) {
                  int irow=stackR[istackR];
                  double gap = rowUpper[irow]-maxR[irow];
                  if (gap>primalTolerance_) {
                    // also see if singletons can go to good objective
                    bool moveSingletons=true;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      int iColumn = column[kk];
                      if (moveSingletons&&j!=iColumn) {
                        if (colUpper[iColumn]>colLower[iColumn]) {
                          if (columnLength[iColumn]!=1) {
                            moveSingletons=false;
                          }
                        }
                      }
                    }
                    if (moveSingletons) {
                      // can fix any with good costs
                      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                           kk++) {
                        int iColumn = column[kk];
                        if (j!=iColumn) {
                          if (colUpper[iColumn]>colLower[iColumn]) {
                            if (columnLength[iColumn]==1) {
                              double value = rowElements[kk];
                              if (direction*objective[iColumn]*value<0.0&&!markC[iColumn]) {
                                // Fix
                                stackC[nstackC]=iColumn;
                                saveL[nstackC]=colLower[iColumn];
                                saveU[nstackC]=colUpper[iColumn];
                                assert (saveU[nstackC]>saveL[nstackC]);
                                if (value>0.0) {
                                  colLower[iColumn]=colUpper[iColumn];
                                } else {
                                  colUpper[iColumn]=colLower[iColumn];
                                }
                                assert (nstackC<nCols);
                                nstackC++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  gap = minR[irow]-rowLower[irow];
                  if (gap>primalTolerance_) {
                    // also see if singletons can go to good objective
                    bool moveSingletons=true;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      int iColumn = column[kk];
                      if (moveSingletons&&j!=iColumn) {
                        if (colUpper[iColumn]>colLower[iColumn]) {
                          if (columnLength[iColumn]!=1) {
                            moveSingletons=false;
                          }
                        }
                      }
                    }
                    if (moveSingletons) {
                      // can fix any with good costs
                      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                        int iColumn = column[kk];
                        if (j!=iColumn) {
                          if (colUpper[iColumn]>colLower[iColumn]) {
                            if (columnLength[iColumn]==1) {
                              double value = rowElements[kk];
                              if (direction*objective[iColumn]*value>0.0&&!markC[iColumn]) {
                                // Fix
                                stackC[nstackC]=iColumn;
                                saveL[nstackC]=colLower[iColumn];
                                saveU[nstackC]=colUpper[iColumn];
                                assert (saveU[nstackC]>saveL[nstackC]);
                                if (value<0.0) {
                                  colLower[iColumn]=colUpper[iColumn];
                                } else {
                                  colUpper[iColumn]=colLower[iColumn];
                                }
                                assert (nstackC<nCols);
                                nstackC++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              /* point back to stack */
              for (istackC=nstackC-1;istackC>=0;istackC--) {
                int icol=stackC[istackC];
                markC[icol]=istackC+1000;
              }
              OsiColCut cc;
              int nTot=0,nFix=0,nInt=0;
              bool ifCut=false;
              for (istackC=1;istackC<nstackC0;istackC++) {
                int icol=stackC0[istackC];
                int istackC1=markC[icol]-1000;
                if (istackC1>=0) {
                  if (CoinMin(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
                    saveL[istackC1]=CoinMin(lo0[istackC],colLower[icol]);
                    if (intVar[icol]||!info.inTree) {
                      element[nFix]=saveL[istackC1];
                      index[nFix++]=icol;
                      nInt++;
                      if (colsol[icol]<saveL[istackC1]-primalTolerance_)
                        ifCut=true;
                    }
                    nfixed++;
                  } 
                }
              }
              if (nFix) {
                nTot=nFix;
                cc.setLbs(nFix,index,element);
                nFix=0;
              }
              for (istackC=1;istackC<nstackC0;istackC++) {
                int icol=stackC0[istackC];
                int istackC1=markC[icol]-1000;
                if (istackC1>=0) {
                  if (CoinMax(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
                    saveU[istackC1]=CoinMax(up0[istackC],colUpper[icol]);
                    if (intVar[icol]||!info.inTree) {
                      element[nFix]=saveU[istackC1];
                      index[nFix++]=icol;
                      nInt++;
                      if (colsol[icol]>saveU[istackC1]+primalTolerance_)
                        ifCut=true;
                    }
                    nfixed++;
                  } else if (!info.inTree&&saveL[0]==0.0&&saveU[0]==1.0) {
                    // See if can do two cut
                    double upperWhenDown = up0[istackC];
                    double lowerWhenDown = lo0[istackC];
                    double upperWhenUp = colUpper[icol];
                    double lowerWhenUp = colLower[icol];
                    double upperOriginal = saveU[istackC1];
                    double lowerOriginal = saveL[istackC1];
                    if (upperWhenDown<lowerOriginal+1.0e-12&&lowerWhenUp>upperOriginal-1.0e-12) {
                      OsiRowCut rc;
                      rc.setLb(lowerOriginal);
                      rc.setUb(lowerOriginal);
                      rc.setEffectiveness(1.0e-5);
                      int index[2];
                      double element[2];
                      index[0]=j;
                      index[1]=icol;
                      element[0]=-(upperOriginal-lowerOriginal);
                      element[1]=1.0;
                      rc.setRow(2,index,element);
                      cs.insert(rc);
                    } else if (upperWhenUp<lowerOriginal+1.0e-12&&lowerWhenDown>upperOriginal-1.0e-12) {
                      OsiRowCut rc;
                      rc.setLb(upperOriginal);
                      rc.setUb(upperOriginal);
                      rc.setEffectiveness(1.0e-5);
                      int index[2];
                      double element[2];
                      index[0]=j;
                      index[1]=icol;
                      element[0]=upperOriginal-lowerOriginal;
                      element[1]=1.0;
                      rc.setRow(2,index,element);
                      cs.insert(rc);
                    } 
                  }
                }
              }
              if (nFix) {
                nTot+=nFix;
                cc.setUbs(nFix,index,element);
              }
              // could tighten continuous as well
              if (nInt) {
                if (ifCut) {
                  cc.setEffectiveness(100.0);
                } else {
                  cc.setEffectiveness(1.0e-5);
                }
#ifdef CGL_DEBUG
                checkBounds(debugger,cc);
#endif
                cs.insert(cc);
              }
            } else {
              goingToTrueBound=0;
            }
            double solMove = up-saveSolval;
            double boundChange;
            /* restore all */
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              if(goingToTrueBound==2&&istackC) {
                // upper disaggregation cut would be
                // xval < upper + (old_upper-upper) (up-jval)
                boundChange = oldU-colUpper[icol];
                if (boundChange>0.0&&oldU<1.0e10&&
                    (colsol[icol]>colUpper[icol]
                     + boundChange*solMove+primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  rc.setLb(-DBL_MAX);
                  rc.setUb(colUpper[icol]+up*boundChange);
                  index[0]=icol;
                  element[0]=1.0;
                  index[1]=j;
                  element[1]= + boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol]-colUpper[icol])/
                    boundChange;
                  assert(newSol>solMove);
                  rc.setEffectiveness(newSol-solMove);
                  if (rc.effectiveness()>disaggEffectiveness) {
                    rc.setRow(2,index,element);
#ifdef CGL_DEBUG
                    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                    rowCut.addCutIfNotDuplicate(rc);
                  }
                }
                // lower disaggregation cut would be
                // xval > lower + (old_lower-lower) (up-jval)
                boundChange = oldL-colLower[icol];
                if (boundChange<0.0&&oldL>-1.0e10&&
                    (colsol[icol]<colLower[icol]
                     + boundChange*solMove-primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  rc.setLb(colLower[icol]+up*boundChange);
                  rc.setUb(DBL_MAX);
                  index[0]=icol;
                  element[0]=1.0;
                  index[1]=j;
                  element[1]= + boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol]-colLower[icol])/
                    boundChange;
                  assert(newSol>solMove);
                  rc.setEffectiveness(newSol-solMove);
                  if (rc.effectiveness()>disaggEffectiveness) {
                    rc.setRow(2,index,element);
#ifdef CGL_DEBUG
                    if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                    rowCut.addCutIfNotDuplicate(rc);
                  }
                }
              }
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              if (oldU>oldL+1.0e-4)
                markC[icol]=0;
              else
                markC[icol]=3;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0&&goingToTrueBound) {
		bool canReplace = info.strengthenRow&&(goingToTrueBound==2);
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    sum += rowElements[kk]*colsol[column[kk]];
                  }
                  if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_||(canReplace&&rowLower[irow]<-1.e20)) {
                    // can be a cut
                    // add gap to integer coefficient
                    // saveU and saveL spare
                    int * index = (int *)saveL;
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]+gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(-DBL_MAX);
                    rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0));
                    // effectiveness is how far j moves
                    rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      //if(canReplace)
                      //printf("c point to row %d\n",irow);
                      rowCut.addCutIfNotDuplicate(rc,(canReplace&&rowLower[irow]<-1.0e20) ? irow : -1);
                    }
                  }
                }
                gap = minR[irow]-rowLower[irow];
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  if (!sum) {
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      sum += rowElements[kk]*colsol[column[kk]];
                    }
                  }
                  if (sum-gap*colsol[j]<rowLower[irow]-primalTolerance_||(canReplace&&rowUpper[irow]>1.0e20)) {
                    // can be a cut
                    // subtract gap from integer coefficient
                    // saveU and saveL spare
                    int * index = (int *)saveL;
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]-gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=-gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(rowLower[irow]-gap*(colUpper[j]-1));
                    rc.setUb(DBL_MAX);
                    // effectiveness is how far j moves
                    rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      //if(canReplace)
                      //printf("d point to row %d\n",irow);
                      rowCut.addCutIfNotDuplicate(rc,(canReplace&&rowUpper[irow]>1.0e20) ? irow : -1);
                    }
                  }
                }
              }
              minR[irow]=saveMin[istackR];
              maxR[irow]=saveMax[istackR];
              markR[irow]=-1;
            }
          }
        }
      }
    }
  }
  if ((!ninfeas&&!rowCut.outOfSpace())&&(info.strengthenRow||
                 !rowCut.numberCuts())) {
    // Try and find ALL big M's
    for (i = 0; i < nRowsSafe; ++i) {
      if ((rowLower[i]>-1.0e20||rowUpper[i]<1.0e20)&&
          (!info.strengthenRow||!info.strengthenRow[i])) {
	int iflagu = 0;
	int iflagl = 0;
	double dmaxup = 0.0;
	double dmaxdown = 0.0;
	int krs = rowStart[i];
	int kre = rowStart[i]+rowLength[i];
        int kInt = -1;
        double valueInteger=0.0;
        // Find largest integer coefficient
        for (k = krs; k < kre; ++k) {
          j = column[k];
          if (intVar[j]) {
            double value=rowElements[k];
            if (colUpper[j]>colLower[j]&&!colLower[j]&&
                fabs(value)>fabs(valueInteger)) {
              kInt=j;
              valueInteger=value;
            }
          }
        }
        if (kInt>=0) {
          double upperBound = CoinMin(colUpper[kInt],(double) INT_MAX);
          for (k = krs; k < kre; ++k) {
            double value=rowElements[k];
            j = column[k];
            if (colUpper[j]==colLower[j]) {
              dmaxup += colUpper[j]*value;
              dmaxdown += colUpper[j]*value;
              continue;
            }
            if (j!=kInt) {
              // continuous
              if (value > 0.0) {
                if (colUpper[j] >= 1e15) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colLower[j] * value;
                }
              } else if (value<0.0) {
                if (colUpper[j] >= 1e15) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colLower[j] * value;
                }
              }
            }
          }
          // end of row
          if (iflagu)
            dmaxup=1.0e31;
          if (iflagl)
            dmaxdown=-1.0e31;
          if (dmaxdown+valueInteger*upperBound>rowLower[i]&&
              dmaxup+valueInteger*upperBound<rowUpper[i]) {
            // check to see if always feasible at 1 but not always at 0
            if (dmaxdown+valueInteger>rowLower[i]&&dmaxup+valueInteger<rowUpper[i]&&
                (dmaxdown<rowLower[i]-primalTolerance_||dmaxup>rowUpper[i]+primalTolerance_)) {
              // can tighten (maybe)
              double saveValue = valueInteger;
              if (valueInteger>0.0) {
                assert (dmaxdown<rowLower[i]);
                valueInteger = rowLower[i]-dmaxdown;
              } else {
                assert (dmaxup>rowUpper[i]);
                valueInteger = rowUpper[i]-dmaxup;
              }
              if (fabs(saveValue-valueInteger)>1.0e-12) {
                // take
                OsiRowCut rc;
                rc.setLb(rowLower[i]);
                rc.setUb(rowUpper[i]);
                int n=0;
                double sum=0.0;
                for (int kk=rowStart[i];kk<rowStart[i]+rowLength[i];kk++) {
                  int j=column[kk];
                  if (j!=kInt) {
                    sum += colsol[j]*rowElements[kk];
                    index[n]=j;
                    element[n++]=rowElements[kk];
                  } else {
                    sum += colsol[j]*valueInteger;
                    assert (rowElements[kk]*valueInteger>=0.0);
#if 0
                    if (fabs(rowElements[kk])>1.01*fabs(valueInteger)) {
                      printf("row %d changing coefficient of %d from %g to %g\n",
                             i,kInt,rowElements[kk],valueInteger);
                    }
#endif
                    if (fabs(valueInteger)>1.0e-12) {
                      index[n]=column[kk];
                      element[n++]=valueInteger;
                    }
                  }
                }
                double gap = 0.0;
                if (sum<rowLower[i])
                  gap=rowLower[i]-sum;
                else if (sum>rowUpper[i])
                  gap=sum-rowUpper[i];
                if (gap>1.0e-4||info.strengthenRow!=NULL) {
                  rc.setEffectiveness(gap);
                  rc.setRow(n,index,element);
                  int returnCode=rowCut.addCutIfNotDuplicate(rc,i);
                  if (returnCode<0)
                    break; // out of space
                }
              }
            }
          }
        }
      }
    }
  }
  delete [] stackC0;
  delete [] lo0;
  delete [] up0;
  delete [] tempL;
  delete [] tempU;
  delete [] markC;
  delete [] stackC;
  delete [] stackR;
  delete [] saveL;
  delete [] saveU;
  delete [] saveMin;
  delete [] saveMax;
  delete [] index;
  delete [] element;
  delete [] djs;
  delete [] colsol;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info.strengthenRow);
  }
  if (!columnCopy_)
    delete columnCopy; // if just created
  return (ninfeas);
}
// Does probing and adding cuts
int CglProbing::probeCliques( const OsiSolverInterface & si, 
                              const OsiRowCutDebugger * debugger, 
                              OsiCuts & cs, 
                              double * colLower, double * colUpper, 
		       CoinPackedMatrix *rowCopy,
		       double * rowLower, double * rowUpper,
		       char * intVar, double * minR, double * maxR, 
		       int * markR, 
		       int * look, int nLook,
                       const CglTreeInfo info) const
{
  // Set up maxes
  int maxStack = info.inTree ? maxStack_ : maxStackRoot_;
  int nRows=rowCopy->getNumRows();
  int nCols=rowCopy->getNumCols();
  double * colsol = new double[nCols];
  double * djs = new double[nCols];
  const double * currentColLower = si.getColLower();
  const double * currentColUpper = si.getColUpper();
  double * tempL = new double [nCols];
  double * tempU = new double [nCols];
  int * markC = new int [nCols];
  int * stackC = new int [2*nCols];
  int * stackR = new int [nRows];
  double * saveL = new double [2*nCols];
  double * saveU = new double [2*nCols];
  double * saveMin = new double [nRows];
  double * saveMax = new double [nRows];
  double * element = new double[nCols];
  int * index = new int[nCols];
  // For trying to extend cliques
  int * cliqueStack=NULL;
  int * cliqueCount=NULL;
  int * to_01=NULL;
  if (!mode_) {
    to_01 = new int[nCols];
    cliqueStack = new int[numberCliques_];
    cliqueCount = new int[numberCliques_];
    int i;
    for (i=0;i<numberCliques_;i++) {
      cliqueCount[i]=cliqueStart_[i+1]-cliqueStart_[i];
    }
    for (i=0;i<nCols;i++) 
      to_01[i]=-1;
    for (i=0;i<number01Integers_;i++) {
      int j=cutVector_[i].sequence;
      to_01[j]=i;
    }
  }
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info.inTree ? nRows/3 : nRows;
  row_cut rowCut(nRowsFake);
  CoinPackedMatrix * columnCopy=NULL;
  if (!rowCopy_) {
    columnCopy = new CoinPackedMatrix(*rowCopy);
    columnCopy->reverseOrdering();
  } else {
    columnCopy = columnCopy_;
  }
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements();
  double movement;
  int i, j, k,kk,jj;
  int jcol,kcol,irow,krow;
  bool anyColumnCuts=false;
  double dbound, value, value2;
  int ninfeas=0;
  int rowCuts;
  double disaggEffectiveness;
  if (mode_) {
    /* clean up djs and solution */
    memcpy(djs,si.getReducedCost(),nCols*sizeof(double));
    memcpy(colsol, si.getColSolution(),nCols*sizeof(double));
    disaggEffectiveness=1.0e-3;
    rowCuts=rowCuts_;
  } else {
    // need to go from a neutral place
    memset(djs,0,nCols*sizeof(double));
    memcpy(colsol, si.getColSolution(),nCols*sizeof(double));
    disaggEffectiveness=-1.0e10;
    if (rowCuts_!=4)
      rowCuts=1;
    else
      rowCuts=4;
  }
  for (i = 0; i < nCols; ++i) {
    /* was if (intVar[i]) */
    if (1) {
      if (colUpper[i]-colLower[i]>1.0e-8) {
	if (colsol[i]<colLower[i]+primalTolerance_) {
	  colsol[i]=colLower[i];
	  djs[i] = CoinMax(0.0,djs[i]);
	} else if (colsol[i]>colUpper[i]-primalTolerance_) {
	  colsol[i]=colUpper[i];
	  djs[i] = CoinMin(0.0,djs[i]);
	} else {
	  djs[i]=0.0;
	}
	/*if (fabs(djs[i])<1.0e-5) 
	  djs[i]=0.0;*/
      }
    }
  }

  int ipass=0,nfixed=-1;

  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff);
  if (!cutoff_available) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  cutoff *= si.getObjSense();
  if (fabs(cutoff)>1.0e30)
    assert (cutoff>1.0e30);
  double current = si.getObjValue();
  // make irrelevant if mode is 0
  if (!mode_)
    cutoff=DBL_MAX;
  /* for both way coding */
  int nstackC0=-1;
  int * stackC0 = new int[maxStack];
  double * lo0 = new double[maxStack];
  double * up0 = new double[maxStack];
  int nstackR,nstackC;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3;
    } else {
      markC[i]=0;
    }
  }
  double tolerance = 1.0e1*primalTolerance_;
  int maxPass = info.inTree ? maxPass_ : maxPassRoot_;
  // If we are going to replace coefficient then we don't need to be effective
  double needEffectiveness = info.strengthenRow ? -1.0e10 : 1.0e-3;
  while (ipass<maxPass&&nfixed) {
    int iLook;
    ipass++;
    nfixed=0;
    for (iLook=0;iLook<nLook;iLook++) {
      double solval;
      double down;
      double up;
      int awayFromBound=1;
      j=look[iLook];
      solval=colsol[j];
      down = floor(solval+tolerance);
      up = ceil(solval-tolerance);
      if(colUpper[j]-colLower[j]<1.0e-8) markC[j]=3;
      if (markC[j]||!intVar[j]) continue;
      double saveSolval = solval;
      if (solval>=colUpper[j]-tolerance||solval<=colLower[j]+tolerance||up==down) {
	awayFromBound=0;
	if (solval<=colLower[j]+2.0*tolerance) {
	  solval = colLower[j]+1.0e-1;
	  down=colLower[j];
	  up=down+1.0;
	} else if (solval>=colUpper[j]-2.0*tolerance) {
	  solval = colUpper[j]-1.0e-1;
	  up=colUpper[j];
	  down=up-1;
	} else {
          // odd
          up=down+1.0;
          solval = down+1.0e-1;
        }
      }
      assert (up<=colUpper[j]);
      assert (down>=colLower[j]);
      assert (up>down);
      if ((solval-down>1.0e-6&&up-solval>1.0e-6)||mode_!=1) {
	int istackC,iway, istackR;
	int way[]={1,2,1};
	int feas[]={1,2,4};
	int feasible=0;
	int notFeasible;
	for (iway=0;iway<3;iway ++) {
	  int fixThis=0;
	  double objVal=current;
	  int goingToTrueBound=0;
	  stackC[0]=j;
	  markC[j]=way[iway];
          double solMovement;
	  if (way[iway]==1) {
	    movement=down-colUpper[j];
            solMovement = down-colsol[j];
	    assert(movement<-0.99999);
	    if (fabs(down-colLower[j])<1.0e-7) {
	      goingToTrueBound=2;
	      down=colLower[j];
	    }
	  } else {
	    movement=up-colLower[j];
            solMovement = up-colsol[j];
	    assert(movement>0.99999);
	    if (fabs(up-colUpper[j])<1.0e-7) {
	      goingToTrueBound=2;
	      up=colUpper[j];
	    }
	  }
	  if (goingToTrueBound&&(colUpper[j]-colLower[j]>1.5||colLower[j]))
	    goingToTrueBound=1;
	  // switch off disaggregation if not wanted
	  if ((rowCuts&1)==0)
	    goingToTrueBound=0;
#ifdef PRINT_DEBUG
	  if (fabs(movement)>1.01) {
	    printf("big %d %g %g %g\n",j,colLower[j],solval,colUpper[j]);
	  }
#endif
	  if (solMovement*djs[j]>0.0)
	    objVal += solMovement*djs[j];
	  nstackC=1;
	  nstackR=0;
	  saveL[0]=colLower[j];
	  saveU[0]=colUpper[j];
          assert (saveU[0]>saveL[0]);
	  notFeasible=0;
	  if (movement<0.0) {
	    colUpper[j] += movement;
	    colUpper[j] = floor(colUpper[j]+0.5);
#ifdef PRINT_DEBUG
	    printf("** Trying %d down to 0\n",j);
#endif
	  } else {
	    colLower[j] += movement;
	    colLower[j] = floor(colLower[j]+0.5);
#ifdef PRINT_DEBUG
	    printf("** Trying %d up to 1\n",j);
#endif
	  }
	  if (fabs(colUpper[j]-colLower[j])<1.0e-6)
	    markC[j]=3; // say fixed
	  istackC=0;
	  /* update immediately */
	  for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
	    irow = row[k];
	    value = columnElements[k];
	    if (markR[irow]==-1) {
	      stackR[nstackR]=irow;
	      markR[irow]=nstackR;
	      saveMin[nstackR]=minR[irow];
	      saveMax[nstackR]=maxR[irow];
	      nstackR++;
	    } else if (markR[irow]==-2) {
	      continue;
	    }
	    /* could check immediately if violation */
	    if (movement>0.0) {
	      /* up */
	      if (value>0.0) {
		/* up does not change - down does */
                if (minR[irow]>-1.0e10)
                  minR[irow] += value;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      } else {
		/* down does not change - up does */
                if (maxR[irow]<1.0e10)
                  maxR[irow] += value;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      }
	    } else {
	      /* down */
	      if (value<0.0) {
		/* up does not change - down does */
                if (minR[irow]>-1.0e10)
                  minR[irow] -= value;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      } else {
		/* down does not change - up does */
                if (maxR[irow]<1.0e10)
                  maxR[irow] -= value;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1;
		  istackC=1;
		  break;
		}
	      }
	    }
	  }
	  while (istackC<nstackC&&nstackC<maxStack) {
	    int jway;
	    jcol=stackC[istackC];
	    jway=markC[jcol];
	    // If not first and fixed then skip
	    if (jway==3&&istackC) {
	      //istackC++;
	      //continue;
              //printf("fixed %d on stack\n",jcol);
	    }
	    // Do cliques
	    if (oneFixStart_&&oneFixStart_[jcol]>=0) {
	      int start;
	      int end;
	      if (colLower[jcol]>saveL[istackC]) {
		// going up
		start = oneFixStart_[jcol];
		end = zeroFixStart_[jcol];
	      } else {
		assert (colUpper[jcol]<saveU[istackC]);
		// going down
		start = zeroFixStart_[jcol];
		end = endFixStart_[jcol];
	      }
	      for (int i=start;i<end;i++) {
		int iClique = whichClique_[i];
		for (int k=cliqueStart_[iClique];k<cliqueStart_[iClique+1];k++) {
		  int kcol = cliqueEntry_[k].sequence;
                  if (jcol==kcol)
                    continue;
		  int kway = cliqueEntry_[k].oneFixes;
                  if (kcol!=jcol) {
                    if (!markC[kcol]) {
                      // not on list yet
                      if (nstackC<2*maxStack) {
                        markC[kcol] = 3; // say fixed
                        fixThis++;
                        stackC[nstackC]=kcol;
                        saveL[nstackC]=colLower[kcol];
                        saveU[nstackC]=colUpper[kcol];
                        assert (saveU[nstackC]>saveL[nstackC]);
                        nstackC++;
                        if (!kway) {
                          // going up
                          double solMovement=1.0-colsol[kcol];
                          if (solMovement>0.0001) {
                            assert (djs[kcol]>=0.0);
                            objVal += djs[kcol]*solMovement;
                          }
                          colLower[kcol]=1.0;
                          /* update immediately */
                          for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                            krow = row[jj];
                            value = columnElements[jj];
                            if (markR[krow]==-1) {
                              stackR[nstackR]=krow;
                              markR[krow]=nstackR;
                              saveMin[nstackR]=minR[krow];
                              saveMax[nstackR]=maxR[krow];
                              nstackR++;
                            } else if (markR[krow]==-2) {
                              continue;
                            }
                            /* could check immediately if violation */
                            /* up */
                            if (value>0.0) {
                              /* up does not change - down does */
                              if (minR[krow]>-1.0e10)
                                minR[krow] += value;
                              if (minR[krow]>rowUpper[krow]+1.0e-5) {
                                colUpper[kcol]=-1.0e50; /* force infeasible */
                                break;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow]<1.0e10)
                                maxR[krow] += value;
                              if (maxR[krow]<rowLower[krow]-1.0e-5) {
                                notFeasible=1;
                                break;
                              }
                            }
                          }
                        } else {
                          // going down
                          double solMovement=0.0-colsol[kcol];
                          if (solMovement<-0.0001) {
                            assert (djs[kcol]<=0.0);
                            objVal += djs[kcol]*solMovement;
                          }
                          colUpper[kcol]=0.0;
                          /* update immediately */
                          for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                            krow = row[jj];
                            value = columnElements[jj];
                            if (markR[krow]==-1) {
                              stackR[nstackR]=krow;
                              markR[krow]=nstackR;
                              saveMin[nstackR]=minR[krow];
                              saveMax[nstackR]=maxR[krow];
                              nstackR++;
                            } else if (markR[krow]==-2) {
                              continue;
                            }
                            /* could check immediately if violation */
                            /* down */
                            if (value<0.0) {
                              /* up does not change - down does */
                              if (minR[krow]>-1.0e10)
                                minR[krow] -= value;
                              if (minR[krow]>rowUpper[krow]+1.0e-5) {
                                notFeasible=1;
                                break;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow]<1.0e10)
                                maxR[krow] -= value;
                              if (maxR[krow]<rowLower[krow]-1.0e-5) {
                                notFeasible=1;
                                break;
                              }
                            }
                          }
                        }
                      }
                    } else if (markC[kcol]==1) {
                      // marked as going to 0
                      assert (!colUpper[kcol]);
                      if (!kway) {
                        // contradiction
                        notFeasible=1;
                        break;
                      }
                    } else if (markC[kcol]==2) {
                      // marked as going to 1
                      assert (colLower[kcol]);
                      if (kway) {
                        // contradiction
                        notFeasible=1;
                        break;
                      }
                    } else {
                      // marked as fixed
                      assert (markC[kcol]==3);
                      int jkway;
                      if (colLower[kcol])
                        jkway=1;
                      else
                        jkway=0;
                      if (kway==jkway) {
                        // contradiction
                        notFeasible=1;
                        break;
                      }
                    }
                  }
		}
		if (notFeasible)
		  break;
	      }
	      if (notFeasible)
		istackC=nstackC+1;
	    }
	    for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
	      // break if found not feasible
	      if (notFeasible)
		break;
	      irow = row[k];
	      /*value = columnElements[k];*/
	      if (markR[irow]!=-2) {
		/* see if anything forced */
		for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
		  double moveUp=0.0;
		  double moveDown=0.0;
		  double newUpper=-1.0,newLower=1.0;
		  kcol=column[kk];
		  bool onList = (markC[kcol]!=0);
		  if (markC[kcol]!=3) {
		    value2=rowElements[kk];
                    int markIt=markC[kcol];
		    if (value2 < 0.0) {
		      if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colUpper[kcol]+
			  (rowUpper[irow]-minR[irow])/value2;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
                            markIt |= 2;
			    newLower = ceil(dbound-primalTolerance_);
			  } else {
			    newLower=dbound;
			    if (newLower+primalTolerance_>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol];
                              markIt |= 2;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
                            }
			  }
			  moveUp = newLower-colLower[kcol];
			}
		      }
		      if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colLower[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 1;
			    newUpper = floor(dbound+primalTolerance_);
			  } else {
			    newUpper=dbound;
			    if (newUpper-primalTolerance_<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol];
                              markIt |= 1;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
                            }
			  }
			  moveDown = newUpper-colUpper[kcol];
			}
		      }
		    } else {
		      /* positive element */
		      if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colUpper[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 2;
			    newLower = ceil(dbound-primalTolerance_);
			  } else {
			    newLower=dbound;
			    if (newLower+primalTolerance_>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol];
                              markIt |= 2;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
			    }
			  }
			  moveUp = newLower-colLower[kcol];
			}
		      }
		      if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colLower[kcol] + 
			  (rowUpper[irow]-minR[irow])/value2;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 1;
			    newUpper = floor(dbound+primalTolerance_);
			  } else {
			    newUpper=dbound;
			    if (newUpper-primalTolerance_<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol];
                              markIt |= 1;
                              //markIt=3;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3;
			    }
			  }
			  moveDown = newUpper-colUpper[kcol];
			}
		      }
		    }
		    if (nstackC<2*maxStack) 
                      markC[kcol] = markIt;
		    if (moveUp&&nstackC<2*maxStack) {
		      fixThis++;
#ifdef PRINT_DEBUG
		      printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]);
		      value=0.0;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]);
			} else {
			  value += rowElements[jj]*colLower[ii];
			}
		      }
		      printf(" - fixed %g\n",value);
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
			}
		      }
		      printf("\n");
#endif
		      if (!onList) {
			stackC[nstackC]=kcol;
			saveL[nstackC]=colLower[kcol];
			saveU[nstackC]=colUpper[kcol];
                        assert (saveU[nstackC]>saveL[nstackC]);
			nstackC++;
			onList=true;
		      }
		      if (newLower>colsol[kcol]) {
			if (djs[kcol]<0.0) {
			  /* should be infeasible */
			  assert (newLower>colUpper[kcol]+primalTolerance_);
			} else {
			  objVal += moveUp*djs[kcol];
			}
		      }
		      if (intVar[kcol])
			newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4));
		      colLower[kcol]=newLower;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
			markC[kcol]=3; // say fixed
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj];
			value = columnElements[jj];
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow;
			  markR[krow]=nstackR;
			  saveMin[nstackR]=minR[krow];
			  saveMax[nstackR]=maxR[krow];
			  nstackR++;
			} else if (markR[krow]==-2) {
			  continue;
			}
			/* could check immediately if violation */
			/* up */
			if (value>0.0) {
			  /* up does not change - down does */
                          if (minR[krow]>-1.0e10)
                            minR[krow] += value*moveUp;
			  if (minR[krow]>rowUpper[krow]+1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			} else {
			  /* down does not change - up does */
                          if (maxR[krow]<1.0e10)
                            maxR[krow] += value*moveUp;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			}
		      }
		    }
		    if (moveDown&&nstackC<2*maxStack) {
		      fixThis++;
#ifdef PRINT_DEBUG
		      printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]);
		      value=0.0;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]);
			} else {
			  value += rowElements[jj]*colLower[ii];
			}
		      }
		      printf(" - fixed %g\n",value);
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj];
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
			}
		      }
		      printf("\n");
#endif
		      if (!onList) {
			stackC[nstackC]=kcol;
			saveL[nstackC]=colLower[kcol];
			saveU[nstackC]=colUpper[kcol];
                        assert (saveU[nstackC]>saveL[nstackC]);
			nstackC++;
			onList=true;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_);
			} else {
			  objVal += moveDown*djs[kcol];
			}
		      }
		      if (intVar[kcol])
			newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4));
		      colUpper[kcol]=newUpper;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
			markC[kcol]=3; // say fixed
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj];
			value = columnElements[jj];
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow;
			  markR[krow]=nstackR;
			  saveMin[nstackR]=minR[krow];
			  saveMax[nstackR]=maxR[krow];
			  nstackR++;
			} else if (markR[krow]==-2) {
			  continue;
			}
			/* could check immediately if violation */
			/* down */
			if (value<0.0) {
			  /* up does not change - down does */
                          if (minR[krow]>-1.0e10)
                            minR[krow] += value*moveDown;
			  if (minR[krow]>rowUpper[krow]+1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			} else {
			  /* down does not change - up does */
                          if (maxR[krow]<1.0e10)
                            maxR[krow] += value*moveDown;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break;
			  }
			}
		      }
		    }
		    if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
		      notFeasible=1;;
		      k=columnStart[jcol]+columnLength[jcol];
		      istackC=nstackC+1;
#ifdef PRINT_DEBUG
		      printf("** not feasible this way\n");
#endif
		      break;
		    }
		  }
		}
	      }
	    }
	    istackC++;
	  }
	  if (!notFeasible) {
	    if (objVal<=cutoff) {
	      feasible |= feas[iway];
	    } else {
#ifdef PRINT_DEBUG
	      printf("not feasible on dj\n");
#endif
	      notFeasible=1;
	      if (iway==1&&feasible==0) {
		/* not feasible at all */
		ninfeas=1;
		j=nCols-1;
		break;
	      }
	    }
	  } else if (iway==1&&feasible==0) {
	    /* not feasible at all */
	    ninfeas=1;
	    j=nCols-1;
            iLook=nLook;
	    ipass=maxPass;
	    break;
	  }
	  if (notFeasible)
	    goingToTrueBound=0;
	  if (iway==2||(iway==1&&feasible==2)) {
	    /* keep */
	    iway=3;
	    nfixed++;
            if (mode_) {
	      OsiColCut cc;
	      int nTot=0,nFix=0,nInt=0;
	      bool ifCut=false;
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC];
		if (intVar[icol]) {
		  if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
		    element[nFix]=colUpper[icol];
		    index[nFix++]=icol;
		    nInt++;
		    if (colsol[icol]>colUpper[icol]+primalTolerance_) {
		      ifCut=true;
		      anyColumnCuts=true;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot=nFix;
		cc.setUbs(nFix,index,element);
		nFix=0;
	      }
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC];
		if (intVar[icol]) {
		  if (colLower[icol]>currentColLower[icol]+1.0e-4) {
		    element[nFix]=colLower[icol];
		    index[nFix++]=icol;
		    nInt++;
		    if (colsol[icol]<colLower[icol]-primalTolerance_) {
		      ifCut=true;
		      anyColumnCuts=true;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot+=nFix;
		cc.setLbs(nFix,index,element);
	      }
	      // could tighten continuous as well
	      if (nInt) {
		if (ifCut) {
		  cc.setEffectiveness(100.0);
		} else {
		  cc.setEffectiveness(1.0e-5);
		}
#ifdef CGL_DEBUG
		checkBounds(debugger,cc);
#endif
		cs.insert(cc);
	      }
	    }
	    for (istackC=0;istackC<nstackC;istackC++) {
	      int icol=stackC[istackC];
	      if (colUpper[icol]-colLower[icol]>primalTolerance_) {
		markC[icol]=0;
	      } else {
		markC[icol]=3;
	      }
	    }
	    for (istackR=0;istackR<nstackR;istackR++) {
	      int irow=stackR[istackR];
	      markR[irow]=-1;
	    }
	  } else {
	    /* is it worth seeing if can increase coefficients
	       or maybe better see if it is a cut */
	    if (iway==0) {
	      nstackC0=CoinMin(nstackC,maxStack);
	      double solMove = saveSolval-down;
	      double boundChange;
	      if (notFeasible) {
		nstackC0=0;
	      } else {
		for (istackC=0;istackC<nstackC0;istackC++) {
		  int icol=stackC[istackC];
		  stackC0[istackC]=icol;
		  lo0[istackC]=colLower[icol];
		  up0[istackC]=colUpper[icol];
		}
	      }
	      /* restore all */
              int nCliquesAffected=0;
              assert (iway==0);
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC];
		double oldU=saveU[istackC];
		double oldL=saveL[istackC];
		if(goingToTrueBound==2&&istackC) {
                  // Work for extending cliques
                  if (!mode_&&numberCliques_) {
                    int i_01 = to_01[icol];
                    if (i_01>=0) {
                      int start;
                      int end;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol];
                        end = endFixStart_[icol];
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol];
                        end = zeroFixStart_[icol];
                      }
                      //if (end>start)
                      //printf("j %d, other %d is in %d cliques\n",
                      //     j,i_01,end-start);
                      for (int i=start;i<end;i++) {
                        int iClique = whichClique_[i];
                        int size = cliqueStart_[iClique+1]-cliqueStart_[iClique];
                        if (cliqueCount[iClique]==size) {
                          // first time
                          cliqueStack[nCliquesAffected++]=iClique;
                        }
                        // decrement counts
                        cliqueCount[iClique]--;
                      }
                    }
                  }
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (jval-down)
		  boundChange = oldU-colUpper[icol];
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(-DBL_MAX);
		    rc.setUb(colUpper[icol]-down*boundChange);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= - boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange;
		    if (mode_) 
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (jval-down)
		  boundChange = oldL-colLower[icol];
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(colLower[icol]-down*boundChange);
		    rc.setUb(DBL_MAX);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]=- boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
#if 0
		      printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,oldL,oldU,colLower[icol],colsol[icol], j, colsol[j]);
		      printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
			     icol,boundChange,j,rc.ub());
#endif
		    }
		  }
		}
		colUpper[icol]=oldU;
		colLower[icol]=oldL;
		markC[icol]=0;
	      }
              if (nCliquesAffected) {
                for (int i=0;i<nCliquesAffected;i++) {
                  int iClique = cliqueStack[i];
                  int size = cliqueCount[iClique];
                  // restore
                  cliqueCount[iClique]= cliqueStart_[iClique+1]-cliqueStart_[iClique];
                  if (!size) {
                    if (logLevel_>1)
                      printf("** could extend clique by adding j!\n");
                  }
                }
              }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR];
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts;
		  double gap = rowUpper[irow]-maxR[irow];
		  double sum=0.0;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]];
		    }
		    if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_||(info.strengthenRow&&rowLower[irow]<-1.0e20)) {
		      // can be a cut
		      // subtract gap from upper and integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        // If strengthenRow point to row
                        //if(info.strengthenRow)
                        //printf("a point to row %d\n",irow);
			rowCut.addCutIfNotDuplicate(rc,rowLower[irow]<-1.0e20 ? irow : -1);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum+gap*colsol[j]<minR[irow]+primalTolerance_||(info.strengthenRow&&rowUpper[irow]>1.0e20)) {
		      // can be a cut
		      // add gap to lower and integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info.strengthenRow)
                        //printf("b point to row %d\n",irow);
			rowCut.addCutIfNotDuplicate(rc,rowUpper[irow]>1.0e20 ? irow : -1);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    } else {
	      if (iway==1&&feasible==3) {
		iway=3;
		/* point back to stack */
		for (istackC=nstackC-1;istackC>=0;istackC--) {
		  int icol=stackC[istackC];
		  markC[icol]=istackC+1000;
		}
		if (mode_) {
		  OsiColCut cc;
		  int nTot=0,nFix=0,nInt=0;
		  bool ifCut=false;
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC];
		    int istackC1=markC[icol]-1000;
		    if (istackC1>=0) {
		      if (CoinMin(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
			saveL[istackC1]=CoinMin(lo0[istackC],colLower[icol]);
			if (intVar[icol]) {
			  element[nFix]=saveL[istackC1];
			  index[nFix++]=icol;
			  nInt++;
			  if (colsol[icol]<saveL[istackC1]-primalTolerance_)
			    ifCut=true;
			}
			nfixed++;
		      }
		    }
		  }
		  if (nFix) {
		    nTot=nFix;
		    cc.setLbs(nFix,index,element);
		    nFix=0;
		  }
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC];
		    int istackC1=markC[icol]-1000;
		    if (istackC1>=0) {
		      if (CoinMax(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
			saveU[istackC1]=CoinMax(up0[istackC],colUpper[icol]);
			if (intVar[icol]) {
			  element[nFix]=saveU[istackC1];
			  index[nFix++]=icol;
			  nInt++;
			  if (colsol[icol]>saveU[istackC1]+primalTolerance_)
			    ifCut=true;
			}
			nfixed++;
		      }
		    }
		  }
		  if (nFix) {
		    nTot+=nFix;
		    cc.setUbs(nFix,index,element);
		  }
		  // could tighten continuous as well
		  if (nInt) {
		    if (ifCut) {
		      cc.setEffectiveness(100.0);
		    } else {
		      cc.setEffectiveness(1.0e-5);
		    }
#ifdef CGL_DEBUG
		    checkBounds(debugger,cc);
#endif
		    cs.insert(cc);
		  }
		}
	      } else {
		goingToTrueBound=0;
	      }
	      double solMove = up-saveSolval;
	      double boundChange;
	      /* restore all */
              int nCliquesAffected=0;
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC];
		double oldU=saveU[istackC];
		double oldL=saveL[istackC];
		if(goingToTrueBound==2&&istackC) {
                  // Work for extending cliques
                  if (!mode_&&numberCliques_&&iway==3) {
                    int i_01 = to_01[icol];
                    if (i_01>=0) {
                      int start;
                      int end;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol];
                        end = endFixStart_[icol];
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol];
                        end = zeroFixStart_[icol];
                      }
                      //if (end>start)
                      //printf("up j %d, other %d is in %d cliques\n",
                      //     j,i_01,end-start);
                      for (int i=start;i<end;i++) {
                        int iClique = whichClique_[i];
                        int size = cliqueStart_[iClique+1]-cliqueStart_[iClique];
                        if (cliqueCount[iClique]==size) {
                          // first time
                          cliqueStack[nCliquesAffected++]=iClique;
                        }
                        // decrement counts
                        cliqueCount[iClique]--;
                      }
                    }
                  }
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (up-jval)
		  boundChange = oldU-colUpper[icol];
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(-DBL_MAX);
		    rc.setUb(colUpper[icol]+up*boundChange);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= + boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (up-jval)
		  boundChange = oldL-colLower[icol];
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc;
		    rc.setLb(colLower[icol]+up*boundChange);
		    rc.setUb(DBL_MAX);
		    index[0]=icol;
		    element[0]=1.0;
		    index[1]=j;
		    element[1]= + boundChange;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange;
		    if (mode_)
		      assert(newSol>solMove);
		    rc.setEffectiveness(newSol-solMove);
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element);
#ifdef CGL_DEBUG
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc);
		    }
		  }
		}
		colUpper[icol]=oldU;
		colLower[icol]=oldL;
                if (oldU>oldL+1.0e-4)
                  markC[icol]=0;
                else
                  markC[icol]=3;
	      }
              if (nCliquesAffected) {
                for (int i=0;i<nCliquesAffected;i++) {
                  int iClique = cliqueStack[i];
                  int size = cliqueCount[iClique];
                  // restore
                  cliqueCount[iClique]= cliqueStart_[iClique+1]-cliqueStart_[iClique];
                  if (!size) {
                    if (logLevel_>1)
                      printf("** could extend clique by adding j!\n");
                  }
                }
              }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR];
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts;
		  double gap = rowUpper[irow]-maxR[irow];
		  double sum=0.0;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]];
		    }
		    if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_||(info.strengthenRow&&rowLower[irow]<-1.0e20)) {
		      // can be a cut
		      // add gap to integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info.strengthenRow)
                        //printf("c point to row %d\n",irow);
			rowCut.addCutIfNotDuplicate(rc,rowLower[irow]<-1.0e20? irow : -1);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum-gap*colsol[j]<rowLower[irow]+primalTolerance_||(info.strengthenRow&&rowUpper[irow]>1.0e20)) {
		      // can be a cut
		      // subtract gap from integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]-gap*(colUpper[j]-1));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info.strengthenRow)
                        //printf("d point to row %d\n",irow);
			rowCut.addCutIfNotDuplicate(rc,rowUpper[irow]>1.0e20 ? irow : -1);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    }
	  }
	}
      }
    }
  }
  delete [] cliqueStack;
  delete [] cliqueCount;
  delete [] to_01;
  delete [] stackC0;
  delete [] lo0;
  delete [] up0;
  delete [] tempL;
  delete [] tempU;
  delete [] markC;
  delete [] stackC;
  delete [] stackR;
  delete [] saveL;
  delete [] saveU;
  delete [] saveMin;
  delete [] saveMax;
  delete [] index;
  delete [] element;
  delete [] djs;
  delete [] colsol;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info.strengthenRow);
  }
  if (!columnCopy_)
    delete columnCopy; // if just created
  return (ninfeas);
}
// Does probing and adding cuts for clique slacks
int 
CglProbing::probeSlacks( const OsiSolverInterface & si, 
                          const OsiRowCutDebugger * debugger, 
                          OsiCuts & cs, 
                          double * colLower, double * colUpper, CoinPackedMatrix *rowCopy,
                          double * rowLower, double * rowUpper,
                          char * intVar, double * minR, double * maxR,int * markR,
                          const CglTreeInfo info) const
{
  if (!numberCliques_)
    return 0;
  // Set up maxes
  int maxProbe = info.inTree ? maxProbe_ : maxProbeRoot_;
  int maxStack = info.inTree ? maxStack_ : maxStackRoot_;
  int nRows=rowCopy->getNumRows();
  int nCols=rowCopy->getNumCols();
  double * colsol = new double[nCols];
  memcpy(colsol, si.getColSolution(),nCols*sizeof(double));
  int rowCuts=rowCuts_;
  double_int_pair * array = new double_int_pair [numberCliques_];
  // look at <= cliques
  int iClique;
  int nLook=0;
  for (iClique=0;iClique<numberCliques_;iClique++) {
    if (!cliqueType_[iClique].equality) {
      double sum=0.0;
      for (int j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
        int iColumn = cliqueEntry_[j].sequence;
        double value = colsol[iColumn];
        if (cliqueEntry_[j].oneFixes)
          sum += value;
        else
          sum -= value;
      }
      double away = fabs(0.5-(sum-floor(sum)));
      if (away<0.49999) {
        array[nLook].infeasibility=away;
        array[nLook++].sequence=iClique;
      }
    }
  }
  std::sort(array,array+nLook,double_int_pair_compare());
  nLook=CoinMin(nLook,maxProbe);
  const double * currentColLower = si.getColLower();
  const double * currentColUpper = si.getColUpper();
  double * tempL = new double [nCols];
  double * tempU = new double [nCols];
  int * markC = new int [nCols];
  int * stackC = new int [2*nCols];
  int * stackR = new int [nRows];
  double * saveL = new double [2*nCols];
  double * saveU = new double [2*nCols];
  double * saveMin = new double [nRows];
  double * saveMax = new double [nRows];
  double * element = new double[nCols];
  int * index = new int[nCols];
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info.inTree ? nRows/3 : nRows;
  row_cut rowCut(nRowsFake);
  CoinPackedMatrix * columnCopy=NULL;
  if (!rowCopy_) {
    columnCopy = new CoinPackedMatrix(*rowCopy);
    columnCopy->reverseOrdering();
  } else {
    columnCopy = columnCopy_;
  }
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements();
  double movement;
  int i, j, k,kk,jj;
  int jcol,kcol,irow,krow;
  bool anyColumnCuts=false;
  double dbound, value, value2;
  int ninfeas=0;
  for (i = 0; i < nCols; ++i) {
    if (colUpper[i]-colLower[i]>1.0e-8) {
      if (colsol[i]<colLower[i]+primalTolerance_) {
        colsol[i]=colLower[i];
      } else if (colsol[i]>colUpper[i]-primalTolerance_) {
        colsol[i]=colUpper[i];
      }
    }
  }

  int ipass=0,nfixed=-1;

  /* for both way coding */
  int nstackC0=-1;
  int * stackC0 = new int[maxStack];
  double * lo0 = new double[maxStack];
  double * up0 = new double[maxStack];
  int nstackR,nstackC;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3;
    } else {
      markC[i]=0;
    }
  }
  double tolerance = 1.0e1*primalTolerance_;
  // If we are going to replace coefficient then we don't need to be effective
  int maxPass = info.inTree ? maxPass_ : maxPassRoot_;
  double needEffectiveness = info.strengthenRow ? -1.0e10 : 1.0e-3;
  while (ipass<maxPass&&nfixed) {
    int iLook;
    ipass++;
    nfixed=0;
    for (iLook=0;iLook<nLook;iLook++) {
      double solval;
      double down;
      double up;
      int iClique=array[iLook].sequence;
      solval=0.0;
      j=0;
      for (j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
        int iColumn = cliqueEntry_[j].sequence;
        double value = colsol[iColumn];
        if (cliqueEntry_[j].oneFixes)
          solval += value;
        else
          solval -= value;
      }
      down = floor(solval+tolerance);
      up = ceil(solval-tolerance);
      int istackC,iway, istackR;
      int way[]={1,2,1};
      int feas[]={1,2,4};
      int feasible=0;
      int notFeasible;
      for (iway=0;iway<3;iway ++) {
        int fixThis=0;
        stackC[0]=j;
        markC[j]=way[iway];
        if (way[iway]==1) {
          movement=down-colUpper[j];
          assert(movement<-0.99999);
          down=colLower[j];
        } else {
          movement=up-colLower[j];
          assert(movement>0.99999);
          up=colUpper[j];
        }
        nstackC=1;
        nstackR=0;
        saveL[0]=colLower[j];
        saveU[0]=colUpper[j];
        assert (saveU[0]>saveL[0]);
        notFeasible=0;
        if (movement<0.0) {
          colUpper[j] += movement;
          colUpper[j] = floor(colUpper[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d down to 0\n",j);
#endif
        } else {
          colLower[j] += movement;
          colLower[j] = floor(colLower[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d up to 1\n",j);
#endif
        }
        if (fabs(colUpper[j]-colLower[j])<1.0e-6)
          markC[j]=3; // say fixed
        istackC=0;
        /* update immediately */
        for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
          irow = row[k];
          value = columnElements[k];
          if (markR[irow]==-1) {
            stackR[nstackR]=irow;
            markR[irow]=nstackR;
            saveMin[nstackR]=minR[irow];
            saveMax[nstackR]=maxR[irow];
            nstackR++;
          } else if (markR[irow]==-2) {
            continue;
          }
          /* could check immediately if violation */
          if (movement>0.0) {
            /* up */
            if (value>0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] += value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] += value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          } else {
            /* down */
            if (value<0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] -= value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] -= value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          }
        }
        while (istackC<nstackC&&nstackC<maxStack) {
          int jway;
          jcol=stackC[istackC];
          jway=markC[jcol];
          // If not first and fixed then skip
          if (jway==3&&istackC) {
            //istackC++;
            //continue;
            //printf("fixed %d on stack\n",jcol);
          }
          // Do cliques
          if (oneFixStart_&&oneFixStart_[jcol]>=0) {
            int start;
            int end;
            if (colLower[jcol]>saveL[istackC]) {
              // going up
              start = oneFixStart_[jcol];
              end = zeroFixStart_[jcol];
            } else {
              assert (colUpper[jcol]<saveU[istackC]);
              // going down
              start = zeroFixStart_[jcol];
              end = endFixStart_[jcol];
            }
            for (int i=start;i<end;i++) {
              int iClique = whichClique_[i];
              for (int k=cliqueStart_[iClique];k<cliqueStart_[iClique+1];k++) {
                int kcol = cliqueEntry_[k].sequence;
                if (jcol==kcol)
                  continue;
                int kway = cliqueEntry_[k].oneFixes;
                if (kcol!=jcol) {
                  if (!markC[kcol]) {
                    // not on list yet
                    if (nstackC<2*maxStack) {
                      markC[kcol] = 3; // say fixed
                      fixThis++;
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      if (!kway) {
                        // going up
                        colLower[kcol]=1.0;
                        /* update immediately */
                        for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                          krow = row[jj];
                          value = columnElements[jj];
                          if (markR[krow]==-1) {
                            stackR[nstackR]=krow;
                            markR[krow]=nstackR;
                            saveMin[nstackR]=minR[krow];
                            saveMax[nstackR]=maxR[krow];
                            nstackR++;
                          } else if (markR[krow]==-2) {
                            continue;
                          }
                          /* could check immediately if violation */
                          /* up */
                          if (value>0.0) {
                            /* up does not change - down does */
                            if (minR[krow]>-1.0e10)
                              minR[krow] += value;
                            if (minR[krow]>rowUpper[krow]+1.0e-5) {
                              colUpper[kcol]=-1.0e50; /* force infeasible */
                              break;
                            }
                          } else {
                            /* down does not change - up does */
                            if (maxR[krow]<1.0e10)
                              maxR[krow] += value;
                            if (maxR[krow]<rowLower[krow]-1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          }
                        }
                      } else {
                        // going down
                        colUpper[kcol]=0.0;
                        /* update immediately */
                        for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                          krow = row[jj];
                          value = columnElements[jj];
                          if (markR[krow]==-1) {
                            stackR[nstackR]=krow;
                            markR[krow]=nstackR;
                            saveMin[nstackR]=minR[krow];
                            saveMax[nstackR]=maxR[krow];
                            nstackR++;
                          } else if (markR[krow]==-2) {
                            continue;
                          }
                          /* could check immediately if violation */
                          /* down */
                          if (value<0.0) {
                            /* up does not change - down does */
                            if (minR[krow]>-1.0e10)
                              minR[krow] -= value;
                            if (minR[krow]>rowUpper[krow]+1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          } else {
                            /* down does not change - up does */
                            if (maxR[krow]<1.0e10)
                              maxR[krow] -= value;
                            if (maxR[krow]<rowLower[krow]-1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          }
                        }
                      }
                    }
                  } else if (markC[kcol]==1) {
                    // marked as going to 0
                    assert (!colUpper[kcol]);
                    if (!kway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  } else if (markC[kcol]==2) {
                    // marked as going to 1
                    assert (colLower[kcol]);
                    if (kway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  } else {
                    // marked as fixed
                    assert (markC[kcol]==3);
                    int jkway;
                    if (colLower[kcol])
                      jkway=1;
                    else
                      jkway=0;
                    if (kway==jkway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  }
                }
              }
              if (notFeasible)
                break;
            }
            if (notFeasible)
              istackC=nstackC+1;
          }
          for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
            // break if found not feasible
            if (notFeasible)
              break;
            irow = row[k];
            /*value = columnElements[k];*/
            if (markR[irow]!=-2) {
              /* see if anything forced */
              for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
                double moveUp=0.0;
                double moveDown=0.0;
                double newUpper=-1.0,newLower=1.0;
                kcol=column[kk];
                bool onList = (markC[kcol]!=0);
                if (markC[kcol]!=3) {
                  value2=rowElements[kk];
                  int markIt=markC[kcol];
                  if (value2 < 0.0) {
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colUpper[kcol]+
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  } else {
                    /* positive element */
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colUpper[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  }
                  if (nstackC<2*maxStack) 
                    markC[kcol] = markIt;
                  if (moveUp&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      onList=true;
                    }
                    if (intVar[kcol])
                      newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4));
                    colLower[kcol]=newLower;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
                      markC[kcol]=3; // say fixed
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* up */
                      if (value>0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveUp;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveUp;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (moveDown&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      onList=true;
                    }
                    if (intVar[kcol])
                      newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4));
                    colUpper[kcol]=newUpper;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6)
                      markC[kcol]=3; // say fixed
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* down */
                      if (value<0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveDown;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveDown;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
                    notFeasible=1;;
                    k=columnStart[jcol]+columnLength[jcol];
                    istackC=nstackC+1;
#ifdef PRINT_DEBUG
                    printf("** not feasible this way\n");
#endif
                    break;
                  }
                }
              }
            }
          }
          istackC++;
        }
        if (!notFeasible) {
          feasible |= feas[iway];
        } else if (iway==1&&feasible==0) {
          /* not feasible at all */
          ninfeas=1;
          j=nCols-1;
          iLook=nLook;
          ipass=maxPass;
          break;
        }
        if (iway==2||(iway==1&&feasible==2)) {
          /* keep */
          iway=3;
          nfixed++;
          if (mode_) {
            OsiColCut cc;
            int nTot=0,nFix=0,nInt=0;
            bool ifCut=false;
            for (istackC=0;istackC<nstackC;istackC++) {
              int icol=stackC[istackC];
              if (intVar[icol]) {
                if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
                  element[nFix]=colUpper[icol];
                  index[nFix++]=icol;
                  nInt++;
                  if (colsol[icol]>colUpper[icol]+primalTolerance_) {
                    ifCut=true;
                    anyColumnCuts=true;
                  }
                }
              }
            }
            if (nFix) {
              nTot=nFix;
              cc.setUbs(nFix,index,element);
              nFix=0;
            }
            for (istackC=0;istackC<nstackC;istackC++) {
              int icol=stackC[istackC];
              if (intVar[icol]) {
                if (colLower[icol]>currentColLower[icol]+1.0e-4) {
                  element[nFix]=colLower[icol];
                  index[nFix++]=icol;
                  nInt++;
                  if (colsol[icol]<colLower[icol]-primalTolerance_) {
                    ifCut=true;
                    anyColumnCuts=true;
                  }
                }
              }
            }
            if (nFix) {
              nTot+=nFix;
              cc.setLbs(nFix,index,element);
            }
            // could tighten continuous as well
            if (nInt) {
              if (ifCut) {
                cc.setEffectiveness(100.0);
              } else {
                cc.setEffectiveness(1.0e-5);
              }
#ifdef CGL_DEBUG
              checkBounds(debugger,cc);
#endif
              cs.insert(cc);
            }
          }
          for (istackC=0;istackC<nstackC;istackC++) {
            int icol=stackC[istackC];
            if (colUpper[icol]-colLower[icol]>primalTolerance_) {
              markC[icol]=0;
            } else {
              markC[icol]=3;
            }
          }
          for (istackR=0;istackR<nstackR;istackR++) {
            int irow=stackR[istackR];
            markR[irow]=-1;
          }
        } else {
          /* is it worth seeing if can increase coefficients
             or maybe better see if it is a cut */
          if (iway==0) {
            nstackC0=CoinMin(nstackC,maxStack);
            if (notFeasible) {
              nstackC0=0;
            } else {
              for (istackC=0;istackC<nstackC0;istackC++) {
                int icol=stackC[istackC];
                stackC0[istackC]=icol;
                lo0[istackC]=colLower[icol];
                up0[istackC]=colUpper[icol];
              }
            }
            /* restore all */
            assert (iway==0);
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              markC[icol]=0;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0) {
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    sum += rowElements[kk]*colsol[column[kk]];
                  }
                  if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_||info.strengthenRow) {
                    // can be a cut
                    // subtract gap from upper and integer coefficient
                    // saveU and saveL spare
                    int * index = (int *)saveL;
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]-gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=-gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(-DBL_MAX);
                    rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0));
                    // effectiveness is how far j moves
                    rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      // If strengthenRow point to row
                      //if(info.strengthenRow)
                      //printf("a point to row %d\n",irow);
                      rowCut.addCutIfNotDuplicate(rc,irow);
                    }
                  }
                }
                gap = minR[irow]-rowLower[irow];
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  if (!sum) {
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      sum += rowElements[kk]*colsol[column[kk]];
                    }
                  }
                  if (sum+gap*colsol[j]<minR[irow]+primalTolerance_||info.strengthenRow) {
                    // can be a cut
                    // add gap to lower and integer coefficient
                    // saveU and saveL spare
                    int * index = (int *)saveL;
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]+gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0));
                    rc.setUb(DBL_MAX);
                    // effectiveness is how far j moves
                    rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                      //if(info.strengthenRow)
                      //printf("b point to row %d\n",irow);
                      rowCut.addCutIfNotDuplicate(rc,irow);
                    }
                  }
                }
              }
              minR[irow]=saveMin[istackR];
              maxR[irow]=saveMax[istackR];
              markR[irow]=-1;
            }
          } else {
            if (iway==1&&feasible==3) {
              iway=3;
              /* point back to stack */
              for (istackC=nstackC-1;istackC>=0;istackC--) {
                int icol=stackC[istackC];
                markC[icol]=istackC+1000;
              }
              if (mode_) {
                OsiColCut cc;
                int nTot=0,nFix=0,nInt=0;
                bool ifCut=false;
                for (istackC=0;istackC<nstackC0;istackC++) {
                  int icol=stackC0[istackC];
                  int istackC1=markC[icol]-1000;
                  if (istackC1>=0) {
                    if (CoinMin(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
                      saveL[istackC1]=CoinMin(lo0[istackC],colLower[icol]);
                      if (intVar[icol]) {
                        element[nFix]=saveL[istackC1];
                        index[nFix++]=icol;
                        nInt++;
                        if (colsol[icol]<saveL[istackC1]-primalTolerance_)
                          ifCut=true;
                      }
                      nfixed++;
                    }
                  }
                }
                if (nFix) {
                  nTot=nFix;
                  cc.setLbs(nFix,index,element);
                  nFix=0;
                }
                for (istackC=0;istackC<nstackC0;istackC++) {
                  int icol=stackC0[istackC];
                  int istackC1=markC[icol]-1000;
                  if (istackC1>=0) {
                    if (CoinMax(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
                      saveU[istackC1]=CoinMax(up0[istackC],colUpper[icol]);
                      if (intVar[icol]) {
                        element[nFix]=saveU[istackC1];
                        index[nFix++]=icol;
                        nInt++;
                        if (colsol[icol]>saveU[istackC1]+primalTolerance_)
                          ifCut=true;
                      }
                      nfixed++;
                    }
                  }
                }
                if (nFix) {
                  nTot+=nFix;
                  cc.setUbs(nFix,index,element);
                }
                // could tighten continuous as well
                if (nInt) {
                  if (ifCut) {
                    cc.setEffectiveness(100.0);
                  } else {
                    cc.setEffectiveness(1.0e-5);
                  }
#ifdef CGL_DEBUG
                  checkBounds(debugger,cc);
#endif
                  cs.insert(cc);
                }
              }
            }
            /* restore all */
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              if (oldU>oldL+1.0e-4)
                markC[icol]=0;
              else
                markC[icol]=3;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0) {
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    sum += rowElements[kk]*colsol[column[kk]];
                  }
                  if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_||info.strengthenRow) {
                    // can be a cut
                    // add gap to integer coefficient
                    // saveU and saveL spare
                    int * index = (int *)saveL;
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-DBL_MAX);
		      rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info.strengthenRow)
                        //printf("c point to row %d\n",irow);
			rowCut.addCutIfNotDuplicate(rc,irow);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum-gap*colsol[j]<rowLower[irow]+primalTolerance_||info.strengthenRow) {
		      // can be a cut
		      // subtract gap from integer coefficient
		      // saveU and saveL spare
		      int * index = (int *)saveL;
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]-gap*(colUpper[j]-1));
		      rc.setUb(DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info.strengthenRow)
                        //printf("d point to row %d\n",irow);
			rowCut.addCutIfNotDuplicate(rc,irow);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    }
	  }
        }
    }
  }
  delete [] stackC0;
  delete [] lo0;
  delete [] up0;
  delete [] tempL;
  delete [] tempU;
  delete [] markC;
  delete [] stackC;
  delete [] stackR;
  delete [] saveL;
  delete [] saveU;
  delete [] saveMin;
  delete [] saveMax;
  delete [] index;
  delete [] element;
  delete [] colsol;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info.strengthenRow);
  }
  if (!columnCopy_)
    delete columnCopy; // if just created
  delete [] array;
  abort();
  return (ninfeas);
}
// Create a copy of matrix which is to be used
// this is to speed up process and to give global cuts
// Can give an array with 1 set to select, 0 to ignore
// column bounds are tightened
// If array given then values of 1 will be set to 0 if redundant
int CglProbing::snapshot ( const OsiSolverInterface & si,
		  char * possible,bool withObjective)
{
  deleteSnapshot();
  // Get basic problem information
  
  numberColumns_=si.getNumCols(); 
  numberRows_=si.getNumRows();
  colLower_ = new double[numberColumns_];
  colUpper_ = new double[numberColumns_];
  memcpy(colLower_,si.getColLower(),numberColumns_*sizeof(double));
  memcpy(colUpper_,si.getColUpper(),numberColumns_*sizeof(double));
  rowLower_= new double [numberRows_+1];
  rowUpper_= new double [numberRows_+1];
  memcpy(rowLower_,si.getRowLower(),numberRows_*sizeof(double));
  memcpy(rowUpper_,si.getRowUpper(),numberRows_*sizeof(double));

  int i;
  if (possible) {
    for (i=0;i<numberRows_;i++) {
      if (!possible[i]) {
	rowLower_[i]=-DBL_MAX;
	rowUpper_[i]=DBL_MAX;
      }
    }
  }
  

  char * intVar = new char[numberColumns_];
  numberIntegers_=0;
  number01Integers_=0;
  for (i=0;i<numberColumns_;i++) {
    if (si.isInteger(i)) {
      intVar[i]=2;  
      numberIntegers_++;
      if (si.isBinary(i)) {
        intVar[i]=1;  
        number01Integers_++;
      }
    } else {
      intVar[i]=0;
    }
  }
    
  rowCopy_ = new CoinPackedMatrix(*si.getMatrixByRow());

  const int * column = rowCopy_->getIndices();
  const CoinBigIndex * rowStart = rowCopy_->getVectorStarts();
  const int * rowLength = rowCopy_->getVectorLengths(); 
  const double * rowElements = rowCopy_->getElements();

  int returnCode=0;
  int ninfeas= 
    tighten(colLower_, colUpper_, column, rowElements,
			 rowStart, rowLength, rowLower_, rowUpper_,
			 numberRows_, numberColumns_, intVar, 5, primalTolerance_);
  if (ninfeas) {
    // let someone else find out
    returnCode = 1;
  }
/*
  QUESTION: If ninfeas > 1 (one or more variables infeasible), shouldn't we
	    bail out here?
*/

  // do integer stuff for mode 0
  cutVector_ = new disaggregation [number01Integers_];
  memset(cutVector_,0,number01Integers_*sizeof(disaggregation));
  number01Integers_=0;
  for (i=0;i<numberColumns_;i++) {
    if (si.isBinary(i))
      cutVector_[number01Integers_++].sequence=i;
  }
  delete [] intVar;

  // now delete rows
  if (possible) {
    for (i=0;i<numberRows_;i++) {
      if (rowLower_[i]<-1.0e30&&rowUpper_[i]>1.0e30) 
	possible[i]=0;
    }
  }
  int * index = new int[numberRows_];
  int nDrop=0,nKeep=0;
  for (i=0;i<numberRows_;i++) {
    if (rowLower_[i]<-1.0e30&&rowUpper_[i]>1.0e30) {
      index[nDrop++]=i;
    } else {
      rowLower_[nKeep]=rowLower_[i];
      rowUpper_[nKeep++]=rowUpper_[i];
    }
  }
  numberRows_=nKeep;
  if (nDrop)
    rowCopy_->deleteRows(nDrop,index);
  delete [] index;
  if (withObjective) {
    // add in objective 
    int * columns = new int[numberColumns_];
    double * elements = new double[numberColumns_];
    int n=0;
    const double * objective = si.getObjCoefficients();
    bool maximize = (si.getObjSense()==-1);
    for (i=0;i<numberColumns_;i++) {
      if (objective[i]) {
        elements[n]= (maximize) ? -objective[i] : objective[i];
        columns[n++]=i;
      }
    }
    rowCopy_->appendRow(n,columns,elements);
    delete [] columns;
    delete [] elements;
    numberRows_++;
  }
  // create column copy
  columnCopy_ = new CoinPackedMatrix(* rowCopy_);
  columnCopy_->reverseOrdering();
  // make sure big enough - in case too many rows dropped
  columnCopy_->setDimensions(numberRows_,numberColumns_);
  rowCopy_->setDimensions(numberRows_,numberColumns_);
  return returnCode;
}
// Delete snapshot
void CglProbing::deleteSnapshot()
{
  delete [] rowLower_;
  delete [] rowUpper_;
  delete [] colLower_;
  delete [] colUpper_;
  delete rowCopy_;
  delete columnCopy_;
  rowCopy_=NULL;
  columnCopy_=NULL;
  rowLower_=NULL;
  rowUpper_=NULL;
  colLower_=NULL;
  colUpper_=NULL;
  int i;
  for (i=0;i<number01Integers_;i++) {
    delete [] cutVector_[i].index;
  }
  delete [] cutVector_;
  numberIntegers_=0;
  number01Integers_=0;
  cutVector_=NULL;
}
// Mode stuff
void CglProbing::setMode(int mode)
{
  if (mode>=0&&mode<3) {
    // take off bottom bit
    mode_ &= ~15;
    mode_ |= mode;
  }
}
int CglProbing::getMode() const
{
  return mode_&15;
}
// Set maximum number of passes per node
void CglProbing::setMaxPass(int value)
{
  if (value>0)
    maxPass_=value;
}
// Get maximum number of passes per node
int CglProbing::getMaxPass() const
{
  return maxPass_;
}
// Set log level
void CglProbing::setLogLevel(int value)
{
  if (value>=0)
    logLevel_=value;
}
// Get log level
int CglProbing::getLogLevel() const
{
  return logLevel_;
}
// Set maximum number of unsatisfied variables to look at
void CglProbing::setMaxProbe(int value)
{
  if (value>=0)
    maxProbe_=value;
}
// Get maximum number of unsatisfied variables to look at
int CglProbing::getMaxProbe() const
{
  return maxProbe_;
}
// Set maximum number of variables to look at in one probe
void CglProbing::setMaxLook(int value)
{
  if (value>0)
    maxStack_=value;
}
// Get maximum number of variables to look at in one probe
int CglProbing::getMaxLook() const
{
  return maxStack_;
}
// Set maximum number of elements in row for scan
void CglProbing::setMaxElements(int value)
{
  if (value>0)
    maxElements_=value;
}
// Get maximum number of elements in row for scan
int CglProbing::getMaxElements() const
{
  return maxElements_;
}
// Set maximum number of passes per node (root node)
void CglProbing::setMaxPassRoot(int value)
{
  if (value>0)
    maxPassRoot_=value;
}
// Get maximum number of passes per node (root node)
int CglProbing::getMaxPassRoot() const
{
  return maxPassRoot_;
}
// Set maximum number of unsatisfied variables to look at (root node)
void CglProbing::setMaxProbeRoot(int value)
{
  if (value>0)
    maxProbeRoot_=value;
}
// Get maximum number of unsatisfied variables to look at (root node)
int CglProbing::getMaxProbeRoot() const
{
  return maxProbeRoot_;
}
// Set maximum number of variables to look at in one probe (root node)
void CglProbing::setMaxLookRoot(int value)
{
  if (value>0)
    maxStackRoot_=value;
}
// Get maximum number of variables to look at in one probe (root node)
int CglProbing::getMaxLookRoot() const
{
  return maxStackRoot_;
}
// Set maximum number of elements in row for scan (root node)
void CglProbing::setMaxElementsRoot(int value)
{
  if (value>0)
    maxElementsRoot_=value;
}
// Get maximum number of elements in row for scan (root node)
int CglProbing::getMaxElementsRoot() const
{
  return maxElementsRoot_;
}
// Set whether to use objective
void CglProbing::setUsingObjective(bool yesNo)
{
  usingObjective_=yesNo;
}
// Get whether objective is being used
int CglProbing::getUsingObjective() const
{
  return usingObjective_;
}
// Decide whether to do row cuts
void CglProbing::setRowCuts(int type)
{
  if (type>-5&&type<5)
    rowCuts_=type;
}
// Returns row cuts generation type
int CglProbing::rowCuts() const
{
  return rowCuts_;
}
// Returns tight lower
const double * CglProbing::tightLower() const
{
  return colLower_;
}
// Returns tight upper
const double * CglProbing::tightUpper() const
{
  return colUpper_;
}
// Returns relaxed Row lower
const double * CglProbing::relaxedRowLower() const
{
  return rowLower_;
}
// Returns relaxed Row upper
const double * CglProbing::relaxedRowUpper() const
{
  return rowUpper_;
}


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing ()
:
CglCutGenerator(),
primalTolerance_(1.0e-07),
mode_(1),
rowCuts_(1),
maxPass_(3),
logLevel_(0),
maxProbe_(100),
maxStack_(50),
maxElements_(100000),
maxPassRoot_(3),
maxProbeRoot_(100),
maxStackRoot_(50),
maxElementsRoot_(1000),
usingObjective_(false)
{

  numberRows_=0;
  numberColumns_=0;
  rowCopy_=NULL;
  columnCopy_=NULL;
  rowLower_=NULL;
  rowUpper_=NULL;
  colLower_=NULL;
  colUpper_=NULL;
  numberIntegers_=0;
  number01Integers_=0;
  cutVector_=NULL;
  numberCliques_=0;
  cliqueType_=NULL;
  cliqueStart_=NULL;
  cliqueEntry_=NULL;
  oneFixStart_=NULL;
  zeroFixStart_=NULL;
  endFixStart_=NULL;
  whichClique_=NULL;
  cliqueRow_=NULL;
  cliqueRowStart_=NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CglProbing::CglProbing (  const CglProbing & rhs)
                                                              :
  CglCutGenerator(rhs),
  primalTolerance_(rhs.primalTolerance_),
  mode_(rhs.mode_),
  rowCuts_(rhs.rowCuts_),
  maxPass_(rhs.maxPass_),
  logLevel_(rhs.logLevel_),
  maxProbe_(rhs.maxProbe_),
  maxStack_(rhs.maxStack_),
  maxElements_(rhs.maxElements_),
  maxPassRoot_(rhs.maxPassRoot_),
  maxProbeRoot_(rhs.maxProbeRoot_),
  maxStackRoot_(rhs.maxStackRoot_),
  maxElementsRoot_(rhs.maxElementsRoot_),
  usingObjective_(rhs.usingObjective_)
{  
  numberRows_=rhs.numberRows_;
  numberColumns_=rhs.numberColumns_;
  numberCliques_=rhs.numberCliques_;
  if (rhs.rowCopy_) {
    rowCopy_= new CoinPackedMatrix(*(rhs.rowCopy_));
    columnCopy_= new CoinPackedMatrix(*(rhs.columnCopy_));
    rowLower_=new double[numberRows_];
    memcpy(rowLower_,rhs.rowLower_,numberRows_*sizeof(double));
    rowUpper_=new double[numberRows_];
    memcpy(rowUpper_,rhs.rowUpper_,numberRows_*sizeof(double));
    colLower_=new double[numberColumns_];
    memcpy(colLower_,rhs.colLower_,numberColumns_*sizeof(double));
    colUpper_=new double[numberColumns_];
    memcpy(colUpper_,rhs.colUpper_,numberColumns_*sizeof(double));
    int i;
    numberIntegers_=rhs.numberIntegers_;
    number01Integers_=rhs.number01Integers_;
    cutVector_=new disaggregation [number01Integers_];
    memcpy(cutVector_,rhs.cutVector_,number01Integers_*sizeof(disaggregation));
    for (i=0;i<number01Integers_;i++) {
      if (cutVector_[i].index) {
	cutVector_[i].index = CoinCopyOfArray(rhs.cutVector_[i].index,cutVector_[i].length);
      }
    }
  } else {
    rowCopy_=NULL;
    columnCopy_=NULL;
    rowLower_=NULL;
    rowUpper_=NULL;
    colLower_=NULL;
    colUpper_=NULL;
    numberIntegers_=0;
    number01Integers_=0;
    cutVector_=NULL;
  }
  if (numberCliques_) {
    cliqueType_ = new cliqueType [numberCliques_];
    memcpy(cliqueType_,rhs.cliqueType_,numberCliques_*sizeof(cliqueType));
    cliqueStart_ = new int [numberCliques_+1];
    memcpy(cliqueStart_,rhs.cliqueStart_,(numberCliques_+1)*sizeof(int));
    int n = cliqueStart_[numberCliques_];
    cliqueEntry_ = new cliqueEntry [n];
    memcpy(cliqueEntry_,rhs.cliqueEntry_,n*sizeof(cliqueEntry));
    oneFixStart_ = new int [numberColumns_];
    memcpy(oneFixStart_,rhs.oneFixStart_,numberColumns_*sizeof(int));
    zeroFixStart_ = new int [numberColumns_];
    memcpy(zeroFixStart_,rhs.zeroFixStart_,numberColumns_*sizeof(int));
    endFixStart_ = new int [numberColumns_];
    memcpy(endFixStart_,rhs.endFixStart_,numberColumns_*sizeof(int));
    int n2=-1;
    for (int i=numberColumns_-1;i>=0;i--) {
      if (oneFixStart_[i]>=0) {
	n2=endFixStart_[i];
	break;
      }
    }
    assert (n==n2);
    whichClique_ = new int [n];
    memcpy(whichClique_,rhs.whichClique_,n*sizeof(int));
    if (rhs.cliqueRowStart_) {
      cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_,numberRows_+1);
      n=cliqueRowStart_[numberRows_];
      cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_,n);
    } else {
      cliqueRow_=NULL;
      cliqueRowStart_=NULL;
    }
  } else {
    cliqueType_=NULL;
    cliqueStart_=NULL;
    cliqueEntry_=NULL;
    oneFixStart_=NULL;
    zeroFixStart_=NULL;
    endFixStart_=NULL;
    cliqueRow_=NULL;
    cliqueRowStart_=NULL;
    whichClique_=NULL;
  }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglProbing::clone() const
{
  return new CglProbing(*this);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CglProbing::~CglProbing ()
{
  // free memory
  delete [] rowLower_;
  delete [] rowUpper_;
  delete [] colLower_;
  delete [] colUpper_;
  delete rowCopy_;
  delete columnCopy_;
  delete [] cliqueType_;
  delete [] cliqueStart_;
  delete [] cliqueEntry_;
  delete [] oneFixStart_;
  delete [] zeroFixStart_;
  delete [] endFixStart_;
  delete [] whichClique_;
  delete [] cliqueRow_;
  delete [] cliqueRowStart_;
  if (cutVector_) {
    for (int i=0;i<number01Integers_;i++) {
      delete [] cutVector_[i].index;
    }
    delete [] cutVector_;
  }
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CglProbing &
CglProbing::operator=(
                                         const CglProbing& rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    primalTolerance_=rhs.primalTolerance_;
    numberRows_=rhs.numberRows_;
    numberColumns_=rhs.numberColumns_;
    delete [] rowLower_;
    delete [] rowUpper_;
    delete [] colLower_;
    delete [] colUpper_;
    delete rowCopy_;
    delete columnCopy_;
    delete [] cliqueType_;
    delete [] cliqueStart_;
    delete [] cliqueEntry_;
    delete [] oneFixStart_;
    delete [] zeroFixStart_;
    delete [] endFixStart_;
    delete [] whichClique_;
    delete [] cliqueRow_;
    delete [] cliqueRowStart_;
    mode_=rhs.mode_;
    rowCuts_=rhs.rowCuts_;
    maxPass_=rhs.maxPass_;
    logLevel_=rhs.logLevel_;
    maxProbe_=rhs.maxProbe_;
    maxStack_=rhs.maxStack_;
    maxElements_ = rhs.maxElements_;
    maxPassRoot_ = rhs.maxPassRoot_;
    maxProbeRoot_ = rhs.maxProbeRoot_;
    maxStackRoot_ = rhs.maxStackRoot_;
    maxElementsRoot_ = rhs.maxElementsRoot_;
    usingObjective_=rhs.usingObjective_;
    numberCliques_=rhs.numberCliques_;
    if (rhs.rowCopy_) {
      rowCopy_= new CoinPackedMatrix(*(rhs.rowCopy_));
      columnCopy_= new CoinPackedMatrix(*(rhs.columnCopy_));
      rowLower_=new double[numberRows_];
      memcpy(rowLower_,rhs.rowLower_,numberRows_*sizeof(double));
      rowUpper_=new double[numberRows_];
      memcpy(rowUpper_,rhs.rowUpper_,numberRows_*sizeof(double));
      colLower_=new double[numberColumns_];
      memcpy(colLower_,rhs.colLower_,numberColumns_*sizeof(double));
      colUpper_=new double[numberColumns_];
      memcpy(colUpper_,rhs.colUpper_,numberColumns_*sizeof(double));
      int i;
      numberIntegers_=rhs.numberIntegers_;
      number01Integers_=rhs.number01Integers_;
      for (i=0;i<number01Integers_;i++) {
        delete [] cutVector_[i].index;
      }
      delete [] cutVector_;
      cutVector_=new disaggregation [number01Integers_];
      memcpy(cutVector_,rhs.cutVector_,number01Integers_*sizeof(disaggregation));
      for (i=0;i<number01Integers_;i++) {
        if (cutVector_[i].index) {
          cutVector_[i].index = CoinCopyOfArray(rhs.cutVector_[i].index,cutVector_[i].length);
        }
      }
    } else {
      rowCopy_=NULL;
      columnCopy_=NULL;
      rowLower_=NULL;
      rowUpper_=NULL;
      colLower_=NULL;
      colUpper_=NULL;
      numberIntegers_=0;
      number01Integers_=0;
      cutVector_=NULL;
    }
    if (numberCliques_) {
      cliqueType_ = new cliqueType [numberCliques_];
      memcpy(cliqueType_,rhs.cliqueType_,numberCliques_*sizeof(cliqueType));
      cliqueStart_ = new int [numberCliques_+1];
      memcpy(cliqueStart_,rhs.cliqueStart_,(numberCliques_+1)*sizeof(int));
      int n = cliqueStart_[numberCliques_];
      cliqueEntry_ = new cliqueEntry [n];
      memcpy(cliqueEntry_,rhs.cliqueEntry_,n*sizeof(cliqueEntry));
      oneFixStart_ = new int [numberColumns_];
      memcpy(oneFixStart_,rhs.oneFixStart_,numberColumns_*sizeof(int));
      zeroFixStart_ = new int [numberColumns_];
      memcpy(zeroFixStart_,rhs.zeroFixStart_,numberColumns_*sizeof(int));
      endFixStart_ = new int [numberColumns_];
      memcpy(endFixStart_,rhs.endFixStart_,numberColumns_*sizeof(int));
      int n2=-1;
      for (int i=numberColumns_-1;i>=0;i--) {
	if (oneFixStart_[i]>=0) {
	  n2=endFixStart_[i];
	  break;
	}
      }
      assert (n==n2);
      whichClique_ = new int [n];
      memcpy(whichClique_,rhs.whichClique_,n*sizeof(int));
      if (rhs.cliqueRowStart_) {
        cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_,numberRows_+1);
        n=cliqueRowStart_[numberRows_];
        cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_,n);
      } else {
        cliqueRow_=NULL;
        cliqueRowStart_=NULL;
      }
    } else {
      cliqueType_=NULL;
      cliqueStart_=NULL;
      cliqueEntry_=NULL;
      oneFixStart_=NULL;
      zeroFixStart_=NULL;
      endFixStart_=NULL;
      whichClique_=NULL;
      cliqueRow_=NULL;
      cliqueRowStart_=NULL;
    }
  }
  return *this;
}

/// This can be used to refresh any inforamtion
void 
CglProbing::refreshSolver(OsiSolverInterface * solver)
{
  if (rowCopy_) {
    // snapshot existed - redo
    snapshot(*solver,NULL);
  }
}
/* Creates cliques for use by probing.
   Can also try and extend cliques as a result of probing (root node).
   Returns number of cliques found.
*/
int 
CglProbing::createCliques( OsiSolverInterface & si, 
			  int minimumSize, int maximumSize, bool extendCliques)
{
  // get rid of what is there
  deleteCliques();
  CoinPackedMatrix matrixByRow(*si.getMatrixByRow());
  int numberRows = si.getNumRows();
  if (!rowCopy_)
    numberRows_=numberRows;
  numberColumns_ = si.getNumCols();

  numberCliques_=0;
  int numberEntries=0;
  int numberIntegers=0;
  int * lookup = new int[numberColumns_];
  int i;
  for (i=0;i<numberColumns_;i++) {
    if (si.isBinary(i))
      lookup[i]=numberIntegers++;
    else
      lookup[i]=-1;
  }

  int * which = new int[numberColumns_];
  int * whichRow = new int[numberRows];
  // Statistics
  int totalP1=0,totalM1=0;
  int numberBig=0,totalBig=0;
  int numberFixed=0;

  // Row copy
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column lengths for slacks
  const int * columnLength = si.getMatrixByCol()->getVectorLengths();

  const double * lower = si.getColLower();
  const double * upper = si.getColUpper();
  const double * rowLower = si.getRowLower();
  const double * rowUpper = si.getRowUpper();
  int iRow;
  for (iRow=0;iRow<numberRows;iRow++) {
    int numberP1=0, numberM1=0;
    int j;
    double upperValue=rowUpper[iRow];
    double lowerValue=rowLower[iRow];
    bool good=true;
    int slack = -1;
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn = column[j];
      int iInteger=lookup[iColumn];
      if (upper[iColumn]-lower[iColumn]<1.0e-8) {
	// fixed
	upperValue -= lower[iColumn]*elementByRow[j];
	lowerValue -= lower[iColumn]*elementByRow[j];
	continue;
      } else if (upper[iColumn]!=1.0||lower[iColumn]!=0.0) {
	good = false;
	break;
      } else if (iInteger<0) {
	good = false;
	break;
      } else {
	if (columnLength[iColumn]==1)
	  slack = iInteger;
      }
      if (fabs(elementByRow[j])!=1.0) {
	good=false;
	break;
      } else if (elementByRow[j]>0.0) {
	which[numberP1++]=iColumn;
      } else {
	numberM1++;
	which[numberIntegers-numberM1]=iColumn;
      }
    }
    int iUpper = (int) floor(upperValue+1.0e-5);
    int iLower = (int) ceil(lowerValue-1.0e-5);
    int state=0;
    if (upperValue<1.0e6) {
      if (iUpper==1-numberM1)
	state=1;
      else if (iUpper==-numberM1)
	state=2;
      else if (iUpper<-numberM1)
	state=3;
    }
    if (!state&&lowerValue>-1.0e6) {
      if (-iLower==1-numberP1)
	state=-1;
      else if (-iLower==-numberP1)
	state=-2;
      else if (-iLower<-numberP1)
	state=-3;
    }
    if (good&&state) {
      if (abs(state)==3) {
	// infeasible
	numberCliques_ = -99999;
	break;
      } else if (abs(state)==2) {
	// we can fix all
	numberFixed += numberP1+numberM1;
	if (state>0) {
	  // fix all +1 at 0, -1 at 1
	  for (i=0;i<numberP1;i++)
	    si.setColUpper(which[i],0.0);
	  for (i=0;i<numberM1;i++)
	    si.setColLower(which[numberIntegers-i-1],
				 1.0);
	} else {
	  // fix all +1 at 1, -1 at 0
	  for (i=0;i<numberP1;i++)
	    si.setColLower(which[i],1.0);
	  for (i=0;i<numberM1;i++)
	    si.setColUpper(which[numberIntegers-i-1],
				 0.0);
	}
      } else {
	int length = numberP1+numberM1;
        totalP1 += numberP1;
        totalM1 += numberM1;
	if (length >= minimumSize&&length<maximumSize) {
	  whichRow[numberCliques_++]=iRow;
	  numberEntries += length;
	} else if (numberP1+numberM1 >= maximumSize) {
	  // too big
	  numberBig++;
	  totalBig += numberP1+numberM1;
	}
      }
    }
  }
  if (numberCliques_<0) {
    if (logLevel_)
      printf("*** Problem infeasible\n");
  } else {
    if (numberCliques_) {
      if (logLevel_)
        printf("%d cliques of average size %g found, %d P1, %d M1\n",
               numberCliques_,
               ((double)(totalP1+totalM1))/((double) numberCliques_),
               totalP1,totalM1);
    } else {
      if (logLevel_>1)
        printf("No cliques found\n");
    }
    if (numberBig) {
      if (logLevel_)
        printf("%d large cliques ( >= %d) found, total %d\n",
	     numberBig,maximumSize,totalBig);
    }
    if (numberFixed) {
      if (logLevel_)
        printf("%d variables fixed\n",numberFixed);
    }
  }
  if (numberCliques_>0) {
    cliqueType_ = new cliqueType [numberCliques_];
    cliqueStart_ = new int [numberCliques_+1];
    cliqueEntry_ = new cliqueEntry [numberEntries];
    oneFixStart_ = new int [numberColumns_];
    zeroFixStart_ = new int [numberColumns_];
    endFixStart_ = new int [numberColumns_];
    whichClique_ = new int [numberEntries];
    numberEntries=0;
    cliqueStart_[0]=0;
    for (i=0;i<numberColumns_;i++) {
      oneFixStart_[i]=-1;
      zeroFixStart_[i]=-1;
      endFixStart_[i]=-1;
    }
    int iClique;
    for (iClique=0;iClique<numberCliques_;iClique++) {
      int iRow=whichRow[iClique];
      int numberP1=0, numberM1=0;
      int j;
      double upperValue=rowUpper[iRow];
      double lowerValue=rowLower[iRow];
      for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	int iColumn = column[j];
	if (upper[iColumn]-lower[iColumn]<1.0e-8) {
	  // fixed
	  upperValue -= lower[iColumn]*elementByRow[j];
	  lowerValue -= lower[iColumn]*elementByRow[j];
	  continue;
	}
	if (elementByRow[j]>0.0) {
	  which[numberP1++]=iColumn;
	} else {
	  numberM1++;
	  which[numberIntegers-numberM1]=iColumn;
	}
      }
      int iUpper = (int) floor(upperValue+1.0e-5);
      int iLower = (int) ceil(lowerValue-1.0e-5);
      int state=0;
      if (upperValue<1.0e6) {
	if (iUpper==1-numberM1)
	  state=1;
      }
      if (!state&&lowerValue>-1.0e6) {
	state=-1;
      }
      assert (abs(state)==1);
      if (iLower==iUpper) {
	cliqueType_[iClique].equality=1;
      } else {
	cliqueType_[iClique].equality=0;
      }
      if (state>0) {
	for (i=0;i<numberP1;i++) {
	  // 1 is strong branch
	  int iColumn = which[i];
	  cliqueEntry_[numberEntries].sequence=iColumn;
	  cliqueEntry_[numberEntries].oneFixes=1;
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
	for (i=0;i<numberM1;i++) {
	  // 0 is strong branch
	  int iColumn = which[numberIntegers-i-1];
	  cliqueEntry_[numberEntries].sequence=iColumn;
	  cliqueEntry_[numberEntries].oneFixes=0;
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
      } else {
	for (i=0;i<numberP1;i++) {
	  // 0 is strong branch
	  int iColumn = which[i];
	  cliqueEntry_[numberEntries].sequence=iColumn;
	  cliqueEntry_[numberEntries].oneFixes=0;
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
	for (i=0;i<numberM1;i++) {
	  // 1 is strong branch
	  int iColumn = which[numberIntegers-i-1];
	  cliqueEntry_[numberEntries].sequence=iColumn;
	  cliqueEntry_[numberEntries].oneFixes=1;
	  numberEntries++;
	  // zero counts
	  oneFixStart_[iColumn]=0;
	  zeroFixStart_[iColumn]=0;
	}
      }
      cliqueStart_[iClique+1]=numberEntries;
    }
    // Now do column lists
    // First do counts
    for (iClique=0;iClique<numberCliques_;iClique++) {
      for (int j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
	int iColumn = cliqueEntry_[j].sequence;
	if (cliqueEntry_[j].oneFixes)
	  oneFixStart_[iColumn]++;
	else
	  zeroFixStart_[iColumn]++;
      }
    }
    // now get starts and use which and end as counters
    numberEntries=0;
    for (int iColumn=0;iColumn<numberColumns_;iColumn++) {
      if (oneFixStart_[iColumn]>=0) {
	int n1=oneFixStart_[iColumn];
	int n2=zeroFixStart_[iColumn];
	oneFixStart_[iColumn]=numberEntries;
	which[iColumn]=numberEntries;
	numberEntries += n1;
	zeroFixStart_[iColumn]=numberEntries;
	endFixStart_[iColumn]=numberEntries;
	numberEntries += n2;
      }
    }
    // now put in
    for (iClique=0;iClique<numberCliques_;iClique++) {
      for (int j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
	int iColumn = cliqueEntry_[j].sequence;
	if (cliqueEntry_[j].oneFixes) {
	  int put = which[iColumn];
	  which[iColumn]++;
	  whichClique_[put]=iClique;
	} else {
	  int put = endFixStart_[iColumn];
	  endFixStart_[iColumn]++;
	  whichClique_[put]=iClique;
	}
      }
    }
  }
  delete [] which;
  delete [] whichRow;
  delete [] lookup;
  return numberCliques_;
}
// Delete all clique information
void 
CglProbing::deleteCliques()
{
  delete [] cliqueType_;
  delete [] cliqueStart_;
  delete [] cliqueEntry_;
  delete [] oneFixStart_;
  delete [] zeroFixStart_;
  delete [] endFixStart_;
  delete [] whichClique_;
  delete [] cliqueRow_;
  delete [] cliqueRowStart_;
  cliqueType_=NULL;
  cliqueStart_=NULL;
  cliqueEntry_=NULL;
  oneFixStart_=NULL;
  zeroFixStart_=NULL;
  endFixStart_=NULL;
  whichClique_=NULL;
  cliqueRow_=NULL;
  cliqueRowStart_=NULL;
  numberCliques_=0;
}
/*
  Returns true if may generate Row cuts in tree (rather than root node).
  Used so know if matrix will change in tree.  Really
  meant so column cut generators can still be active
  without worrying code.
  Default is true
*/
bool 
CglProbing::mayGenerateRowCutsInTree() const
{
  return rowCuts_>0;
}
// Sets up clique information for each row
void 
CglProbing::setupRowCliqueInformation(const OsiSolverInterface & si) 
{
  if (!numberCliques_)
    return;
  CoinPackedMatrix * rowCopy;
  if (!rowCopy_) {
    // create from current
    numberRows_=si.getNumRows(); 
    numberColumns_=si.getNumCols(); 
    rowCopy = new CoinPackedMatrix(*si.getMatrixByRow());
  } else {
    rowCopy = rowCopy_;
    assert(numberRows_<=si.getNumRows()); 
    assert(numberColumns_==si.getNumCols()); 
  }
  assert(numberRows_&&numberColumns_);
  cliqueRowStart_ = new int [numberRows_+1];
  cliqueRowStart_[0]=0;
  // Temporary array while building list
  cliqueEntry ** array = new cliqueEntry * [numberRows_];
  // Which cliques in use
  int * which = new int[numberCliques_];
  int * count = new int[numberCliques_];
  int * back =new int[numberColumns_];
  CoinZeroN(count,numberCliques_);
  CoinFillN(back,numberColumns_,-1);
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * lower = si.getColLower();
  const double * upper = si.getColUpper();
  int iRow;
  for (iRow=0;iRow<numberRows_;iRow++) {
    int j;
    int numberFree=0;
    int numberUsed=0;
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn=column[j];
      if (upper[iColumn]>lower[iColumn]) {
        back[iColumn]=j-rowStart[iRow];
        numberFree++;
        for (int k=oneFixStart_[iColumn];k<endFixStart_[iColumn];k++) {
          int iClique = whichClique_[k];
          if (!count[iClique]) {
            which[numberUsed++]=iClique;
          }
          count[iClique]++;
        }
      }
    }
    // find largest cliques
    bool finished=false;
    int numberInThis=0;
    cliqueEntry * entries = NULL;
    array[iRow]=entries;
    while (!finished) {
      int largest=1;
      int whichClique=-1;
      for (int i=0;i<numberUsed;i++) {
        int iClique = which[i];
        if (count[iClique]>largest) {
          largest=count[iClique];
          whichClique=iClique;
        }
      }
      // Add in if >1 (but not if all as that means clique==row)
      if (whichClique>=0&&largest<numberFree) {
        if (!numberInThis) {
          int length=rowLength[iRow];
          entries = new cliqueEntry [length];
          array[iRow]=entries;
          for (int i=0;i<length;i++) {
            entries[i].oneFixes=0;
            entries[i].sequence=numberColumns_+1;
          }
        }
        // put in (and take out all counts)
        for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
          int iColumn=column[j];
          if (upper[iColumn]>lower[iColumn]) {
            bool found=false;
            int k;
            for ( k=oneFixStart_[iColumn];k<endFixStart_[iColumn];k++) {
              int iClique = whichClique_[k];
              if (iClique==whichClique) {
                found=true;
                break;
              }
            }
            if (found) {
              for ( k=oneFixStart_[iColumn];k<endFixStart_[iColumn];k++) {
                int iClique = whichClique_[k];
                count[iClique]--;
              }
              for (k=cliqueStart_[whichClique];k<cliqueStart_[whichClique+1];k++) {
                if (cliqueEntry_[k].sequence==(unsigned int) iColumn) {
                  int iback=back[iColumn];
                  entries[iback].sequence=numberInThis;
                  entries[iback].oneFixes=cliqueEntry_[k].oneFixes;
                  break;
                }
              }
            }
          }
        }
        numberInThis++;
      } else {
        finished=true;
      }
    }
    if (numberInThis) 
      cliqueRowStart_[iRow+1]=cliqueRowStart_[iRow]+rowLength[iRow];
    else
      cliqueRowStart_[iRow+1]=cliqueRowStart_[iRow];
    for (int i=0;i<numberUsed;i++) {
      int iClique = which[i];
      count[iClique]=0;
    }
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn=column[j];
      back[iColumn]=-1;
    }
  }
  delete [] which;
  delete [] count;
  delete [] back;
  // Now put info in one array
  cliqueRow_ = new cliqueEntry [cliqueRowStart_[numberRows_]];
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (array[iRow]) {
      int start = cliqueRowStart_[iRow];
      memcpy(cliqueRow_+start,array[iRow],rowLength[iRow]*sizeof(cliqueEntry));
      delete [] array[iRow];
    }
  }
  delete [] array;
  if (rowCopy!=rowCopy_)
    delete rowCopy;
}
// Create C++ lines to get to current state
std::string
CglProbing::generateCpp( FILE * fp) 
{
  CglProbing other;
  fprintf(fp,"0#include \"CglProbing.hpp\"\n");
  fprintf(fp,"3  CglProbing probing;\n");
  if (getMode()!=other.getMode())
    fprintf(fp,"3  probing.setMode(%d);\n",getMode());
  else
    fprintf(fp,"4  probing.setMode(%d);\n",getMode());
  if (getMaxPass()!=other.getMaxPass())
    fprintf(fp,"3  probing.setMaxPass(%d);\n",getMaxPass());
  else
    fprintf(fp,"4  probing.setMaxPass(%d);\n",getMaxPass());
  if (getLogLevel()!=other.getLogLevel())
    fprintf(fp,"3  probing.setLogLevel(%d);\n",getLogLevel());
  else
    fprintf(fp,"4  probing.setLogLevel(%d);\n",getLogLevel());
  if (getMaxProbe()!=other.getMaxProbe())
    fprintf(fp,"3  probing.setMaxProbe(%d);\n",getMaxProbe());
  else
    fprintf(fp,"4  probing.setMaxProbe(%d);\n",getMaxProbe());
  if (getMaxLook()!=other.getMaxLook())
    fprintf(fp,"3  probing.setMaxLook(%d);\n",getMaxLook());
  else
    fprintf(fp,"4  probing.setMaxLook(%d);\n",getMaxLook());
  if (getMaxElements()!=other.getMaxElements())
    fprintf(fp,"3  probing.setMaxElements(%d);\n",getMaxElements());
  else
    fprintf(fp,"4  probing.setMaxElements(%d);\n",getMaxElements());
  if (getMaxPassRoot()!=other.getMaxPassRoot())
    fprintf(fp,"3  probing.setMaxPassRoot(%d);\n",getMaxPassRoot());
  else
    fprintf(fp,"4  probing.setMaxPassRoot(%d);\n",getMaxPassRoot());
  if (getMaxProbeRoot()!=other.getMaxProbeRoot())
    fprintf(fp,"3  probing.setMaxProbeRoot(%d);\n",getMaxProbeRoot());
  else
    fprintf(fp,"4  probing.setMaxProbeRoot(%d);\n",getMaxProbeRoot());
  if (getMaxLookRoot()!=other.getMaxLookRoot())
    fprintf(fp,"3  probing.setMaxLookRoot(%d);\n",getMaxLookRoot());
  else
    fprintf(fp,"4  probing.setMaxLookRoot(%d);\n",getMaxLookRoot());
  if (getMaxElementsRoot()!=other.getMaxElementsRoot())
    fprintf(fp,"3  probing.setMaxElementsRoot(%d);\n",getMaxElementsRoot());
  else
    fprintf(fp,"4  probing.setMaxElementsRoot(%d);\n",getMaxElementsRoot());
  if (rowCuts()!=other.rowCuts())
    fprintf(fp,"3  probing.setRowCuts(%d);\n",rowCuts());
  else
    fprintf(fp,"4  probing.setRowCuts(%d);\n",rowCuts());
  if (getUsingObjective()!=other.getUsingObjective())
    fprintf(fp,"3  probing.setUsingObjective(%s);\n",getUsingObjective() ? "true" : "false");
  else
    fprintf(fp,"4  probing.setUsingObjective(%s);\n",getUsingObjective() ? "true" : "false");
  if (getAggressiveness()!=other.getAggressiveness())
    fprintf(fp,"3  probing.setAggressiveness(%d);\n",getAggressiveness());
  else
    fprintf(fp,"4  probing.setAggressiveness(%d);\n",getAggressiveness());
  return "probing";
}
