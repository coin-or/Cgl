
/*
  Copyright (C) 2010, International Business Machines Corporation and others.
  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

  $Id$
*/

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinPackedMatrix.hpp"
#include "CglProbing.hpp"
#include "CglProbingRowCut.hpp"

// To enable debugging, set CGL_DEBUG in CglProbingDebug.hpp
#include "CglProbingDebug.hpp"


/*
  Parameters useObj, useCutoff, and cutoff are unused as of 110208. Added to
  probe(), added here for symmetry.

  An assert down in mode 0 (snapshot) processing in gutsOfGenerateCuts
  says that any cut returned from probeCliques has two coefficients. So
  when I go through this, I should be limited to disaggregation cuts. No
  coefficient strengthening.  But a quick look through the cases at the
  end of the probe loop says that isn't true. I can clearly see the code
  block for coefficient strengthening.  -- lh, 110211 --

  My guess, after a few more days of browsing, is that the mode 0 code in
  gutsOfGenerateCuts is obsolete (unused since April 2008 r629) and implements
  discovery of cuts that really should be found by the current cut generation
  subroutines. So the goal in revising probeCliques is to figure out just
  what's affected by the presence of cliques and (in so far as possible)
  implement the control structure and functionality now in probe().
  -- lh, 110215 --
*/

int CglProbing::probeCliques (const OsiSolverInterface & si, 
			      OsiCuts & cs, 
			      double *colLower, double *colUpper, 
			      CoinPackedMatrix *rowCopy,
			      CoinPackedMatrix *columnCopy,
			      const int *realRows,
			      double *rowLower, double *rowUpper,
			      char *intVar,
			      double *minR, double *maxR, int *markR, 
			      CglTreeInfo * info,
			      bool useObj, bool useCutoff, double cutoff
			     ) const
{
# if CGL_DEBUG > 0
  std::cout << "Entering CglProbing::probeCliques." << std::endl ;
  int numberRowCutsBefore = cs.sizeRowCuts() ;
  int numberColCutsBefore = cs.sizeColCuts() ;

  const OsiRowCutDebugger *debugger = si.getRowCutDebugger() ;
# endif

  assert(false) ;

/*
  PREPARATION

  All the code down through `PASS LOOP: HEAD' is preparation. Do all of the
  setup that will not change over the nested loops that do the work of probing.
  Note that rowCopy may have considerably fewer rows than the original system.
*/

  int nRows = rowCopy->getNumRows() ;
  int nCols = rowCopy->getNumCols() ;
  int maxStack = info->inTree ? maxStack_ : maxStackRoot_ ;
/*
  Allocate working arrays.

  TODO: Note that there's no equivalent to the `one big block' allocation in
        probe().   -- lh, 110215 --
*/
  double * colsol = new double[nCols] ;
  double * djs = new double[nCols] ;
  const double * currentColLower = si.getColLower() ;
  const double * currentColUpper = si.getColUpper() ;
  double * tempL = new double [nCols] ;
  double * tempU = new double [nCols] ;
  int * markC = new int [nCols] ;
  int * stackC = new int [2*nCols] ;
  int * stackR = new int [nRows] ;
  double * saveL = new double [2*nCols] ;
  double * saveU = new double [2*nCols] ;
  double * saveMin = new double [nRows] ;
  double * saveMax = new double [nRows] ;
  double * element = new double[nCols] ;
  int * index = new int[nCols] ;

  // For trying to extend cliques
  int * cliqueStack=NULL ;
  int * cliqueCount=NULL ;
  int * to_01=NULL ;
  if (!mode_) {
    to_01 = new int[nCols] ;
    cliqueStack = new int[numberCliques_] ;
    cliqueCount = new int[numberCliques_] ;
    int i ;
    for (i=0;i<numberCliques_;i++) {
      cliqueCount[i]=cliqueStart_[i+1]-cliqueStart_[i] ;
    }
    for (i=0;i<nCols;i++) 
      to_01[i]=-1 ;
    for (i=0;i<number01Integers_;i++) {
      int j=cutVector_[i].sequence ;
      to_01[j]=i ;
    }
  }
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info->inTree ? nRows/3 : nRows ;
  CglProbingRowCut rowCut(nRowsFake, !info->inTree) ;
  const int * column = rowCopy->getIndices() ;
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts() ;
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements() ;
  const int * row = columnCopy->getIndices() ;
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts() ;
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements() ;
  double movement ;
  int i, j, k,kk,jj ;
  int kcol,krow ;
  bool anyColumnCuts=false ;
  double dbound, value, value2 ;
  int ninfeas=0 ;
  int rowCuts ;
  double disaggEffectiveness ;
  if (mode_) {
    /* clean up djs and solution */
    CoinMemcpyN(si.getReducedCost(),nCols,djs) ;
    CoinMemcpyN( si.getColSolution(),nCols,colsol) ;
    disaggEffectiveness=1.0e-3 ;
    rowCuts=rowCuts_ ;
  } else {
    // need to go from a neutral place
    memset(djs,0,nCols*sizeof(double)) ;
    CoinMemcpyN( si.getColSolution(),nCols,colsol) ;
    disaggEffectiveness=-1.0e10 ;
    if (rowCuts_!=4)
      rowCuts=1 ;
    else
      rowCuts=4 ;
  }
  for (i = 0; i < nCols; ++i) {
    /* was if (intVar[i]) */
    if (1) {
      if (colUpper[i]-colLower[i]>1.0e-8) {
	if (colsol[i]<colLower[i]+primalTolerance_) {
	  colsol[i]=colLower[i] ;
	  djs[i] = CoinMax(0.0,djs[i]) ;
	} else if (colsol[i]>colUpper[i]-primalTolerance_) {
	  colsol[i]=colUpper[i] ;
	  djs[i] = CoinMin(0.0,djs[i]) ;
	} else {
	  djs[i]=0.0 ;
	}
	/*if (fabs(djs[i])<1.0e-5) 
	  djs[i]=0.0;*/
      }
    }
  }

  int ipass=0,nfixed=-1 ;
/*
  Chop now that useCutoff and cutoff come in as parameters.
  double cutoff ;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit,cutoff) ;
  if (!cutoff_available||usingObjective_<0) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity() ;
  }
  cutoff *= si.getObjSense() ;
  if (fabs(cutoff)>1.0e30)
    assert (cutoff>1.0e30) ;
*/
  double current = si.getObjValue() ;
  // make irrelevant if mode is 0
  if (!mode_)
    cutoff=DBL_MAX ;
  /* for both way coding */
  int nstackC0=-1 ;
  int * stackC0 = new int[maxStack] ;
  double * lo0 = new double[maxStack] ;
  double * up0 = new double[maxStack] ;
  int nstackR,nstackC ;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3 ;
    } else {
      markC[i]=0 ;
    }
  }
  double tolerance = 1.0e1*primalTolerance_ ;
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_ ;
  // If we are going to replace coefficient then we don't need to be effective
  double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3 ;
/*
  PASS LOOP: HEAD
*/
  while (ipass<maxPass&&nfixed) {
    int iLook ;
    ipass++ ;
    nfixed=0 ;
    for (iLook=0;iLook<numberThisTime_;iLook++) {
      double solval ;
      double down ;
      double up ;
      int awayFromBound=1 ;
      j=lookedAt_[iLook] ;
      solval=colsol[j] ;
      down = floor(solval+tolerance) ;
      up = ceil(solval-tolerance) ;
      if(colUpper[j]-colLower[j]<1.0e-8) markC[j]=3 ;
      if (markC[j]||!intVar[j]) continue ;
      double saveSolval = solval ;
      if (solval>=colUpper[j]-tolerance||solval<=colLower[j]+tolerance||up==down) {
	awayFromBound=0 ;
	if (solval<=colLower[j]+2.0*tolerance) {
	  solval = colLower[j]+1.0e-1 ;
	  down=colLower[j] ;
	  up=down+1.0 ;
	} else if (solval>=colUpper[j]-2.0*tolerance) {
	  solval = colUpper[j]-1.0e-1 ;
	  up=colUpper[j] ;
	  down=up-1 ;
	} else {
          // odd
          up=down+1.0 ;
          solval = down+1.0e-1 ;
        }
      }
      assert (up<=colUpper[j]) ;
      assert (down>=colLower[j]) ;
      assert (up>down) ;
      if ((solval-down>1.0e-6&&up-solval>1.0e-6)||mode_!=1) {
	int istackC,iway, istackR ;
	int way[]={1,2,1} ;
	int feas[]={1,2,4} ;
	int feasible=0 ;
	int notFeasible ;
/*
  PROBE LOOP: HEAD
*/
	for (iway=0;iway<3;iway ++) {
	  int fixThis=0 ;
	  double objVal=current ;
	  int goingToTrueBound=0 ;
	  stackC[0]=j ;
	  markC[j]=way[iway] ;
          double solMovement ;
	  if (way[iway]==1) {
	    movement=down-colUpper[j] ;
            solMovement = down-colsol[j] ;
	    assert(movement<-0.99999) ;
	    if (fabs(down-colLower[j])<1.0e-7) {
	      goingToTrueBound=2 ;
	      down=colLower[j] ;
	    }
	  } else {
	    movement=up-colLower[j] ;
            solMovement = up-colsol[j] ;
	    assert(movement>0.99999) ;
	    if (fabs(up-colUpper[j])<1.0e-7) {
	      goingToTrueBound=2 ;
	      up=colUpper[j] ;
	    }
	  }
	  if (goingToTrueBound&&(colUpper[j]-colLower[j]>1.5||colLower[j]))
	    goingToTrueBound=1 ;
	  // switch off disaggregation if not wanted
	  if ((rowCuts&1)==0)
	    goingToTrueBound=0 ;
#ifdef PRINT_DEBUG
	  if (fabs(movement)>1.01) {
	    printf("big %d %g %g %g\n",j,colLower[j],solval,colUpper[j]) ;
	  }
#endif
	  if (solMovement*djs[j]>0.0)
	    objVal += solMovement*djs[j] ;
	  nstackC=1 ;
	  nstackR=0 ;
	  saveL[0]=colLower[j] ;
	  saveU[0]=colUpper[j] ;
          assert (saveU[0]>saveL[0]) ;
	  notFeasible=0 ;
	  if (movement<0.0) {
	    colUpper[j] += movement ;
	    colUpper[j] = floor(colUpper[j]+0.5) ;
#ifdef PRINT_DEBUG
	    printf("** Trying %d down to 0\n",j) ;
#endif
	  } else {
	    colLower[j] += movement ;
	    colLower[j] = floor(colLower[j]+0.5) ;
#ifdef PRINT_DEBUG
	    printf("** Trying %d up to 1\n",j) ;
#endif
	  }
	  if (fabs(colUpper[j]-colLower[j])<1.0e-6) 
	    markC[j]=3; // say fixed
	  istackC=0 ;
	  /* update immediately */
/*
  Propagate effect of probe. This loop will be replaced with updateRowBounds.
*/
	  for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
	    int irow = row[k] ;
	    value = columnElements[k] ;
	    assert (markR[irow]!=-2) ;
	    if (markR[irow]==-1) {
	      stackR[nstackR]=irow ;
	      markR[irow]=nstackR ;
	      saveMin[nstackR]=minR[irow] ;
	      saveMax[nstackR]=maxR[irow] ;
	      nstackR++ ;
#if 0
	    } else if (markR[irow]==-2) {
	      continue ;
#endif
	    }
	    /* could check immediately if violation */
	    if (movement>0.0) {
	      /* up */
	      if (value>0.0) {
		/* up does not change - down does */
                if (minR[irow]>-1.0e10)
                  minR[irow] += value ;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1 ;
		  istackC=1 ;
		  break ;
		}
	      } else {
		/* down does not change - up does */
                if (maxR[irow]<1.0e10)
                  maxR[irow] += value ;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1 ;
		  istackC=1 ;
		  break ;
		}
	      }
	    } else {
	      /* down */
	      if (value<0.0) {
		/* up does not change - down does */
                if (minR[irow]>-1.0e10)
                  minR[irow] -= value ;
		if (minR[irow]>rowUpper[irow]+1.0e-5) {
		  notFeasible=1 ;
		  istackC=1 ;
		  break ;
		}
	      } else {
		/* down does not change - up does */
                if (maxR[irow]<1.0e10)
                  maxR[irow] -= value ;
		if (maxR[irow]<rowLower[irow]-1.0e-5) {
		  notFeasible=1 ;
		  istackC=1 ;
		  break ;
		}
	      }
	    }
	  }
/*
  PROBE LOOP: BEGIN PROPAGATION
*/
	  while (istackC<nstackC&&nstackC<maxStack) {
	    int jway ;
	    int jcol =stackC[istackC] ;
	    jway=markC[jcol] ;
	    // If not first and fixed then skip
	    if (jway==3&&istackC) {
	      //istackC++ ;
	      //continue ;
              //printf("fixed %d on stack\n",jcol) ;
	    }
	    // Do cliques
	    if (oneFixStart_&&oneFixStart_[jcol]>=0) {
	      int start ;
	      int end ;
	      if (colLower[jcol]>saveL[istackC]) {
		// going up
		start = oneFixStart_[jcol] ;
		end = zeroFixStart_[jcol] ;
	      } else {
		assert (colUpper[jcol]<saveU[istackC]) ;
		// going down
		start = zeroFixStart_[jcol] ;
		end = endFixStart_[jcol] ;
	      }
	      for (int i=start;i<end;i++) {
		int iClique = whichClique_[i] ;
		for (int k=cliqueStart_[iClique];k<cliqueStart_[iClique+1];k++) {
		  int kcol = sequenceInCliqueEntry(cliqueEntry_[k]) ;
                  if (jcol==kcol)
                    continue ;
		  int kway = oneFixesInCliqueEntry(cliqueEntry_[k]) ;
                  if (kcol!=jcol) {
                    if (!markC[kcol]) {
                      // not on list yet
                      if (nstackC<2*maxStack) {
                        markC[kcol] = 3; // say fixed
                        fixThis++ ;
                        stackC[nstackC]=kcol ;
                        saveL[nstackC]=colLower[kcol] ;
                        saveU[nstackC]=colUpper[kcol] ;
                        assert (saveU[nstackC]>saveL[nstackC]) ;
                        nstackC++ ;
                        if (!kway) {
                          // going up
                          double solMovement=1.0-colsol[kcol] ;
                          if (solMovement>0.0001) {
                            assert (djs[kcol]>=0.0) ;
                            objVal += djs[kcol]*solMovement ;
                          }
                          colLower[kcol]=1.0 ;
                          /* update immediately */
                          for (int jj =columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                            krow = row[jj] ;
                            value = columnElements[jj] ;
			    assert (markR[krow]!=-2) ;
                            if (markR[krow]==-1) {
                              stackR[nstackR]=krow ;
                              markR[krow]=nstackR ;
                              saveMin[nstackR]=minR[krow] ;
                              saveMax[nstackR]=maxR[krow] ;
                              nstackR++ ;
#if 0
                            } else if (markR[krow]==-2) {
                              continue ;
#endif
                            }
                            /* could check immediately if violation */
                            /* up */
                            if (value>0.0) {
                              /* up does not change - down does */
                              if (minR[krow]>-1.0e10)
                                minR[krow] += value ;
                              if (minR[krow]>rowUpper[krow]+1.0e-5) {
                                colUpper[kcol]=-1.0e50; /* force infeasible */
                                break ;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow]<1.0e10)
                                maxR[krow] += value ;
                              if (maxR[krow]<rowLower[krow]-1.0e-5) {
                                notFeasible=1 ;
                                break ;
                              }
                            }
                          }
                        } else {
                          // going down
                          double solMovement=0.0-colsol[kcol] ;
                          if (solMovement<-0.0001) {
                            assert (djs[kcol]<=0.0) ;
                            objVal += djs[kcol]*solMovement ;
                          }
                          colUpper[kcol]=0.0 ;
                          /* update immediately */
                          for (int jj =columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                            krow = row[jj] ;
                            value = columnElements[jj] ;
			    assert (markR[krow]!=-2) ;
                            if (markR[krow]==-1) {
                              stackR[nstackR]=krow ;
                              markR[krow]=nstackR ;
                              saveMin[nstackR]=minR[krow] ;
                              saveMax[nstackR]=maxR[krow] ;
                              nstackR++ ;
#if 0
                            } else if (markR[krow]==-2) {
                              continue ;
#endif
                            }
                            /* could check immediately if violation */
                            /* down */
                            if (value<0.0) {
                              /* up does not change - down does */
                              if (minR[krow]>-1.0e10)
                                minR[krow] -= value ;
                              if (minR[krow]>rowUpper[krow]+1.0e-5) {
                                notFeasible=1 ;
                                break ;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow]<1.0e10)
                                maxR[krow] -= value ;
                              if (maxR[krow]<rowLower[krow]-1.0e-5) {
                                notFeasible=1 ;
                                break ;
                              }
                            }
                          }
                        }
                      }
                    } else if (markC[kcol]==1) {
                      // marked as going to 0
                      assert (!colUpper[kcol]) ;
                      if (!kway) {
                        // contradiction
                        notFeasible=1 ;
                        break ;
                      }
                    } else if (markC[kcol]==2) {
                      // marked as going to 1
                      assert (colLower[kcol]) ;
                      if (kway) {
                        // contradiction
                        notFeasible=1 ;
                        break ;
                      }
                    } else {
                      // marked as fixed
                      assert (markC[kcol]==3) ;
                      int jkway ;
                      if (colLower[kcol])
                        jkway=1 ;
                      else
                        jkway=0 ;
                      if (kway==jkway) {
                        // contradiction
                        notFeasible=1 ;
                        break ;
                      }
                    }
                  }
		}
		if (notFeasible)
		  break ;
	      }
	      if (notFeasible)
		istackC=nstackC+1 ;
	    }
	    for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
	      // break if found not feasible
	      if (notFeasible)
		break ;
	      int irow = row[k] ;
	      /*value = columnElements[k];*/
	      assert (markR[irow]!=-2) ;
#if 0
	      if (markR[irow]!=-2) {
#endif
		/* see if anything forced */
		for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
		  double moveUp=0.0 ;
		  double moveDown=0.0 ;
		  double newUpper=-1.0,newLower=1.0 ;
		  kcol=column[kk] ;
		  bool onList = (markC[kcol]!=0) ;
		  if (markC[kcol]!=3) {
		    value2=rowElements[kk] ;
                    int markIt=markC[kcol] ;
		    if (value2 < 0.0) {
		      if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colUpper[kcol]+
			  (rowUpper[irow]-minR[irow])/value2 ;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
                            markIt |= 2 ;
			    newLower = ceil(dbound-primalTolerance_) ;
			  } else {
			    newLower=dbound ;
			    if (newLower+primalTolerance_>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol] ;
                              markIt |= 2 ;
                              //markIt=3 ;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3 ;
                            }
			  }
			  moveUp = newLower-colLower[kcol] ;
			}
		      }
		      if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colLower[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2 ;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 1 ;
			    newUpper = floor(dbound+primalTolerance_) ;
			  } else {
			    newUpper=dbound ;
			    if (newUpper-primalTolerance_<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol] ;
                              markIt |= 1 ;
                              //markIt=3 ;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3 ;
                            }
			  }
			  moveDown = newUpper-colUpper[kcol] ;
			}
		      }
		    } else {
		      /* positive element */
		      if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
			  rowLower[irow]>-1.0e10) {
			dbound = colUpper[kcol] + 
			  (rowLower[irow]-maxR[irow])/value2 ;
			if (dbound > colLower[kcol] + primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 2 ;
			    newLower = ceil(dbound-primalTolerance_) ;
			  } else {
			    newLower=dbound ;
			    if (newLower+primalTolerance_>colUpper[kcol]&&
				newLower-primalTolerance_<=colUpper[kcol]) {
			      newLower=colUpper[kcol] ;
                              markIt |= 2 ;
                              //markIt=3 ;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3 ;
			    }
			  }
			  moveUp = newLower-colLower[kcol] ;
			}
		      }
		      if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
			  rowUpper[irow]<1.0e10) {
			dbound = colLower[kcol] + 
			  (rowUpper[irow]-minR[irow])/value2 ;
			if (dbound < colUpper[kcol] - primalTolerance_) {
			  if (intVar[kcol]) {
			    markIt |= 1 ;
			    newUpper = floor(dbound+primalTolerance_) ;
			  } else {
			    newUpper=dbound ;
			    if (newUpper-primalTolerance_<colLower[kcol]&&
				newUpper+primalTolerance_>=colLower[kcol]) {
			      newUpper=colLower[kcol] ;
                              markIt |= 1 ;
                              //markIt=3 ;
			    } else {
                              // avoid problems - fix later ?
                              markIt=3 ;
			    }
			  }
			  moveDown = newUpper-colUpper[kcol] ;
			}
		      }
		    }
		    if (nstackC<2*maxStack) {
                      markC[kcol] = markIt ;
		    }
		    if (moveUp&&nstackC<2*maxStack) {
		      fixThis++ ;
#ifdef PRINT_DEBUG
		      printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]) ;
		      value=0.0 ;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj] ;
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]) ;
			} else {
			  value += rowElements[jj]*colLower[ii] ;
			}
		      }
		      printf(" - fixed %g\n",value) ;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj] ;
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]) ;
			}
		      }
		      printf("\n") ;
#endif
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
                        assert (saveU[nstackC]>saveL[nstackC]) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newLower>colsol[kcol]) {
			if (djs[kcol]<0.0) {
			  /* should be infeasible */
			  assert (newLower>colUpper[kcol]+primalTolerance_) ;
			} else {
			  objVal += moveUp*djs[kcol] ;
			}
		      }
		      if (intVar[kcol])
			newLower = CoinMax(colLower[kcol],ceil(newLower-1.0e-4)) ;
		      colLower[kcol]=newLower ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]=3; // say fixed
		      }
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj] ;
			value = columnElements[jj] ;
			assert (markR[krow]!=-2) ;
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow ;
			  markR[krow]=nstackR ;
			  saveMin[nstackR]=minR[krow] ;
			  saveMax[nstackR]=maxR[krow] ;
			  nstackR++ ;
#if 0
			} else if (markR[krow]==-2) {
			  continue ;
#endif
			}
			/* could check immediately if violation */
			/* up */
			if (value>0.0) {
			  /* up does not change - down does */
                          if (minR[krow]>-1.0e10)
                            minR[krow] += value*moveUp ;
			  if (minR[krow]>rowUpper[krow]+1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break ;
			  }
			} else {
			  /* down does not change - up does */
                          if (maxR[krow]<1.0e10)
                            maxR[krow] += value*moveUp ;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break ;
			  }
			}
		      }
		    }
		    if (moveDown&&nstackC<2*maxStack) {
		      fixThis++ ;
#ifdef PRINT_DEBUG
		      printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]) ;
		      value=0.0 ;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj] ;
			if (colUpper[ii]-colLower[ii]>primalTolerance_) {
			  printf("(%d, %g) ",ii,rowElements[jj]) ;
			} else {
			  value += rowElements[jj]*colLower[ii] ;
			}
		      }
		      printf(" - fixed %g\n",value) ;
		      for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
			int ii=column[jj] ;
			if (colUpper[ii]-colLower[ii]<primalTolerance_) {
			  printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]) ;
			}
		      }
		      printf("\n") ;
#endif
		      if (!onList) {
			stackC[nstackC]=kcol ;
			saveL[nstackC]=colLower[kcol] ;
			saveU[nstackC]=colUpper[kcol] ;
                        assert (saveU[nstackC]>saveL[nstackC]) ;
			nstackC++ ;
			onList=true ;
		      }
		      if (newUpper<colsol[kcol]) {
			if (djs[kcol]>0.0) {
			  /* should be infeasible */
			  assert (colLower[kcol]>newUpper+primalTolerance_) ;
			} else {
			  objVal += moveDown*djs[kcol] ;
			}
		      }
		      if (intVar[kcol])
			newUpper = CoinMin(colUpper[kcol],floor(newUpper+1.0e-4)) ;
		      colUpper[kcol]=newUpper ;
		      if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
			markC[kcol]=3; // say fixed
		      }
		      /* update immediately */
		      for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
			krow = row[jj] ;
			value = columnElements[jj] ;
			assert (markR[krow]!=-2) ;
			if (markR[krow]==-1) {
			  stackR[nstackR]=krow ;
			  markR[krow]=nstackR ;
			  saveMin[nstackR]=minR[krow] ;
			  saveMax[nstackR]=maxR[krow] ;
			  nstackR++ ;
#if 0
			} else if (markR[krow]==-2) {
#endif
			  continue ;
			}
			/* could check immediately if violation */
			/* down */
			if (value<0.0) {
			  /* up does not change - down does */
                          if (minR[krow]>-1.0e10)
                            minR[krow] += value*moveDown ;
			  if (minR[krow]>rowUpper[krow]+1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break ;
			  }
			} else {
			  /* down does not change - up does */
                          if (maxR[krow]<1.0e10)
                            maxR[krow] += value*moveDown ;
			  if (maxR[krow]<rowLower[krow]-1.0e-5) {
			    colUpper[kcol]=-1.0e50; /* force infeasible */
			    break ;
			  }
			}
		      }
		    }
		    if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
		      notFeasible=1; ;
		      k=columnStart[jcol]+columnLength[jcol] ;
		      istackC=nstackC+1 ;
#ifdef PRINT_DEBUG
		      printf("** not feasible this way\n") ;
#endif
		      break ;
		    }
		  }
		}
#if 0
	      }
#endif
	    }
	    istackC++ ;
	  }
/*
  PROBE LOOP: END PROPAGATION

  PROBE LOOP: INFEASIBILITY
*/
	  if (!notFeasible) {
	    if (objVal<=cutoff) {
	      feasible |= feas[iway] ;
	    } else {
#ifdef PRINT_DEBUG
	      printf("not feasible on dj\n") ;
#endif
	      notFeasible=1 ;
	      if (iway==1&&feasible==0) {
		/* not feasible at all */
		ninfeas=1 ;
		j=nCols-1 ;
		break ;
	      }
	    }
	  } else if (iway==1&&feasible==0) {
	    /* not feasible at all */
	    ninfeas=1 ;
	    j=nCols-1 ;
            iLook=numberThisTime_ ;
	    ipass=maxPass ;
	    break ;
	  }
	  if (notFeasible)
	    goingToTrueBound=0 ;
/*
  PROBE LOOP: MONOTONE (BREAK AND KEEP)
*/
	  if (iway==2||(iway==1&&feasible==2)) {
	    /* keep */
	    iway=3 ;
	    nfixed++ ;
            if (mode_) {
	      OsiColCut cc ;
	      int nTot=0,nFix=0,nInt=0 ;
	      bool ifCut=false ;
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC] ;
		if (intVar[icol]) {
		  if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
		    element[nFix]=colUpper[icol] ;
		    index[nFix++]=icol ;
		    nInt++ ;
		    if (colsol[icol]>colUpper[icol]+primalTolerance_) {
		      ifCut=true ;
		      anyColumnCuts=true ;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot=nFix ;
		cc.setUbs(nFix,index,element) ;
		nFix=0 ;
	      }
	      for (istackC=0;istackC<nstackC;istackC++) {
		int icol=stackC[istackC] ;
		if (intVar[icol]) {
		  if (colLower[icol]>currentColLower[icol]+1.0e-4) {
		    element[nFix]=colLower[icol] ;
		    index[nFix++]=icol ;
		    nInt++ ;
		    if (colsol[icol]<colLower[icol]-primalTolerance_) {
		      ifCut=true ;
		      anyColumnCuts=true ;
		    }
		  }
		}
	      }
	      if (nFix) {
		nTot+=nFix ;
		cc.setLbs(nFix,index,element) ;
	      }
	      // could tighten continuous as well
	      if (nInt) {
		if (ifCut) {
		  cc.setEffectiveness(100.0) ;
		} else {
		  cc.setEffectiveness(1.0e-5) ;
		}
#if CGL_DEBUG > 0
		CglProbingDebug::checkBounds(si,cc) ;
#endif
		cs.insert(cc) ;
	      }
	    }
	    for (istackC=0;istackC<nstackC;istackC++) {
	      int icol=stackC[istackC] ;
	      if (colUpper[icol]-colLower[icol]>primalTolerance_) {
		markC[icol]=0 ;
	      } else {
		markC[icol]=3 ;
	      }
	    }
	    for (istackR=0;istackR<nstackR;istackR++) {
	      int irow=stackR[istackR] ;
	      markR[irow]=-1 ;
	    }
	  }
/*
  PROBE LOOP: ITERATE

  First case is down feasible, we want to go on to the up probe.
*/
	  else {
	    /* is it worth seeing if can increase coefficients
	       or maybe better see if it is a cut */
	    if (iway==0) {
	      nstackC0=CoinMin(nstackC,maxStack) ;
	      double solMove = saveSolval-down ;
	      double boundChange ;
	      if (notFeasible) {
		nstackC0=0 ;
	      } else {
		for (istackC=0;istackC<nstackC0;istackC++) {
		  int icol=stackC[istackC] ;
		  stackC0[istackC]=icol ;
		  lo0[istackC]=colLower[icol] ;
		  up0[istackC]=colUpper[icol] ;
		}
	      }
	      /* restore all */
              int nCliquesAffected=0 ;
              assert (iway==0) ;
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC] ;
		double oldU=saveU[istackC] ;
		double oldL=saveL[istackC] ;
		if(goingToTrueBound==2&&istackC) {
                  // Work for extending cliques
                  if (!mode_&&numberCliques_) {
                    int i_01 = to_01[icol] ;
                    if (i_01>=0) {
                      int start ;
                      int end ;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol] ;
                        end = endFixStart_[icol] ;
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol] ;
                        end = zeroFixStart_[icol] ;
                      }
                      //if (end>start)
                      //printf("j %d, other %d is in %d cliques\n",
                      //     j,i_01,end-start) ;
                      for (int i=start;i<end;i++) {
                        int iClique = whichClique_[i] ;
                        int size = cliqueStart_[iClique+1]-cliqueStart_[iClique] ;
                        if (cliqueCount[iClique]==size) {
                          // first time
                          cliqueStack[nCliquesAffected++]=iClique ;
                        }
                        // decrement counts
                        cliqueCount[iClique]-- ;
                      }
                    }
                  }
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (jval-down)
		  boundChange = oldU-colUpper[icol] ;
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc ;
		    rc.setLb(-DBL_MAX) ;
		    rc.setUb(colUpper[icol]-down*boundChange) ;
		    index[0]=icol ;
		    element[0]=1.0 ;
		    index[1]=j ;
		    element[1]= - boundChange ;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange ;
		    if (mode_) 
		      assert(newSol>solMove) ;
		    rc.setEffectiveness(newSol-solMove) ;
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false) ;
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc) ;
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (jval-down)
		  boundChange = oldL-colLower[icol] ;
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc ;
		    rc.setLb(colLower[icol]-down*boundChange) ;
		    rc.setUb(DBL_MAX) ;
		    index[0]=icol ;
		    element[0]=1.0 ;
		    index[1]=j ;
		    element[1]=- boundChange ;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange ;
		    if (mode_)
		      assert(newSol>solMove) ;
		    rc.setEffectiveness(newSol-solMove) ;
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false) ;
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc) ;
#if 0
		      printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,oldL,oldU,colLower[icol],colsol[icol], j, colsol[j]) ;
		      printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
			     icol,boundChange,j,rc.ub()) ;
#endif
		    }
		  }
		}
		colUpper[icol]=oldU ;
		colLower[icol]=oldL ;
		markC[icol]=0 ;
	      }
              if (nCliquesAffected) {
                for (int i=0;i<nCliquesAffected;i++) {
                  int iClique = cliqueStack[i] ;
                  int size = cliqueCount[iClique] ;
                  // restore
                  cliqueCount[iClique]= cliqueStart_[iClique+1]-cliqueStart_[iClique] ;
                  if (!size) {
                    if (logLevel_>1)
                      printf("** could extend clique by adding j!\n") ;
                  }
                }
              }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR] ;
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts ;
		  double gap = rowUpper[irow]-maxR[irow] ;
		  double sum=0.0 ;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]] ;
		    }
		    if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_||(info->strengthenRow&&rowLower[irow]<-1.0e20)) {
		      // can be a cut
		      // subtract gap from upper and integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL) ;
		      double * element = saveU ;
		      int n=0 ;
                      bool coefficientExists=false ;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk] ;
			  element[n++]=rowElements[kk] ;
			} else {
			  double value=rowElements[kk]-gap ;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk] ;
			    element[n++]=value ;
			  }
			  coefficientExists=true ;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j ;
			element[n++]=-gap ;
		      }
		      OsiRowCut rc ;
		      rc.setLb(-DBL_MAX) ;
		      rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0)) ;
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap) ;
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        // If strengthenRow point to row
                        //if(info->strengthenRow)
                        //printf("a point to row %d\n",irow) ;
#ifdef STRENGTHEN_PRINT
		      if (rowLower[irow]<-1.0e20) {
			printf("5Cut %g <= ",rc.lb()) ;
			int k ;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k] ;
			  printf("%g*",element[k]) ;
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn) ;
			  else
			    printf("x%d ",iColumn) ;
			}
			printf("<= %g\n",rc.ub()) ;
			printf("Row %g <= ",rowLower[irow]) ;
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k] ;
			  printf("%g*",rowElements[k]) ;
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn) ;
			  else
			    printf("x%d ",iColumn) ;
			}
			printf("<= %g\n",rowUpper[irow]) ;
		      }
#endif
		      int realRow = (rowLower[irow]<-1.0e20) ? irow : -1 ;
		      if (realRows&&realRow>0)
			realRow=realRows[realRow] ;
		      rc.setWhichRow(realRow) ;
		      rowCut.addCutIfNotDuplicate(rc) ;
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow] ;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]] ;
		      }
		    }
		    if (sum+gap*colsol[j]<minR[irow]+primalTolerance_||(info->strengthenRow&&rowUpper[irow]>1.0e20)) {
		      // can be a cut
		      // add gap to lower and integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL) ;
		      double * element = saveU ;
		      int n=0 ;
                      bool coefficientExists=false ;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk] ;
			  element[n++]=rowElements[kk] ;
			} else {
			  double value=rowElements[kk]+gap ;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk] ;
			    element[n++]=value ;
			  }
			  coefficientExists=true ;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j ;
			element[n++]=gap ;
		      }
		      OsiRowCut rc ;
		      rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0)) ;
		      rc.setUb(DBL_MAX) ;
		      // effectiveness is how far j moves
		      rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap) ;
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("b point to row %d\n",irow) ;
#ifdef STRENGTHEN_PRINT
		      if (rowUpper[irow]>1.0e20) {
			printf("6Cut %g <= ",rc.lb()) ;
			int k ;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k] ;
			  printf("%g*",element[k]) ;
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn) ;
			  else
			    printf("x%d ",iColumn) ;
			}
			printf("<= %g\n",rc.ub()) ;
			printf("Row %g <= ",rowLower[irow]) ;
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k] ;
			  printf("%g*",rowElements[k]) ;
			  if (si.isInteger(iColumn))
			    printf("i%d ",iColumn) ;
			  else
			    printf("x%d ",iColumn) ;
			}
			printf("<= %g\n",rowUpper[irow]) ;
		      }
#endif
		      int realRow = (rowUpper[irow]>1.0e20) ? irow : -1 ;
		      if (realRows&&realRow>0)
			realRow=realRows[realRow] ;
		      rc.setWhichRow(realRow) ;
		      rowCut.addCutIfNotDuplicate(rc) ;
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR] ;
		maxR[irow]=saveMax[istackR] ;
		markR[irow]=-1 ;
	      }
	    }
/*
  PROBE LOOP: DOWN AND UP FEASIBLE

  Down and up feasible. We'll do cuts and be done.
*/
	    else {
	      if (iway==1&&feasible==3) {
		iway=3 ;
		/* point back to stack */
		for (istackC=nstackC-1;istackC>=0;istackC--) {
		  int icol=stackC[istackC] ;
		  markC[icol]=istackC+1000 ;
		}
		if (mode_) {
		  OsiColCut cc ;
		  int nTot=0,nFix=0,nInt=0 ;
		  bool ifCut=false ;
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC] ;
		    int istackC1=markC[icol]-1000 ;
		    if (istackC1>=0) {
		      if (CoinMin(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
			saveL[istackC1]=CoinMin(lo0[istackC],colLower[icol]) ;
			if (intVar[icol]) {
			  element[nFix]=saveL[istackC1] ;
			  index[nFix++]=icol ;
			  nInt++ ;
			  if (colsol[icol]<saveL[istackC1]-primalTolerance_)
			    ifCut=true ;
			}
			nfixed++ ;
		      }
		    }
		  }
		  if (nFix) {
		    nTot=nFix ;
		    cc.setLbs(nFix,index,element) ;
		    nFix=0 ;
		  }
		  for (istackC=0;istackC<nstackC0;istackC++) {
		    int icol=stackC0[istackC] ;
		    int istackC1=markC[icol]-1000 ;
		    if (istackC1>=0) {
		      if (CoinMax(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
			saveU[istackC1]=CoinMax(up0[istackC],colUpper[icol]) ;
			if (intVar[icol]) {
			  element[nFix]=saveU[istackC1] ;
			  index[nFix++]=icol ;
			  nInt++ ;
			  if (colsol[icol]>saveU[istackC1]+primalTolerance_)
			    ifCut=true ;
			}
			nfixed++ ;
		      }
		    }
		  }
		  if (nFix) {
		    nTot+=nFix ;
		    cc.setUbs(nFix,index,element) ;
		  }
		  // could tighten continuous as well
		  if (nInt) {
		    if (ifCut) {
		      cc.setEffectiveness(100.0) ;
		    } else {
		      cc.setEffectiveness(1.0e-5) ;
		    }
#if CGL_DEBUG > 0
		    CglProbingDebug::checkBounds(si,cc) ;
#endif
		    cs.insert(cc) ;
		  }
		}
	      } else {
		goingToTrueBound=0 ;
	      }
	      double solMove = up-saveSolval ;
	      double boundChange ;
	      /* restore all */
              int nCliquesAffected=0 ;
	      for (istackC=nstackC-1;istackC>=0;istackC--) {
		int icol=stackC[istackC] ;
		double oldU=saveU[istackC] ;
		double oldL=saveL[istackC] ;
		if(goingToTrueBound==2&&istackC) {
                  // Work for extending cliques
                  if (!mode_&&numberCliques_&&iway==3) {
                    int i_01 = to_01[icol] ;
                    if (i_01>=0) {
                      int start ;
                      int end ;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol] ;
                        end = endFixStart_[icol] ;
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol] ;
                        end = zeroFixStart_[icol] ;
                      }
                      //if (end>start)
                      //printf("up j %d, other %d is in %d cliques\n",
                      //     j,i_01,end-start) ;
                      for (int i=start;i<end;i++) {
                        int iClique = whichClique_[i] ;
                        int size = cliqueStart_[iClique+1]-cliqueStart_[iClique] ;
                        if (cliqueCount[iClique]==size) {
                          // first time
                          cliqueStack[nCliquesAffected++]=iClique ;
                        }
                        // decrement counts
                        cliqueCount[iClique]-- ;
                      }
                    }
                  }
		  // upper disaggregation cut would be
		  // xval < upper + (old_upper-upper) (up-jval)
		  boundChange = oldU-colUpper[icol] ;
		  if (boundChange>0.0&&oldU<1.0e10&&
		      (!mode_||colsol[icol]>colUpper[icol]
		      + boundChange*solMove+primalTolerance_)) {
		    // create cut
		    OsiRowCut rc ;
		    rc.setLb(-DBL_MAX) ;
		    rc.setUb(colUpper[icol]+up*boundChange) ;
		    index[0]=icol ;
		    element[0]=1.0 ;
		    index[1]=j ;
		    element[1]= + boundChange ;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colUpper[icol])/
		      boundChange ;
		    if (mode_)
		      assert(newSol>solMove) ;
		    rc.setEffectiveness(newSol-solMove) ;
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false) ;
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc) ;
		    }
		  }
		  // lower disaggregation cut would be
		  // xval > lower + (old_lower-lower) (up-jval)
		  boundChange = oldL-colLower[icol] ;
		  if (boundChange<0.0&&oldL>-1.0e10&&
		      (!mode_||colsol[icol]<colLower[icol]
		      + boundChange*solMove-primalTolerance_)) {
		    // create cut
		    OsiRowCut rc ;
		    rc.setLb(colLower[icol]+up*boundChange) ;
		    rc.setUb(DBL_MAX) ;
		    index[0]=icol ;
		    element[0]=1.0 ;
		    index[1]=j ;
		    element[1]= + boundChange ;
		    // effectiveness is how far j moves
		    double newSol = (colsol[icol]-colLower[icol])/
		      boundChange ;
		    if (mode_)
		      assert(newSol>solMove) ;
		    rc.setEffectiveness(newSol-solMove) ;
		    if (rc.effectiveness()>disaggEffectiveness) {
		      rc.setRow(2,index,element,false) ;
#if CGL_DEBUG > 0
		      if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
		      rowCut.addCutIfNotDuplicate(rc) ;
		    }
		  }
		}
		colUpper[icol]=oldU ;
		colLower[icol]=oldL ;
                if (oldU>oldL+1.0e-4)
                  markC[icol]=0 ;
                else
                  markC[icol]=3 ;
	      }
              if (nCliquesAffected) {
                for (int i=0;i<nCliquesAffected;i++) {
                  int iClique = cliqueStack[i] ;
                  int size = cliqueCount[iClique] ;
                  // restore
                  cliqueCount[iClique]= cliqueStart_[iClique+1]-cliqueStart_[iClique] ;
                  if (!size) {
                    if (logLevel_>1)
                      printf("** could extend clique by adding j!\n") ;
                  }
                }
              }
	      for (istackR=0;istackR<nstackR;istackR++) {
		int irow=stackR[istackR] ;
		// switch off strengthening if not wanted
		if ((rowCuts&2)!=0&&goingToTrueBound) {
		  bool ifCut=anyColumnCuts ;
		  double gap = rowUpper[irow]-maxR[irow] ;
		  double sum=0.0 ;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			 kk++) {
		      sum += rowElements[kk]*colsol[column[kk]] ;
		    }
		    if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_||(info->strengthenRow&&rowLower[irow]<-1.0e20)) {
		      // can be a cut
		      // add gap to integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL) ;
		      double * element = saveU ;
		      int n=0 ;
                      bool coefficientExists=false ;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk] ;
			  element[n++]=rowElements[kk] ;
			} else {
			  double value=rowElements[kk]+gap ;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk] ;
			    element[n++]=value ;
			  }
			  coefficientExists=true ;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j ;
			element[n++]=gap ;
		      }
		      OsiRowCut rc ;
		      rc.setLb(-DBL_MAX) ;
		      rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0)) ;
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap) ;
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("c point to row %d\n",irow) ;
#ifdef STRENGTHEN_PRINT
			if (rowLower[irow]<-1.0e20) {
			  printf("7Cut %g <= ",rc.lb()) ;
			  int k ;
			  for ( k=0;k<n;k++) {
			    int iColumn = index[k] ;
			    printf("%g*",element[k]) ;
			    if (si.isInteger(iColumn))
			      printf("i%d ",iColumn) ;
			    else
			      printf("x%d ",iColumn) ;
			  }
			  printf("<= %g\n",rc.ub()) ;
			  printf("Row %g <= ",rowLower[irow]) ;
			  for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			    int iColumn = column[k] ;
			    printf("%g*",rowElements[k]) ;
			    if (si.isInteger(iColumn))
			      printf("i%d ",iColumn) ;
			    else
			      printf("x%d ",iColumn) ;
			  }
			  printf("<= %g\n",rowUpper[irow]) ;
			}
  #endif
			int realRow = (rowLower[irow]<-1.0e20) ? irow : -1 ;
			if (realRows&&realRow>0)
			  realRow=realRows[realRow] ;
			rc.setWhichRow(realRow) ;
			rowCut.addCutIfNotDuplicate(rc) ;
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow] ;
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]] ;
		      }
		    }
		    if (sum-gap*colsol[j]<rowLower[irow]+primalTolerance_||(info->strengthenRow&&rowUpper[irow]>1.0e20)) {
		      // can be a cut
		      // subtract gap from integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL) ;
		      double * element = saveU ;
		      int n=0 ;
                      bool coefficientExists=false ;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow] ;
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk] ;
			  element[n++]=rowElements[kk] ;
			} else {
			  double value=rowElements[kk]-gap ;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk] ;
			    element[n++]=value ;
			  }
			  coefficientExists=true ;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j ;
			element[n++]=-gap ;
		      }
		      OsiRowCut rc ;
		      rc.setLb(rowLower[irow]-gap*(colUpper[j]-1)) ;
		      rc.setUb(DBL_MAX) ;
		      // effectiveness is how far j moves
		      rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap) ;
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false) ;
#if CGL_DEBUG > 0
			if (debugger) assert(!debugger->invalidCut(rc)); 
#endif
                        //if(info->strengthenRow)
                        //printf("d point to row %d\n",irow) ;
#ifdef STRENGTHEN_PRINT
			if (rowUpper[irow]>1.0e20) {
			  printf("8Cut %g <= ",rc.lb()) ;
			  int k ;
			  for ( k=0;k<n;k++) {
			    int iColumn = index[k] ;
			    printf("%g*",element[k]) ;
			    if (si.isInteger(iColumn))
			      printf("i%d ",iColumn) ;
			    else
			      printf("x%d ",iColumn) ;
			  }
			  printf("<= %g\n",rc.ub()) ;
			  printf("Row %g <= ",rowLower[irow]) ;
			  for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			    int iColumn = column[k] ;
			    printf("%g*",rowElements[k]) ;
			    if (si.isInteger(iColumn))
			      printf("i%d ",iColumn) ;
			    else
			      printf("x%d ",iColumn) ;
			  }
			  printf("<= %g\n",rowUpper[irow]) ;
			}
  #endif
			int realRow = (rowUpper[irow]>1.0e20) ? irow : -1 ;
			if (realRows&&realRow>0)
			  realRow=realRows[realRow] ;
			rc.setWhichRow(realRow) ;
			rowCut.addCutIfNotDuplicate(rc) ;
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR] ;
		maxR[irow]=saveMax[istackR] ;
		markR[irow]=-1 ;
	      }
	    }      // end of PROBE: DOWN AND UP FEASIBLE
	  }      // end of PROBE LOOP: ITERATE      
	}      // PROBE LOOP: END
      }
    }      // LOOKEDAT LOOP: END
  }      // PASS LOOP: END
  delete [] cliqueStack ;
  delete [] cliqueCount ;
  delete [] to_01 ;
  delete [] stackC0 ;
  delete [] lo0 ;
  delete [] up0 ;
  delete [] tempL ;
  delete [] tempU ;
  delete [] markC ;
  delete [] stackC ;
  delete [] stackR ;
  delete [] saveL ;
  delete [] saveU ;
  delete [] saveMin ;
  delete [] saveMax ;
  delete [] index ;
  delete [] element ;
  delete [] djs ;
  delete [] colsol ;

  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info->strengthenRow,0) ;
  }

# if CGL_DEBUG > 0
  int numberRowCutsAfter = cs.sizeRowCuts() ;
  int numberColCutsAfter = cs.sizeColCuts() ;

  std::cout
    << "Leaving CglProbing::probeCliques, ninfeas " << ninfeas << ", "
    << (numberRowCutsAfter-numberRowCutsBefore) << " row cuts, "
    << (numberColCutsAfter-numberColCutsBefore) << " col cuts."
    << std::endl ;
 
  if ((numberRowCutsAfter-numberRowCutsBefore) > 0) {
    std::cout << "  row cuts:" << std::endl ;
    for (int k = numberRowCutsBefore ; k < numberRowCutsAfter ; k++) {
      OsiRowCut thisCut = cs.rowCut(k) ;
      thisCut.print() ;
    }
  }
  if ((numberColCutsAfter-numberColCutsBefore) > 0) {
    std::cout << "  col cuts:" << std::endl ;
    for (int k = numberColCutsBefore ; k < numberColCutsAfter ; k++) {
      OsiColCut thisCut = cs.colCut(k) ;
      thisCut.print() ;
    }
  }
# endif
  return (ninfeas) ;
}
