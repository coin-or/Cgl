// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <string>
#include <cassert>
#include <cmath>
#include <cfloat>

#include "CglProcess.hpp"
//#include "CglMessage.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"
#include "CglCutGenerator.hpp"
#include "CoinTime.hpp"


OsiSolverInterface *
CglProcess::preProcess(OsiSolverInterface & model, 
                       bool makeEquality, int numberPasses)
{
  originalModel_ = & model;
  numberSolvers_ = numberPasses;
  model_ = new OsiSolverInterface * [numberSolvers_];
  modifiedModel_ = new OsiSolverInterface * [numberSolvers_];
  presolve_ = new OsiPresolve * [numberSolvers_];
  for (int i=0;i<numberSolvers_;i++) {
    model_[i]=NULL;
    modifiedModel_[i]=NULL;
    presolve_[i]=NULL;
  }
  startModel_=&model;
  CoinPackedMatrix matrixByRow(*originalModel_->getMatrixByRow());
  int numberRows = originalModel_->getNumRows();
  int numberColumns = originalModel_->getNumCols();
  
  // We want to add columns
  int numberSlacks=0;
  int * rows = new int[numberRows];
  double * element =new double[numberRows];
  
  int iRow;
  
  int numberCliques=0;
  int * which = new int[numberColumns];
  
  // Statistics
  int totalP1=0,totalM1=0;
  int numberFixed=0;
  // May just find it is infeasible
  bool feasible=true;
  
  // Row copy
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();
  
  const double * lower = originalModel_->getColLower();
  const double * upper = originalModel_->getColUpper();
  const double * rowLower = originalModel_->getRowLower();
  const double * rowUpper = originalModel_->getRowUpper();
  
  for (iRow=0;iRow<numberRows;iRow++) {
    int numberP1=0, numberM1=0;
    int j;
    double upperValue=rowUpper[iRow];
    double lowerValue=rowLower[iRow];
    bool good=true;
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn = column[j];
      if (upper[iColumn]-lower[iColumn]<1.0e-8) {
        // fixed
        upperValue -= lower[iColumn]*elementByRow[j];
        lowerValue -= lower[iColumn]*elementByRow[j];
        continue;
      } else if (!originalModel_->isBinary(iColumn)) {
        good = false;
        break;
      }
      if (fabs(elementByRow[j])!=1.0) {
        good=false;
        break;
      } else if (elementByRow[j]>0.0) {
        which[numberP1++]=iColumn;
      } else {
        numberM1++;
        which[numberColumns-numberM1]=iColumn;
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
        feasible=false;
        break;
      } else if (abs(state)==2) {
        // we can fix all
        numberFixed += numberP1+numberM1;
        int i;
        if (state>0) {
          // fix all +1 at 0, -1 at 1
          for (i=0;i<numberP1;i++)
            originalModel_->setColUpper(which[i],0.0);
          for (i=0;i<numberM1;i++)
            originalModel_->setColLower(which[numberColumns-i-1],1.0);
        } else {
          // fix all +1 at 1, -1 at 0
          for (i=0;i<numberP1;i++)
            originalModel_->setColLower(which[i],1.0);
          for (i=0;i<numberM1;i++)
            originalModel_->setColUpper(which[numberColumns-i-1],0.0);
        }
      } else {
        numberCliques++;
        if (iLower!=iUpper) {
          element[numberSlacks]=state;
          rows[numberSlacks++]=iRow;
        }
        if (state>0) {
          totalP1 += numberP1;
          totalM1 += numberM1;
        } else {
          totalP1 += numberM1;
          totalM1 += numberP1;
        }
      }
    }
  }
  delete [] which;
  if (!feasible) {
    printf("*** Problem infeasible\n");
    return NULL;
  } else {
    if (numberCliques) {
      printf("%d cliques of average size %g found, %d P1, %d M1\n",
	     numberCliques,
	     ((double)(totalP1+totalM1))/((double) numberCliques),
	     totalP1,totalM1);
      printf("%d of these could be converted to equality constraints\n",
             numberSlacks);
    }
    if (numberFixed)
      printf("%d variables fixed\n",numberFixed);
  }
  if (numberSlacks&&makeEquality) {
    // add variables to make equality rows
    // Get new model
    startModel_ = originalModel_->clone();
    for (int i=0;i<numberSlacks;i++) {
      int iRow = rows[i];
      double value = element[i];
      double lowerValue = 0.0;
      double upperValue = 1.0;
      double objValue  = 0.0;
      CoinPackedVector column(1,&iRow,&value);
      startModel_->addCol(column,lowerValue,upperValue,objValue);
      // set integer
      startModel_->setInteger(numberColumns+i);
      if (value >0)
	startModel_->setRowLower(iRow,rowUpper[iRow]);
      else
	startModel_->setRowUpper(iRow,rowLower[iRow]);
    }
  }
  delete [] rows;
  delete [] element;
   
  OsiSolverInterface * returnModel;
  int numberChanges;
  if (!numberSolvers_) {
    // just fix
    startModel_->initialSolve();
    OsiSolverInterface * newModel = modified(startModel_,false,numberChanges);
    if (startModel_!=originalModel_)
      delete startModel_;
    startModel_=newModel;
    returnModel=startModel_;
  } else {
    OsiSolverInterface * presolvedModel;
    OsiSolverInterface * oldModel = startModel_;
    for (int iPass=0;iPass<numberSolvers_;iPass++) {
      OsiPresolve * pinfo = new OsiPresolve();
      presolvedModel = pinfo->presolvedModel(*oldModel,1.0e-8,true,5);
      if (!presolvedModel) {
        returnModel=NULL;
        break;
      }
      model_[iPass]=presolvedModel;
      presolve_[iPass]=pinfo;
      bool constraints = iPass<numberPasses-1;
      presolvedModel->initialSolve();
      OsiSolverInterface * newModel = modified(presolvedModel,constraints,numberChanges);
      returnModel=newModel;
      if (!newModel) {
        break;
      }
      modifiedModel_[iPass]=newModel;
      if (!numberChanges) {
        numberSolvers_=iPass+1;
        break;
      }
    }
  }
  if (!returnModel) 
    printf("** found to be infeasible\n");
  return returnModel;
}

void
CglProcess::postProcess(OsiSolverInterface & modelIn)
{
  OsiSolverInterface * modelM = &modelIn;
  for (int iPass=numberSolvers_-1;iPass>=0;iPass--) {
    OsiSolverInterface * model = model_[iPass];
    int numberColumns = modelM->getNumCols();
    const double * solutionM = modelM->getColSolution();
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (modelM->isInteger(iColumn)) {
        double value = solutionM[iColumn];
        double value2 = floor(value+0.5);
        assert (fabs(value-value2)<1.0e-3);
        model->setColLower(iColumn,value2);
        model->setColUpper(iColumn,value2);
      }
    }
    model->resolve();
    presolve_[iPass]->postsolve(false);
    delete modifiedModel_[iPass];;
    delete model_[iPass];;
    delete presolve_[iPass];
    modifiedModel_[iPass]=NULL;
    model_[iPass]=NULL;
    presolve_[iPass]=NULL;
    if (iPass)
      modelM = modifiedModel_[iPass-1];
    else
      modelM = startModel_;
  }
  // should be back to startModel_;
  OsiSolverInterface * model = originalModel_;
  // Use number of columns in original
  int numberColumns = model->getNumCols();
  const double * solutionM = modelM->getColSolution();
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (modelM->isInteger(iColumn)) {
      double value = solutionM[iColumn];
      double value2 = floor(value+0.5);
      assert (fabs(value-value2)<1.0e-3);
      model->setColLower(iColumn,value2);
      model->setColUpper(iColumn,value2);
    }
  }
  model->resolve();
}
/* Return model with useful modifications.  
   If constraints true then returns any x+y=1 or x-y=0 constraints
*/
OsiSolverInterface * 
CglProcess::modified(OsiSolverInterface * model,
                     bool constraints,
                     int & numberChanges)
{
  return model->clone();
}

/* Default Constructor

*/
CglProcess::CglProcess() 

:
  originalModel_(NULL),
  startModel_(NULL),
  numberSolvers_(0),
  model_(NULL),
  modifiedModel_(NULL),
  presolve_(NULL),
  handler_(NULL),
  defaultHandler_(true),
  appData_(NULL),
  numberCutGenerators_(0),
  generator_(NULL)
{
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  //messages_ = CglMessage();
}

// Copy constructor.

CglProcess::CglProcess(const CglProcess & rhs)
:
  numberSolvers_(rhs.numberSolvers_),
  defaultHandler_(rhs.defaultHandler_),
  appData_(rhs.appData_),
  numberCutGenerators_(rhs.numberCutGenerators_)
{
  if (defaultHandler_) {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(rhs.handler_->logLevel());
  } else {
    handler_ = rhs.handler_;
  }
  if (numberCutGenerators_) {
    generator_ = new CglCutGenerator * [numberCutGenerators_];
    for (int i=0;i<numberCutGenerators_;i++) {
      generator_[i]=rhs.generator_[i]->clone();
    }
  } else {
    generator_=NULL;
  }
  if (rhs.originalModel_) {
    originalModel_ = rhs.originalModel_;
    // If no make equality then solvers are same
    if (rhs.originalModel_!=rhs.startModel_) {
      startModel_=rhs.startModel_->clone();
    } else {
      startModel_=originalModel_;
    }
  } else {
    originalModel_=NULL;
    startModel_=NULL;
  }
  if (numberSolvers_) {
    model_ = new OsiSolverInterface * [numberSolvers_];
    modifiedModel_ = new OsiSolverInterface * [numberSolvers_];
    presolve_ = new OsiPresolve * [numberSolvers_];
    for (int i=0;i<numberSolvers_;i++) {
      model_[i]=rhs.model_[i]->clone();
      modifiedModel_[i]=rhs.modifiedModel_[i]->clone();
      presolve_[i]=new OsiPresolve(*rhs.presolve_[i]);
    }
  } else {
    model_=NULL;
    presolve_=NULL;
  }
}
  
// Assignment operator 
CglProcess & 
CglProcess::operator=(const CglProcess& rhs)
{
  if (this!=&rhs) {
    gutsOfDestructor();
    numberSolvers_=rhs.numberSolvers_;
    defaultHandler_=rhs.defaultHandler_;
    appData_=rhs.appData_;
    numberCutGenerators_=rhs.numberCutGenerators_;
    if (defaultHandler_) {
      handler_ = new CoinMessageHandler();
      handler_->setLogLevel(rhs.handler_->logLevel());
    } else {
      handler_ = rhs.handler_;
    }
    if (numberCutGenerators_) {
      generator_ = new CglCutGenerator * [numberCutGenerators_];
      for (int i=0;i<numberCutGenerators_;i++) {
        generator_[i]=rhs.generator_[i]->clone();
      }
    }
    if (rhs.originalModel_) {
      originalModel_ = rhs.originalModel_;
      // If no make equality then solvers are same
      if (rhs.originalModel_!=rhs.startModel_) {
        startModel_=rhs.startModel_->clone();
      } else {
        startModel_=originalModel_;
      }
    } else {
      originalModel_=NULL;
      startModel_=NULL;
    }
    if (numberSolvers_) {
      model_ = new OsiSolverInterface * [numberSolvers_];
      modifiedModel_ = new OsiSolverInterface * [numberSolvers_];
      presolve_ = new OsiPresolve * [numberSolvers_];
      for (int i=0;i<numberSolvers_;i++) {
        model_[i]=rhs.model_[i]->clone();
        modifiedModel_[i]=rhs.modifiedModel_[i]->clone();
        presolve_[i]=new OsiPresolve(*rhs.presolve_[i]);
      }
    } else {
      model_=NULL;
      presolve_=NULL;
    }
  }
  return *this;
}
  
// Destructor 
CglProcess::~CglProcess ()
{
  gutsOfDestructor();
}
// Clears out as much as possible (except solver)
void 
CglProcess::gutsOfDestructor()
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  if (startModel_!=originalModel_) 
    delete startModel_;
  startModel_=NULL;
  //delete originalModel_;
  originalModel_=NULL;
  int i;
  for (i=0;i<numberCutGenerators_;i++) {
    delete generator_[i];
  }
  delete [] generator_;
  generator_=NULL;
  for (i=0;i<numberSolvers_;i++) {
    delete model_[i];
    delete modifiedModel_[i];
    delete presolve_[i];
  }
  delete [] model_;
  delete [] modifiedModel_;
  delete [] presolve_;
  model_=NULL;
  presolve_=NULL;
}
// Add one generator
void 
CglProcess::addCutGenerator(CglCutGenerator * generator)
{
  CglCutGenerator ** temp = generator_;
  generator_ = new CglCutGenerator * [numberCutGenerators_+1];
  memcpy(generator_,temp,numberCutGenerators_*sizeof(CglCutGenerator *));
  delete[] temp ;
  generator_[numberCutGenerators_]=generator->clone(); 
}
//#############################################################################
// Set/Get Application Data
// This is a pointer that the application can store into and retrieve
// This field is the application to optionally define and use.
//#############################################################################

void CglProcess::setApplicationData(void * appData)
{
  appData_ = appData;
}
//-----------------------------------------------------------------------------
void * CglProcess::getApplicationData() const
{
  return appData_;
}
/* Set cutoff bound on the objective function.
   
When using strict comparison, the bound is adjusted by a tolerance to
avoid accidentally cutting off the optimal solution.
*/
void 
CglProcess::setCutoff(double value) 
{
  // Solvers know about direction
  double direction = originalModel_->getObjSense();
  originalModel_->setDblParam(OsiDualObjectiveLimit,value*direction); 
}

// Get the cutoff bound on the objective function - always as minimize
double 
CglProcess::getCutoff() const
{ 
  double value ;
  originalModel_->getDblParam(OsiDualObjectiveLimit,value) ;
  return value * originalModel_->getObjSense() ;
}
