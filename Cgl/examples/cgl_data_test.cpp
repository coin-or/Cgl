#include <cstdio>
#include<stdlib.h>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "OsiSolverParameters.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CglRedSplit.hpp"

int main(int argc, char **argv) 
{
  char *f_name_lp, *last_dot_pos, f_name[256], *f_name_pos;
  int i, ncol;

  if((argc < 2) || (argc > 2)) {
    printf("### ERROR: main(): Usage: cgl_data_test input_file_name.mps\n");
    exit(1);
  }

  f_name_lp = strdup(argv[1]);
  f_name_pos = strrchr(f_name_lp, '/');
  if(f_name_pos != NULL) {
    strcpy(f_name, &(f_name_pos[1]));
  }
  else {
    strcpy(f_name, f_name_lp);
  }
  last_dot_pos = strrchr(f_name, '.');
  last_dot_pos[0] = '\0';

  OsiClpSolverInterface *clp = new OsiClpSolverInterface;
  clp->messageHandler()->setLogLevel(0);
  if(strcmp(&(f_name_lp[strlen(f_name_lp)-3]), ".lp") == 0) {
    clp->readLp(f_name_lp);    
  }
  else {
    if(strcmp(&(f_name_lp[strlen(f_name_lp)-4]), ".mps") == 0) {
      clp->readMps(f_name_lp);    
    }
    else {
      printf("### ERROR: unrecognized file type\n");
      exit(1);
    }
  }
  ncol = clp->getNumCols();
  clp->initialSolve();

  printf("LP value: %12.2f\n", clp->getObjValue());

  OsiCuts cuts;

  // Define parameters for CglRedSplit generator
  CglParam cpar;
  cpar.setMAX_SUPPORT(ncol+1);
  CglRedSplitParam rspar(cpar);

  // Create a cut generator with the given parameters
  CglRedSplit cutGen(rspar);

  // Create an empty object to hold data for the CglRedSplit generator
  CglRedSplitData rsdat;

  char *colType = new char[ncol];
  for(i=0; i<ncol; i++) {
    if(clp->isContinuous(i)) {
      colType[i] = 'C';
    }
    else {
      colType[i] = 'I';
    }
  }

  // set data that will not change
  rsdat.setNcol(clp->getNumCols()); 
  rsdat.setColType(colType);

  int round, max_rounds = 10;
  for(round=0; round<max_rounds; round++) {

    // Setup the data; only members used by the CglRedSPlitGenerator
    // see CglRedSplit::generateCuts(const CglRedSplitData rsData, OsiCuts &cs)

    rsdat.setNrow(clp->getNumRows()); 
    const double *curr_colLower = clp->getColLower();
    const double *curr_colUpper = clp->getColUpper();
    const double *curr_rowLower = clp->getRowLower();
    const double *curr_rowUpper = clp->getRowUpper();
    const double *curr_rowRhs = clp->getRightHandSide();
    const double *xlp = clp->getColSolution();
    const double *row_act = clp->getRowActivity();
    const CoinPackedMatrix *byRow = clp->getMatrixByRow();

    rsdat.setSolverPtr(clp);
    rsdat.setOptimalBasisIsAvailable(clp->optimalBasisIsAvailable());
    rsdat.setColLower(curr_colLower);
    rsdat.setColUpper(curr_colUpper);
    rsdat.setRowLower(curr_rowLower);
    rsdat.setRowUpper(curr_rowUpper);
    rsdat.setRowRhs(curr_rowRhs);
    
    rsdat.setSeparateThis(xlp);
    rsdat.setRowActivity(row_act);
    rsdat.setMatrixByRow(byRow);

    // Use the following for calling the cut generator with the data
    cutGen.generateCuts(rsdat, cuts);

    // Use the following for standard way to call the cut generator
     cutGen.generateCuts(*clp, cuts);

    int ncuts = cuts.sizeRowCuts();

    const OsiRowCut **newRowCuts = new const OsiRowCut * [ncuts];
    for(i=0; i<ncuts; i++) {
      newRowCuts[i] = &cuts.rowCut(i); 
    }
    clp->applyRowCuts(ncuts, newRowCuts);
    delete[] newRowCuts;

    printf("round %4d: %4d generated cuts  new objective value: %12.2f\n", 
	   round, ncuts, clp->getObjValue());

    clp->resolve();  

    if(clp->isAbandoned()) {
      printf("###ERROR: Numerical difficulties in Solver\n");
      exit(1);
    }
  
    if(clp->isProvenPrimalInfeasible()) {
      printf("### WARNING: Problem is infeasible\n");
    }
  }

  delete clp;
  free(f_name_lp);
  delete[] colType;

  return(0);
}
