// Copyright (C) 2005, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     21/07/05
//---------------------------------------------------------------------------
#ifndef CglLandPSimplex_H
#define CglLandPSimplex_H

#include <iostream>

#include "CglLandP.hpp"

#include "OsiSolverInterface.hpp"
#include "CoinMessage.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinPackedMatrix.hpp"
namespace LAP {
/** Types of messages for lift-and-project simplex.*/
enum LAPS_messages
{
  Separating,
  FoundImprovingRow,
  FoundBestImprovingCol,
  WarnFailedBestImprovingCol,
  LogHead,
  PivotLog,
  FinishedOptimal,
  HitLimit,
  WarnBadSigmaComputation,
  WarnBadRowComputation,
  WarnGiveUpRow,
  PivotFailedSigmaUnchanged,
  PivotFailedSigmaIncreased,
  FailedSigmaIncreased,
  WarnBadRhsComputation,
  WarnFailedPivotTol,
  WarnFailedPivotIIf,
  DUMMY_END
};
/** Message handler for lift-and-project simplex. */
class LandPMessages : public CoinMessages
{
public:

  /** Constructor */
  LandPMessages();
};



class CglLandPSimplex
{
public:
    struct TabRow
  {
    TabRow():num(-1), row(NULL), rhs(0), n_(0) {}
    TabRow(const TabRow & source):num(source.num), row(NULL), rhs(source.rhs), n_(source.n_)
    {
      row = new double[n_];
      CoinCopyN(source.row, n_, row);
    }
    void size(int n){if(row != NULL) delete [] row; row = new double[n]; n_=n;}
    ~TabRow(){
      if (row != NULL)
	delete [] row;
    }
    void print(std::ostream & os, int width = 9, const int * nonBasics = NULL,
 int m = 0)
    {
      os.width(3);
      os.precision(4);
      os.setf(std::ios_base::right, std::ios_base::adjustfield);
      os<< num <<": ";
      for(int j = 0 ; j < m ; j++)
	{
	  os.width(width);
	  os.precision(3);
	  //      os.setf(std::ios_base::fixed, std::ios_base::floatfield);
	  os.setf(std::ios_base::right, std::ios_base::adjustfield);
	  os<<row[nonBasics[j]]<<" ";
	}
      
      os.width(width);
      os.precision(4);
      //    os.setf(std::ios_base::fixed, std::ios_base::floatfield);
      os.setf(std::ios_base::right, std::ios_base::adjustfield);
      os<<rhs;
      
      os<<std::endl;

    }
    int num;
    double *row;
    double rhs;
    int n_;
  };
  /** Usefull onstructor */
  CglLandPSimplex(const OsiSolverInterface &si, 
		  const CglLandP::CachedData &cached,
                  bool reducedSpace, int numPivots);
  /** Destructor */
  ~CglLandPSimplex();
  /**Update cached information in case of basis change in a round*/
  void cacheUpdate(const CglLandP::CachedData &cached, bool reducedSpace = 0);
  /** reset the solver to optimal basis */
  bool resetSolver(const CoinWarmStartBasis * basis);
  /** Perfom pivots to find the best cuts */
  bool findBestCut(int var, OsiRowCut & cut, const CglLandP::CachedData &cached, const CglLandP::Parameters & params);
  /** Find Gomory cut (i.e. don't do extra setup required for pivots).*/
  bool generateMig(int row, OsiRowCut &cut, const CglLandP::CachedData &cached, const CglLandP::Parameters & params) const;
  
  /// No default constructor
  CglLandPSimplex()
#ifdef LandP_DEBUG
    :debug_(0,0)
#endif
  {
    std::cerr<<"Class CglLandPSimplex has no default constructor"<<std::endl;
    throw CoinError("No default constructor","CglLandPSimplex","CglLandPSimplex()");
  }
    void setLogLevel(int level)
  { handler_->setLogLevel(level);}


  void setSi(OsiSolverInterface *si)
  {si_ = si;}
  void freeSi()
  {delete si_;}

  struct extraInfo {
    extraInfo():
#if LandP_DEBUG > 1
      Nnz(0), Nneg(0),
#endif
      AngleToObj(0),
      numPivots(0), depth(0.),
      nNegativeRcRows(0), bestRow(0),
      maxBestRow(0), bestRc(0.),
      maxRc(-DBL_MAX)
    {}
#if LandP_DEBUG > 1
    int Nnz;
    int Nneg;
#endif
    double AngleToObj;
    int numPivots;
    double depth;
    int nNegativeRcRows;
    int bestRow;
    int maxBestRow;
    double bestRc;
    double maxRc;
  };
  mutable extraInfo extra;

protected:
  /** Perform a change in the basis (direction is 1 if leaving variable is going to ub, 0 otherwise)*/
  bool changeBasis(int incoming, int leaving, int direction);
  /** Find a row which can be used to perform an improving pivot return index of the cut or -1 if none exists
    * (i.e., find the leaving variable).*/
  int findCutImprovingPivotRow( int &direction, int &gammaSign, double tolerance);
  /** Find a row which can be used to perform an improving pivot the fast way
    * (i.e., find the leaving variable).
    \return index of the cut or -1 if none exists. */
  int fastFindCutImprovingPivotRow( int &direction, int &gammaSign, double tolerance);
  /** Rescan reduced costs tables */
  int rescanReducedCosts( int &direction, int &gammaSign, double tolerance);
  /** Find the column which leads to the best cut (i.e., find incoming variable).*/
  int fastFindBestPivotColumn(int direction, int gammaSign,
                          double pivotTol, double rhsTol, 
                          bool reducedSpace,  
                          bool allowNonImproving, 
                          double &bestSigma);
  /** Find the column which leads to the best cut (i.e., find incoming variable).*/
  int findBestPivotColumn(int direction, 
                          double pivotTol, bool reducedSpace, bool allowDegeneratePivot,
                          bool modularize);
#if 0
int plotCGLPobj(int direction, double gammaTolerance,
                double pivotTol, bool reducedSpace, bool allowDegenerate, bool modularize);
#endif

  /** Find incoming and leaving variables which lead to the most violated
      adjacent normalized lift-and-project cut. 
      \remark At this point reduced costs should be already computed.
      \return incoming variable variable,
      \param leaving variable
      \param direction leaving direction
      \param numTryRows number rows tried
      \param pivotTol pivot tolerance
      \param reducedSpace separaration space (reduced or full)
      \param allowNonStrictlyImproving wether or not to allow non stricly improving pivots.      
    */
  int findBestPivot(int &leaving, int & direction, 
                    const CglLandP::Parameters & params);
  /** Compute the objective value of the Cglp with linear combintation of the two rows by gamma */
  double computeCglpObjective(double gamma, bool strengthen);
  /** Compute the objective value of the Cglp for given row and rhs (if strengthening shall be applied
    	row should have been modularized).*/
  double computeCglpObjective(TabRow & row) const;
  /** return the coefficients of the intersection cut */
   inline double intersectionCutCoef(double alpha_i, double beta) const;
  /** return the coefficients of the strengthened intersection cut 
   * takes one extra argument seens needs to consider variable type.
   */
  inline double strengthenedIntersectionCutCoef(int i, double alpha_i, double beta) const;
  /** return the coefficient of the new row (combining row_k + gamma row_i).
   */
  inline double newRowCoefficient(int j, double gamma) const;
  /** compute the modularized row coefficient for an integer variable (does not check integrity
   \bug not the same as in Ionut's pseudo code, it is not the same*/
  inline double modularizedCoef(double alpha, double beta) const;
  /** Compute the reduced cost of Cglp */
  double computeCglpRedCost(int direction, int gammaSign);
  /** Compute the value of sigma and thau (which are constants for a row i as defined in Mike Perregaard thesis */
  void computeRedCostConstantsInRow();
  /** Create the intersection cut of row k*/
  void createIntersectionCut(const TabRow & row, OsiRowCut &cut) const;
  /** Create strenghtened row */
	//  void createIntersectionCut(double * row);
  /** Create MIG cut from row k*/
  void createMIG( TabRow & row, OsiRowCut &cut) const;
  /** Compute normalization coefficient corresponding to sum of multipliers
      for cuts derived from row */
  double normCoef(TabRow &row);
  /** scale the cut passed as argument using provided normalization factor*/
  void scale(OsiRowCut &cut, double norma);
  /** scale the cut passed as argument*/
  void scale(OsiRowCut &cut);
  /** Modularize row.*/
  void modularizeRow(TabRow & row);
  /** Get the row i of the tableau */
  void pullTableauRow(TabRow & row) const;
  /** Adjust the row of the tableau to reflect leaving variable direction */
  void adjustTableauRow(int var, TabRow & row, int direction);
  /** reset the tableau row after a call to adjustTableauRow */
  void resetOriginalTableauRow(int var, TabRow & row, int direction);
  /**Get lower bound for variable or constraint */
  inline double getLoBound(int index){return loBounds_[index];}
  /**Get upper bound for variable or constraint */
  inline double getUpBound(int index){return upBounds_[index];}
  /** print the tableau of current basis. */
  void printTableau(std::ostream & os);
private:

  enum lpSolver {clp, cplex, xpress};
  lpSolver solver_;
  void updateM1_M2_M3(TabRow & row, double tolerance, bool recucedSpace, bool alwaysComputeCheap);

  /// @name Work infos
  /// @{
  /** Source row for cut */
  mutable CglLandPSimplex::TabRow row_k_;
  /** Perturbed source row */
  CglLandPSimplex::TabRow perturbed_row_k_;
  /** Row of leaving candidate*/
  CglLandPSimplex::TabRow row_i_;
  /** Stores the new row. */
  CglLandPSimplex::TabRow newRow_;
  /**vector to sort the gammas*/
  CoinPackedVector gammas_;
  /**first work vector in row space.*/
  double * rWk1_;
  /**scond work vector in row space.*/
  double * rWk2_;
  /**third work vector in row space.*/
  double * rWk3_;
  /**fourth work vector in row space.*/
  double * rWk4_;
  /** integer valued work vector on the rows */
  int * rIntWork_;
  /** Flag rows which we don't want to try anymore */
  bool * rowFlags_;
  /** Flag columns for which contribution has to be computed in reduced cost(if in reduced space remove nonbasic structurals) */
  bool *colHasToComputeContrib_;
  /** Flag columns which have to be considered for leaving the basis */
  bool *colCandidateToLeave_;
  /** Store the basics variable */
  int * basics_;
  /** Stores the nonBasicVariables */
  int * nonBasics_;
  /** Stores the variables which are always in M1 for a given k*/
  int * inM1_;
  /** Stores the variables which are always in M2 for a given k*/
  int * inM2_;
  /** Stores the variables which could be either in M1 or M2 */
  int * inM3_;
  /** temporarily stores the value of thau for the current row */
  double tau_;
  /** stores the cglp value of the normalized cut obtained from row k_ */
  double sigma_;
  /** Keep track of basis status */
  CoinWarmStartBasis basis_;
  /** Pointer to the solution to cut (need to be modified after each pivot because we are only considering slacks).*/
   double * colsolToCut_;
 /** Pointer to the current basic solution.*/
   double * colsol_;
  /// cached numcols
  int numcols_;
  ///cached numrows
  int numrows_;
  // for fast access to lower bounds (both cols and rows)
  double * loBounds_;
  // for fast access to upper bounds (both cols and rows)
  double * upBounds_;
  /// Say if we are in a sequence of degenerate pivots
  bool inDegenerateSequence_;
  /// Value for the reduced cost chosen for pivoting
  double chosenReducedCostVal_;
  /// pointer to array of integer info for both structural and slacks
  const bool * integers_;
  /// @}
  /// @name Interfaces to the solver
  /// @{
  /** Pointer to the solver interface */
  OsiSolverInterface * si_;
  ///@}
  /// Own the data or not?
  bool own_;


  /// number of negative reduced costs in current iteration
  int nNegativeRc_;
  /// number of rows with a <0 rc in current iteration
  int nNegativeRcRows_; 
#ifdef LandP_DEBUG

#if LandP_DEBUG > 1
/** Create MIG cut from row k and express it in the original non-basic space*/
void put_in_non_basic_init_space( OsiRowCut &cut);
#endif
  friend class DebugData;
  /** Data for debugging */
  class DebugData
  {
  public:
    DebugData(int n, int m):
      bestNewRow_(NULL),
      req(1e-05),
      eq(1e-5)
#if LandP+DEBUG> 1
      , initialTableau_(), initialBasics_(NULL), initialBasis_(NULL)
#endif
      
    { 
      bestNewRow_ = new double[n + m];
#if LandP_DEBUG > 1
      initialBasics_ = new int[m];
      initialNonBasics_ = new int[n];
      initialColsol_ = new double[n + m];
      trueInitialSol_ = new double[n + m];
#endif
    }
    ~DebugData()
    {
      delete [] bestNewRow_;
#if LandP_DEBUG > 1
      delete [] initialBasics_;
      delete [] initialNonBasics_;
      delete [] initialColsol_;
      delete [] trueInitialSol_;
      if(initialBasis_)
	delete initialBasis_;
#endif
    }
    /** stores the new row as computed when looking for the best entering column.*/
    double * bestNewRow_;
    /** stores the new rhs as computed when looking for the best entering column.*/
    double bestNewRhs_;
    /** Stores the next Cglp objective function */
    double newSigma_;
    /** a relative equality checker */
    CoinRelFltEq req;
    /** an absolute equality checker */
    CoinAbsFltEq eq;

#if LandP_DEBUG > 1
    /** Get the tableau corresponding to the current basis.*/
    void getCurrentTableau(OsiSolverInterface &si, CglLandPSimplex &lap);
      /** Stores the initial tableau */
    CoinPackedMatrix initialTableau_;
    /** Stores the initial basic variables */
    int * initialBasics_;
    /** Stores the initial non basics*/
    int *initialNonBasics_;
    /** Stores the initial basis */
    CoinWarmStartBasis * initialBasis_;
    /** stores the initial solution */
    double * initialColsol_;
    /** stores the initial solution with unmodified non basics*/
    double * trueInitialSol_;
#endif
  };
  DebugData debug_;
#endif
  /// @}
  /** Message handler. */
  CoinMessageHandler * handler_;
  /** Messages. */
  CoinMessages messages_;
};

/** return the coefficients of the intersection cut */
double 
CglLandPSimplex::intersectionCutCoef(double alpha_i, double beta) const
{
  if(alpha_i>0) return alpha_i* (1 - beta);
  else return -alpha_i * beta;// (1 - beta);
}

/** return the coefficients of the strengthened intersection cut */
double 
CglLandPSimplex::strengthenedIntersectionCutCoef(int i, double alpha_i, double beta) const
{
  //  double ratio = beta/(1-beta);
  if( (!integers_[i]))
    return intersectionCutCoef(alpha_i, beta);
  else
  {
    double f_i = alpha_i - floor(alpha_i);
    if(f_i < beta)
      return f_i*(1- beta);
    else
      return (1 - f_i)*beta;//(1-beta);
  }
}

/** return the coefficient of the new row (combining row_k + gamma row_i).
*/
double 
CglLandPSimplex::newRowCoefficient(int j, double gamma) const
{
  return row_k_.row[j] + gamma * row_i_.row[j];
}

/** compute the modularized row coefficient for an integer variable*/
double 
CglLandPSimplex::modularizedCoef(double alpha, double beta) const
{
  double f_i = alpha - floor(alpha);
  if(f_i <= beta)
    return f_i;
  else
    return f_i - 1;
}




}
#endif
