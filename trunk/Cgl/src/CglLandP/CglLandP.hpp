// Copyright (C) 2005, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami                                                  
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     07/21/05
//---------------------------------------------------------------------------
#ifndef CglLandP_H
#define CglLandP_H

//Several level of Debug
// 1 - A few simple sanity check
// 2 - do extra computations on cut quality
// 3 - Output cut in initial non-basic space (only if logLevel >= 3 )
// 4 - Double check (compute in two different ways) reduced cost and f+ f-
//#define LandP_DEBUG 1


#include "CglLandPValidator.hpp"
#include "CglCutGenerator.hpp"
#include "CglParam.hpp"

#ifdef DO_STAT
#include "CglLandPStats.hpp"
#endif
#include <iostream>
#include <fstream>
class CoinWarmStartBasis;
/** Performs one round of Lift & Project using CglLandPSimplex
    to build cuts
*/

namespace LAP {
enum LapMessagesTypes
{
  BEGIN_ROUND,
  END_ROUND,
  DURING_SEP,
  CUT_REJECTED,
  CUT_FAILED,
  LAP_CUT_FAILED_DO_MIG,
  LAP_MESSAGES_DUMMY_END
};
/** Output messages for Cgl */
class LapMessages : public CoinMessages
{
public:
  /** Constructor */
  LapMessages( );
};
class CglLandPSimplex;
}

class CglLandP : public CglCutGenerator
{
  friend void CglLandPUnitTest(OsiSolverInterface *si, const std::string & mpsDir);

  friend class LAP::CglLandPSimplex;

public:

  enum SelectionRules {
    mostNegativeRc /** select most negative reduced cost */,
    bestPivot /** select best possible pivot.*/
  };
  /** Class storing parameters.
      \remark I take all parameters from Ionut's code */
  class Parameters : public CglParam {
  public:
    /** Default constructor (with default values)*/
    Parameters();
    /** Copy constructor */
    Parameters(const Parameters &other);
    /** Assignment opertator */
    Parameters & operator=(const Parameters &other);
    /// @name integer parameters
    ///@{
    
    /** Max number of pivots before we generate the cut 
      \default 20 */
    int pivotLimit;
    /** Max number of pivots at regular nodes. Put a value if you want it lower than the global pivot limit.
     \default 100.*/
    int pivotLimitInTree;
    /** Maximum number of cuts generated at a given round*/
    int maxCutPerRound;
    /** Maximum number of failed pivots before aborting */
    int failedPivotLimit;
    /** maximum number of consecutive degenerate pivots
      \default 0 */
    int degeneratePivotLimit;
    ///@}
    /// @name double parameters
    ///@{
    /** Tolerance for small pivots values (should be the same as the solver */
    double pivotTol;
    /** A variable have to be at least away from integrity to be generated */
    double away;
    /** Total time limit for cut generation.*/
    mutable double timeLimit;
    /** Time limit for generating a single cut.*/
    double singleCutTimeLimit;
    ///@}

    /// @name Flags
    ///@{
    /** Do we use tableau row or the disjunction (I don't really get that there should be a way to always use the tableau)*/
    bool useTableauRow;
    /** Do we apply Egon Balas's Heuristic for modularized cuts */
    bool modularize;
    /** Do we strengthen the final cut (always do if modularize is 1) */
    bool strengthen;
    /** Work in the reduced space (only non-structurals enter the basis) */
    int reducedSpace;
    /** Apply perturbation procedure. */
    bool perturb;
    /** Scale the cuts */
    bool scaleCuts;
    /** Which rule to apply for choosing entering and leaving variables.*/
    SelectionRules pivotSelection;
    ///@}
  };


  /** Constructor for the class*/
   CglLandP(const CglLandP::Parameters &params = CglLandP::Parameters(), 
	    const CglValidator &validator = CglValidator());
  /** Destructor */
  ~CglLandP();
  /** Copy constructor */
   CglLandP(const CglLandP &source);
  /** Assignment operator */
  CglLandP& operator=(const CglLandP &rhs);
  /** Clone function */
  CglCutGenerator * clone() const;

#ifdef DO_STAT
  /** display averages on last rounds*/
  void displayStats(){roundsStats_.displaySumUp(*logStream);}

  void setIdString(const std::string &id){roundsStats_.setIdString(id);}
#endif

    /**@name Generate Cuts */
  //@{

  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
			    const CglTreeInfo info = CglTreeInfo()) const;

  //@}

virtual bool needsOptimalBasis() const
{
  return true;
}

	CglValidator & validator(){return validator_;}
	/** set level of log for cut generation procedure :
    <ol start=0 >
		<li> for none </li>
		<li> for log at begin and end of procedure + at some time interval </li>
		<li> for log at every cut generated </li>
		</ol>
		*/
  void setLogLevel(int level)
  { handler_->setLogLevel(level);
  }
#ifdef DO_STAT
  void setLogStream(std::ostream & os)
  {
    releaseLogStream();
    logStream = &os;
  }

void openLogStream(std::string name){
  releaseLogStream();
  logStream = new std::ofstream(name.c_str(),std::ios::app);
  logAssignCount = new int;
  (*logAssignCount) = 1;
}

inline void releaseLogStream(){
  if(logAssignCount && *logAssignCount){
    (*logAssignCount)--;
    if(*logAssignCount == 0){
      delete logStream;
      delete logAssignCount;
      logAssignCount=NULL;
    }
    else logAssignCount = NULL;
}
}

void lookupProblem(const std::string &s);

#endif
  class NoBasisError : public CoinError
  {
  public:
    NoBasisError(): CoinError("No basis available","LandP",""){}
  };

  class SimplexInterfaceError : public CoinError
  {
  public:
    SimplexInterfaceError(): CoinError("Invalid conversion to simplex interface", "CglLandP","CglLandP"){}
  };
  Parameters & parameter() {return params_;}
private:
  Parameters params_;

  /** Some informations that will be changed by the pivots and that we want to keep*/
  struct CachedData
  {
    CachedData(int nBasics = 0 , int nNonBasics = 0);
    CachedData(const CachedData & source);

    CachedData& operator=(const CachedData &source);
    /** Get the data from a problem */
    void getData(const OsiSolverInterface &si);
    ~CachedData();
    /** Indices of basic variables in starting basis (ordered if variable basics_[i] s basic in row i)*/
    int * basics_;
    /** Indices of non-basic variables */
    int *nonBasics_;
    /** number of basics variables */
    int nBasics_;
    /** number of non-basics */
    int nNonBasics_;
    /** Optimal basis */
    CoinWarmStartBasis * basis_;
    /** Stores the value of the solution to cut */
    double * colsol_;
    /** Stores the values of the slacks */
    double * slacks_;
    /** Stores wheter slacks are integer constrained */
    bool * integers_;
  };
  mutable CachedData cached_;
  /** message handler */
  CoinMessageHandler * handler_;
  /** messages */
  CoinMessages messages_;
  /** cut validator */
  CglValidator validator_;
#ifdef DO_STAT
  /** store statistics on separation */
  mutable roundsStatistics roundsStats_;
  /** log file output */
  mutable std::ostream * logStream;
  /**counter of number time logStream used.*/
  mutable int * logAssignCount;
  /** set to 1 if generator is used (and print statistics).*/
  mutable bool used;
#endif
};
void CglLandPUnitTest(OsiSolverInterface *si, const std::string & mpsDir);

#endif

