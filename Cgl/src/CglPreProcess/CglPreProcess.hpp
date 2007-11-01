// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglPreProcess_H
#define CglPreProcess_H

#include <string>
#include <vector>

#include "CoinFinite.hpp"
#include "CoinMessageHandler.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiPresolve.hpp"
#include "CglCutGenerator.hpp"

//#############################################################################

/** Class for preProcessing and postProcessing.

    While cuts can be added at any time in the tree, some cuts are actually just
    stronger versions of existing constraints.  In this case they can replace those
    constraints rather than being added as new constraints.  This is awkward in the
    tree but reasonable at the root node.

    This is a general process class which uses other cut generators to strengthen
    constraints, establish that constraints are redundant, fix variables and
    find relationships such as x + y == 1.

    Presolve will also be done.

    If row names existed they may be replaced by R0000000 etc

*/

class CglPreProcess  {
  
public:

  ///@name Main methods 
  //@{
  /** preProcess problem - returning new problem.
      If makeEquality true then <= cliques converted to ==.
      Presolve will be done numberPasses times.

      Returns NULL if infeasible

      This version uses default strategy.  For more control copy and edit
      code from this function i.e. call preProcessNonDefault
  */
  OsiSolverInterface * preProcess(OsiSolverInterface & model, 
                                  bool makeEquality=false, int numberPasses=5);
  /** preProcess problem - returning new problem.
      If makeEquality true then <= cliques converted to ==.
      Presolve will be done numberPasses times.

      Returns NULL if infeasible

      This version assumes user has added cut generators to CglPreProcess object
      before calling it.  As an example use coding in preProcess
      If makeEquality is 1 add slacks to get cliques,
      if 2 add slacks to get sos (but only if looks plausible) and keep sos info
  */
  OsiSolverInterface * preProcessNonDefault(OsiSolverInterface & model, 
                                  int makeEquality=0, int numberPasses=5,
					    int tuning=0);
  /// Creates solution in original model
  void postProcess(OsiSolverInterface &model);
  /** Tightens primal bounds to make dual and branch and cutfaster.  Unless
      fixed or integral, bounds are slightly looser than they could be.
      Returns non-zero if problem infeasible
      Fudge for branch and bound - put bounds on columns of factor *
      largest value (at continuous) - should improve stability
      in branch and bound on infeasible branches (0.0 is off)
  */
  int tightenPrimalBounds(OsiSolverInterface & model,double factor=0.0);
  /** Fix some of problem - returning new problem.
      Uses reduced costs.
      Optional signed character array
      1 always keep, -1 always discard, 0 use djs

  */
  OsiSolverInterface * someFixed(OsiSolverInterface & model, 
                                 double fractionToKeep=0.25,
                                 bool fixContinuousAsWell=false,
                                 char * keep=NULL) const;
  /// If we have a cutoff - fix variables
  int reducedCostFix(OsiSolverInterface & model);
  //@}

  //---------------------------------------------------------------------------

  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false if the value of the parameter is out of range.

     The get methods return the value of the parameter.

  */
  //@{
  /** Set cutoff bound on the objective function.

    When using strict comparison, the bound is adjusted by a tolerance to
    avoid accidentally cutting off the optimal solution.
  */
  void setCutoff(double value) ;

  /// Get the cutoff bound on the objective function - always as minimize
  double getCutoff() const;
  /// The original solver associated with this model.
  inline OsiSolverInterface * originalModel() const
  { return originalModel_;}
  /// Solver after making clique equalities (may == original)
  inline OsiSolverInterface * startModel() const
  { return startModel_;}
  /// Copies of solver at various stages after presolve
  inline OsiSolverInterface * modelAtPass(int iPass) const
  { if (iPass>=0&&iPass<numberSolvers_) return model_[iPass]; else return NULL;}
  /// Copies of solver at various stages after presolve after modifications
  inline OsiSolverInterface * modifiedModel(int iPass) const
  { if (iPass>=0&&iPass<numberSolvers_) return modifiedModel_[iPass]; else return NULL;}
  /// Matching presolve information
  inline OsiPresolve * presolve(int iPass) const
  { if (iPass>=0&&iPass<numberSolvers_) return presolve_[iPass]; else return NULL;}
  /** Return a pointer to the original columns (with possible  clique slacks)
      MUST be called before postProcess otherwise you just get 0,1,2.. */
  const int * originalColumns() const;
  /** Return a pointer to the original rows
      MUST be called before postProcess otherwise you just get 0,1,2.. */
  const int * originalRows() const;
  /// Number of SOS if found
  inline int numberSOS() const
  { return numberSOS_;}
  /// Type of each SOS
  inline const int * typeSOS() const
  { return typeSOS_;}
  /// Start of each SOS
  inline const int * startSOS() const
  { return startSOS_;}
  /// Columns in SOS
  inline const int * whichSOS() const
  { return whichSOS_;}
  /// Weights for each SOS column
  inline const double * weightSOS() const
  { return weightSOS_;}
  /// Pass in prohibited columns 
  void passInProhibited(const char * prohibited,int numberColumns);
  /// Updated prohibited columns
  inline const char * prohibited()
  { return prohibited_;}
  //@}

  ///@name Cut generator methods 
  //@{
  /// Get the number of cut generators
  inline int numberCutGenerators() const
  { return numberCutGenerators_;}
  /// Get the list of cut generators
  inline CglCutGenerator ** cutGenerators() const
  { return generator_;}
  ///Get the specified cut generator
  inline CglCutGenerator * cutGenerator(int i) const
  { return generator_[i];}
  /** Add one generator - up to user to delete generators.
  */
  void addCutGenerator(CglCutGenerator * generator);
//@}
    
  /**@name Setting/Accessing application data */
  //@{
    /** Set application data.

	This is a pointer that the application can store into and
	retrieve.
	This field is available for the application to optionally
	define and use.
    */
    void setApplicationData (void * appData);

    /// Get application data
    void * getApplicationData() const;
  //@}
  
  //---------------------------------------------------------------------------

  /**@name Message handling */
  //@{
  /// Pass in Message handler (not deleted at end)
  void passInMessageHandler(CoinMessageHandler * handler);
  /// Set language
  void newLanguage(CoinMessages::Language language);
  inline void setLanguage(CoinMessages::Language language)
  {newLanguage(language);}
  /// Return handler
  inline CoinMessageHandler * messageHandler() const
  {return handler_;}
  /// Return messages
  inline CoinMessages messages() 
  {return messages_;}
  /// Return pointer to messages
  inline CoinMessages * messagesPointer() 
  {return &messages_;}
  //@}
  //---------------------------------------------------------------------------


  ///@name Constructors and destructors etc
  //@{
  /// Constructor
  CglPreProcess(); 
  
  /// Copy constructor .
  CglPreProcess(const CglPreProcess & rhs);
  
  /// Assignment operator 
  CglPreProcess & operator=(const CglPreProcess& rhs);

  /// Destructor 
  ~CglPreProcess ();
  
  /// Clears out as much as possible
  void gutsOfDestructor();
  //@}
private:

  ///@name private methods
  //@{
  /** Return model with useful modifications.  
      If constraints true then adds any x+y=1 or x-y=0 constraints
      If NULL infeasible
  */
  OsiSolverInterface * modified(OsiSolverInterface * model,
                                bool constraints,
                                int & numberChanges,
                                int iBigPass,
				int numberPasses);
  /// create original columns and rows
  void createOriginalIndices() const;
  /// Make continuous variables integer
  void makeInteger();
  //@}

//---------------------------------------------------------------------------

private:
  ///@name Private member data 
  //@{

  /// The original solver associated with this model.
  OsiSolverInterface * originalModel_;
  /// Solver after making clique equalities (may == original)
  OsiSolverInterface * startModel_;
  /// Number of solvers at various stages
  int numberSolvers_;
  /// Copies of solver at various stages after presolve
  OsiSolverInterface ** model_;
  /// Copies of solver at various stages after presolve after modifications
  OsiSolverInterface ** modifiedModel_;
  /// Matching presolve information
  OsiPresolve ** presolve_;

   /// Message handler
  CoinMessageHandler * handler_;

  /** Flag to say if handler_ is the default handler.
  
    The default handler is deleted when the model is deleted. Other
    handlers (supplied by the client) will not be deleted.
  */
  bool defaultHandler_;

  /// Cgl messages
  CoinMessages messages_;

  /// Pointer to user-defined data structure
  void * appData_;
  /// Original column numbers
  mutable int * originalColumn_;
  /// Original row numbers
  mutable int * originalRow_;
  /// Number of cut generators
  int numberCutGenerators_;
  /// Cut generators
  CglCutGenerator ** generator_;
  /// Number of SOS if found
  int numberSOS_;
  /// Type of each SOS
  int * typeSOS_;
  /// Start of each SOS
  int * startSOS_;
  /// Columns in SOS
  int * whichSOS_;
  /// Weights for each SOS column
  double * weightSOS_;
  /// Number of columns in original prohibition set
  int numberProhibited_;
  /// Columns which should not be presolved e.g. SOS
  char * prohibited_;
 //@}
};

#endif
