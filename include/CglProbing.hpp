// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CglProbing_H
#define CglProbing_H

#include <string>

#include "CglCutGenerator.hpp"

/** Probing Cut Generator Class */
class CglProbing : public CglCutGenerator {
   friend void CglProbingUnitTest(const OsiSolverInterface * siP,
				  const std::string mpdDir );
 
public:
    
  
  /**@name Generate Cuts */
  //@{
  /** Generate probing/disaggregation cuts for the model of the 
      solver interface, si.

      Insert the generated cuts into OsiCut, cs.

      If a "snapshot" of a matrix exists then this will be used.
      Presumably this will give global cuts and will be faster.
      No check is done to see if cuts will be global.

      Otherwise use current matrix.

      Both row cuts and column cuts may be returned

      The mode options are:
      0) Only unsatisfied integer variables will be looked at.
         If no information exists for that variable then
	 probing will be done so as a by-product you "may" get a fixing
	 or infeasibility.  This will be fast and is only available
         if a snapshot exists (otherwise as 1).  
	 The bounds in the snapshot are the ones used.
      1) Look at unsatisfied integer variables, using current bounds.
         Probing will be done on all looked at.
      2) Look at all integer variables, using current bounds.
         Probing will be done on all

  */
  virtual void generateCuts( const OsiSolverInterface & si, 
			     OsiCuts & cs) const;
  //@}

  /**@name snapshot */
  //@{
  /// Create a copy of matrix which is to be used
  /// this is to speed up process and to give global cuts
  /// Can give an array with 1 set to select, 0 to ignore
  /// column bounds are tightened
  /// If array given then values of 1 will be set to 0 if redundant
  void snapshot ( const OsiSolverInterface & si,
		  char * possible=NULL);
  //@}

  /**@name deleteSnapshot */
  //@{
  /// Deletes snapshot
  void deleteSnapshot ( );
  //@}

  /**@name Get tighter column bounds */
  //@{
  /// Lower
  const double * tightLower() const;
  /// Upper
  const double * tightUpper() const;
  //@}

  /**@name Change mode */
  //@{
  /// Set
  void setMode(int mode);
  /// Get
  int getMode() const;
  //@}

  /**@name Change maxima */
  //@{
  /// Set maximum number of passes per node
  void setMaxPass(int value);
  /// Get maximum number of passes per node
  int getMaxPass() const;
  /// Set maximum number of unsatisfied variables to look at
  void setMaxProbe(int value);
  /// Get maximum number of unsatisfied variables to look at
  int getMaxProbe() const;
  /// Set maximum number of variables to look at in one probe
  void setMaxLook(int value);
  /// Get maximum number of variables to look at in one probe
  int getMaxLook() const;
  //@}

  /**@name Stop or restart row cuts (otherwise just fixing from probing) */
  //@{
  /// Set
  /// 0 no cuts, 1 just disaggregation type, 2 coefficient ( 3 both)
  void setRowCuts(int type);
  /// Get
  int rowCuts() const;
  //@}

  /**@name Whether use objective as constraint */
  //@{
  /// Set
  void setUsingObjective(bool yesNo);
  /// Get
  int getUsingObjective() const;
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglProbing ();
 
  /// Copy constructor 
  CglProbing (
    const CglProbing &);

  /// Assignment operator 
  CglProbing &
    operator=(
    const CglProbing& rhs);
  
  /// Destructor 
  virtual
    ~CglProbing ();
  //@}
      
private:
  
 // Private member methods
  /**@name probe */
  //@{
  /// Does probing and adding cuts
  int probe( const OsiSolverInterface & si, 
	     const OsiRowCutDebugger * debugger, 
	     OsiCuts & cs, 
	     double * colLower, double * colUpper, OsiPackedMatrix *rowCopy,
	     double * rowLower, double * rowUpper,
	     char * intVar, double * minR, double * maxR, int * markR, 
	     double * movement, int * look, int nlook) const;
  //@}

  // Private member data

  /**@name Private member data */
  //@{
  /// Row copy
  OsiPackedMatrix * rowCopy_;
  /// Lower bounds on rows
  double * rowLower_;
  /// Upper bounds on rows
  double * rowUpper_;
  /// Lower bounds on columns
  double * colLower_;
  /// Upper bounds on columns
  double * colUpper_;
  /// Number of rows in snapshot (0 if no snapshot)
  int numberRows_;
  /// Number of columns in snapshot (must == current)
  int numberColumns_;
  /// Tolerance to see if infeasible
  double primalTolerance_;
  /// Mode - 0 lazy using snapshot, 1 just unsatisfied, 2 all
  int mode_;
  /// Row cuts flag
  /// 0 no cuts, 1 just disaggregation type, 2 coefficient ( 3 both)
  int rowCuts_;
  /// Maximum number of passes to do in probing
  int maxPass_;
  /// Maximum number of unsatisfied variables to probe
  int maxProbe_;
  /// Maximum number of variables to look at in one probe
  int maxStack_;
  /// Whether to include objective as constraint
  bool usingObjective_;
  /// Number of integer variables
  int numberIntegers_;
  /// Disaggregation cuts
  typedef struct{
    int sequence; // integer variable
    // newValue will be NULL if no probing done yet
    // lastLBWhenAt1 gives length of newValue;
    int lastUBWhenAt0; // last UB changed when integer at lb 
    int lastLBWhenAt0; // last LB changed when integer at lb
    int lastUBWhenAt1; // last UB changed when integer at ub
    int lastLBWhenAt1; // last LB changed when integer at ub
    int * index; // columns whose bounds will be changed
    double * newValue; // new values 
  } disaggregation;
  disaggregation * cutVector_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglProbing class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void CglProbingUnitTest(const OsiSolverInterface * siP,
			const std::string mpdDir );
  
#endif
