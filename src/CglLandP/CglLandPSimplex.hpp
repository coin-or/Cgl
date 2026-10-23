// Copyright (C) 2005-2009 Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     21/07/05
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------
#ifndef CglLandPSimplex_H
#define CglLandPSimplex_H

#include <iostream>
#include <vector>

#include "CglConfig.h"
#include "CglLandP.hpp"

#include "OsiSolverInterface.hpp"
#include "CoinMessage.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinPackedMatrix.hpp"

#ifdef CGL_HAS_OSICLP
#include "OsiClpSolverInterface.hpp"
#endif

#include "CglLandPTabRow.hpp"
#include "CglLandPUtils.hpp"
#include "CglLandPMessages.hpp"

//#define APPEND_ROW
#define OLD_COMPUTATION

namespace LAP
{
/** Forward declaration of class to store extra debug data.*/
class DebugData;

class CglLandPSimplex
{
public:
    /** Usefull onstructor */
    CglLandPSimplex(const OsiSolverInterface &si,
                    const CglLandP::CachedData &cached,
                    const CglLandP::Parameters &params,
                    Validator &validator);
    /** Destructor */
    ~CglLandPSimplex();
    /**Update cached information in case of basis change in a round*/
    void cacheUpdate(const CglLandP::CachedData &cached, bool reducedSpace = 0);
    /** reset the solver to optimal basis */
    bool resetSolver(const CoinWarmStartBasis * basis);
    /** Perfom pivots to find the best cuts */
    bool optimize(int var, OsiRowCut & cut, const CglLandP::CachedData &cached, const CglLandP::Parameters & params);
    /** Find Gomory cut (i.e. don't do extra setup required for pivots).*/
    bool generateMig(int row, OsiRowCut &cut, const CglLandP::Parameters & params);

    /** Find extra constraints in current tableau.*/
    int generateExtraCuts(const CglLandP::CachedData &cached, const CglLandP::Parameters & params);
    /** Generate a constrainte for a row of the tableau different from the source row.*/
    int generateExtraCut(int i, const CglLandP::CachedData & cached,
                         const CglLandP::Parameters& params);

    void genThisBasisMigs(const CglLandP::CachedData &cached,
                          const CglLandP::Parameters & params) ;

    /** insert all extra cuts in cs.*/
    int insertAllExtr(OsiCuts & cs, CoinRelFltEq eq);

    void setLogLevel(int level)
    {
        handler_->setLogLevel(level);
    }


    /** Hand the simplex a solver it should take ownership of.
        Kept for API compatibility; generateCuts() no longer uses it, because
        optimize() installs its own clone of cached.solver_ and threw this one
        away unread. */
    void setSi(OsiSolverInterface *si)
    {
        si_ = si;
        ownSi_ = true;
#ifdef CGL_HAS_OSICLP
#ifndef CBC_OTHER_SOLVER
        clp_ = getClpSolver(si);
#else
        OsiClpSolverInterface * clpSi = getClpSolver(si_);
        if (clpSi)
        {
            clp_ = clpSi;
        }
#endif
#endif
    }
    void freeSi()
    {
        /* No assert on si_ != NULL: the second freeSi() of the retry path in
           CglLandP::generateCuts reaches here after the first already nulled it,
           and a borrowed si_ must be dropped, not deleted. */
        if (ownSi_)
            delete si_;
        si_ = NULL;
        ownSi_ = false;
#ifdef CGL_HAS_OSICLP
        clp_ = NULL;
#endif
    }

    Cuts& extraCuts()
    {
        return cuts_;
    }
    void loadBasis(const OsiSolverInterface &si,
                   std::vector<int> &M1, std::vector<int> &M2,
                   int k);


    int getNumCols() const
    {
        return ncols_;
    }

    int getNumRows() const
    {
        return nrows_;
    }

    const CoinWarmStartBasis  * getBasis() const
    {
        return basis_;
    }
    const int * getNonBasics() const
    {
        return nonBasics_;
    }

    const int * getBasics() const
    {
        return basics_;
    }

    void outPivInfo(int ncuts)
    {
        handler_->message(RoundStats, messages_)<<ncuts<<numPivots_
        <<numSourceRowEntered_
        <<numIncreased_<<CoinMessageEol;
    }
#ifdef APPEND_ROW
    /** Append source row to tableau.*/
    void append_row(int row_num, bool modularize) ;
    /** Update appended row after a pivot.*/
    void update_row(TabRow &row);

    void check_mod_row(TabRow &row);
#endif
protected:
    /** Perform a change in the basis (direction is 1 if leaving variable is going to ub, 0 otherwise)*/
    bool changeBasis(int incoming, int leaving, int direction, 
#ifndef OLD_COMPUTATION
                     bool recompute_source_row,
#endif
                     bool modularize);
    /** Find a row which can be used to perform an improving pivot the fast way
      * (i.e., find the leaving variable).
      \return index of the row on which to pivot or -1 if none exists. */
    int fastFindCutImprovingPivotRow( int &direction, int &gammaSign, double tolerance, bool flagPositiveRows);
    /** Rescan reduced costs tables */
    int rescanReducedCosts( int &direction, int &gammaSign, double tolerance);
    /** Find the column which leads to the best cut (i.e., find incoming variable).*/
    int fastFindBestPivotColumn(int direction, int gammaSign,
                                double pivotTol, double rhsTol,
                                bool reducedSpace,
                                bool allowNonImproving,
                                double &bestSigma, bool modularize);

    /** Find incoming and leaving variables which lead to the most violated
      adjacent normalized lift-and-project cut.
      \remark At this point reduced costs should be already computed.
      \return incoming variable variable,
      \param leaving variable
      \param direction leaving direction
      */
    int findBestPivot(int &leaving, int & direction,
                      const CglLandP::Parameters & params);


    /** Compute the objective value of the Cglp for given row and rhs (if strengthening shall be applied
      row should have been modularized).*/
    double computeCglpObjective(const TabRow & row, bool modularize = false) const;
    /** return the coefficients of the strengthened intersection cut
      * takes one extra argument seens needs to consider variable type.
      */
    inline double strengthenedIntersectionCutCoef(int i, double alpha_i, double beta) const;
    /** return the coefficient of the new row (combining row_k + gamma row_i).
      */
    inline double newRowCoefficient(int j, double gamma) const;
    /** Create the intersection cut of row k*/
    void createIntersectionCut(TabRow & row, OsiRowCut &cut) const;
    /** Compute the normalization factor of the cut.*/
    double normalizationFactor(const TabRow & row) const;
    /** Scale the cut by factor.*/
    void scaleCut(OsiRowCut & cut, double factor) const;
    /** Create strenghtened row */
    //  void createIntersectionCut(double * row);
    /** Create MIG cut from row k*/
    void createMIG( TabRow & row, OsiRowCut &cut) const;
    /** Get the row i of the tableau */
    void pullTableauRow(TabRow & row) const;
    /** pullTableauRow, but reading a solver other than si_.
        The only members pullTableauRow reads of the solver are si_ and clp_, so
        pointing those at \p src for the duration is the whole of it. Used by
        optimize() to test the row-length gate against cached.solver_ before it
        commits to a clone. */
    void pullTableauRowFrom(OsiSolverInterface * src, TabRow & row);
    /** Adjust the row of the tableau to reflect leaving variable direction */
    void adjustTableauRow(int var, TabRow & row, int direction);
    /** reset the tableau row after a call to adjustTableauRow */
    void resetOriginalTableauRow(int var, TabRow & row, int direction);
    /**Get lower bound for variable or constraint */
    inline double getLoBound(int index) const
    {
        return lo_bounds_[original_index_[index]];
    }
    /**Get upper bound for variable or constraint */
    inline double getUpBound(int index) const
    {
        return up_bounds_[original_index_[index]];
    }
    /** Access to value in solution to cut (indexed in reduced problem) */
    inline double getColsolToCut(int index) const
    {
        return colsolToCut_[original_index_[index]];
    }
    bool isGtConst(int index) const
    {
        return (index >= ncols_ && lo_bounds_[original_index_[index]] < -1e-10 && up_bounds_[original_index_[index]] <= 1e-09);
    }
    /** Access to value in solution to cut (indexed in reduced problem) */
    inline void setColsolToCut(int index, double value)
    {
        colsolToCut_[original_index_[index]] = value;
    }
    /** Get the basic status of a variable (structural or slack).*/
    inline CoinWarmStartBasis::Status getStatus(int index) const
    {
        if (index < ncols_) return basis_->getStructStatus(index);
        return basis_->getArtifStatus(index - ncols_);
    }
    /** Say if variable index by i in current tableau is integer.*/
    inline bool isInteger(int index) const
    {
        return integers_[original_index_[index]];
    }
    /** Compute normalization weights.*/
    void computeWeights(CglLandP::LHSnorm norm, CglLandP::Normalization type,
                        CglLandP::RhsWeightType rhs);
    /** Evenutaly multiply a by w if normed_weights_ is not empty.*/
    double normedCoef(double a, int ii) const
    {
        if (norm_weights_.empty())
        {
            return a;
        }
        else
        {
            return a*norm_weights_[ii];
        }
    }
    /** print the tableau of current basis. */
    void printTableau(std::ostream & os);

  /** Print everything .*/
  void printEverything();
    /** print the tableau of current basis. */
    void printTableauLateX(std::ostream & os);
    void printRowLateX(std::ostream & os, int i);
    void printCutLateX(std::ostream & os, int i);

    /** Print CGLP basis corresponding to current tableau and source row.*/
    void printCglpBasis(std::ostream& os = std::cout);
    /** Put variables in M1 M2 and M3 according to their sign.*/
    void get_M1_M2_M3(const TabRow & row,
                      std::vector<int> &M1,
                      std::vector<int> &M2,
                      std::vector<int> &M3);
    /** Put a vector in structural sapce.*/
    void eliminate_slacks(double * vec) const;
private:
    /// No default constructor
    CglLandPSimplex();
    /// No copy constructor
    CglLandPSimplex(const CglLandPSimplex&);
    /// No assignment operator
    CglLandPSimplex& operator=(const CglLandPSimplex&);
#ifdef CGL_HAS_OSICLP
    /** Pointer to OsiClpSolverInterface if used.*/
    OsiClpSolverInterface * clp_;
#endif

    /** Update values in M1 M2 and M3 before an iteration.*/
    void updateM1_M2_M3(TabRow & row, double tolerance, bool alwaysComputeCheap);
    /** Remove rows from current tableau.*/
    void removeRows(int nDelete, const int * rowsIdx);


    void compute_p_q_r_s(double gamma, int gammaSign, double &p, double & q, double & r , double &s);


    /// @name Work infos
    /// @{
    /** Source row for cut */
    TabRow row_k_;
    /** Original version of source row (without modularization).*/
    TabRow original_row_k_;
    /** Row of leaving candidate*/
    TabRow row_i_;
#ifndef NDEBUG
    TabRow new_row_;
#endif
    /**vector to sort the gammas*/
    CoinPackedVector gammas_;
    /**first work vector in row space.*/
    std::vector<double> rWk1_;
    /**scond work vector in row space.*/
    std::vector<double> rWk2_;
    /**third work vector in row space.*/
    std::vector<double> rWk3_;
    /**fourth work vector in row space.*/
    std::vector<double> rWk4_;
    /** integer valued work vector on the rows */
    std::vector<int> rIntWork_;
    /** Flag rows which we don't want to try anymore */
    bool * rowFlags_;
    /** The p and r of the CGLP reduced cost gammaSign * (q*r - p*s) / r.  Both
        depend only on the source row, the point to cut and the subspace, none of
        which changes while a pivot is being chosen, so they are constants of the
        pivot iteration rather than of the candidate.  Filled by
        computeInvariantPR, read by exactRowReducedCosts. */
    double pInvariant_;
    double rInvariant_;
    /** The fractional part of the source row's basic variable, also invariant. */
    double fzeroInvariant_;
    /** Numerator q*r - p*s that exactRowReducedCosts predicted for the candidate
        it last returned, and the slot it returned, or -1 for none.  Declared
        unconditionally so that the profile and non-profile builds agree on
        sizeof(*this); only a profile build reads them, to check the prediction
        against the value fastFindBestPivotColumn recomputes for itself. */
    double predictedRcNum_;
    int predictedSlot_;
    /** Flag columns which are in the subspace (usualy remove nonbasic structurals in subspace) */
    std::vector<bool> col_in_subspace;
    /** Flag columns which have to be considered for leaving the basis */
    bool *colCandidateToLeave_;
    /** Store the basics variable */
    int * basics_;
    /** Stores the nonBasicVariables */
    int * nonBasics_;
    /** Stores the variables which are always in M1 for a given k*/
    std::vector<int> M1_;
    /** Stores the variables which are always in M2 for a given k*/
    std::vector<int> M2_;
    /** Stores the variables which could be either in M1 or M2 */
    std::vector<int> M3_;
    /** stores the cglp value of the normalized cut obtained from row k_ */
    double sigma_;
    /** Keep track of basis status */
    CoinWarmStartBasis * basis_;
    /** Pointer to the solution to cut (need to be modified after each pivot because we are only considering slacks).*/
    double * colsolToCut_;
    /** Pointer to the current basic solution.*/
    double * colsol_;
    /// cached numcols in original problem
    int ncols_orig_;
    ///cached numrows in original problem
    int nrows_orig_;
    /// cached number of columns in reduced size problem
    int ncols_;
    /// Cached number of rows in reduced size problem
    int nrows_;
    // for fast access to lower bounds (both cols and rows)
    std::vector<double> lo_bounds_;
    // for fast access to upper bounds (both cols and rows)
    std::vector<double> up_bounds_;
    /// Say if we are in a sequence of degenerate pivots
    bool inDegenerateSequence_;
    /// Value for the reduced cost chosen for pivoting
    double chosenReducedCostVal_;
    /// pointer to array of integer info for both structural and slacks
    const bool * integers_;
    /// Original index of variable before deletions.
    std::vector<int> original_index_;
    /// Stores extra cuts which are generated along the procedure
    Cuts cuts_;
    /// @}
    /// @name Interfaces to the solver
    /// @{
    /** Pointer to the solver interface */
    OsiSolverInterface * si_;
    /** Whether si_ is ours to delete.
        The constructor borrows the caller's solver, setSi() and optimize() both
        install a clone, and optimize() deletes whatever it finds before cloning.
        Without this flag the first of those deletes would free the *caller's*
        solver, so this is what lets generateCuts() stop cloning defensively. */
    bool ownSi_;
    ///@}
    /// Own the data or not?
    bool own_;
    /// A pointer to a cut validator
    Validator & validator_;
    /// Weights for the normalization constraint
    std::vector<double> norm_weights_;
    /// Weight for rhs of normalization constraint.*/
    double rhs_weight_;

    /// number of rows with a <0 rc in current iteration
    int nNegativeRcRows_;
    /** Check that the basis is correct.*/
    bool checkBasis();

    /** Record the number of pivots.*/
    int numPivots_;
    /** Record the number of times the source row entered the basis.*/
    int numSourceRowEntered_;
    /** Record the number of times that sigma increased.*/
    int numIncreased_;

    /** Message handler. */
    CoinMessageHandler * handler_;
    /** Messages. */
    CoinMessages messages_;
#ifndef NDEBUG
    double bestSigma_;
#endif

protected:
    /** \name Slow versions of the function (old versions do not work).*/
    /** @{ */
    /** Compute the reduced cost of Cglp */
    double computeCglpRedCost(int direction, int gammaSign, double tau);
    /** Compute the value of sigma and thau (which are constants for a row i as defined in Mike Perregaard thesis */
    double computeRedCostConstantsInRow();
    /** Compute the objective value of the Cglp with linear combintation of the two rows by gamma */
    double computeCglpObjective(double gamma, bool strengthen, TabRow &row);
    /** Compute the objective value of the Cglp with linear combintation of the row_k_ and gamma row_i_ */
    double computeCglpObjective(double gamma, bool strengthen);
    /** Find a row which can be used to perform an improving pivot return index of the cut or -1 if none exists
      * (i.e., find the leaving variable).*/
    int findCutImprovingPivotRow( int &direction, int &gammaSign, double tolerance);
    /** Index of a (direction, gammaSign) candidate in the ul/vl/uu/vu order the
        four reduced-cost tables are read in, i.e. rWk1_, rWk3_, rWk2_, rWk4_. */
    static int rcSlot(int direction, int gammaSign)
    {
        return (direction > 0 ? 2 : 0) + (gammaSign > 0 ? 1 : 0);
    }
    /** The reduced-cost table holding candidate slot \p s. */
    double * rcSlotTable(int s)
    {
        double * const t[4] = {&rWk1_[0], &rWk3_[0], &rWk2_[0], &rWk4_[0]};
        return t[s];
    }
    /** Recompute pInvariant_, rInvariant_ and fzeroInvariant_ for row_k_. */
    void computeInvariantPR(bool reducedSpace);
    /** Exact reduced cost of every candidate of the row currently in row_i_.

        Considers every candidate whose direction is available, not just the ones
        the tabulated costs believed in, and retires -- writes the 10. "not
        improving" sentinel into -- each one it proves nonimproving.  \p skipMask
        excludes candidates already tried, so a caller looping on this call always
        makes progress.  \return the slot of the best genuinely improving candidate,
        or -1 if the row has none; \p bestRc receives its reduced cost. */
    int exactRowReducedCosts(bool reducedSpace, int skipMask, double &bestRc);
    /** Find the column which leads to the best cut (i.e., find incoming variable).*/
    int findBestPivotColumn(int direction,
                            double pivotTol, bool reducedSpace, bool allowDegeneratePivot,
                            bool modularize);
#if 1
    int plotCGLPobj(int direction, double gammaTolerance,
                    double pivotTol, bool reducedSpace, bool allowDegenerate, bool modularize);
#endif

    /** @} */
};


/** return the coefficients of the strengthened intersection cut */
double CglLandPSimplex::strengthenedIntersectionCutCoef(int i, double alpha_i, double beta) const
{
    //  double ratio = beta/(1-beta);
    if ( (!integers_[i]))
        return intersectionCutCoef(alpha_i, beta);
    else
    {
        double f_i = alpha_i - floor(alpha_i);
        if (f_i < beta)
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
    return row_k_[j] + gamma * row_i_[j];
}

}
#endif
