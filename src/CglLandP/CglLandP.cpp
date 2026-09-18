// Copyright (C) 2005-2009, Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     07/21/05
//
//
// This code is licensed under the terms of the Eclipse Public License (EPL).
//---------------------------------------------------------------------------
#include "CglLandP.hpp"
#include "CglLandPSimplex.hpp"
#include "OsiRowCutDebugger.hpp"

#define INT_INFEAS(value) fabs(value - floor(value+0.5))

#include "CglConfig.h"

#ifdef CGL_HAS_OSICLP
#include "OsiClpSolverInterface.hpp"
#endif

#define CLONE_SI //Solver is cloned between two cuts

#include "CoinTime.hpp"
#include "CglGomory.hpp"
#include "CoinFactorization.hpp"
#include <fstream>

#ifdef CGL_LANDP_PROFILE
/*
 * Per-stage attribution for landp-bench, compiled in only under
 * -DCGL_LANDP_PROFILE.
 *
 * It exists because the obvious question about this generator -- how much of the
 * time is the per-candidate solver cloning, as against the pivoting it exists to
 * enable -- cannot be answered by the control flags alone. --pivot-limit=0 removes
 * the clones AND the pivots together, so it gives their sum and not the split.
 *
 * Two levels: the clones and the candidate loop here, and the pivot loop's stages
 * further down (see lpProfStage). The second level was left out originally on the
 * grounds that a clock call inside the pivot loop would distort what it measures,
 * and that was right at the time -- see the note there for why it no longer is.
 *
 * THE ACCUMULATORS ARE FILE-STATIC, so a profile run must be serial -- which is
 * already the rule for any run whose times get quoted. And the counters are the
 * part to trust: a clone count is exact, while every timer here pays the clock's
 * own overhead.
 */
/* Kept, and now expected to read 0: the two clones in generateCuts that these two
   counted were dead -- optimize() deleted them unread -- and have been removed.
   A nonzero value here means a clone has come back. */
static double lpProfCloneTime = 0.0;
static double lpProfLoopTime = 0.0;    /* the whole candidate loop */
static double lpProfSetupTime = 0.0;   /* getData + simplex ctor + candidate sort */
static double lpProfRangeTime = 0.0;   /* the range-row clone-and-augment path */
static double lpProfOptCloneTime = 0.0; /* the clone inside optimize() */
static int lpProfClones = 0;           /* see lpProfCloneTime: now always 0 */
static int lpProfOptClones = 0;        /* clones made inside optimize() */
static int lpProfCandidates = 0;
static int lpProfRetries = 0;          /* candidates that took the :750 retry */

/* Called from CglLandPSimplex::optimize, which is where the clones that the
   algorithm actually pivots on are made. Split out from lpProfCloneTime on
   purpose: those two are the removable ones and these two are not, so lumping
   them would hide exactly the number this profile exists to produce. */
void cglLandPProfileAddOptClone(double secs)
{
  lpProfOptCloneTime += secs;
  ++lpProfOptClones;
}

/* optimize()'s six exits, in the order they appear in CglLandPSimplex.cpp:
   0 entry length gate (BEFORE any pivot -- its clone was pure waste)
   1 cglp objective went nonnegative after a pivot
   2 pivot failed and sigma increased
   3 pivot failed and sigma unchanged
   4 exit length gate (after pivoting)
   5 a cut was built. */
static int lpProfOptExit[6] = {0, 0, 0, 0, 0, 0};
static int lpProfPivots = 0;

void cglLandPProfileOptExit(int which, int numPivots)
{
  if (which >= 0 && which < 6) ++lpProfOptExit[which];
  lpProfPivots += numPivots;
}

/* si_->pivot()'s return code is NOT a boolean. Clp documents (and
   ClpSimplex.cpp:9181 implements) 0 = okay, 1 = "inaccuracy forced
   re-factorization", which is a pivot that WAS taken, and -1 = "would be
   singular", the only code under which nothing happened. changeBasis treats
   every nonzero alike, so splitting them decides between two very different
   diagnoses: "the algorithm cannot pivot from here" and "the algorithm threw
   away a pivot it had already made". */
static int lpProfPivCode[3] = {0, 0, 0}; /* ok / refactorized(+) / singular(-) */

void cglLandPProfilePivotCode(int code)
{
  if (code == 0) ++lpProfPivCode[0];
  else if (code > 0) ++lpProfPivCode[1];
  else ++lpProfPivCode[2];
}

/* Per-stage timers INSIDE the pivot loop. The comment at the top of this block
   used to say nesting this deep would distort what it measures, and that was true
   of the code as it stood: LandP completed zero pivots, so the loop body ran once
   per candidate and a clock call was a visible fraction of it. With pivoting
   repaired the loop body is ~1.4 ms and the eight CoinGetTimeOfDay calls below are
   ~200 ns of it, so the distortion is ~0.01% and the attribution is worth having.

   The stages are the loop body in source order. LP_STAGE_ROW is the one to watch:
   fastFindCutImprovingPivotRow walks every column of M1_ and M2_ through the
   column-major matrix, i.e. all of nnz(A), on every single pivot. */
enum { LP_ST_UPD = 0,   /* updateM1_M2_M3 + computeCglpObjective */
       LP_ST_ROW,       /* fastFindCutImprovingPivotRow */
       LP_ST_COL,       /* fastFindBestPivotColumn */
       LP_ST_RESCAN,    /* rescanReducedCosts retry loop */
       LP_ST_PIVOT,     /* si_->pivot() alone */
       LP_ST_CBREST,    /* rest of changeBasis: bookkeeping + row_k_ update */
       LP_ST_FACT,      /* the every-40-pivots explicit factorize() */
       LP_ST_EXACT,     /* exactRowReducedCosts */
       LP_ST_N };
static double lpProfStage[LP_ST_N] = {0., 0., 0., 0., 0., 0., 0., 0.};
/* Call counts, free to collect and worth more than the times on their own: with
   them, "st_rescan is 38% of the loop" separates into "each call is slow" and
   "it is called many times per pivot", which point at different fixes. */
static int lpProfStageN[LP_ST_N] = {0, 0, 0, 0, 0, 0, 0, 0};
static const char *const lpProfStageName[LP_ST_N] =
  {"st_upd", "st_row", "st_col", "st_rescan", "st_pivot", "st_cbrest", "st_fact",
   "st_exact"};

void cglLandPProfileAddStage(int which, double secs)
{
  if (which >= 0 && which < LP_ST_N)
    {
      lpProfStage[which] += secs;
      ++lpProfStageN[which];
    }
}

/* Why fastFindBestPivotColumn returned no column. Order matches the LP_COL_*
   macros in CglLandPSimplex.cpp. Only nopivtol and degen return -1 and so charge
   the maxTryRow budget; the other three return -2 and are free, which is what
   lets the retry loop run until it runs out of unflagged rows. */
enum { LP_CF_N = 5 };
static int lpProfColFail[LP_CF_N] = {0, 0, 0, 0, 0};
static const char *const lpProfColFailName[LP_CF_N] =
  {"cf_nogamma", "cf_mistakenrc", "cf_nopivtol", "cf_tinypiv", "cf_degen"};

void cglLandPProfileColFail(int which)
{
  if (which >= 0 && which < LP_CF_N) ++lpProfColFail[which];
}

/* Bins the exact reduced cost of a row the table had promised was improving.
   "tiny" is |exact| below the pivot tolerance CglLandP itself uses (1e-4 by
   default, but the bins are absolute so they can be read against any of them);
   "sign" means the exact value is positive by more than that, i.e. the table did
   not merely lose precision, it pointed the search at the wrong row. */
static int lpProfMistakenTiny = 0;
static int lpProfMistakenSign = 0;
static double lpProfMistakenWorst = 0.0;   /* largest positive exact rc seen */
static double lpProfMistakenTableSum = 0.0;/* sum of the table's promises */
static double lpProfMistakenExactSum = 0.0;/* sum of the exact recomputations */

/* Which of the four reduced-cost tables the row came from. rescanReducedCosts
   reads ul_i, vl_i, uu_i, vu_i and each maps to exactly one (direction, gammaSign)
   pair, so binning by that pair says whether the disagreement is a property of ONE
   table -- a localized sign error -- or spread over all four. Order:
     0 ul (dir -1, sign -1)   1 vl (dir -1, sign +1)
     2 uu (dir +1, sign -1)   3 vu (dir +1, sign +1)  */
static int lpProfMisBin[4] = {0, 0, 0, 0};
static int lpProfAgreeBin[4] = {0, 0, 0, 0};
/* Mistaken cases where the exact value is the table's value NEGATED, to within
   relative 1e-9. A large count here means the table is not imprecise, it is one
   sign flip away from correct, which is a very different repair. */
static int lpProfMisNegated = 0;
/* Ambiguous-column census, split the same way. mis_amb0 / agree_amb0 are the
   decisive cells: the perturbation explanation predicts mis_amb0 == 0 (no mistaken
   call is free of ambiguous columns) and agree_amb0 == the whole agree population
   minus the lucky draws. A nonzero mis_amb0 falsifies it outright. */
static int lpProfMisAmb0 = 0, lpProfAgreeAmb0 = 0;
static double lpProfMisAmbSum = 0., lpProfAgreeAmbSum = 0.;

static int lpRcBin(int direction, int gammaSign)
{
  return (direction > 0 ? 2 : 0) + (gammaSign > 0 ? 1 : 0);
}

void cglLandPProfileMistaken(double exact, double table, int direction,
                             int gammaSign, int nAmb)
{
  if (nAmb == 0) ++lpProfMisAmb0;
  lpProfMisAmbSum += nAmb;
  if (exact > 1e-6) { ++lpProfMistakenSign;
                      if (exact > lpProfMistakenWorst) lpProfMistakenWorst = exact; }
  else ++lpProfMistakenTiny;
  lpProfMistakenTableSum += table;
  lpProfMistakenExactSum += exact;
  ++lpProfMisBin[lpRcBin(direction, gammaSign)];
  const double scale = fabs(table) > 1. ? fabs(table) : 1.;
  if (fabs(exact + table) <= 1e-9 * scale) ++lpProfMisNegated;
}

/* The SAME pair on the calls where the exact test AGREED with the table. This is
   the control for the mistaken bin above, and it is what decides between the two
   readings of that data:

     - if table and exact are close here and only differ on the mistaken calls,
       the table is a good approximation with a tolerance boundary problem;
     - if they differ by a consistent factor here TOO, then the two formulas do
       not compute the same quantity at all, they merely usually share a sign,
       and the "mistake" is a units mismatch rather than a numerical one.

   Sums rather than a max, because a ratio of means is what separates those. */
static int lpProfAgreeN = 0;
static double lpProfAgreeExactSum = 0.0;
static double lpProfAgreeTableSum = 0.0;
static double lpProfAgreeWorstRatio = 0.0;

void cglLandPProfileAgree(double exact, double table, int direction,
                          int gammaSign, int nAmb)
{
  if (nAmb == 0) ++lpProfAgreeAmb0;
  lpProfAgreeAmbSum += nAmb;
  ++lpProfAgreeN;
  ++lpProfAgreeBin[lpRcBin(direction, gammaSign)];
  lpProfAgreeExactSum += exact;
  lpProfAgreeTableSum += table;
  if (table < -1e-12)
    {
      const double ratio = exact / table;   /* both negative => positive ratio */
      if (ratio > lpProfAgreeWorstRatio) lpProfAgreeWorstRatio = ratio;
    }
}

/* Loop-iteration census. COUNTED, not timed, so a parallel sweep collects it
   correctly -- which matters because the question it answers is "is a sparse
   rewrite worth writing?" and that must be settled before spending the effort,
   not after.

   w_col_dense is what the p/q/r/s accumulation loop in fastFindBestPivotColumn
   costs today: one iteration per nonbasic per call, and it is called ~5x per
   pivot. w_col_sparse is what it would cost visiting only columns that can
   contribute -- a column with row_k==0 AND row_i==0 adds exactly zero to all
   four of p,q,r,s and inserts no gamma, so the ratio of the two is the ceiling
   on what a sparse rewrite can buy.

   w_gamma / w_gammadrop are the same question for gammas_.sortIncrElement():
   dropped means the mistaken-rc test then threw the sorted array away.

   w_slot_* are the exact-retry accounting.  A call to exactRowReducedCosts
   examines the (up to) four candidates of the row already in hand: w_slot_retired
   counts the ones its exact reduced cost proves non-improving, which fastFind-
   BestPivotColumn would otherwise have been called on and lost, and w_slot_live
   counts the calls that found a genuine sibling and so skipped a rescan plus a
   tableau row solve.  w_pred_bad / w_pred_sign are the verification: the value
   predicted here must equal the one the column search then computes itself, so
   w_pred_bad must be 0 (DblEqAssert is compiled out at -O2, hence a counter). */
enum { LP_W_COLDENSE = 0, LP_W_COLSPARSE, LP_W_GAMMA, LP_W_GAMMADROP,
       LP_W_PULLDENSE, LP_W_PULLSPARSE, LP_W_M3COL, LP_W_M3ITER,
       LP_W_SLOTCALL, LP_W_SLOTRETIRE, LP_W_SLOTLIVE, LP_W_PREDBAD,
       LP_W_PREDSIGN, LP_W_N };
static double lpProfWork[LP_W_N] = {0., 0., 0., 0., 0., 0., 0., 0.,
                                    0., 0., 0., 0., 0.};
static const char *const lpProfWorkName[LP_W_N] =
  {"w_col_dense", "w_col_sparse", "w_gamma", "w_gammadrop",
   "w_pull_dense", "w_pull_sparse", "w_m3col", "w_m3iter",
   "w_slot_calls", "w_slot_retired", "w_slot_live", "w_pred_bad",
   "w_pred_sign"};

/* double, not long: these reach 1e9 on one fixture and a double counts exactly
   up to 2^53, which is well past anything a single instance can produce. */
void cglLandPProfileWork(int which, double n)
{
  if (which >= 0 && which < LP_W_N) lpProfWork[which] += n;
}

void cglLandPProfileReset()
{
  for (int i = 0; i < LP_W_N; ++i) lpProfWork[i] = 0.;
  for (int i = 0; i < 6; ++i) lpProfOptExit[i] = 0;
  for (int i = 0; i < 3; ++i) lpProfPivCode[i] = 0;
  for (int i = 0; i < LP_ST_N; ++i) { lpProfStage[i] = 0.; lpProfStageN[i] = 0; }
  for (int i = 0; i < LP_CF_N; ++i) lpProfColFail[i] = 0;
  lpProfMistakenTiny = lpProfMistakenSign = 0;
  lpProfMistakenWorst = lpProfMistakenTableSum = lpProfMistakenExactSum = 0.0;
  lpProfAgreeN = lpProfMisNegated = 0;
  lpProfMisAmb0 = lpProfAgreeAmb0 = 0;
  lpProfMisAmbSum = lpProfAgreeAmbSum = 0.;
  for (int i = 0; i < 4; ++i) lpProfMisBin[i] = lpProfAgreeBin[i] = 0;
  lpProfAgreeExactSum = lpProfAgreeTableSum = lpProfAgreeWorstRatio = 0.0;
  lpProfPivots = 0;
  lpProfCloneTime = lpProfLoopTime = lpProfSetupTime = lpProfRangeTime = 0.0;
  lpProfOptCloneTime = 0.0;
  lpProfClones = lpProfCandidates = lpProfRetries = lpProfOptClones = 0;
}

void cglLandPProfilePrint(const char *tag)
{
  /* One line, tab-separated, tagged, so a serial sweep reduces with awk. The
     residual is loop time not attributed to the generateCuts clones, i.e.
     optimize() itself including ITS clone. */
  printf("[landp-prof]\t%s\tsetup\t%.6f\tloop\t%.6f\tclone\t%.6f\toptclone\t%.6f"
         "\trange\t%.6f\tclones\t%d\toptclones\t%d\tcand\t%d\tretries\t%d"
         "\tpivots\t%d\texit_len0\t%d\texit_sigma\t%d\texit_pfup\t%d"
         "\texit_pfsame\t%d\texit_len1\t%d\texit_cut\t%d"
         "\tpiv_ok\t%d\tpiv_refact\t%d\tpiv_singular\t%d",
    tag, lpProfSetupTime, lpProfLoopTime, lpProfCloneTime, lpProfOptCloneTime,
    lpProfRangeTime, lpProfClones, lpProfOptClones, lpProfCandidates,
    lpProfRetries, lpProfPivots, lpProfOptExit[0], lpProfOptExit[1],
    lpProfOptExit[2], lpProfOptExit[3], lpProfOptExit[4], lpProfOptExit[5],
    lpProfPivCode[0], lpProfPivCode[1], lpProfPivCode[2]);
  for (int i = 0; i < LP_ST_N; ++i)
    printf("\t%s\t%.6f", lpProfStageName[i], lpProfStage[i]);
  for (int i = 0; i < LP_ST_N; ++i)
    printf("\tn_%s\t%d", lpProfStageName[i], lpProfStageN[i]);
  for (int i = 0; i < LP_CF_N; ++i)
    printf("\t%s\t%d", lpProfColFailName[i], lpProfColFail[i]);
  printf("\tmis_tiny\t%d\tmis_sign\t%d\tmis_worst\t%.6g\tmis_tablesum\t%.6g"
         "\tmis_exactsum\t%.6g\tagree_n\t%d\tagree_exactsum\t%.6g"
         "\tagree_tablesum\t%.6g\tagree_worstratio\t%.6g",
         lpProfMistakenTiny, lpProfMistakenSign, lpProfMistakenWorst,
         lpProfMistakenTableSum, lpProfMistakenExactSum, lpProfAgreeN,
         lpProfAgreeExactSum, lpProfAgreeTableSum, lpProfAgreeWorstRatio);
  printf("\tmis_negated\t%d\tmis_amb0\t%d\tagree_amb0\t%d"
         "\tmis_ambsum\t%.0f\tagree_ambsum\t%.0f",
         lpProfMisNegated, lpProfMisAmb0, lpProfAgreeAmb0,
         lpProfMisAmbSum, lpProfAgreeAmbSum);
  {
    static const char *const bn[4] = {"ul", "vl", "uu", "vu"};
    for (int i = 0; i < 4; ++i)
      printf("\tmis_%s\t%d\tagree_%s\t%d", bn[i], lpProfMisBin[i],
             bn[i], lpProfAgreeBin[i]);
  }
  for (int i = 0; i < LP_W_N; ++i)
    printf("\t%s\t%.0f", lpProfWorkName[i], lpProfWork[i]);
  printf("\n");
  fflush(stdout);
}

#define LP_PROF_T0(v) const double v = CoinGetTimeOfDay()
#define LP_PROF_ADD(acc, t0) (acc) += CoinGetTimeOfDay() - (t0)
#define LP_PROF_INC(c) ++(c)
#else
#define LP_PROF_T0(v)
#define LP_PROF_ADD(acc, t0)
#define LP_PROF_INC(c)
#endif

namespace LAP
{
//Setup output messages
LapMessages::LapMessages( )
        :CoinMessages(LAP_MESSAGES_DUMMY_END)
{
    strcpy(source_,"Lap");
    addMessage(BEGIN_ROUND,CoinOneMessage( 1, 2,"Starting %s round %d variable considered for separation."));
    addMessage(END_ROUND,CoinOneMessage(2, 2,"End ouf %s round %d cut generated in %g seconds."));
    addMessage(DURING_SEP,CoinOneMessage(3,1,"After %g seconds, separated %d cuts."));
    addMessage(CUT_REJECTED, CoinOneMessage(4,1,"Cut rejected for %s."));
    addMessage(CUT_FAILED,CoinOneMessage(5,1,"Generation failed."));
    addMessage(CUT_GAP, CoinOneMessage(7,1,"CUTGAP after %i pass objective is %g"));
    addMessage(LAP_CUT_FAILED_DO_MIG, CoinOneMessage(3006,1,"Failed to generate a cut generate a Gomory cut instead"));
}
}
using namespace LAP;
CglLandP::Parameters::Parameters():
        CglParam(),
        pivotLimit(20),
        pivotLimitInTree(10),
        maxCutPerRound(5000),
        failedPivotLimit(1),
        degeneratePivotLimit(0),
        extraCutsLimit(5),
	maximumCandidates(1000000),
	maximumCutLength(10000),
        pivotTol(1e-4),
        away(5e-4),
        timeLimit(COIN_DBL_MAX),
        singleCutTimeLimit(COIN_DBL_MAX),
        rhsWeight(1.),
        useTableauRow(true),
        modularize(false),
        strengthen(true),
        countMistakenRc(false),
        sepSpace(Fractional),
        perturb(true),
        exactRetry(true),
        exactBest(false),
        preLengthGate(true),
        normalization(Unweighted),
        rhsWeightType(Fixed),
        lhs_norm(L1),
        generateExtraCuts(none),
        pivotSelection(mostNegativeRc)
{
    EPS = 1e-08;
}

CglLandP::Parameters::Parameters(const Parameters &other):
        CglParam(other),
        pivotLimit(other.pivotLimit),
        pivotLimitInTree(other.pivotLimitInTree),
        maxCutPerRound(other.maxCutPerRound),
        failedPivotLimit(other.failedPivotLimit),
        degeneratePivotLimit(other.degeneratePivotLimit),
        extraCutsLimit(other.extraCutsLimit),
	maximumCandidates(other.maximumCandidates),
	maximumCutLength(other.maximumCutLength),
        pivotTol(other.pivotTol),
        away(other.away),
        timeLimit(other.timeLimit),
        singleCutTimeLimit(other.singleCutTimeLimit),
        rhsWeight(other.rhsWeight),
        useTableauRow(other.useTableauRow),
        modularize(other.modularize),
        strengthen(other.strengthen),
        countMistakenRc(other.countMistakenRc),
        sepSpace(other.sepSpace),
        perturb(other.perturb),
        exactRetry(other.exactRetry),
        exactBest(other.exactBest),
        preLengthGate(other.preLengthGate),
        normalization(other.normalization),
        rhsWeightType(other.rhsWeightType),
        lhs_norm(other.lhs_norm),
        generateExtraCuts(other.generateExtraCuts),
        pivotSelection(other.pivotSelection)
{}

CglLandP::Parameters & CglLandP::Parameters::operator=(const Parameters &other)
{
    if (this != &other)
    {
        CglParam::operator=(other);
        pivotLimit = other.pivotLimit;
        pivotLimitInTree = other.pivotLimitInTree;
        maxCutPerRound = other.maxCutPerRound;
        failedPivotLimit = other.failedPivotLimit;
        /* Was other.failedPivotLimit -- a copy-paste, and not a harmless one: the
           two members have different defaults (1 and 0), so assigning a
           default-constructed Parameters used to turn degenerate pivots ON. The
           copy constructor above always had this right, which is why it survived:
           clone() uses the copy constructor and only CglLandP::operator= reaches
           here. */
        degeneratePivotLimit = other.degeneratePivotLimit;
        extraCutsLimit = other.extraCutsLimit;
	maximumCandidates = other.maximumCandidates;
	maximumCutLength = other.maximumCutLength;
        pivotTol = other.pivotTol;
        away = other.away;
        timeLimit = other.timeLimit;
        singleCutTimeLimit = other.singleCutTimeLimit;
        rhsWeight = other.rhsWeight;
        useTableauRow = other.useTableauRow;
        modularize = other.modularize;
        strengthen = other.strengthen;
        countMistakenRc = other.countMistakenRc;
        sepSpace = other.sepSpace;
        perturb = other.perturb;
        exactRetry = other.exactRetry;
        exactBest = other.exactBest;
        preLengthGate = other.preLengthGate;
        normalization = other.normalization;
        rhsWeightType = other.rhsWeightType;
        lhs_norm = other.lhs_norm;
        generateExtraCuts = other.generateExtraCuts;
        pivotSelection = other.pivotSelection;
    }
    return *this;
}

CglLandP::CachedData::CachedData(int nBasics, int nNonBasics):
        basics_(NULL), nonBasics_(NULL), nBasics_(nBasics),
        nNonBasics_(nNonBasics), basis_(NULL), colsol_(NULL),
        slacks_(NULL), integers_(NULL), solver_(NULL)
{
    if (nBasics_>0)
    {
        basics_ = new int[nBasics_];
        integers_ = new bool [nNonBasics_ + nBasics_];
    }
    if (nNonBasics_>0)
        nonBasics_ = new int[nNonBasics_];
    if (nBasics_ + nNonBasics_ > 0)
    {
        colsol_ = new double[nBasics_ + nNonBasics_];
        slacks_ = &colsol_[nNonBasics_];
    }
}

CglLandP::CachedData::CachedData(const CachedData &source):
  basics_(NULL), nonBasics_(NULL), nBasics_(source.nBasics_),
        nNonBasics_(source.nNonBasics_), basis_(NULL),
        colsol_(NULL), slacks_(NULL), integers_(NULL), solver_(NULL)
{
    if (nBasics_>0)
    {
        basics_ = new int[nBasics_];
        CoinCopyN(source.basics_, nBasics_, basics_);
        integers_ = new bool [nNonBasics_ + nBasics_];
        CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);
    }
    if (nNonBasics_>0)
    {
        nonBasics_ = new int[nNonBasics_];
        CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
    }
    if (nBasics_ + nNonBasics_ > 0)
    {
        colsol_ = new double[nBasics_ + nNonBasics_];
        slacks_ = &colsol_[nNonBasics_];
        CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
    }
    if (source.basis_!=NULL)
        basis_ = new CoinWarmStartBasis(*source.basis_);
    if (source.solver_!=NULL)
      solver_ = source.solver_->clone();
}

CglLandP::CachedData& CglLandP::CachedData::operator=(const CachedData &source)
{
    if (this != &source)
    {
        nBasics_ = source.nBasics_;
        nNonBasics_ = source.nNonBasics_;
        delete [] basics_;
        basics_ = NULL;
        delete [] nonBasics_;
        nonBasics_ = NULL;
        delete [] basis_;
        basis_ = NULL;
        delete [] colsol_;
        colsol_ = NULL;
        delete [] slacks_;
        slacks_ = NULL;
        delete [] integers_;
        integers_ = NULL;
        if (nBasics_>0)
        {
            basics_ = new int[nBasics_];
            CoinCopyN(source.basics_, nBasics_, basics_);
            integers_ = new bool [nBasics_ + nNonBasics_];
            CoinCopyN(source.integers_, nBasics_ + nNonBasics_, integers_);
        }
        if (nNonBasics_>0)
        {
            nonBasics_ = new int[nNonBasics_];
            CoinCopyN(source.nonBasics_, nBasics_, nonBasics_);
        }
        if (nBasics_ + nNonBasics_ > 0)
        {
            colsol_ = new double[nBasics_ + nNonBasics_];
            slacks_ = &colsol_[nNonBasics_];
            CoinCopyN(source.colsol_, nBasics_ + nNonBasics_, colsol_);
        }
        if (source.basis_!=NULL)
            basis_ = new CoinWarmStartBasis(*source.basis_);
        delete solver_;
	if (source.solver_)
	  solver_ = source.solver_->clone();
    }
    return *this;
}

void
CglLandP::CachedData::getData(const OsiSolverInterface &si)
{
    int nBasics = si.getNumRows();
    int nNonBasics = si.getNumCols();
    if (basis_ != NULL)
        delete basis_;
    basis_ = dynamic_cast<CoinWarmStartBasis *> (si.getWarmStart());
    if (!basis_)
        throw NoBasisError();

    if (nBasics_ > 0 || nBasics != nBasics_)
    {
        delete [] basics_;
        basics_ = NULL;
    }
    if (basics_ == NULL)
    {
        basics_ = new int[nBasics];
        nBasics_ = nBasics;
    }

    if (nNonBasics_ > 0 || nNonBasics != nNonBasics_)
    {
        delete [] nonBasics_;
        nonBasics_ = NULL;
    }
    if (nonBasics_ == NULL)
    {
        nonBasics_ = new int[nNonBasics];
        nNonBasics_ = nNonBasics;
    }
    int n = nBasics + nNonBasics;
    if ( nBasics_ + nNonBasics_ > 0 || nBasics_ + nNonBasics_ != n)
    {
        delete [] colsol_;
        delete [] integers_;
        integers_ = NULL;
        colsol_ = NULL;
        slacks_ = NULL;
    }
    if (colsol_ == NULL)
    {
        colsol_ = new double[n];
        slacks_ = &colsol_[nNonBasics];
    }

    if (integers_ == NULL)
    {
        integers_ = new bool[n];
    }

    
    const double * rowLower = si.getRowLower();
    const double * rowUpper = si.getRowUpper();
    //determine which slacks are integer
    const CoinPackedMatrix * m = si.getMatrixByCol();
    const double * elems = m->getElements();
    const int * inds = m->getIndices();
    const CoinBigIndex * starts = m->getVectorStarts();
    const int * lengths = m->getVectorLengths();
    //    int numElems = m->getNumElements();
    int numCols = m->getNumCols();
    assert(numCols == nNonBasics_);
    //   int numRows = m->getNumRows();
    CoinFillN(integers_ ,n, true);
    for (int i = 0 ;  i < numCols ; i++)
    {
        if (si.isContinuous(i))
            integers_[i] = false;
    }
    bool * integerSlacks = integers_ + numCols;
    for (int i = 0 ; i < nBasics ; i++)
    {
        if (rowLower[i] > -1e50 && INT_INFEAS(rowLower[i]) > 1e-15)
            integerSlacks[i] = false;
        if (rowUpper[i] < 1e50 && INT_INFEAS(rowUpper[i]) > 1e-15)
            integerSlacks[i] = false;
    }
    for (int i = 0 ;  i < numCols ; i++)
    {
        CoinBigIndex end = starts[i] + lengths[i];
        if (integers_[i])
        {
            for (CoinBigIndex k=starts[i] ; k < end; k++)
            {
                if (integerSlacks[inds[k]] && INT_INFEAS(elems[k])>1e-15 )
                    integerSlacks[inds[k]] = false;
            }
        }
        else
        {
            for (CoinBigIndex k=starts[i] ; k < end; k++)
            {
                if (integerSlacks[inds[k]])
                    integerSlacks[inds[k]] = false;
            }
        }
    }

    CoinCopyN(si.getColSolution(), si.getNumCols(), colsol_);
    CoinCopyN(si.getRowActivity(), si.getNumRows(), slacks_);
    for (int i = 0 ; i < si.getNumRows() ; i++)
    {
        slacks_[i]*=-1;
        if (rowLower[i]>-1e50)
        {
            slacks_[i] += rowLower[i];
        }
        else
        {
            slacks_[i] += rowUpper[i];
        }
    }
    //Now fill the arrays;
    nNonBasics = 0;
    nBasics = 0;



    //For having the index variables correctly ordered we need to access to OsiSimplexInterface
    {
        OsiSolverInterface * ncSi = (const_cast<OsiSolverInterface *>(&si));
        ncSi->enableSimplexInterface(0);
        ncSi->getBasics(basics_);
	// Save enabled solver
	solver_ = si.clone();
#ifdef CGL_HAS_OSICLP
	OsiClpSolverInterface * clpSi = getClpSolver(solver_);
	const OsiClpSolverInterface * clpSiRhs = getConstClpSolver(&si);
	if (CBC_SKIP_CLP_TEST||clpSi)
	  clpSi->getModelPtr()->copyEnabledStuff(clpSiRhs->getModelPtr());;
#endif
        ncSi->disableSimplexInterface();
    }

    int numStructural = basis_->getNumStructural();
    for (int i = 0 ; i < numStructural ; i++)
    {
        if (basis_->getStructStatus(i)== CoinWarmStartBasis::basic)
        {
	    nBasics++;
            //Basically do nothing
        }
        else
        {
            nonBasics_[nNonBasics++] = i;
        }
    }

    int numArtificial = basis_->getNumArtificial();
    int numStruct = basis_->getNumStructural();
    for (int i = 0 ; i < numArtificial ; i++)
    {
        if (basis_->getArtifStatus(i)== CoinWarmStartBasis::basic)
        {
            //Just check number of basics
            nBasics++;
        }
        else
        {
            nonBasics_[nNonBasics++] = i + numStruct;
        }
    }
}
void
CglLandP::CachedData::clean(){
    if (basics_!=NULL)
        delete [] basics_;
    basics_ = NULL;
    if (nonBasics_!=NULL)
        delete [] nonBasics_;
    nonBasics_ = NULL;
    if (colsol_ != NULL)
        delete [] colsol_;
    colsol_ = NULL;
    delete basis_;
    basis_ = NULL;
    if (integers_)
        delete [] integers_;
    integers_ = NULL;

   nBasics_ = 0;
   nNonBasics_ = 0;
   delete solver_;
   solver_ = NULL;
}
CglLandP::CachedData::~CachedData()
{
    if (basics_!=NULL)
        delete [] basics_;
    if (nonBasics_!=NULL)
        delete [] nonBasics_;
    if (colsol_ != NULL)
        delete [] colsol_;
    delete basis_;
    if (integers_)
        delete [] integers_;
    delete solver_;
}

CglLandP::CglLandP(const CglLandP::Parameters &params,
                   const LAP::Validator &validator):
        params_(params), cached_(), validator_(validator), numrows_(-1),
        numcols_(-1),originalColLower_(NULL), originalColUpper_(NULL),
        canLift_(false),
        extraCuts_()
{
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(0);
    messages_ = LapMessages();
}


CglLandP::~CglLandP()
{
    delete handler_;
    if (originalColLower_ != NULL)
        delete [] originalColLower_;
    if (originalColUpper_ != NULL)
        delete [] originalColUpper_;
}

CglLandP::CglLandP(const CglLandP & source):
        CglCutGenerator(source),
        params_(source.params_), cached_(source.cached_),
        validator_(source.validator_), numrows_(source.numrows_),numcols_(source.numcols_),
        originalColLower_(NULL), originalColUpper_(NULL),
        canLift_(source.canLift_),
        extraCuts_(source.extraCuts_)
{
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(source.handler_->logLevel());
    messages_ = LapMessages();
    if (numcols_ != -1)
    {
        assert(numcols_ > 0);
        assert(originalColLower_!=NULL);
        assert(originalColUpper_!=NULL);
        originalColLower_ = new double[numcols_];
        originalColUpper_ = new double[numcols_];
        CoinCopyN(source.originalColLower_,numcols_,originalColLower_);
        CoinCopyN(source.originalColUpper_,numcols_,originalColUpper_);
    }
}

/** Assignment operator */
CglLandP& CglLandP::operator=(const CglLandP &rhs)
{
    if (this != &rhs)
    {
        params_ = rhs.params_;
        cached_ = rhs.cached_;
        validator_ = rhs.validator_;
        extraCuts_ = rhs.extraCuts_;
    }
    return *this;
}


CglCutGenerator *
CglLandP::clone() const
{
    return new CglLandP(*this);
}

extern double restaurationTime;

struct cutsCos
{
    int i;
    int j;
    double angle;
    cutsCos(int  i_, int j_ , double angle_):i(i_), j(j_), angle(angle_)
    {
    }
    bool operator<(const cutsCos&other)const
    {
        return angle > other.angle;
    }
};


void
CglLandP::scanExtraCuts(OsiCuts& cs, const double * colsol) const
{
    //int numAdded = 0;
    for (int i = extraCuts_.sizeRowCuts() - 1; i > -1 ; i--)
    {
        double violation = extraCuts_.rowCut(i).violated(colsol);
        if (violation > 0.)
        {
            cs.insert(extraCuts_.rowCut(i));
            //numAdded++;
            //      std::cout<<"A cut computed in a previous iteration is violated by "<<violation<<"."<<std::endl;
            //extraCuts_.eraseRowCut(i);
        }
    }
    //  std::cout<<"Added "<<numAdded<<" previously generated cuts."<<std::endl;
}

void
CglLandP::generateCuts(const OsiSolverInterface & si, OsiCuts & cs,
                       const CglTreeInfo info )
{
    int numberRanges = 0;
    if (numrows_<0)
    {
      numrows_ = si.getNumRows();
      // but switch off? if ranges
      const double * rowLower = si.getRowLower();
      const double * rowUpper = si.getRowUpper();
      for (int i=0;i<numrows_;i++) {
	if (rowLower[i]<rowUpper[i]) {
	  if (rowLower[i]> -1.0e50 && rowUpper[i] < 1.0e50) {
	    numberRanges++;
	  }
	}
      }
      if (numberRanges && false) {
	params_.maximumCutLength =-1;
	return;
      }
    } else if (params_.maximumCutLength < 0) {
      return;
    }
// scanExtraCuts(cs, si.getColSolution());
    Parameters params = params_;
    params.rhsWeight = numrows_ + 2;

    handler_->message(CUT_GAP, messages_)<<info.pass<<si.getObjValue() <<CoinMessageEol;

    if (info.inTree)   //put lower pivot limit
    {
        params.pivotLimit = std::min(params.pivotLimit, params.pivotLimitInTree);
        params.countMistakenRc = true;
    }
    if (params.timeLimit < 0)
    {
        params.pivotLimit = 0;
    }

    assert(si.basisIsAvailable());


#ifdef APPEND_ROW
    OsiSolverInterface * t_si = si.clone();
    if (params.modularize)
    {
      if (numberRanges) {
	// modify
      }
      int new_idx = si.getNumCols();
      int v_idx[1] = {new_idx};
      double v_val[1] = {-1};
      CoinPackedVector v(1, v_idx, v_val, false);
      t_si->addCol(CoinPackedVector(), 0, 1, 0);
      t_si->setInteger(new_idx);
      t_si->addRow(v,0, 0);
      t_si->resolve();
    }
#else
    const OsiSolverInterface * t_si = &si;
    int numberCuts = cs.sizeRowCuts(); 
    double * upper = NULL;
    CoinBigIndex * starts = NULL;
    int * row = NULL;
    int nThrownAway=0;
    if (numberRanges) {
      // modify by adding slacks
      LP_PROF_T0(tRange);
      OsiSolverInterface *tt_si = si.clone();
      const double * rowLower = si.getRowLower();
      const double * rowUpper = si.getRowUpper();
      const double * rowActivity = si.getRowActivity();
      upper = new double[4*numberRanges];
      double * lower = upper+numberRanges;
      double * obj = lower+numberRanges;
      double * element = obj+numberRanges;
      memset(lower,0,2*numberRanges*sizeof(double));
      starts = new CoinBigIndex[numberRanges+1];
      row = new int[numberRanges];
      // First add columns
      numberRanges = 0;
      starts[0] = 0;
      for (int i=0;i<numrows_;i++) {
	if (rowLower[i]<rowUpper[i]) {
	  if (rowLower[i]> -1.0e50 && rowUpper[i] < 1.0e50) {
	    double value = rowActivity[i];
	    if (value-rowLower[i] < rowUpper[i]-value) {
	      // keep rowLower
	      //tt_si->setRowUpper(i,rowLower[i]);
	      element[numberRanges] = -1.0;
	    } else {
	      // keep rowUpper
	      //tt_si->setRowLower(i,rowUpper[i]);
	      element[numberRanges] = 1.0;
	    }
	    row[numberRanges] = i;
	    upper[numberRanges++] = rowUpper[i]-rowLower[i];
	    starts[numberRanges] = numberRanges;
	  }
	}
      }
      tt_si->addCols(numberRanges,starts,row,element,lower,upper,obj);
      // Now basis and values
      int nCol = si.getNumCols();
      double * solution = CoinCopyOfArray(tt_si->getColSolution(),
					  nCol+numberRanges);
      CoinWarmStartBasis * basis
	= dynamic_cast<CoinWarmStartBasis *>(tt_si->getWarmStart());
      numberRanges = 0;
      for (int i=0;i<numrows_;i++) {
	if (rowLower[i]<rowUpper[i]) {
	  if (rowLower[i]> -1.0e50 && rowUpper[i] < 1.0e50) {
	    double value = rowActivity[i];
	    if (value-rowLower[i] < rowUpper[i]-value) {
	      // keep rowLower
	      tt_si->setRowUpper(i,rowLower[i]);
	      value -= rowLower[i];
	    } else {
	      // keep rowUpper
	      tt_si->setRowLower(i,rowUpper[i]);
	      value = rowUpper[i]-value;
	    }
	    solution[nCol+numberRanges] = value;
	    if (basis->getArtifStatus(i) ==
		CoinWarmStartBasis::basic) {
	      // set basic
	      basis->setStructStatus(nCol+numberRanges,
				     CoinWarmStartBasis::basic);
	      basis->setArtifStatus(i,
				    CoinWarmStartBasis::atLowerBound);
	    }
	    numberRanges++;
	  }
	}
      }
      // update solution and basis
      tt_si->setColSolution(solution);
      tt_si->setWarmStart(basis);
      delete basis;
      delete [] solution;
      tt_si->resolve();
      t_si = tt_si;
      LP_PROF_ADD(lpProfRangeTime, tRange);
    }
#endif

    LP_PROF_T0(tSetup);
    cached_.getData(*t_si);
    CglLandPSimplex landpSi(*t_si, cached_, params, validator_);
    if (params.generateExtraCuts == CglLandP::AllViolatedMigs)
    {
        landpSi.genThisBasisMigs(cached_, params);
    }
    landpSi.setLogLevel(handler_->logLevel());
    int nCut = 0;

    std::vector<int> indices;
    getSortedFractionalIndices(indices,cached_, params);
    if (indices.size()>params.maximumCandidates)
      indices.resize(params.maximumCandidates);
    LP_PROF_ADD(lpProfSetupTime, tSetup);

#ifndef NDEBUG
    int numrows = si.getNumRows();
#endif

#ifdef DO_STAT
    //Get informations on current optimum
    {
        OsiSolverInterface * gapTester = si.clone();
        gapTester->resolve();

        roundsStats_.analyseOptimalBasis(gapTester,info.pass, numrows_);
        delete gapTester;
    }
#endif

    params_.timeLimit += CoinCpuTime();
    CoinRelFltEq eq(1e-04);

    LP_PROF_T0(tLoop);
    for (unsigned int i = 0; i < indices.size() && nCut < params.maxCutPerRound &&
            nCut < cached_.nBasics_ ; i++)
    {

        //Check for time limit
        int iRow = indices[i];
        LP_PROF_INC(lpProfCandidates);
        assert(iRow < numrows);
        OsiRowCut cut;
        int code=1;

        /* No clone here. optimize() opens by doing `delete si_; si_ =
           cached.solver_->clone()` unconditionally, so the clone that used to be
           made at this point was destroyed before a single line of it was read --
           one full LP copy per candidate, plus a second one on the retry path
           below, for nothing. cached_.solver_ is itself a clone of *t_si taken by
           cached_.getData(*t_si) above, so optimize() starts from the same matrix
           and basis either way.

           Two setters went with it: setDblParam(OsiDualObjectiveLimit,
           COIN_DBL_MAX) and setLogLevel(0). They were applied to the discarded
           clone, so they never reached the solver that pivots; they are dead
           today and removing them changes nothing. Whether the *live* clone ought
           to have them is a separate question -- answering it yes would be a
           behaviour change, so it is not folded in here. */

        int generated = 0;
        if (params.pivotLimit == 0)
        {
            generated = landpSi.generateMig(iRow, cut, params);
        }
        else
        {
            generated = landpSi.optimize(iRow, cut, cached_, params);
            if (params.generateExtraCuts == CglLandP::AllViolatedMigs)
            {
                landpSi.genThisBasisMigs(cached_, params);
            }
            landpSi.resetSolver(cached_.basis_);
        }
        code = 0;
        if (generated)
            code = validator_(cut, cached_.colsol_, si, params, originalColLower_, originalColUpper_);
        if (!generated || code)
        {
            if (params.pivotLimit !=0)
            {
                handler_->message(LAP_CUT_FAILED_DO_MIG, messages_)<<validator_.failureString(code)<<CoinMessageEol;
                landpSi.freeSi();
                LP_PROF_INC(lpProfRetries);
                /* Same dead clone as above: optimize() re-clones cached.solver_
                   on entry whatever pivotLimit is, so the retry needs nothing
                   from here but the freeSi() that precedes it. */
                params.pivotLimit = 0;
                if (landpSi.optimize(iRow, cut, cached_, params))
                {
                    code = validator_(cut, cached_.colsol_, si, params, originalColLower_, originalColUpper_);
                }
                else if (!code)
                {
                    // Neither attempt produced anything, so `cut` is still
                    // default constructed and there is nothing to validate. That
                    // only happens when the first optimize() failed too, which
                    // leaves code at 0 -- and code 0 means "accepted", so the
                    // empty cut was inserted into the cut set as a free row
                    // (-inf <= 0 <= +inf) and charged against the nCut budget.
                    // Reachable through any of optimize()'s failure exits --
                    // all of them return before the cut is built -- most easily
                    // the pivot-failure ones in CglLandPSimplex.cpp. Not observed
                    // firing on this corpus, so this is a latent defect closed by
                    // inspection rather than one caught in the act.
                    code = Validator::EmptyCut;
                }
                params.pivotLimit = params_.pivotLimit;
            }
        }

        if (params.pivotLimit != 0)
        {
            landpSi.freeSi();
        }
        if (code)
        {
            handler_->message(CUT_REJECTED, messages_)<<
            validator_.failureString(code)<<CoinMessageEol;
        }
        else
        {
  	    if (numberRanges) {
	      int numberColumns = si.getNumCols();
	      int n = cut.row().getNumElements();
	      const int *column = cut.row().getIndices();
	      const double *element = cut.row().getElements();
	      bool oddSlack = false;
	      for (int i=0;i<n;i++) {
		if (column[i]>=numberColumns) {
		  oddSlack=true;
		}
	      }
	      if (oddSlack) {
		nThrownAway++;
#if 0
		// for now throw away
		printf("odd cut\n %d entries ncol=%d, nrange %d %g<=%g\n",
		       n,numberColumns,numberRanges,cut.lb(),cut.ub());
		for (int i=0;i<n;i++) {
		  printf("(%g*x%d) ",column[i],element[i]);
		  if ((i%5)==4)
		    printf("\n");
		}
		if ((n%5))
		  printf("\n");
#endif
		// The cut references a slack column of a range row, so it cannot
		// be expressed in the original column space at all. Skip it,
		// rather than replacing it with a default OsiRowCut: that got
		// inserted into the cut set as a free row (-inf <= 0 <= +inf) and
		// counted against the nCut budget. Every other rejection in this
		// loop already neither inserts nor counts.
		continue;
	      }
	    }
            if (canLift_)
            {
                cut.setGloballyValid(true);
            }
#ifdef CHECK_KNOWN_SOLUTION
	    const OsiRowCutDebugger *debugger = si.getRowCutDebugger();
	    if (debugger) {
	      if (debugger->invalidCut(cut)) {
		printf("BAD cut\n");
		exit(0);
	      }
	    }
	    //CoinAssert (!debugger->invalidCut(*cut));
#endif
            cs.insertIfNotDuplicate(cut, eq);
            //cs.insertIfNotDuplicate(cut);
            {
                //std::cout<<"Violation "<<cut.violated(cached_.colsol_)<<std::endl;
                nCut++;
            }
        }
    }

    LP_PROF_ADD(lpProfLoopTime, tLoop);

    Cuts& extra = landpSi.extraCuts();
    for (int i = 0 ; i < cached_.nNonBasics_; i++)
    {
        OsiRowCut * cut = extra.rowCut(i);
        if (cut == NULL) continue;
        int code = validator_(*cut, cached_.colsol_, si, params,
                              originalColLower_, originalColUpper_);
        if (code)
        {
            handler_->message(LAP_CUT_FAILED_DO_MIG, messages_)
            <<validator_.failureString(code)<<CoinMessageEol;
        }
        else
        {
            cs.insertIfNotDuplicate(*cut, eq);
            {
                nCut++;
            }
        }
        delete cut;
    }

    landpSi.outPivInfo(nCut);
    params_.timeLimit -= CoinCpuTime();

    cached_.clean();
#ifdef APPEND_ROW
    assert(t_si != &si);
    delete t_si;
#else
    if (t_si != &si) {
      delete [] upper;
      delete [] starts;
      delete [] row;
      delete t_si;
      assert (numberRanges);
      //printf("Ranges (%d) - switching off CglLandP after this %d cuts, %d thrown\n",
      //     numberRanges,
      //     cs.sizeRowCuts()-numberCuts, 
      //     nThrownAway);
      params_.maximumCutLength = -params_.maximumCutLength;
    }
#endif
    return;
}


template < class S, class T, class U >
class StableCompare
{
public:
    inline bool operator()(const CoinTriple<S,T,U>& t1,
                           const CoinTriple<S,T,U>& t2) const
    {
        return (t1.third < t2.third) ||
               ((t1.third == t2.third) && (t1.second < t2.second));
    }

};

template <class T1,class T2>
struct StableExternalComp
{
    const std::vector<T1> &vec_1_;
    const std::vector<T2> &vec_2_;
    StableExternalComp(const std::vector<T1> &vec_1,
                       const std::vector<T2> &vec_2):
            vec_1_(vec_1),
            vec_2_(vec_2)
    {
    }
    CoinRelFltEq eq;
    bool operator()(int i, int j)
    {
        bool result = (vec_1_[i] < vec_1_[j]) ||
                      ( ((vec_1_[i]== vec_1_[j]))
                        && (vec_2_[i] < vec_2_[j]));
        return result;
    }

};
void
CglLandP::getSortedFractionalIndices(std::vector<int> &frac_indices,
                                     const CachedData &data,
                                     const CglLandP::Parameters & params) const
{
    std::vector<int> colIndices;
    std::vector<double> values;
    std::vector<int> indices;
    for (int i = 0 ; i < data.nBasics_ ; i++)
    {
        const int& iCol = data.basics_[i];
        if (iCol >= data.nNonBasics_ ||
                !data.integers_[iCol] ||
                INT_INFEAS(data.colsol_[iCol]) <= params.away)
            continue;
        const double value = INT_INFEAS(data.colsol_[iCol]);

        frac_indices.push_back(i);
        indices.push_back(static_cast<int>(values.size()));
        values.push_back(- value);
        colIndices.push_back(iCol);
    }
    std::sort(indices.begin(), indices.end(),StableExternalComp<double, int>(values,colIndices));
    colIndices = frac_indices;
    for (unsigned int i = 0; i < indices.size() ; i++)
    {
        frac_indices[i] = colIndices[indices[i]];
    }

}


