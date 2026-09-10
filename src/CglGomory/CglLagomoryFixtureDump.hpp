// Copyright (C) 2026, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/**
 * CglLagomoryFixtureDump -- serialize one *Lagrangean* Gomory call so it can be
 * replayed offline.
 *
 * WHY THIS EXISTS SEPARATELY FROM CbcGomoryFixtureDump. The 330 fixtures written
 * by that header cannot exercise lagomory at all, and this was checked
 * empirically rather than argued: running gomory-bench over them with
 * `--orig-solver --gomory-type=11` and `=12` gives 0 cuts and a separation time
 * of ~4 microseconds, and with `=21` / `=22` gives cut statistics
 * *byte-identical* to plain `--gomory-type=0` (50v-10: 28 row cuts, totalViol
 * 218.2481753, objImprove 197.701 in all three). The reason is one line: the
 * Cbc-side dump fires at `currentPassNumber_ == 1` (CbcModel.cpp:9494), i.e.
 * before any cut has been added to the LP, so `numberRows == numberOriginalRows`
 * and there are *no cut rows to dualize*. The `.meta` of every existing fixture
 * says the same thing out loud -- `rows 491` next to `infoFormulationRows 491`.
 *
 * So a lagomory fixture has to be captured later in the cut loop, and from a
 * different place.
 *
 * WHY THE CAPTURE IS INSIDE THE GENERATOR AND NOT INSIDE CbcModel. Three reasons,
 * in increasing order of how much time they save:
 *
 *  1. `originalSolver_` is a private member of CglGomory with no getter, and it
 *     is what the Lagrangean pass actually generates cuts from (`useSolver`).
 *     CbcModel has no access to it.
 *  2. It cannot be faithfully reconstructed from `si` either. Deleting rows
 *     [numberOriginalRows, numberRows) from `si` looks equivalent, but the
 *     Lagrangean pass never syncs *row* bounds -- only column bounds
 *     (CglGomory.cpp:128-131) -- so if CBC tightened a formulation row's rhs
 *     during the cut loop, `si`'s copy of that row and `originalSolver_`'s copy
 *     differ, and the reconstruction would silently be the wrong LP.
 *  3. Dumping from inside the `whenToDo` gate lands automatically on the first
 *     call that does Lagrangean work, with no guessing about which CBC pass
 *     number that is. The gate has three arms (CglGomory.cpp:121-126) and two of
 *     them do not mention the pass at all.
 *
 * WHAT A REPLAY NEEDS, AND WHERE EACH PIECE COMES FROM. The Lagrangean pass reads
 * exactly five things from outside the generator:
 *
 *   si's augmented matrix + row bounds   -> `.mps.gz` (rows 0..numberOriginalRows-1
 *                                          are the formulation, the rest are the
 *                                          cuts to be dualized or copied)
 *   si's column bounds                   -> same `.mps.gz` (+ `.ctype`, see below)
 *   si's basis                           -> `.bas`, truncated to numberOriginalRows
 *                                          by the algorithm itself at :224
 *   si's row duals, getRowPrice()        -> `.pi`, and this is the objective
 *                                          perturbation: obj[j] -= pi[i]*a[i][j]
 *   si's column solution                 -> `.sol`, used *only* by the final
 *                                          violation filter at :293-323
 *   the original formulation             -> `.orig.mps.gz` + `.orig.ctype`
 *
 * `.pi` and `.sol` are dumped even though a replay that reads `.bas` and resolves
 * recomputes both. They are the fidelity check: a fixture whose resolve lands on a
 * different vertex is a fixture that measures a different experiment, and without
 * a stored dual vector that failure is invisible. (Four of eight CglOddWheel
 * fixtures stopped reproducing when presolve moved, and only a stored-vs-recomputed
 * comparison would have caught it early.) The bench compares and reports; it does
 * not inject.
 *
 * WHAT IS *NOT* CAPTURED, AND WHY THAT IS SOUND. `originalSolver_`'s column bounds
 * and objective are provenance only. The algorithm overwrites the bounds from `si`
 * on entry (:128-131) and restores the objective from its own saved copy on exit
 * (:293), so on the second and later calls the values sitting in `originalSolver_`
 * are the *previous* node's -- stale, and immediately overwritten. A replay that
 * overwrites them from `si` in the same way is therefore faithful regardless of
 * what was stored. The `.orig.mps.gz` matters for its *matrix and row bounds*.
 *
 * Not captured either: ClpSimplex's factorization, scaling and perturbation state.
 * A replay starts from a fresh factorization of the same basis. This is the one
 * place a replay can legitimately diverge from the dumping run, so the exactness
 * gate for any optimization is before-vs-after *within* the bench, never
 * bench-vs-CBC.
 *
 * THE .ctype SIDECAR, AGAIN. MPS cannot express "fixed and integer": lb == ub
 * takes writeMps's " FX " branch, which has no integer form, and the column reads
 * back continuous. For the Lagrangean pass this bites harder than for plain
 * Gomory, because `intVar[]` is classified from `si`'s *node* bounds (:82-99),
 * and inside the tree many integer columns are fixed -- exactly the ones MPS
 * loses. A lost marker moves a column from `intVar=4` to `intVar=0`, which changes
 * its coefficient contribution in the cut derivation and bumps `numberNonInteger`,
 * which feeds the rhs relaxation ladder. Both solvers get their own sidecar.
 *
 * Entirely behind CGL_DUMP_LAGOMORY_FIXTURE and off by default. Header-only static
 * functions, so no Makefile.am/Makefile.in changes are needed.
 *
 * Environment:
 *   CGL_LAGOMORY_FIXTURE_DIR    output directory
 *                               (default ~/instances/miplib/2017+spp/lagomoryFixtures)
 *   CGL_LAGOMORY_FIXTURE_NAME   instance base name, when the MPS problem name is
 *                               unhelpful (drivers normally set this)
 *   CGL_LAGOMORY_FIXTURE_MAX    stop after this many dumps in one process
 *                               (default 1 -- the first qualifying call)
 **/

#ifndef CglLagomoryFixtureDump_H
#define CglLagomoryFixtureDump_H

#ifdef CGL_DUMP_LAGOMORY_FIXTURE

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>
#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"
#include "OsiSolverInterface.hpp"

/// mkdir -p, so a driver need not pre-create the tree.
static void cglLagomoryFixtureMkdirP(const std::string &path)
{
  for (size_t i = 1; i <= path.size(); ++i) {
    if (i < path.size() && path[i] != '/')
      continue;
    const std::string part = path.substr(0, i);
#ifdef _WIN32
    _mkdir(part.c_str());
#else
    mkdir(part.c_str(), 0775);
#endif
  }
}

/// Output directory, from the environment or the default fixture location.
static std::string cglLagomoryFixtureDir()
{
  const char *env = getenv("CGL_LAGOMORY_FIXTURE_DIR");
  if (env && *env)
    return std::string(env);
  const char *home = getenv("HOME");
  std::string base = home && *home ? std::string(home) : std::string(".");
  return base + "/instances/miplib/2017+spp/lagomoryFixtures";
}

/// Instance base name. The environment wins, because by the time a cut generator
/// runs the solver's problem name has usually been replaced by CglPreProcess.
static std::string cglLagomoryFixtureName(const OsiSolverInterface &si)
{
  const char *env = getenv("CGL_LAGOMORY_FIXTURE_NAME");
  if (env && *env)
    return std::string(env);
  std::string name;
  si.getStrParam(OsiProbName, name);
  if (name.empty())
    name = "unnamed";
  for (size_t i = 0; i < name.size(); ++i) {
    if (name[i] == '/' || name[i] == ' ')
      name[i] = '_';
  }
  return name;
}

/**
 * Write an MPS whose column indices survive the round trip.
 *
 * writeMps drops a column with no matrix entry, which shifts every later column
 * index and so invalidates the `.ctype`, the `.bas` and the `.sol` all at once. A
 * single redundant final row covering those columns pins them at their original
 * index; `paddedColumns` in the `.meta` says how many there were and the loader
 * deletes the extra row. Note where that row lands: *after* the cut rows, so a
 * loader that forgets to delete it would hand the generator one extra "cut row"
 * to dualize with a garbage dual.
 *
 * Names are regenerated so that the `.bas` written next matches them. A `.bas`
 * carrying names the `.mps` does not have is skipped by readBasis *and returns
 * success*, leaving the replay to start from a cold LP -- which for a
 * needsOptimalBasis generator means no cuts at all rather than merely a different
 * vertex.
 */
static int cglLagomoryFixtureWriteMps(const OsiSolverInterface &si,
  const std::string &path, int &paddedColumns)
{
  OsiSolverInterface *clone = si.clone();
  paddedColumns = 0;

  const CoinPackedMatrix *byCol = clone->getMatrixByCol();
  std::vector< int > empties;
  if (byCol) {
    const int *len = byCol->getVectorLengths();
    for (int j = 0; j < clone->getNumCols(); ++j) {
      if (len[j] == 0)
        empties.push_back(j);
    }
  }
  paddedColumns = (int)empties.size();

  if (!empties.empty()) {
    // A row that cannot cut anything off: its bounds are the interval its own
    // activity can take, so it is redundant by construction. Deliberately not
    // truly free -- a doubly-infinite row would be dropped by writeMps and would
    // take the column indices with it.
    double padLb = 0.0, padUb = 0.0;
    const double *cl = clone->getColLower();
    const double *cu = clone->getColUpper();
    bool finite = true;
    for (size_t k = 0; k < empties.size(); ++k) {
      const double lo = cl[empties[k]], up = cu[empties[k]];
      if (lo <= -1.0e30 || up >= 1.0e30) {
        finite = false;
        break;
      }
      padLb += lo < 0.0 ? lo : 0.0;
      padUb += up > 0.0 ? up : 0.0;
    }
    if (!finite) {
      padLb = -1.0e29;
      padUb = 1.0e29;
    }
    const std::vector< double > coefs(empties.size(), 1.0);
    clone->addRow((int)empties.size(), &empties[0], &coefs[0], padLb, padUb);
  }

  const int nr = clone->getNumRows(), nc = clone->getNumCols();
  std::vector< std::string > rowStore, colStore;
  std::vector< const char * > rowPtrs, colPtrs;
  rowStore.reserve(nr);
  colStore.reserve(nc);
  for (int i = 0; i < nr; ++i)
    rowStore.push_back(clone->getRowName(i));
  for (int j = 0; j < nc; ++j)
    colStore.push_back(clone->getColName(j));
  rowPtrs.reserve(nr);
  colPtrs.reserve(nc);
  for (int i = 0; i < nr; ++i)
    rowPtrs.push_back(rowStore[i].c_str());
  for (int j = 0; j < nc; ++j)
    colPtrs.push_back(colStore[j].c_str());

  const int rc = clone->writeMpsNative(path.c_str(),
    rowPtrs.empty() ? NULL : &rowPtrs[0],
    colPtrs.empty() ? NULL : &colPtrs[0], 2, 1);
  delete clone;
  return rc;
}

/**
 * Write the basis. OsiClp's writeBasisNative() emits FREEIEEE, which round-trips
 * doubles exactly and matches the `.bas` files of preProcessedInstances so the two
 * sets stay interchangeable. The base-class implementation is a no-op that still
 * returns success, so success is confirmed by the file existing.
 */
static bool cglLagomoryFixtureWriteBasis(const OsiSolverInterface &si,
  const std::string &path)
{
  OsiSolverInterface *clone = si.clone();
  clone->writeBasisNative(path.c_str());
  delete clone;
  FILE *probe = fopen(path.c_str(), "r");
  if (probe) {
    fclose(probe);
    return true;
  }
  return false;
}

/// LP solution, in the shape of preProcessedInstances/*.sol: an "=obj=" line then
/// one line per nonzero as "<index> <name> <value>".
static bool cglLagomoryFixtureWriteSol(const OsiSolverInterface &si,
  const std::string &path)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;
  fprintf(fp, "=obj= %.15g\n", si.getObjValue());
  const double *x = si.getColSolution();
  if (x) {
    const int n = si.getNumCols();
    for (int j = 0; j < n; ++j) {
      if (x[j] == 0.0)
        continue;
      fprintf(fp, "%5d %-24s %.15g\n", j, si.getColName(j).c_str(), x[j]);
    }
  }
  fclose(fp);
  return true;
}

/**
 * Row duals. Every row is written, zeros included, because the file is indexed by
 * row and a sparse form would need the row count to be trusted from elsewhere.
 * `formulation` is recorded in the header line so a reader can tell at a glance
 * how much of the vector is the part that actually perturbs the objective: only
 * rows >= numberOriginalRows are dualized.
 */
static bool cglLagomoryFixtureWritePi(const OsiSolverInterface &si,
  int numberOriginalRows, const std::string &path)
{
  const double *pi = si.getRowPrice();
  if (!pi)
    return false;
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;
  const int m = si.getNumRows();
  fprintf(fp, "pi %d formulation %d\n", m, numberOriginalRows);
  for (int i = 0; i < m; ++i)
    fprintf(fp, "%d %.15g\n", i, pi[i]);
  fclose(fp);
  return true;
}

/**
 * Column types, which the MPS cannot carry for fixed integers.
 * `isContinuous()` is the source of truth rather than `getColType()`, which
 * derives Binary/GeneralInteger from the current bounds.
 */
static bool cglLagomoryFixtureWriteColTypes(const OsiSolverInterface &si,
  const std::string &path, int &integerColumns)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;
  const int n = si.getNumCols();
  const char *ct = si.getColType(true);
  integerColumns = 0;
  fprintf(fp, "cols %d\n", n);
  for (int j = 0; j < n; ++j) {
    if (si.isContinuous(j))
      continue;
    ++integerColumns;
    fprintf(fp, "%d %d\n", j, (int)ct[j]);
  }
  fclose(fp);
  return true;
}

/// CglGomory's own fractionality test (CglGomory.cpp:363-369), reproduced so the
/// candidate count is the generator's and not an approximation of it. The 1e-9
/// relative snap matters: a column at 3.9999999997 is *not* a candidate.
static double cglLagomoryFixtureAboveInteger(double value)
{
  const double value2 = floor(value);
  const double value3 = floor(value + 0.5);
  if (fabs(value3 - value) < 1.0e-9 * (fabs(value3) + 1.0))
    return 0.0;
  return value - value2;
}

/**
 * Characterise the Lagrangean work this call represents, so fixtures can be ranked
 * without replaying them. Every counter here is a *cause* of cost or of cut yield
 * in the wrapper, not in the shared Gomory core:
 *
 *   cutRows          rows to dualize or copy = numberRows - numberOriginalRows.
 *                    Zero means the wrapper re-solves an unchanged objective and
 *                    reproduces plain Gomory exactly (which is why whenToDo==2
 *                    with no cut rows is duplicated work, and why the existing
 *                    fixtures are useless here).
 *   cutRowNz         elements in those rows: the exact cost of the dualization
 *                    loop at :182-202.
 *   dualizedRows     cut rows with a nonzero dual, i.e. the ones that actually
 *                    move the objective. A cut row at its bound with pi == 0
 *                    contributes nothing and is pure loop overhead.
 *   integralCutRows  cut rows whose every coefficient is integral -- the
 *                    gomoryType==2 "copy instead of dualize" population
 *                    (:149-177). Also the addRows/deleteRows churn per call.
 *   piNorm           sum |pi| over cut rows: how far the perturbed objective is
 *                    from the real one, hence how different a vertex the resolve
 *                    is likely to land on.
 *   fractionalInts   fractional integer columns in si's solution. Not the Gomory
 *                    candidate count (that is taken at the *perturbed* vertex,
 *                    which does not exist yet), but it bounds how much there is
 *                    to cut at all.
 */
static void cglLagomoryFixtureCount(const OsiSolverInterface &si,
  int numberOriginalRows, int &cutRows, double &cutRowNz, int &dualizedRows,
  int &integralCutRows, double &piNorm, int &fractionalInts)
{
  const int m = si.getNumRows();
  cutRows = m - numberOriginalRows;
  cutRowNz = 0.0;
  dualizedRows = 0;
  integralCutRows = 0;
  piNorm = 0.0;
  fractionalInts = 0;

  const double *pi = si.getRowPrice();
  const CoinPackedMatrix *byRow = si.getMatrixByRow();
  if (byRow && cutRows > 0) {
    const CoinBigIndex *rowStart = byRow->getVectorStarts();
    const int *rowLength = byRow->getVectorLengths();
    const double *element = byRow->getElements();
    for (int i = numberOriginalRows; i < m; ++i) {
      cutRowNz += (double)rowLength[i];
      if (pi && pi[i] != 0.0) {
        ++dualizedRows;
        piNorm += fabs(pi[i]);
      }
      bool allIntegral = true;
      for (CoinBigIndex k = rowStart[i]; k < rowStart[i] + rowLength[i]; ++k) {
        const double v = element[k];
        if (v != floor(v + 0.5)) {
          allIntegral = false;
          break;
        }
      }
      if (allIntegral)
        ++integralCutRows;
    }
  }

  const double *x = si.getColSolution();
  const char *intInfo = si.getColType();
  const double *cl = si.getColLower();
  const double *cu = si.getColUpper();
  if (x && intInfo) {
    const int n = si.getNumCols();
    for (int j = 0; j < n; ++j) {
      if (!intInfo[j] || cu[j] <= cl[j] + 0.5)
        continue;
      if (cglLagomoryFixtureAboveInteger(x[j]) > 0.0)
        ++fractionalInts;
    }
  }
}

/**
 * The `.meta`. Every field is either an input to the replay or a ranking key, and
 * the ones whose derivation is not obvious are commented at their write site.
 */
static bool cglLagomoryFixtureWriteMeta(const OsiSolverInterface &si,
  const OsiSolverInterface &orig, const std::string &path, int numberOriginalRows,
  int gomoryType, int pass, int options, int inTree, int numberTimesStalled,
  int limit, int limitAtRoot, int dynamicLimitInTree, double away,
  double awayAtRoot, double conditionNumberMultiplier,
  double largestFactorMultiplier, int alternateFactorization, int paddedColumns,
  int integerColumns, int origPaddedColumns, int origIntegerColumns,
  int callNumber)
{
  FILE *fp = fopen(path.c_str(), "w");
  if (!fp)
    return false;

  fprintf(fp, "tag lagomory\n");
  fprintf(fp, "rows %d\n", si.getNumRows());
  fprintf(fp, "cols %d\n", si.getNumCols());
  fprintf(fp, "elements %d\n", si.getNumElements());
  fprintf(fp, "lpOptimal %d\n", (int)si.isProvenOptimal());
  fprintf(fp, "objSense %g\n", si.getObjSense());
  fprintf(fp, "objValue %.15g\n", si.getObjValue());
  fprintf(fp, "paddedColumns %d\n", paddedColumns);
  fprintf(fp, "integerColumns %d\n", integerColumns);

  // The formulation/cut split. This is THE number that makes the fixture a
  // lagomory fixture rather than a Gomory one: rows [0,formulationRows) are
  // dualized-from, rows [formulationRows,rows) are dualized-away.
  fprintf(fp, "formulationRows %d\n", numberOriginalRows);

  fprintf(fp, "origRows %d\n", orig.getNumRows());
  fprintf(fp, "origCols %d\n", orig.getNumCols());
  fprintf(fp, "origElements %d\n", orig.getNumElements());
  fprintf(fp, "origPaddedColumns %d\n", origPaddedColumns);
  fprintf(fp, "origIntegerColumns %d\n", origIntegerColumns);

  // Whether si's formulation rows still agree with orig's. If this is nonzero the
  // reconstruct-by-deleting-rows shortcut would have produced the wrong LP, which
  // is reason 2 in the header for dumping orig separately. Cheap to compute and
  // the only evidence that the decision was necessary rather than defensive.
  int rowBoundDrift = 0;
  {
    const int mm = numberOriginalRows < orig.getNumRows()
      ? numberOriginalRows : orig.getNumRows();
    const double *sl = si.getRowLower(), *su = si.getRowUpper();
    const double *ol = orig.getRowLower(), *ou = orig.getRowUpper();
    for (int i = 0; i < mm; ++i) {
      if (fabs(sl[i] - ol[i]) > 1.0e-12 || fabs(su[i] - ou[i]) > 1.0e-12)
        ++rowBoundDrift;
    }
  }
  fprintf(fp, "rowBoundDrift %d\n", rowBoundDrift);

  int cutRows, dualizedRows, integralCutRows, fractionalInts;
  double cutRowNz, piNorm;
  cglLagomoryFixtureCount(si, numberOriginalRows, cutRows, cutRowNz, dualizedRows,
    integralCutRows, piNorm, fractionalInts);
  fprintf(fp, "cutRows %d\n", cutRows);
  fprintf(fp, "cutRowNz %.0f\n", cutRowNz);
  fprintf(fp, "dualizedRows %d\n", dualizedRows);
  fprintf(fp, "integralCutRows %d\n", integralCutRows);
  fprintf(fp, "piNorm %.15g\n", piNorm);
  fprintf(fp, "fractionalInts %d\n", fractionalInts);

  // The call's own configuration, read off the generator rather than re-derived
  // from CBC's setup code -- these are the values the dumping run actually used.
  // gomoryType is the composite: gomoryType%10 is 1 (dualize all) or 2 (copy the
  // all-integral cut rows in as real rows), gomoryType/10 is 1 (only when cut rows
  // exist) or 2 (always).
  fprintf(fp, "gomoryType %d\n", gomoryType);
  fprintf(fp, "numberTimesStalled %d\n", numberTimesStalled);
  fprintf(fp, "limit %d\n", limit);
  fprintf(fp, "limitAtRoot %d\n", limitAtRoot);
  fprintf(fp, "dynamicLimitInTree %d\n", dynamicLimitInTree);
  fprintf(fp, "away %.15g\n", away);
  fprintf(fp, "awayAtRoot %.15g\n", awayAtRoot);
  fprintf(fp, "conditionNumberMultiplier %.15g\n", conditionNumberMultiplier);
  fprintf(fp, "largestFactorMultiplier %.15g\n", largestFactorMultiplier);
  fprintf(fp, "alternateFactorization %d\n", alternateFactorization);

  // The CglTreeInfo the call was made with. Not decoration: bit 512 (parentModel)
  // and bit 1024 (must-call-again) appear in the whenToDo gate itself, bit 16
  // (globalCuts) drives the setGloballyValid arm at :335-344 that has no inTree
  // guard, and infoPass selects an absolute-vs-relative accuracy test in the core.
  fprintf(fp, "infoPass %d\n", pass);
  fprintf(fp, "infoOptions %d\n", options);
  fprintf(fp, "infoInTree %d\n", inTree);

  fprintf(fp, "callNumber %d\n", callNumber);
  fclose(fp);
  return true;
}

/**
 * Dump one Lagrangean Gomory call. Call from inside the `whenToDo` gate, *before*
 * the column bounds and objective of `orig` are overwritten.
 *
 * Returns true if a fixture was written. Refuses (returns false) when
 * `numberRows == numberOriginalRows`: such a call has nothing to dualize and
 * reproduces plain Gomory exactly, which is the very degeneracy that made the
 * existing 330 fixtures useless -- writing more of them would be writing the same
 * mistake to a new directory.
 */
static bool cglDumpLagomoryFixture(const OsiSolverInterface &si,
  const OsiSolverInterface &orig, int numberOriginalRows, int gomoryType,
  int pass, int options, int inTree, int numberTimesStalled, int limit,
  int limitAtRoot, int dynamicLimitInTree, double away, double awayAtRoot,
  double conditionNumberMultiplier, double largestFactorMultiplier,
  int alternateFactorization)
{
  static int callNumber = 0;
  static int written = 0;
  ++callNumber;

  int budget = 1;
  {
    const char *env = getenv("CGL_LAGOMORY_FIXTURE_MAX");
    if (env && *env)
      budget = atoi(env);
  }
  if (written >= budget)
    return false;
  if (si.getNumRows() <= numberOriginalRows)
    return false;

  const std::string dir = cglLagomoryFixtureDir();
  cglLagomoryFixtureMkdirP(dir);
  const std::string stem = dir + "/" + cglLagomoryFixtureName(si) + ".lagomory";

  int paddedColumns = 0, origPaddedColumns = 0;
  // The ".mps.gz" suffix is what lands on disk: writeMpsNative always gzips and
  // appends the extension itself when told to compress.
  if (cglLagomoryFixtureWriteMps(si, stem + ".mps", paddedColumns) != 0) {
    fprintf(stderr, "[lagomory-fixture] writeMps failed for %s\n", stem.c_str());
    return false;
  }
  if (cglLagomoryFixtureWriteMps(orig, stem + ".orig.mps", origPaddedColumns) != 0) {
    fprintf(stderr, "[lagomory-fixture] writeMps(orig) failed for %s\n", stem.c_str());
    return false;
  }
  if (!cglLagomoryFixtureWriteBasis(si, stem + ".bas")) {
    fprintf(stderr, "[lagomory-fixture] no basis written for %s -- a replay of a\n"
                    "  needsOptimalBasis generator without one produces no cuts, so\n"
                    "  the fixture is discarded rather than left misleading.\n",
      stem.c_str());
    return false;
  }
  cglLagomoryFixtureWriteSol(si, stem + ".sol");
  cglLagomoryFixtureWritePi(si, numberOriginalRows, stem + ".pi");

  int integerColumns = 0, origIntegerColumns = 0;
  cglLagomoryFixtureWriteColTypes(si, stem + ".ctype", integerColumns);
  cglLagomoryFixtureWriteColTypes(orig, stem + ".orig.ctype", origIntegerColumns);

  cglLagomoryFixtureWriteMeta(si, orig, stem + ".meta", numberOriginalRows,
    gomoryType, pass, options, inTree, numberTimesStalled, limit, limitAtRoot,
    dynamicLimitInTree, away, awayAtRoot, conditionNumberMultiplier,
    largestFactorMultiplier, alternateFactorization, paddedColumns,
    integerColumns, origPaddedColumns, origIntegerColumns, callNumber);

  ++written;
  fprintf(stderr, "[lagomory-fixture] wrote %s (call %d: %d rows, %d formulation,"
                  " %d cut rows, gomoryType %d)\n",
    stem.c_str(), callNumber, si.getNumRows(), numberOriginalRows,
    si.getNumRows() - numberOriginalRows, gomoryType);
  return true;
}

#endif /* CGL_DUMP_LAGOMORY_FIXTURE */
#endif /* CglLagomoryFixtureDump_H */
