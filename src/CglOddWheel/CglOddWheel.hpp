/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class for separating violated odd-cycles. It contains
 * a lifting module that tries to transform the odd-cycles
 * into odd-wheels.
 *
 * @file CglOddWheel.hpp
 * @brief Odd-wheel cut separator
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo.santos@gmail.com
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef _CglOddWheel_h_
#define _CglOddWheel_h_

#include "CglCutGenerator.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinOddWheelSeparator.hpp"

class CoinConflictGraph;

/**
 * Class for separating violated odd-cycles. It contains
 * a lifting module that tries to transform the odd-cycles
 * into odd-wheels.
 **/
class CGLLIB_EXPORT CglOddWheel : public CglCutGenerator
{
public:
  /**
   * Number of cuts separated.
   **/
  static size_t sepCuts;

  /**
   * Execution time spent for the clique
   * cut separator.
   **/
  static double sepTime;

  /**
   * Default constructor
   *
   * @param extMethod strategy that will be used to lift odd cycles,
   * transforming them into odd wheels: 0 = no lifting, 1 = only one
   * variable as wheel center, 2 = a clique as wheel center.
   **/
  CglOddWheel(size_t extMethod = 2);

  /**
   * Copy constructor
   **/
  CglOddWheel(const CglOddWheel& rhs);

  /**
   * Clone
   **/
  virtual CglCutGenerator * clone() const;

  /**
   * Generate clique cuts for the model data contained
   * in si. The generated cuts are inserted into and returned
   * in the collection of cuts cs.
   **/
  virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info = CglTreeInfo() );

  /**
   * Destructor
   **/
  virtual ~CglOddWheel();

  /**
   * Refresh the conflict graph if necessary.
   **/
  virtual void refreshSolver(OsiSolverInterface *solver);

  /**
   * Set the strategy that will be used to lift odd cycles,
   * transforming them into odd wheels: 0 = no lifting, 1 = only one
   * variable as wheel center, 2 = a clique as wheel center.
   **/
  void setExtendingMethod(size_t extMethod);

  /**
   * Return the strategy used to lift odd cycles.
   **/
  size_t getExtendingMethod() const { return extMethod_; }

  /**
   * Have the separator cross-check its two ways of building the auxiliary
   * graph's arcs, reporting the result in Stats::sep.prepareMismatches. Off by
   * default; see CoinOddWheelSeparator::setVerifyPrepare(). Diagnostic only --
   * it roughly doubles graph preparation time.
   **/
  void setVerifyPrepare(bool verify) { verifyPrepare_ = verify; }

  /**
   * Enable (default) or disable the separator's futility gate, which proves per
   * active node that no shortest-path call from it can yield a cut and skips the
   * call. See CoinOddWheelSeparator::buildFutilityGate(). Turning it off exists
   * so the "no cut lost" claim can be *checked* by replaying a fixture both ways
   * and comparing every reported field, not just asserted from the proof.
   **/
  void setUseGate(bool use) { useGate_ = use; }

  /**
   * Certify every odd wheel against the conflict graph before it is translated
   * into a row cut, reporting the outcome in the Stats::cert* counters. Off by
   * default.
   *
   * This is a *self-contained* validity proof, which is what makes it worth
   * having next to Osi's row-cut debugger: the debugger needs a known feasible
   * solution, and on this fixture set most reference files do not supply one --
   * 9 of 32 hold an LP relaxation point, which every valid cut cuts off. The
   * certificate needs nothing but the graph.
   *
   * What it proves. The node-space inequality is
   *   sum_{v in C} z_v + alpha * sum_{w in W} z_w <= k,   k = floor(|C|/2),
   * and it is valid for every 0/1 point respecting the conflict graph as soon
   * as: |C| is odd, consecutive nodes of C conflict (so C is an odd cycle and
   * at most k of it can be 1), the nodes of C are distinct, every centre in W
   * conflicts with every node of C, W is a clique, and alpha <= k. The two
   * cases are then immediate -- some w = 1 forces all of C and the rest of W to
   * 0, giving alpha <= k; otherwise the cycle bound applies. The final check
   * re-derives the column-space cut from (C, W, alpha) independently and
   * compares it with the one actually emitted, so the complement translation is
   * covered too.
   *
   * Cost is O(|C|^2 + |W|^2 + |C||W|) per wheel, no enumeration.
   **/
  void setCheckValidity(bool check) { checkValidity_ = check; }

  /**
   * Counters and per-stage times of the last generateCuts() call.
   * Unlike the static sepCuts/sepTime totals above these are per call,
   * which is what a profiling harness needs.
   **/
  struct Stats {
    CoinOddWheelSeparator::Stats sep; /**< the separator's own counters and stage times */
    size_t cutsBeforePool;            /**< odd wheels handed to the cut pool */
    size_t cutsDuplicatedIdx;         /**< odd wheels where a column appeared twice and was merged */
    size_t cutsZeroCoefs;             /**< coefficients that cancelled to zero while merging */
    size_t cutsEmpty;                 /**< odd wheels dropped: every coefficient cancelled */
    size_t cutsAfterPool;             /**< survivors of the cut pool's dominance filter */
    size_t rowCutsAdded;              /**< row cuts actually inserted into cs */
    double tSetup;                    /**< the doubled x_/rc_ arrays */
    double tSeparator;                /**< separator construction plus searchOddWheels() */
    double tCutPool;                  /**< index translation, cut pool, insertion into cs */

    /* setCheckValidity(): the first is how many wheels were certified, the
     * next five must all be 0, and the last three measure how much of the win
     * comes from complemented nodes -- the reason the doubled graph exists. */
    size_t certChecked;               /**< odd wheels put through the certificate */
    size_t certBadCycle;              /**< FAILS: |C| not odd, a node repeated, or consecutive nodes not in conflict */
    size_t certBadCenterAdj;          /**< FAILS: a wheel centre does not conflict with every node of C */
    size_t certBadCenterClq;          /**< FAILS: the wheel centres are not pairwise in conflict */
    size_t certBadAlpha;              /**< FAILS: the centre coefficient exceeds floor(|C|/2) */
    size_t certBadTranslate;          /**< FAILS: the emitted cut is not the node-space one re-derived */
    size_t certComplCycle;            /**< wheels whose cycle uses at least one complemented node */
    size_t certComplCenter;           /**< wheels whose centre uses at least one complemented node */
    size_t certComplPair;             /**< wheels where some column appears both plain and complemented */

    /* How often the choice of alpha can matter at all. alpha must be the
     * cycle's own floor(|C|/2); deriving it from the already-translated rhs
     * instead subtracts one per complemented cycle node, which is only
     * observable on a wheel that has a centre *and* a complemented cycle. The
     * second counter is where that arithmetic went non-positive, i.e. where the
     * lifting silently stopped strengthening anything. */
    size_t certCenterOnComplCycle;    /**< wheels with a centre whose cycle traverses a complemented node */
    size_t certComplAtLeastK;         /**< wheels whose cycle has at least floor(|C|/2) complemented nodes */
  };

  /**
   * Statistics of the last generateCuts() call.
   **/
  inline const Stats &stats() const { return stats_; }

private:
  /**
   * Check if it is necessary realloc the memory
   * for the data structures.
   **/
  void checkMemory(const size_t newNumCols);

  /**
   * Prove one odd wheel valid against the conflict graph and check that the cut
   * about to be emitted is the one it implies. See setCheckValidity() for what
   * is proven; the outcome lands in the Stats::cert* counters.
   *
   * @param cgraph the doubled conflict graph the wheel was found in
   * @param numCols columns of the model (so a node >= numCols is a complement)
   * @param cycle nodes of the odd cycle C, in cycle order
   * @param cycleSize |C|
   * @param center nodes of the wheel centre W (may be empty)
   * @param centerSize |W|
   * @param alpha coefficient actually given to the centres, 0 if none were used
   * @param idxs columns of the emitted cut
   * @param coefs coefficients of the emitted cut
   * @param nz length of idxs/coefs
   * @param rhs right-hand side of the emitted cut
   **/
  void certifyOddWheel(const CoinConflictGraph *cgraph, size_t numCols,
    const size_t *cycle, size_t cycleSize,
    const size_t *center, size_t centerSize, double alpha,
    const int *idxs, const double *coefs, int nz, double rhs);

  /**
   * Capacity of storage of the data structures.
   **/
  size_t cap_;

  /**
   * Auxiliary arrays used to store the indexes
   * of a clique.
   **/
  int *idxs_, *idxMap_;

  /**
   * Auxiliary array used to store the coefficients
   * of a clique.
   **/
  double *coefs_;

  /**
   * Current solution of the LP relaxation
   **/
  double *x_;

  /**
   * Current reduced costs of the variables
   **/
  double *rc_;

  /**
   * Auxiliary structure used to temporary
   * store a cut.
   **/
  OsiRowCut osrc_;

  /**
   * Lifting strategy: 0 = no lifting,
   * 1 = only one variable as wheel center,
   * 2 = a clique as wheel center
   **/
  size_t extMethod_;

  /**
   * Whether the separator cross-checks its arc construction,
   * see setVerifyPrepare().
   **/
  bool verifyPrepare_;

  /**
   * Whether the separator's futility gate runs, see setUseGate().
   **/
  bool useGate_;

  /**
   * Whether every odd wheel is certified against the conflict graph
   * before being emitted, see setCheckValidity().
   **/
  bool checkValidity_;

  /**
   * Counters and per-stage times, see stats().
   **/
  Stats stats_;
};

#endif
