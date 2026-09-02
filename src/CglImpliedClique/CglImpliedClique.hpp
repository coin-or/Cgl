/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class responsible for separating "implied clique" cuts: it roots a
 * clique search at the complement of every binary variable y (i.e. the
 * literal "y = 0") and greedily grows it with conflict-graph neighbours,
 * disaggregating and strengthening rows shaped like
 *
 *     x1 + x2 + ... + xk <= M*y      (all variables binary)
 *
 * -- normally used to model "x1 OR x2 OR ... -> y" -- into the much
 * stronger clique-style inequality
 *
 *     x1 + x2 + ... + xk - y <= k - 1     (i.e. every xi <= y disaggregated,
 *                                          plus any further conflicting
 *                                          literals found via the conflict
 *                                          graph, including complemented
 *                                          ones, e.g. (1 - xj) <= y).
 *
 * Unlike CglBKClique (which runs a general Bron-Kerbosch search over the
 * whole fractional-vertex induced subgraph once per round), this generator
 * is deliberately narrow and cheap: one hub-rooted greedy extension per
 * binary variable, using only the already-built conflict graph -- no row
 * parsing is needed, since the "x1+...+xk <= M*y" structure is exactly
 * what the conflict graph already encodes as pairwise (xi, y=0) conflicts.
 *
 * Benchmarking (Cbc/test/impclique-bench.cpp and
 * Cbc/test/impclique-vs-clique-bench.cpp, replayed against 761
 * post-preprocessing fixtures from mip-sanity-data and MIPLIB 2017+spp)
 * found that most of the same cliques are eventually rediscovered by
 * CglBKClique given enough rounds -- but (a) a small number of instances
 * (e.g. bab6, irish-electricity) show a large, durable bound improvement
 * that CglBKClique's round-based search does not reach even after many
 * rounds, and (b) several more instances (e.g. unitcal_7) reach the same
 * final bound but reach it strictly sooner when this generator runs
 * alongside CglBKClique -- valuable in itself since CBC's cut loop and
 * node budget both reward faster bound convergence. Separation cost is
 * negligible in all tested cases (sub-millisecond per LP even on the
 * largest post-preprocessing fixtures), so it is designed to run on every
 * round without an opt-in switch, mirroring CglBKClique's own defaults.
 *
 * @file CglImpliedClique.hpp
 * @brief Implied-clique disaggregation cut separator
 *
 **/

#ifndef _CglImpliedClique_h_
#define _CglImpliedClique_h_

#include <CglCutGenerator.hpp>
#include <atomic>
#include <vector>

class CoinConflictGraph;

/**
 * Class responsible for separating "implied clique" cuts, rooted at the
 * complement of a binary variable and grown greedily via the conflict
 * graph. See the file header for the underlying idea and rationale.
 **/
class CGLLIB_EXPORT CglImpliedClique : public CglCutGenerator {
public:
  /**
   * Default constructor
   **/
  CglImpliedClique();

  /**
   * Copy constructor
   **/
  CglImpliedClique(const CglImpliedClique &rhs);

  /**
   * Clone
   **/
  virtual CglCutGenerator *clone() const;

  /**
   * Generate implied-clique cuts for the model data contained in si. The
   * generated cuts are inserted into and returned in the collection of
   * cuts cs.
   **/
  virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info = CglTreeInfo());

  /**
   * Destructor
   **/
  virtual ~CglImpliedClique();

  /**
   * Refresh the conflict graph if necessary.
   **/
  virtual void refreshSolver(OsiSolverInterface *solver);

  /**
   * Set whether complemented-literal candidates (i.e. "col = 0", meaning
   * (1 - col)) may join a clique alongside plain "col = 1" literals.
   * Benchmarking found this never hurts and sometimes finds cuts that are
   * otherwise entirely missed (e.g. gesa2, trdcrooms), so it defaults on.
   **/
  void setMixedLiterals(bool mixedLiterals) { mixedLiterals_ = mixedLiterals; }
  bool getMixedLiterals() const { return mixedLiterals_; }

  /**
   * Disable clique growth: only the k=1 disaggregation (x_i <= y for each
   * conflicting candidate) is separated, no greedy extension. Provided as
   * a diagnostic control; benchmarking found extension is sometimes
   * strictly necessary for the cut to appear at all (e.g. momentum1), so
   * this defaults off.
   **/
  void setSingletonOnly(bool singletonOnly) { singletonOnly_ = singletonOnly; }
  bool getSingletonOnly() const { return singletonOnly_; }

  /**
   * Minimum violation (lhs - rhs, in the normalized "hub" scale, i.e.
   * comparable to a fraction of one unit of y) for a cut to be kept.
   **/
  void setMinViol(double minViol) { minViol_ = minViol; }
  double getMinViol() const { return minViol_; }

  /**
   * Cap on the number of conflict-graph node visits (across all hubs) in
   * a single generateCuts() call, to bound worst-case cost on pathological
   * dense graphs. 0 = unlimited (default).
   **/
  void setMaxNodeVisits(size_t maxNodeVisits) { maxNodeVisits_ = maxNodeVisits; }
  size_t getMaxNodeVisits() const { return maxNodeVisits_; }

  /**
   * Number of cuts separated.
   **/
  static std::atomic< size_t > sepCuts_;

  /**
   * Execution time spent for the implied-clique cut separator.
   **/
  static std::atomic< double > sepTime_;

private:
  /**
   * Whether to consider complemented-literal candidates.
   **/
  bool mixedLiterals_;

  /**
   * Whether to disable clique growth (k=1 disaggregation only).
   **/
  bool singletonOnly_;

  /**
   * Minimum violation of a cut to be kept.
   **/
  double minViol_;

  /**
   * Node-visit budget for a single generateCuts() call. 0 = unlimited.
   **/
  size_t maxNodeVisits_;
};

#endif // CglImpliedClique_HPP
