/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * @file CglImpliedClique.cpp
 * @brief Implied-clique disaggregation cut separator
 *
 **/

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>
#include <OsiCuts.hpp>
#include <OsiRowCut.hpp>
#include <CoinTime.hpp>

#include "CglImpliedClique.hpp"
#include "CoinConflictGraph.hpp"
#include "CoinCutPool.hpp"

#define ICLQ_EPS 1e-6

std::atomic< size_t > CglImpliedClique::sepCuts_(0);
std::atomic< double > CglImpliedClique::sepTime_(0.0);

CglImpliedClique::CglImpliedClique()
  : mixedLiterals_(true)
  , singletonOnly_(false)
  , minViol_(0.02)
  , maxNodeVisits_(0)
{
}

CglImpliedClique::CglImpliedClique(const CglImpliedClique &rhs)
  : CglCutGenerator(rhs)
  , mixedLiterals_(rhs.mixedLiterals_)
  , singletonOnly_(rhs.singletonOnly_)
  , minViol_(rhs.minViol_)
  , maxNodeVisits_(rhs.maxNodeVisits_)
{
}

CglImpliedClique::~CglImpliedClique()
{
}

CglCutGenerator *CglImpliedClique::clone() const
{
  return new CglImpliedClique(*this);
}

void CglImpliedClique::refreshSolver(OsiSolverInterface *solver)
{
  solver->checkCGraph();
  solver->getColType(true);
}

namespace {

/// A clique-cut member is a *literal*, not just a column: node is a raw
/// conflict-graph node id (col for the "col=1" literal, col+numCols for
/// the "col=0" literal, i.e. the complement of col).
struct Lit {
  int col;
  bool neg;
};

/// Greedily grow a clique among candidates (already known to all conflict
/// with the hub's complement literal), ordered by descending LP weight.
/// O(|candidates|^2) conflict checks, same routine used for exploration
/// and benchmarking (Cbc/test/impclique-bench.cpp, impclique-explore.cpp).
std::vector< int > greedyClique(const CoinConflictGraph *cg,
  const std::vector< std::pair< double, int > > &sortedCandidates)
{
  std::vector< int > clique;
  clique.reserve(sortedCandidates.size());
  for (const auto &wc : sortedCandidates) {
    const int cand = wc.second;
    bool ok = true;
    for (int m : clique) {
      if (!cg->conflicting(static_cast< size_t >(cand), static_cast< size_t >(m))) {
        ok = false;
        break;
      }
    }
    if (ok)
      clique.push_back(cand);
  }
  return clique;
}

} // namespace

void CglImpliedClique::generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info)
{
  const int numCols = si.getNumCols();
  if (numCols == 0 || si.getNumRows() == 0) {
    return;
  }

  const double startSep = CoinCpuTime();
  const CoinConflictGraph *cg = si.getCGraph();
  if (!cg || cg->size() == 0) {
    return;
  }

  if (static_cast< size_t >(numCols) != cg->size() / 2) {
    fprintf(stderr, "Invalid conflict graph! Number of columns %d ... in graph %lu\n",
      numCols, (unsigned long)(cg->size() / 2));
    abort();
  }

  const double *sol = si.getColSolution();
  const char *colType = si.getColType(true);

  std::vector< size_t > temp(cg->size());
  std::vector< char > iv(cg->size(), 0);
  size_t nodeVisits = 0;

  // Reused per-hub scratch buffers for cut assembly, sized to the widest
  // possible clique (every column, positive or negated literal).
  std::vector< int > idxs(numCols);
  std::vector< double > coefs(numCols);
  std::vector< int > idxMap(numCols, -1);

  CoinCutPool cutpool(sol, numCols);

  for (int y = 0; y < numCols; y++) {
    if (colType[y] == 0)
      continue;
    const double yVal = sol[y];
    if (yVal > 1.0 - ICLQ_EPS)
      continue;

    const size_t yComp = static_cast< size_t >(y) + numCols;
    if (yComp >= cg->size())
      continue;

    if (maxNodeVisits_ > 0 && nodeVisits >= maxNodeVisits_)
      break;

    const auto conf = cg->conflictingNodes(yComp, temp.data(), iv.data());
    nodeVisits += conf.first;

    // Weighted-by-LP-value candidates: positive-literal candidates are the
    // column node itself; a complemented-literal candidate ("col=0", i.e.
    // (1-col)) is the col+numCols node -- only considered when
    // mixedLiterals_ is set.
    std::vector< std::pair< double, int > > candidates;
    for (size_t k = 0; k < conf.first; k++) {
      const size_t node = conf.second[k];
      int col;
      bool neg;
      if (node < static_cast< size_t >(numCols)) {
        col = static_cast< int >(node);
        neg = false;
      } else if (mixedLiterals_ && node < static_cast< size_t >(2 * numCols)) {
        col = static_cast< int >(node - static_cast< size_t >(numCols));
        neg = true;
      } else {
        continue;
      }
      if (col == y)
        continue;
      if (colType[col] == 0)
        continue;
      const double litVal = neg ? (1.0 - sol[col]) : sol[col];
      if (litVal > ICLQ_EPS)
        candidates.push_back(std::make_pair(litVal, static_cast< int >(node)));
    }
    if (candidates.empty())
      continue;

    std::sort(candidates.begin(), candidates.end(),
      [](const std::pair< double, int > &a, const std::pair< double, int > &b) {
        return a.first > b.first;
      });

    std::vector< int > cliqueNodes;
    if (singletonOnly_) {
      cliqueNodes.push_back(candidates.front().second);
    } else {
      cliqueNodes = greedyClique(cg, candidates);
    }

    double lhs = 0.0;
    for (int node : cliqueNodes) {
      const bool neg = node >= numCols;
      const int col = neg ? node - numCols : node;
      lhs += neg ? (1.0 - sol[col]) : sol[col];
    }

    if (lhs - yVal <= minViol_)
      continue;

    // Assemble the row cut: sum_pos x_i + sum_neg (1-x_j) <= y  <=>
    // sum_pos x_i - sum_neg x_j - y <= -negCount. A column can in
    // principle appear via both its positive and negated literal (they
    // trivially conflict with each other too), in which case their
    // coefficients cancel; mirror CglBKClique::insertCuts's de-duplication
    // so that degenerate case is handled the same way everywhere.
    int cutSize = 0;
    int rhsAdj = 0;
    int duplicated = 0;
    for (int node : cliqueNodes) {
      const bool neg = node >= numCols;
      const int col = neg ? node - numCols : node;
      if (idxMap[col] == -1) {
        idxMap[col] = cutSize;
        idxs[cutSize] = col;
        coefs[cutSize] = neg ? -1.0 : 1.0;
        cutSize++;
      } else {
        coefs[idxMap[col]] += neg ? -1.0 : 1.0;
        duplicated++;
      }
      if (neg)
        rhsAdj--;
    }

    if (duplicated > 0) {
      int last = 0;
      rhsAdj = 0;
      for (int k = 0; k < cutSize; k++) {
        if (coefs[k] == -1.0 || coefs[k] == 1.0) {
          idxs[last] = idxs[k];
          coefs[last] = coefs[k];
          if (coefs[k] == -1.0)
            rhsAdj--;
          last++;
        }
      }
      cutSize = last;
    }

    // add the hub itself with coefficient -1
    idxs[cutSize] = y;
    coefs[cutSize] = -1.0;
    cutSize++;

    cutpool.add(idxs.data(), coefs.data(), cutSize, static_cast< double >(rhsAdj));

    for (int node : cliqueNodes) {
      const int col = node >= numCols ? node - numCols : node;
      idxMap[col] = -1;
    }
  }

  cutpool.removeNullCuts();

  OsiRowCut rc;
  const size_t numberRowCutsBefore = cs.sizeRowCuts();
  for (size_t i = 0; i < cutpool.numCuts(); i++) {
    rc.setRow(cutpool.cutSize(i), cutpool.cutIdxs(i), cutpool.cutCoefs(i));
    rc.setUb(cutpool.cutRHS(i));
    cs.insertIfNotDuplicate(rc);
  }

  const size_t numberRowCutsAfter = cs.sizeRowCuts();
  CglImpliedClique::sepCuts_.fetch_add(numberRowCutsAfter - numberRowCutsBefore, std::memory_order_relaxed);

  if (!info.inTree && ((info.options & 4) == 4 || ((info.options & 8) && !info.pass))) {
    for (size_t i = numberRowCutsBefore; i < numberRowCutsAfter; i++) {
      cs.rowCutPtr(i)->setGloballyValid();
    }
  }

  {
    const double delta = CoinCpuTime() - startSep;
    double old = CglImpliedClique::sepTime_.load(std::memory_order_relaxed);
    while (!CglImpliedClique::sepTime_.compare_exchange_weak(old, old + delta, std::memory_order_relaxed))
      ;
  }
}
