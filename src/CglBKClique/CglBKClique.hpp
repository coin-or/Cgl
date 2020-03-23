#ifndef _CglBKClique_h_
#define _CglBKClique_h_

#include <CglCutGenerator.hpp>

class CoinCliqueList;
class CoinConflictGraph;

class CglBKClique : public CglCutGenerator {
public:
    CglBKClique();
    CglBKClique(const CglBKClique& rhs);
    virtual CglCutGenerator * clone() const;
    virtual void generateCuts(const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info = CglTreeInfo());
    virtual ~CglBKClique();
    virtual void refreshSolver(OsiSolverInterface *solver);

    void setMaxCallsBK(size_t maxCallsBK);
    size_t getMaxCallsBK() const { return maxCallsBK_; }
    void setExtendingMethod(size_t extMethod);
    size_t getExtendingMethod() const { return extMethod_; }
    size_t getNumCallsBK() const { return callsBK_; }

    const double getMinFrac() const { return minFrac_; }
    const double getMinViol() const { return minViol_; }
    const double getMinWeight() const { return minWeight_; }
    void setMinFrac(const double minFrac);
    void setMinViol(const double minViol);
    void setPivotingStrategy(const size_t pivotingStrategy);

    static size_t sepCuts_;
    static double sepTime_;

private:
    void checkMemory(const size_t newNumCols);
    CoinCliqueList* separateCliques(const OsiSolverInterface &si);
    CoinCliqueList* extendCliques(const OsiSolverInterface &si, const CoinCliqueList *initialCliques);
    void insertCuts(const OsiSolverInterface &si, const CglTreeInfo &info, const CoinCliqueList *cliques, OsiCuts &cs);

    size_t cap_;
    double *rc_;
    int *idxs_, *idxMap_;
    double *coefs_;
    size_t *inducedVert_, *currClq_;
    double *vertexWeight_;

    double minFrac_;
    double minViol_;
    double minWeight_;

    size_t pivotingStrategy_;
    size_t extMethod_;
    size_t maxCallsBK_;
    size_t callsBK_;
    bool completeBK_;

    OsiRowCut osrc_;
};

#endif // CglBKClique_HPP
