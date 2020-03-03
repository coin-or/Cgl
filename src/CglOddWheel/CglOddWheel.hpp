#ifndef _CglOddWheel_h_
#define _CglOddWheel_h_

#include "CglCutGenerator.hpp"
#include "CoinConflictGraph.hpp"

class CoinConflictGraph;

class CglOddWheel : public CglCutGenerator
{
public:

    static size_t sepCuts;
    static double sepTime;

    CglOddWheel(size_t extMethod = 2);
    CglOddWheel(const CglOddWheel& rhs);
    virtual CglCutGenerator * clone() const;
    virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info = CglTreeInfo() );
    virtual ~CglOddWheel();

    void setExtendingMethod(size_t extMethod);
    size_t getExtendingMethod() const { return extMethod_; }

private:
    void checkMemory(const size_t newNumCols);

    size_t cap_;
    int *idxs_, *idxMap_;
    double *coefs_;

    double *x_, *rc_;

    OsiRowCut osrc_;

    /*Extending method:
     * 0 - no extension
     * 1 - wheel center with only one variable
     * 2 - wheel center formed by a clique
    */
    size_t extMethod_;
};

#endif
