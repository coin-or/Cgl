#ifndef _CglOddHoleWC_h_
#define _CglOddHoleWC_h_

#include "CglCutGenerator.hpp"

class CGLLIB_EXPORT CglOddHoleWC : public CglCutGenerator
{
public:

    static size_t sepCuts;
    static double sepTime;

    CglOddHoleWC();
    CglOddHoleWC(const CglOddHoleWC& rhs);
    virtual CglCutGenerator * clone() const;
    virtual void generateCuts( const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info = CglTreeInfo() );
    virtual ~CglOddHoleWC();
};

#endif
