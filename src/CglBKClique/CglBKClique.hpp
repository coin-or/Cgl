#ifndef _CglBKClique_h_
#define _CglBKClique_h_

#include "CglCutGenerator.hpp"

class CGLLIB_EXPORT CglBKClique : public CglCutGenerator {
public:

    static size_t sepCuts;
    static double sepTime;

    CglBKClique();
    CglBKClique(const CglBKClique& rhs);
    virtual CglCutGenerator * clone() const;
    virtual void generateCuts( const OsiSolverInterface &si, OsiCuts &cs, const CglTreeInfo info = CglTreeInfo() );
    virtual ~CglBKClique();

    void setMaxItBK(size_t _maxItBK);
    size_t getMaxItBK() { return maxItBK; }
    void setMaxItBKExt(size_t _maxItBKExt);
    size_t getMaxItBKExt() { return maxItBKExt; }
    void setExtendingMethod(size_t _extMethod);
    size_t getExtendingMethod() { return extMethod; }

private:

    size_t maxItBK, maxItBKExt, extMethod;
};

#endif // CglBKClique_HPP
