#ifndef CGLCLIQUESTRENGTHENING_HPP
#define CGLCLIQUESTRENGTHENING_HPP

#include "CoinMessageHandler.hpp"

class OsiSolverInterface;
class CoinConflictGraph;

class CglCliqueStrengthening {
public:
    CglCliqueStrengthening();
    CglCliqueStrengthening(const CglCliqueStrengthening &rhs);
    CglCliqueStrengthening &operator=(const CglCliqueStrengthening &rhs);
    ~CglCliqueStrengthening();
    void gutsOfDestructor();

    void strengthenCliques(OsiSolverInterface &model, const CoinConflictGraph *cgraph, size_t extMethod = 4);
    int constraintsExtended() const { return nExtended_; }
    int constraintsDominated() const { return nDominated_; }

    /**@name Message handling */
    //@{
    /// Pass in Message handler (not deleted at end)
    void passInMessageHandler(CoinMessageHandler * handler);
    /// Set language
    void newLanguage(CoinMessages::Language language);
    inline void setLanguage(CoinMessages::Language language)
    {newLanguage(language);}
    /// Return handler
    inline CoinMessageHandler * messageHandler() const
    {return handler_;}
    /// Return messages
    inline CoinMessages messages()
    {return messages_;}
    /// Return pointer to messages
    inline CoinMessages * messagesPointer()
    {return &messages_;}
    //@}

private:
    int nExtended_, nDominated_;

    CoinMessageHandler * handler_;
    bool defaultHandler_;
    CoinMessages messages_;
};


#endif //CGLCLIQUESTRENGTHENING_HPP
