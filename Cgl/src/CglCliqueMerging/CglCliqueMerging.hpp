#ifndef CglCliqueMerging_H
#define CglCliqueMerging_H

#include "CoinMessageHandler.hpp"

/** @brief merge cliques in a MIP
 *
 *  Merge cliques to produce a stronger
 *  mixed integer programming formulation
 *
 *  @param osi an OsiSolverInterface object containing a Mixed-Integer Program (MIP)
 *  @param maxExtensions maximum number of larger cliques generated from a single clique
 *  @param maxItBk maximum number of iteration in the Bron-kerbosch algorithm to extend a clique
 *  @param nExtended fills number of cliques that were extended
 *  @param nDominated fills number of cliques that were dominated
 **/

class CglCliqueMerging {
public:
	CglCliqueMerging();
	CglCliqueMerging(const CglCliqueMerging &rhs);
	CglCliqueMerging &operator=(const CglCliqueMerging &rhs);
	~CglCliqueMerging();
	void gutsOfDestructor();

    void mergeCliques(OsiSolverInterface &model);

	inline size_t maxExtensions() const {
		return maxExtensions_;
	}

	inline void setMaxExtension(size_t maxExtensions) {
		if (maxExtensions > 0) {
			maxExtensions_ = maxExtensions;
		}
	}

	inline size_t maxItBK() const {
		return maxItBK_;
	}

	inline void setMaxItBK(size_t maxItBK) {
		if (maxItBK > 0) {
			maxItBK_ = maxItBK;
		}
	}

	inline int constraintsExtended() const {
		return nExtended_;
	}

	inline int constraintsDominated() const {
		return nDominated_;
	}

	//---------------------------------------------------------------------------

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
  //---------------------------------------------------------------------------

private:
	size_t maxExtensions_;
    size_t maxItBK_;
    int nExtended_;
    int nDominated_;

  CoinMessageHandler * handler_;
  bool defaultHandler_;
  CoinMessages messages_;
};

#endif

