/**
 *
 * This file is part of the COIN-OR CBC MIP Solver
 *
 * Class that implements a conflict-based preprocessing.
 * It tries to extend set packing constraints considering
 * the conflict graph and using a greedy strategy.
 *
 * @file CglCliqueStrengthening.hpp
 * @brief Conflict-based preprocessing
 * @author Samuel Souza Brito and Haroldo Gambini Santos
 * Contact: samuelbrito@ufop.edu.br and haroldo@ufop.edu.br
 * @date 03/27/2020
 *
 * \copyright{Copyright 2020 Brito, S.S. and Santos, H.G.}
 * \license{This This code is licensed under the terms of the Eclipse Public License (EPL).}
 *
 **/

#ifndef CGLCLIQUESTRENGTHENING_HPP
#define CGLCLIQUESTRENGTHENING_HPP

#include "CoinMessageHandler.hpp"
#include "CglConfig.h"

class OsiSolverInterface;
class CoinConflictGraph;

/**
 * Class that implements a conflict-based preprocessing.
 * It tries to extend set packing constraints considering
 * the conflict graph and using a greedy strategy.
 **/
class CGLLIB_EXPORT CglCliqueStrengthening {
public:
  /**
   * Default constructor
   **/
  CglCliqueStrengthening();

  /**
   * Copy constructor
   **/
  CglCliqueStrengthening(const CglCliqueStrengthening &rhs);

  /**
   * Assignment operator
   **/
  CglCliqueStrengthening &operator=(const CglCliqueStrengthening &rhs);

  /**
   * Destructor
   **/
  ~CglCliqueStrengthening();

  /**
   * Clears out as much as possible
   **/
  void gutsOfDestructor();

  /**
   * Tries to strengthen set packing constraints of a model.
   * After strengthening (extending), dominated constraints
   * are removed (clique merging).
   *
   *
   * Extension method: 0 = no extension;1 = random;
   * 2 = max degree; 3 = max modified degree;
   * 4 = reduced cost (inversely proportional);
   * 5 = reduced cost (inversely proportional) + modified degree.
   **/
  void strengthenCliques(OsiSolverInterface &model, size_t extMethod = 4);

  /**
   * Return the number of set packing constraints extended.
   **/
  int constraintsExtended() const { return nExtended_; }

  /**
   * Return the number of set packing constraints dominated
   * by the extended constraints.
   **/
  int constraintsDominated() const { return nDominated_; }

  /**
   * Pass in Message handler (not deleted at end)
   **/
  void passInMessageHandler(CoinMessageHandler * handler);
  
  /**
   * Set language
   **/
  void newLanguage(CoinMessages::Language language);
  
  /**
   * New language
   **/
  inline void setLanguage(CoinMessages::Language language)
  {newLanguage(language);}

  /**
   * Return message handler
   **/
  inline CoinMessageHandler * messageHandler() const
  {return handler_;}

  /**
   * Return messages
   **/
  inline CoinMessages messages()
  {return messages_;}

  /**
   * Return a pointer to messages
   **/
  inline CoinMessages * messagesPointer()
  {return &messages_;}

private:
  /**
   * Number of extended constraints.
   **/
  int nExtended_;

  /**
   * Number of dominated constraints.
   **/
  int nDominated_;

  /**
   * Message handler
   **/
  CoinMessageHandler * handler_;

  /**
   * Flag to say if handler_ is the default handler.
   * The default handler is deleted when the model
   * is deleted. Other handlers (supplied by the client)
   * will not be deleted.
   **/
  bool defaultHandler_;

  /**
   * Messages
   **/
  CoinMessages messages_;
};


#endif //CGLCLIQUESTRENGTHENING_HPP
