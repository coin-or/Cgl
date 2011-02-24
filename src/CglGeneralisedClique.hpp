
// $Id$
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglGeneralisedClique_H
#define CglGeneralisedClique_H

/*
  The content of this file doesn't really justify its existence. It's more of
  a placeholder. Eventually, generalised cliques need to be reimplemented
  properly. When that happens, the associated structures should be defined
  here.
  
  At the time this file was created, only CglKnapsackCover has active code
  that uses generalised cliques. All of the code in CglProbing associated
  with generalised cliques has been removed. All of the code in CglPreProcess
  associated with generalised cliques is disabled. In general, look for the
  symbol CLIQUE_ANALYSIS to find related code in other projects.

  -- lh, 110224 --
*/

/*! \brief Utility structure to encode clique information.

  The only purpose is to hide the details of the encoding. The msb of #fixes
  is set to 1 if setting the variable to 1 fixes all other variables in the
  clique, 0 if setting the variable to 0 fixes all other variables in the
  clique. The remaining 31 bits are the index of the variable.
*/
typedef struct {
  //unsigned int oneFixed:1; //  nonzero if variable to 1 fixes all
  //unsigned int sequence:31; //  variable (in matrix) (but also see cliqueRow_)
  /// Encoded: msb indicates strong-0/strong-1, remaining bits are index.
  unsigned int fixes;
} cliqueEntry;

/*! \brief Extract the index from the clique entry
    \related cliqueEntry
*/
inline int sequenceInCliqueEntry(const cliqueEntry & cEntry)
{ return cEntry.fixes&0x7fffffff;}

/*! \brief Set the index in the clique entry
    \related cliqueEntry
*/
inline void setSequenceInCliqueEntry(cliqueEntry & cEntry,int sequence)
{ cEntry.fixes = sequence|(cEntry.fixes&0x80000000);}

/*! \brief Mark the entry as strong-1
    \related cliqueEntry
*/
inline bool oneFixesInCliqueEntry(const cliqueEntry & cEntry)
{ return (cEntry.fixes&0x80000000)!=0;}

/*! \brief Mark the entry as strong-0
    \related cliqueEntry
*/
inline void setOneFixesInCliqueEntry(cliqueEntry & cEntry,bool oneFixes)
{ cEntry.fixes = (oneFixes ? 0x80000000 : 0)|(cEntry.fixes&0x7fffffff);}

#endif
