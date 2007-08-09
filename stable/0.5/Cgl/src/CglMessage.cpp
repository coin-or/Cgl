// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "CoinPragma.hpp"
#include "CglMessage.hpp"
/// Structure for use by CglMessage.cpp
typedef struct {
  CGL_Message internalNumber;
  int externalNumber; // or continuation
  char detail;
  const char * message;
} Cgl_message;
static Cgl_message us_english[]=
{
  {CGL_INFEASIBLE,0,1,"Cut generators found to be infeasible!"},
  {CGL_CLIQUES,1,2,"%d cliques of average size %g"},
  {CGL_FIXED,2,1,"%d variables fixed"},
  {CGL_PROCESS_STATS,3,1,"%d fixed, %d tightened bounds, %d strengthened rows, %d substitutions"},
  {CGL_SLACKS,4,2,"%d converted to equality constraints"},
  {CGL_PROCESS_STATS2,4,1,"processed model has %d rows, %d columns (%d integer) and %d elements"},
  {CGL_PROCESS_SOS1,5,1,"%d SOS with %d members"},
  {CGL_PROCESS_SOS2,6,1,"%d SOS (%d members out of %d) with %d overlaps - too much overlap or too many others"},
  {CGL_UNBOUNDED,0,7,"Continuous relaxation is unbounded!"},
  {CGL_ELEMENTS_CHANGED1,8,2,"%d elements changed"},
  {CGL_ELEMENTS_CHANGED2,9,3,"element in row %d for column %d changed from %g to %g"},
  {CGL_DUMMY_END,999999,0,""}
};
/* Constructor */
CglMessage::CglMessage(Language language) :
  CoinMessages(sizeof(us_english)/sizeof(Cgl_message))
{
  language_=language;
  strcpy(source_,"Cgl");
  class_ = 3; // Cuts
  Cgl_message * message = us_english;

  while (message->internalNumber!=CGL_DUMMY_END) {
     CoinOneMessage oneMessage(message->externalNumber,message->detail,
			       message->message);
     addMessage(message->internalNumber,oneMessage);
     message ++;
}

}
