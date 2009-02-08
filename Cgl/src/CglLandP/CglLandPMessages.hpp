// Copyright (C) 2005, 2007 Pierre Bonami and others.  All Rights Reserved.
// Author:   Pierre Bonami
//           Tepper School of Business
//           Carnegie Mellon University, Pittsburgh, PA 15213
// Date:     21/07/05
//---------------------------------------------------------------------------
#ifndef CglLandPMessages_H
#define CglLandPMessages_H

#include "CoinMessage.hpp"
#include "CoinMessageHandler.hpp"

namespace LAP
{
/** Forward declaration of class to store extra debug data.*/
class DebugData;
/** Types of messages for lift-and-project simplex.*/
enum LAP_messages {
    Separating,
    FoundImprovingRow,
    FoundBestImprovingCol,
    WarnFailedBestImprovingCol,
    LogHead,
    PivotLog,
    FinishedOptimal,
    HitLimit,
    NumberNegRc,
    NumberZeroRc,
    NumberPosRc,
    WeightsStats,
    WarnBadSigmaComputation,
    WarnBadRowComputation,
    WarnGiveUpRow,
    PivotFailedSigmaUnchanged,
    PivotFailedSigmaIncreased,
    FailedSigmaIncreased,
    WarnBadRhsComputation,
    WarnFailedPivotTol,
    WarnFailedPivotIIf,
    DUMMY_END
};
/** Message handler for lift-and-project simplex. */
class LandPMessages : public CoinMessages
{
public:

    /** Constructor */
    LandPMessages();
};
}
#endif
