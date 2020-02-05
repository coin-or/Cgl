
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_cgl_default.h */
#ifndef CGLLIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
#define CGLLIB_EXPORT __declspec(dllexport)
#else
#define CGLLIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_cgl_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/


/* Define to 1 if the CoinUtils package is used */
#define COIN_HAS_COINUTILS 1

/* Define to 1 if the Osi package is used */
#define COIN_HAS_OSI 1

/* Define to 1 if the Osi package is used */
#define COIN_HAS_CLP 1

/* Define to 1 if the Clp package is used */
#define COIN_HAS_OSICLP 1

