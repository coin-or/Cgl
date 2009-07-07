# Copyright (C) 2006, 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                         unitTest for Coin                            #
########################################################################

noinst_PROGRAMS = unitTest

unitTest_SOURCES = unitTest.cpp

# List libraries of COIN projects
unitTest_LDADD = ../src/libCgl.la

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglFlowCover` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglGomory` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglKnapsackCover` \
        -I`$(CYGPATH_W) $(srcdir)/../src/CglLandP` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglOddHole` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglProbing` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglMixedIntegerRounding2` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglResidualCapacity` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglRedSplit` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglTwomir` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglClique` \
	-I`$(CYGPATH_W) $(srcdir)/../src/CglSimpleRounding`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# Depending of what solvers are available, we add the corresponding files,
# libraries and include dirs

# Depending of what solvers are available, we add the corresponding files,
# libraries and include dirs

if COIN_HAS_CLP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiClp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la
endif
if COIN_HAS_CLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(CLPINSTDIR)/lib` -lClp \
	`cat $(CLPDOCDIR)/clp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  unitTest_LDADD += $(CLPOBJDIR)/src/libClp.la
  LIBS += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif
endif

if COIN_HAS_CBC
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiCbc
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCbc`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiCbc/libOsiCbc.la
endif
if COIN_HAS_CBC_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CBCINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(CBCINSTDIR)/lib` -lCbc \
	`cat $(CBCDOCDIR)/cbc_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CBCSRCDIR)/src`
  unitTest_LDADD += $(CBCOBJDIR)/src/libCbc.la
  LIBS += `cat $(CBCOBJDIR)/cbc_addlibs.txt`
endif
endif

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiCpx
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif
endif

if COIN_HAS_DYLP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiDylp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiDylp`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiDylp/libOsiDylp.la
endif
if COIN_HAS_DYLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(DYLPINSTDIR)/lib` -lDylp \
	`cat $(DYLPDOCDIR)/dylp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
	-I`$(CYGPATH_W) $(DYLPOBJDIR)/inc`
  unitTest_LDADD += $(DYLPOBJDIR)/src/libDylp.la
  LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
endif
endif

if COIN_HAS_FMP
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiFmp
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiFmp`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiFmp/libOsiFmp.la
endif
endif

if COIN_HAS_GLPK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiGlpk
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la
endif
endif

if COIN_HAS_MSK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiMsk
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiMsk`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiMsk/libOsiMsk.la
endif
endif

if COIN_HAS_OSL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiOsl
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
endif
endif

if COIN_HAS_SPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiSpx
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
endif
endif

if COIN_HAS_SYMPHONY
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiSym
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la
endif
if COIN_HAS_SYMPHONY_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(SYMPHONYINSTDIR)/lib` -lSym \
	`cat $(SYMPHONYDOCDIR)/sym_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYSRCDIR)/src`
  unitTest_LDADD += $(SYMPHONYOBJDIR)/src/libSym.la
  LIBS += `cat $(SYMPHONYOBJDIR)/sym_addlibs.txt`
endif
endif

if COIN_HAS_VOL
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiVol
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiVol`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiVol.la
endif
if COIN_HAS_VOL_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(VOLINSTDIR)/lib` -lVol \
	`cat $(VOLDOCDIR)/vol_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLSRCDIR)/src`
  unitTest_LDADD += $(VOLOBJDIR)/src/libVol.la
  LIBS += `cat $(VOLOBJDIR)/vol_addlibs.txt`
endif
endif

if COIN_HAS_XPR
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)`
if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsiXpr
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
endif
endif

if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(OSIINSTDIR)/lib` -lOsi \
	`cat $(COINUTILSDOCDIR)/coinutils_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc`
  unitTest_LDADD += $(OSIOBJDIR)/src/libOsi.la
  LIBS += `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`
endif

if COIN_HAS_COINUTILS_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(COINUTILSINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(COINUTILSOBJDIR)/lib` -lCoinUtils \
	`cat $(OSIDOCDIR)/osi_addlibs.txt`
else
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`
  unitTest_LDADD += $(COINUTILSOBJDIR)/src/libCoinUtils.la
  LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt`
endif

# List additional libraries
LIBS += $(ADDLIBS)

test: unitTest$(EXEEXT)
	./unitTest$(EXEEXT)

.PHONY: test

########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

DISTCLEANFILES =
