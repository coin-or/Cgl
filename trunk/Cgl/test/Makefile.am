# Copyright (C) 2006 International Business Machines and others.
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
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(srcdir)/../src` \
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
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(CLPOBJDIR)/src/libClp.la
  LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif

if COIN_HAS_CBC
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCbc`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiCbc/libOsiCbc.la \
	$(CBCOBJDIR)/src/libCbc.la
  LIBS += `cat $(CBCOBJDIR)/cbc_addlibs.txt`
endif

if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
endif

if COIN_HAS_DYLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
	-I`$(CYGPATH_W) $(DYLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiDylp`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiDylp/libOsiDylp.la \
	$(DYLPOBJDIR)/src/Dylp/libDylp.la
  LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
endif

if COIN_HAS_FMP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiFmp`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiFmp/libOsiFmp.la
endif

if COIN_HAS_GLPK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiGlpk`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiGlpk/libOsiGlpk.la
endif

if COIN_HAS_MSK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiMsk`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiMsk/libOsiMsk.la
endif

if COIN_HAS_OSL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiOsl`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiOsl/libOsiOsl.la
endif

if COIN_HAS_SPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSpx`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiSpx/libOsiSpx.la
endif

if COIN_HAS_SYMPHONY
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiSym`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiSym/libOsiSym.la \
	$(SYMPHONYOBJDIR)/src/libSym.la
  LIBS += `cat $(SYMPHONYOBJDIR)/sym_addlibs.txt`
endif

if COIN_HAS_VOL
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiVol`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiVol/libOsiVol.la \
	$(VOLOBJDIR)/src/libVol.la
  LIBS += `cat $(VOLOBJDIR)/vol_addlibs.txt`
endif

if COIN_HAS_XPR
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiXpr`
  unitTest_LDADD += $(OSIOBJDIR)/src/OsiXpr/libOsiXpr.la
endif

# Finally the Coin library (many of the above need it)
unitTest_LDADD += \
	$(OSIOBJDIR)/src/libOsi.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

# List additional libraries
LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt`\
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt` \
	$(ADDLIBS)

test: unitTest$(EXEEXT)
	./unitTest$(EXEEXT)

.PHONY: test

########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

DISTCLEANFILES =
