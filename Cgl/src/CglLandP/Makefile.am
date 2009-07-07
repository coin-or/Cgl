# Copyright (C) 2006, 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 243 2006-06-01 20:11:02Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                           libCglLandP                                #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglLandP.la

# List all source files for this library, including headers
libCglLandP_la_SOURCES = \
	CglLandP.cpp  CglLandP.hpp \
        CglLandPUtils.cpp  CglLandPUtils.hpp \
        CglLandPTabRow.cpp CglLandPTabRow.hpp \
        CglLandPMessages.hpp CglLandPMessages.cpp \
        CglLandPSimplex.cpp CglLandPSimplex.hpp \
        CglLandPValidator.cpp CglLandPValidator.hpp \
        CglLandPTest.cpp


# This is for libtool (on Windows)
libCglLandP_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../CglGomory`

if COIN_HAS_COINUTILS_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(COINUTILSINSTDIR)/include/coin`
else
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`
endif

if COIN_HAS_OSI_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSIINSTDIR)/include/coin`
else
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp`
endif

if COIN_HAS_CLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPINSTDIR)/include/coin`
else
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` 
endif


# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	CglLandP.hpp CglLandPSimplex.hpp CglLandPValidator.hpp 

