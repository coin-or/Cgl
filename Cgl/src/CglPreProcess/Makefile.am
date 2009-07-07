# Copyright (C) 2006, 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                         libCglPreProcess                             #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglPreProcess.la

# List all source files for this library, including headers
libCglPreProcess_la_SOURCES = \
	CglPreProcess.cpp  CglPreProcess.hpp

# This is for libtool (on Windows)
libCglPreProcess_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../CglProbing` \
	-I`$(CYGPATH_W) $(srcdir)/../CglDuplicateRow`

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
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	CglPreProcess.hpp
