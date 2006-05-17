# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# We need to recurse into all Cgl subproject subdirectories
SUBDIRS = $(CGL_SUBDIRS)

########################################################################
#                               libCgl                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCgl.la

# List all source files for this library, including headers
libCgl_la_SOURCES = \
	CglCutGenerator.cpp CglCutGenerator.hpp\
	CglMessage.cpp CglMessage.hpp \
	CglStored.cpp CglStored.hpp

# This is for libtool (on Windows)
libCgl_la_LDFLAGS = # -no-undefined

# We want to have all the sublibraries from the Cgl subprojects collected into
# this library
libCgl_la_LIBADD = $(CGL_SUBLIBS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
CXXFLAGS += \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	CglCutGenerator.hpp \
	CglMessage.hpp \
	CglStored.hpp

########################################################################
#     Generating the configure header that is to be distributed        #
########################################################################

# We don't do anything here yet, since so far nothing important is in the
# config.h file
