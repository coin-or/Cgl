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
# in the 'lib/coin' directory
libcoindir = $(libdir)/coin
libcoin_LTLIBRARIES = libCgl.la

# List all source files for this library, including headers
libCgl_la_SOURCES = \
	CglConfig.h \
	CglCutGenerator.cpp CglCutGenerator.hpp\
	CglMessage.cpp CglMessage.hpp \
	CglStored.cpp CglStored.hpp \
	CglParam.cpp CglParam.hpp \
	CglTreeInfo.cpp CglTreeInfo.hpp

# This is for libtool (on Windows)
libCgl_la_LDFLAGS = $(LT_LDFLAGS)

# We want to have all the sublibraries from the Cgl subprojects collected into
# this library
libCgl_la_LIBADD = $(CGL_SUBLIBS)

# Since automake is not doing this on its own, we need to declare the
# dependencies to the subdirectory libraries here
libCgl_la_DEPENDENCIES = $(libCgl_la_LIBADD)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS = $(COINUTILS_CFLAGS) $(OSI_CFLAGS)

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../inc/config_cgl.h \
	CglCutGenerator.hpp \
	CglMessage.hpp \
	CglStored.hpp \
	CglParam.hpp \
	CglTreeInfo.hpp

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = CglConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)
