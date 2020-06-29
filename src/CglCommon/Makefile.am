# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                               libCgl                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the $libdir directory
lib_LTLIBRARIES = libCgl.la

# List all source files for this library, including headers
libCgl_la_SOURCES = \
	CglConfig.h \
	CglCutGenerator.cpp CglCutGenerator.hpp\
	CglMessage.cpp CglMessage.hpp \
	CglStored.cpp CglStored.hpp \
	CglParam.cpp CglParam.hpp \
	CglTreeInfo.cpp CglTreeInfo.hpp

# We want to have all the sublibraries from the Cgl subprojects collected into
# this library
libCgl_la_LIBADD = $(CGL_SUBLIBS) $(CGLLIB_LFLAGS)
libCgl_la_DEPENDENCIES = $(CGL_SUBLIBS)

# This is for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = $(CGLLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in $(includedir)/coin-or.
includecoindir = $(includedir)/coin-or
includecoin_HEADERS = \
	CglCutGenerator.hpp \
	CglMessage.hpp \
	CglStored.hpp \
	CglParam.hpp \
	CglTreeInfo.hpp

install-exec-local:
	$(install_sh_DATA) config_cgl.h $(DESTDIR)$(includecoindir)/CglConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/CglConfig.h
