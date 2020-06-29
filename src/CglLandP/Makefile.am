# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                           libCglLandP                                #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglLandP.la

# List all source files for this library, including headers
libCglLandP_la_SOURCES = CglLandPTabRow.hpp CglLandPUtils.cpp CglLandPUtils.hpp

if COIN_HAS_OSICLP
libCglLandP_la_SOURCES += \
	CglLandP.cpp  CglLandP.hpp \
	CglLandPTabRow.cpp \
	CglLandPMessages.h CglLandPMessages.cpp \
	CglLandPSimplex.cpp CglLandPSimplex.hpp \
	CglLandPValidator.cpp CglLandPValidator.hpp \
	CglLandPTest.cpp
endif

# This is for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = $(CGLLIB_CFLAGS) -I$(srcdir)/../CglGomory

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in $(includedir)/coin-or.
includecoindir = $(includedir)/coin-or
if COIN_HAS_OSICLP
includecoin_HEADERS = CglLandP.hpp CglLandPValidator.hpp
endif
