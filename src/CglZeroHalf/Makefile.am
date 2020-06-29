# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                       libCglZeroHalf                           #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglZeroHalf.la

# List all source files for this library, including headers
libCglZeroHalf_la_SOURCES = \
	CglZeroHalf.cpp  CglZeroHalf.hpp \
	Cgl012cut.cpp  Cgl012cut.hpp \
	CglZeroHalfTest.cpp

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
includecoin_HEADERS = CglZeroHalf.hpp Cgl012cut.hpp
