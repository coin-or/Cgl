# Copyright (C) 2020 COIN-OR
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

########################################################################
#                           libCglBKClique                     #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglBKClique.la

# List all source files for this library, including headers
libCglBKClique_la_SOURCES = CglBKClique.cpp CglBKClique.hpp

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
includecoin_HEADERS = CglBKClique.hpp
