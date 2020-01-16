
AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                           libCglBKClique                     #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglBKClique.la

# List all source files for this library, including headers
libCglBKClique_la_SOURCES = \
	CglBKClique.cpp CglBKClique.hpp

# This is for libtool (on Windows)
libCglBKClique_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -DCGL_BUILD $(CGLLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in pkgincludedir.

includecoindir = $(pkgincludedir)
includecoin_HEADERS = CglBKClique.hpp

