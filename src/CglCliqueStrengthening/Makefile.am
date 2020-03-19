AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                           libCglCliqueStrengthening                     #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglCliqueStrengthening.la

# List all source files for this library, including headers
libCglCliqueStrengthening_la_SOURCES = \
	CglCliqueStrengthening.cpp CglCliqueStrengthening.hpp

# This is for libtool (on Windows)
libCglCliqueStrengthening_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -DCGL_BUILD $(CGLLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in pkgincludedir.

includecoindir = $(pkgincludedir)
includecoin_HEADERS = CglCliqueStrengthening.hpp

