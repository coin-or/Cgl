AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                           libCglOddHoleWC                            #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglOddHoleWC.la

# List all source files for this library, including headers
libCglOddHoleWC_la_SOURCES = \
	CglOddHoleWC.cpp CglOddHoleWC.hpp \
	oddhs.cpp oddhs.h \
	spaths.cpp spaths.h

# This is for libtool (on Windows)
libCglOddHoleWC_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -DCGL_BUILD $(CGLLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in pkgincludedir.

includecoindir = $(pkgincludedir)
includecoin_HEADERS = CglOddHoleWC.hpp
