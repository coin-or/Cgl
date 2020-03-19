# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                          unitTest for Cgl                            #
########################################################################

noinst_PROGRAMS = unitTest

unitTest_SOURCES = unitTest.cpp

# List libraries we need to link in

unitTest_LDADD = ../src/CglCommon/libCgl.la $(CGLLIB_LIBS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -DCGL_BUILD
AM_CPPFLAGS += -I$(srcdir)/../src/CglSimpleRounding
AM_CPPFLAGS += -I$(srcdir)/../src/CglKnapsackCover
AM_CPPFLAGS += -I$(srcdir)/../src/CglOddHole
AM_CPPFLAGS += -I$(srcdir)/../src/CglProbing
AM_CPPFLAGS += -I$(srcdir)/../src/CglGomory
AM_CPPFLAGS += -I$(srcdir)/../src/CglLandP
AM_CPPFLAGS += -I$(srcdir)/../src/CglMixedIntegerRounding
AM_CPPFLAGS += -I$(srcdir)/../src/CglMixedIntegerRounding2
AM_CPPFLAGS += -I$(srcdir)/../src/CglResidualCapacity
AM_CPPFLAGS += -I$(srcdir)/../src/CglRedSplit
AM_CPPFLAGS += -I$(srcdir)/../src/CglRedSplit2
AM_CPPFLAGS += -I$(srcdir)/../src/CglTwomir
AM_CPPFLAGS += -I$(srcdir)/../src/CglClique
AM_CPPFLAGS += -I$(srcdir)/../src/CglFlowCover
AM_CPPFLAGS += -I$(srcdir)/../src/CglZeroHalf
AM_CPPFLAGS += $(CGLLIB_CFLAGS) $(CGLUNITTEST_CFLAGS)

if COIN_HAS_SAMPLE
AM_CPPFLAGS += -DSAMPLEDIR=\"`$(CYGPATH_W) @SAMPLE_DATA@ | sed -e 's/\\\\/\\\\\\\\/g'`\"
endif

AM_CPPFLAGS += -DTESTDIR=\"`$(CYGPATH_W) $(srcdir)/CglTestData | sed -e 's/\\\\/\\\\\\\\/g'`\" 

unitTest_LDADD += $(CGLUNITTEST_LFLAGS)

test: unitTest$(EXEEXT)
	./unitTest$(EXEEXT)

.PHONY: test