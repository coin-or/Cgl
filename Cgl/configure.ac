# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cgl],[0.1],[coin-cgl@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/CglMessage.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Cgl)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(Cbc)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(DyLP)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(Sym)
AC_COIN_HAS_PROJECT(Vol)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Those are not really required for the compilation of Cgl, but we
# need to know which are available in order to compile the unitTest
# correctly.  Therefore, we test those again, even though those tests
# are done already for Osi (and the libraries are already listed in
# osi_addlibs.txt)

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],[glp_lpx_simplex])
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosekdl.h],[MSK_openmosek])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

#############################################################################
#                   Determine list of all Cgl subprojects                   #
#############################################################################

# In order to avoid that one has to update this configure.as and the
# Makefile.am whenever a new cut generator is added, we automatically
# determine the list of all subdirectories of src with Makefile
# (CGL_SUBDIRS), and also a list of all sublibraries relative to src
# (CGL_SUBLIBS) that can be used in the src/Makefile.am.  (We still need
# to manually add the Makefile in the AC_CONFIG_FILES list below)1

AC_SUBST(CGL_SUBDIRS)
AC_SUBST(CGL_SUBLIBS)

curr_dir=`pwd`
cd $srcdir/src
for file in `find . -name Makefile.in`; do
  if test $file != ./Makefile.in; then
    newdir=`echo $file | sed -e s%/Makefile.in%% | sed -e s%^./%% `
    CGL_SUBDIRS="$CGL_SUBDIRS $newdir"
    CGL_SUBLIBS="$CGL_SUBLIBS $newdir/lib${newdir}.la"
  fi
done
cd $curr_dir

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
		 src/CglAllDifferent/Makefile
		 src/CglClique/Makefile
		 src/CglDuplicateRow/Makefile
		 src/CglFlowCover/Makefile
		 src/CglGomory/Makefile
		 src/CglKnapsackCover/Makefile
		 src/CglLandP/Makefile
		 src/CglLiftAndProject/Makefile
		 src/CglMixedIntegerRounding/Makefile
		 src/CglMixedIntegerRounding2/Makefile
		 src/CglOddHole/Makefile
		 src/CglPreProcess/Makefile
		 src/CglProbing/Makefile
		 src/CglRedSplit/Makefile
		 src/CglSimpleRounding/Makefile
		 src/CglTwomir/Makefile
                 test/Makefile
		 cgl_addlibs.txt])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_cgl.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
