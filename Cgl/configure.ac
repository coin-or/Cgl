# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Cgl],[split],[cgl@list.coin-or.org])

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
#                            COIN-OR components                             #
#############################################################################

AC_COIN_HAS_MODULE(CoinUtils,          [coinutils = split])
if test $coin_has_coinutils != yes ; then
  #if we have to use the fallback routines, then obviously it will not find CoinUtils dependencies
  #thus, we check whether ThirdParty/Blas and ThirdParty/Lapack have been setup and
  #assume that if present then they are dependencies of CoinUtils
  AC_COIN_HAS_MODULE_FALLBACK(Blas,   [coinblas],   [../ThirdParty/Blas],   [required=0 dodefine=0 doconditional=0])
  AC_COIN_HAS_MODULE_FALLBACK(Lapack, [coinlapack], [../ThirdParty/Lapack], [required=0 dodefine=0 doconditional=0])
fi
AC_COIN_HAS_MODULE_FALLBACK(CoinUtils, [coinutils], [../CoinUtils])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
#so if we had to do find Blas and Lapack separately, then add their linker flags to the CoinUtils flags
#actually the libs may already be included in libCoinUtils, but we may still need the fortran library flags
COINUTILS_LIBS="$COINUTILS_LIBS $BLAS_LIBS $LAPACK_LIBS"

AC_COIN_HAS_MODULE(Osi,               [osi = split])
AC_COIN_HAS_MODULE_FALLBACK(Osi,      [osi], [../Osi])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi not available.])
fi

AC_COIN_HAS_MODULE(Sample,          [coindatasample],                   [required=0])
AC_COIN_HAS_MODULE_FALLBACK(Sample, [coindatasample], [../Data/Sample], [required=0])

AC_COIN_HAS_MODULE(Netlib,          [coindatanetlib],                   [required=0])
AC_COIN_HAS_MODULE_FALLBACK(Netlib, [coindatanetlib], [../Data/Netlib], [required=0])

# Clp can be used by CglLandP
AC_COIN_HAS_MODULE(OsiClp, [clp osi-clp])
if test $coin_has_osiclp != yes ; then
  # in the fallback, we first check for OsiClp
  # if it is there, then we also check for Clp and add the Clp flags to the OsiClp flags
  AC_COIN_HAS_MODULE_FALLBACK(OsiClp, [osi-clp], [../Clp])
  if test $coin_has_osiclp = yes ; then
    AC_COIN_HAS_MODULE_FALLBACK(Clp, [clp], [../Clp], [dodefine=0 doconditional=0])
    OSICLP_CFLAGS="$OSICLP_CFLAGS $CLP_CFLAGS"
    OSICLP_LIBS="$OSICLP_LIBS $CLP_LIBS"
  else
    coin_has_osiclp=no
  fi
fi

#############################################################################
#                       Solvers for unittesting                             #
#############################################################################

# Those are not really required for the compilation of Cgl, but we
# need to know which are available in order to compile the unitTest
# correctly.

AC_COIN_HAS_MODULE(OsiCpx,  [osi-cplex],    [required=0])
AC_COIN_HAS_MODULE(OsiGlpk, [osi-glpk],     [required=0])
#AC_COIN_HAS_MODULE(OsiGrb,  [osi-gurobi],   [required=0])
AC_COIN_HAS_MODULE(OsiMsk,  [osi-mosek],    [required=0])
AC_COIN_HAS_MODULE(OsiXpr,  [osi-xpress],   [required=0])
AC_COIN_HAS_MODULE(OsiVol,  [osi-vol],      [required=0])
AC_COIN_HAS_MODULE(OsiDyLP, [osi-dylp],     [required=0])

#############################################################################
#                   Determine list of all Cgl subprojects                   #
#############################################################################

# In order to avoid that one has to update this configure.ac and the
# Makefile.am whenever a new cut generator is added, we automatically
# determine the list of all subdirectories of src with Makefile
# (CGL_SUBDIRS), and also a list of all sublibraries relative to src
# (CGL_SUBLIBS) that can be used in the src/Makefile.am.  (We still need
# to manually add the Makefile in the AC_CONFIG_FILES list below)

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

# Install files for unitTest
AC_COIN_VPATH_LINK([test/CglTestData/capPlan1.mps test/CglTestData/l152lav.mps test/CglTestData/egout.mps])

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
		 src/CglResidualCapacity/Makefile
		 src/CglSimpleRounding/Makefile
		 src/CglTwomir/Makefile
                 test/Makefile
		 cgl_addlibs.txt
		 cgl.pc
		 cgl-uninstalled.pc])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_cgl.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
