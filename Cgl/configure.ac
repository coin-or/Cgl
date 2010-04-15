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
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_MODULE(CoinUtils,  [coinutils = split])
AC_COIN_HAS_MODULE(Osi,        [osi = split])
AC_COIN_HAS_MODULE(Sample,     [coindatasample], [required=0])
AC_COIN_HAS_MODULE(Netlib,     [coindatanetlib], [required=0])
# Clp can be used by CglLandP
AC_COIN_HAS_MODULE(OsiClp,        [osi-clp])

# fallback if no pkg-config is available
if test "$coin_has_coinutils" != yes ; then
  AC_COIN_HAS_PROJECT(CoinUtils)
  if test $coin_has_coinutils != notGiven && test $coin_has_coinutils != skipping ; then
    COINUTILS_CFLAGS="-I`${CYGPATH_W} ${COINUTILSSRCDIR}/src` -I`${CYGPATH_W} ${COINUTILSOBJDIR}/inc`"
    COINUTILS_LIBS="${COINUTILSOBJDIR}/src/libCoinUtils.la `cat ${COINUTILSOBJDIR}/coinutils_addlibs.txt`"
  else
    AC_MSG_ERROR("Required package CoinUtils unavailable.")
  fi
fi

if test "$coin_has_osi" != yes ; then
  AC_COIN_HAS_PROJECT(Osi)
  if test $coin_has_osi != notGiven && test $coin_has_osi != skipping ; then
    OSI_CFLAGS="-I`${CYGPATH_W} ${OSISRCDIR}/src` -I`${CYGPATH_W} ${OSIOBJDIR}/inc`"
    OSI_LIBS="${OSIOBJDIR}/src/libOsi.la `cat ${OSIOBJDIR}/osi_addlibs.txt`"
  else
    AC_MSG_ERROR("Required package Osi unavailable.")
  fi
fi

if test "$coin_has_clp" != yes ; then
  AC_COIN_HAS_PROJECT(Clp)
  if test $coin_has_clp != notGiven && test $coin_has_clp != skipping ; then
    CLP_CFLAGS="-I`${CYGPATH_W} ${CLPSRCDIR}/src` -I`${CYGPATH_W} ${CLPSRCDIR}/src/OsiClp` -I`${CYGPATH_W} ${CLPOBJDIR}/inc`"
    CLP_LIBS="${CLPOBJDIR}/src/libClp.la ${CLPOBJDIR}/src/OsiClp/libOsiClp.la `cat ${CLPOBJDIR}/clp_addlibs.txt`"
  fi    
  coin_has_osiclp=yes
  AC_DEFINE(COIN_HAS_CLP,[1],[Define to 1 if the OsiClp module is available])
else 
  # to please configure
  AM_CONDITIONAL(COIN_HAS_CLP, [test 0 = 1])
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
