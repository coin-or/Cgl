# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Cgl],[devel],[https://github.com/coin-or/Cgl/issues/new],[coin-or-cgl],[https://github.com/coin-or/Cgl])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test whether
# the package is actually there.
AC_CONFIG_SRCDIR(src/CglCommon/CglMessage.hpp)

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.
AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Initialize libtool.
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                   Math headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h
AC_COIN_CHECK_MATH_HDRS

#############################################################################
#                   Determine list of all Cgl subprojects                   #
#############################################################################

# To avoid updating the main Makefile.am whenever a new cut generator is
# added, automatically determine all subdirectories of src with a Makefile.am
# and create variables CGL_SUBDIRS (all subdirectories) and CGL_SUBLIBS
# (all sublibraries) for use in Makefile.am. Exclude CglCommon
# because that's where we do final assembly of libCgl. It is still necessary
# to manually add an entry to write the Makefile in the AC_CONFIG_FILES
# list below.

AC_SUBST(CGL_SUBDIRS)
AC_SUBST(CGL_SUBLIBS)

curr_dir=`pwd`
cd $srcdir/src
for file in `ls */Makefile.am`; do
  newdir=`echo $file | sed -e s'%\(.*\)/Makefile.am%\1%'`
  if test "$newdir" != CglCommon ; then
    CGL_SUBDIRS="$CGL_SUBDIRS src/$newdir"
    CGL_SUBLIBS="$CGL_SUBLIBS ../$newdir/lib${newdir}.la"
  fi
done
CGL_SUBDIRS="$CGL_SUBDIRS src/CglCommon"
cd $curr_dir

#############################################################################
#                            COIN-OR components                             #
#############################################################################

AC_COIN_CHK_PKG(CoinUtils,[CglLib CglUnitTest])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi

AC_COIN_CHK_PKG(Osi,[CglLib CglUnitTest])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi not available.])
fi

# Clp can be used by CglLandP, and also by the unitTest.
AC_COIN_CHK_PKG(Clp,[CglLib CglUnitTest],[clp])
AC_COIN_CHK_PKG(OsiClp,[CglLib CglUnitTest],[osi-clp])

# Finalize libs and cflags for building CglLib components, then add one
# last component to access header files in CglCommon. This should be absolute
# from top_srcdir so it works for both the individual cut generators and the
# unit test.

AC_COIN_FINALIZE_FLAGS(CglLib)
CGLLIB_CFLAGS="-I\$(top_srcdir)/src/CglCommon $CGLLIB_CFLAGS"

#############################################################################
#                       Solvers for unittesting                             #
#############################################################################

# These are not required for the compilation of Cgl, but we need to
# know which are available in order to compile the unitTest correctly. The
# CglUnitTest unit is simply a convenient handle, allowing us to accumulate
# linker and compiler flags in a single variable.

AC_COIN_CHK_PKG(OsiCpx,[CglUnitTest],[osi-cplex])
AC_COIN_CHK_PKG(OsiGlpk,[CglUnitTest],[osi-glpk])
AC_COIN_CHK_PKG(OsiGrb,[CglUnitTest],[osi-gurobi])
AC_COIN_CHK_PKG(OsiMsk,[CglUnitTest],[osi-mosek])
AC_COIN_CHK_PKG(OsiXpr,[CglUnitTest],[osi-xpress])
AC_COIN_CHK_PKG(OsiVol,[CglUnitTest],[osi-vol])
AC_COIN_CHK_PKG(OsiDyLP,[CglUnitTest],[osi-dylp])
AC_COIN_CHK_PKG(OsiSpx,[CglUnitTest],[osi-soplex])

AC_COIN_FINALIZE_FLAGS(CglUnitTest)

# We just need the COIN_HAS_SAMPLE and SAMPLE_DATA variables for Sample.

AC_COIN_CHK_PKG(Sample,,[coindatasample],,dataonly)


#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN([CoinUtils Osi Clp])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)

AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/CglCommon/Makefile
                 src/CglAllDifferent/Makefile
                 src/CglBKClique/Makefile
                 src/CglClique/Makefile
                 src/CglCliqueStrengthening/Makefile
                 src/CglDuplicateRow/Makefile
                 src/CglFlowCover/Makefile
                 src/CglGMI/Makefile
                 src/CglGomory/Makefile
                 src/CglKnapsackCover/Makefile
                 src/CglLandP/Makefile
                 src/CglLiftAndProject/Makefile
                 src/CglMixedIntegerRounding/Makefile
                 src/CglMixedIntegerRounding2/Makefile
                 src/CglOddHole/Makefile
                 src/CglOddWheel/Makefile
                 src/CglPreProcess/Makefile
                 src/CglProbing/Makefile
                 src/CglRedSplit/Makefile
                 src/CglRedSplit2/Makefile
                 src/CglResidualCapacity/Makefile
                 src/CglSimpleRounding/Makefile
                 src/CglTwomir/Makefile
                 src/CglZeroHalf/Makefile
                 test/Makefile
                 cgl.pc])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADERS([src/CglCommon/config.h src/CglCommon/config_cgl.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
