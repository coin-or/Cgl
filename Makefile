# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -g

LIBNAME := Cgl

LIBSRC := 
LIBSRC += CglCutGenerator.cpp
LIBSRC += CglGomory.cpp
LIBSRC += CglKnapsackCover.cpp
LIBSRC += CglLiftAndProject.cpp
LIBSRC += CglOddHole.cpp
LIBSRC += CglProbing.cpp
LIBSRC += CglSimpleRounding.cpp

TESTSRC :=
TESTSRC += CglGomoryTest.cpp
TESTSRC += CglKnapsackCoverTest.cpp
TESTSRC += CglOddHoleTest.cpp
TESTSRC += CglProbingTest.cpp
TESTSRC += CglSimpleRoundingTest.cpp

export LibType OptLevel LIBNAME LIBSRC

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export CoinDir = $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Common/make
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

export ExtraIncDir := ${OsiIncDir}
export ExtraLibDir := ${OsiLibDir}
export ExtraLibName := ${OsiLibName}
export ExtraDefine := 

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install libosi library clean doc

default: install

install library : % :
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*

clean doc : % :
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*
