# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Osi 
# - Osi*/Makefile for the libs you have specified above

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install library clean unitTest

default: install

unitTest : install
	(cd Test && ${MAKE} unitTest)

install clean library : % :
	(cd ../Osi && ${MAKE} inst-libOsi)
	${MAKE} -f Makefile.Cgl $*
