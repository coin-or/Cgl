# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Cgl

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install clean library unitTest libdepend doc

default: install

libdepend:
	(cd ../Osi && ${MAKE} inst-libOsi)

install library: libdepend
	${MAKE} -f Makefile.Cgl $@

doc:
	doxygen $(MakefileDir)/doxygen.conf

unitTest: install
	(cd Test && ${MAKE} unitTest)

clean: 
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep
