# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -O

# The path to the root of the install directory structure.
# include files will be installed in ${InstallDir}/include, libraries in 
# ${InstallDir}/lib, executables in ${InstallDir}/bin.
InstallDir := ..

# Detect packages. If a path is set to PROBE then we'll check a couple of
# reasonable places. If ...IncDir is set to NONE then the package will not be
# included. Specify the exact path if you know the package exists but we can't
# find it. 
# List all the libraries in the package to be linked against.

# Cpx
# Never use PROBE for cplex. The directory structure is way too convoluted for
# us to autodetect.
# (A contribution to do autodetect is more than welcome!)
CpxIncDir := NONE
#CpxIncDir := /usr/local/ilog/cplex71/include
CpxLibDir := /usr/local/ilog/cplex71/lib/i86_linux2_glibc2.1_egcs1.1/static_mt
CpxLibName := libcplex.a

# Osl
#OslIncDir := NONE
OslIncDir := PROBE
OslLibDir := PROBE
OslLibName := libosl.so

# Vol
#VolIncDir := NONE
VolIncDir := PROBE
VolLibDir := PROBE
VolLibName := libvol.so

# Xpr
XprIncDir := NONE
#XprIncDir := PROBE
XprLibDir := PROBE
XprLibName := libxosl.so libmp-opt.so

# Osi
OsiIncDir := PROBE
OsiLibDir := PROBE
OsiLibName := libosi.so

##############################################################################
# You should not need to edit below this line.
##############################################################################

# The location of the customized Makefiles
MakefileDir := ../Common/make

include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.detectSolver
include ${MakefileDir}/Makefile.detectOsi

##############################################################################

CXXFLAGS += $(OPTFLAG)
CXXFLAGS += -Iinclude -I../Common/include
ifeq ($(OptLevel),-g)
    CXXFLAGS += -DCGL_DEBUG
endif

###############################################################################

space:= $(empty) $(empty)
OptVersion := $(subst $(space),_,$(OptLevel))

TARGETDIR := $(UNAME)$(OptVersion)
DEPDIR := dep

VPATH := . : include : Junk : ${TARGETDIR} : ${DEPDIR}

#########################################################################

CXXFLAGS += $(addprefix -I,$(OsiIncDir))

CXXFLAGS += $(addprefix -I,$(DETECTINCDIRS))
CXXFLAGS += $(addprefix -D,$(SOLVERDEFINES))

LIBDIRS := $(DETECTLIBDIRS) $(InstallDir)/lib
LIBS    := libcgl.so $(OsiLibName) $(DETECTLIBNAMES)

LDFLAGS := $(addprefix -L,$(LIBDIRS))
LDFLAGS += $(call ifprefix,$(SHLINKPREFIX),$(LIBDIRS))
LDFLAGS += $(patsubst lib%,-l%,$(basename $(LIBS)))

###############################################################################

LIBSRC := $(filter-out %Test.cpp, $(shell /bin/ls Cgl*.cpp))
LIBOBJ := $(addprefix $(TARGETDIR)/, $(LIBSRC:.cpp=.o))
LIBDEP := $(addprefix $(DEPDIR)/, $(LIBSRC:.cpp=.d))

###############################################################################

TESTSRC := $(wildcard *Test.cpp)
TESTOBJ := $(addprefix $(TARGETDIR)/, $(TESTSRC:.cpp=.o))
TESTDEP := $(addprefix $(DEPDIR)/, $(TESTSRC:.cpp=.d))

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install all libcgl clean doc unitTest

default: install
all:     unitTest

install: libcgl
	@echo "Installing include files..."
	@mkdir -p ${InstallDir}/include
# *FIXME* : do we _really_ need all the headers?
	@${CP} include/*.hpp ../Common/include/*.hpp ${InstallDir}/include
	@echo "Installing libraries..."
	@mkdir -p ${InstallDir}/lib
	@${CP} $(TARGETDIR)/libcgl$(OptVersion)$(LIBEXT) ${InstallDir}/lib
	@rm -f ${InstallDir}/lib/libcgl$(LIBEXT)
	@cd ${InstallDir}/lib; \
		ln -s libcgl$(OptVersion)$(LIBEXT) libcgl$(LIBEXT)

###############################################################################

include ${MakefileDir}/Makefile.rules

###############################################################################

libcgl : $(TARGETDIR)/libcgl$(OptVersion)$(LIBEXT)

$(TARGETDIR)/libcgl$(OptVersion)$(LIBEXT): $(LIBOBJ)
	@rm -rf Junk
	@echo ""
	@echo Creating library $(notdir $@)
	@mkdir -p $(TARGETDIR)
	@rm -f $@
	$(LD) $@ $(LIBLDFLAGS) $(LIBOBJ)

###############################################################################

unitTest: $(TARGETDIR)/unitTest

$(TARGETDIR)/unitTest : install $(TESTOBJ)
	@rm -rf Junk
	@echo ""
	@echo Creating unitTest
	@mkdir -p $(TARGETDIR)
	@rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(TESTOBJ) $(LDFLAGS) $(SYSLD) -lm

###############################################################################

doc:
	doxygen $(MakefileDir)/doxygen.conf

clean :
	rm -rf Junk
	@rm -rf $(DEPDIR)
	@rm -rf $(TARGETDIR)

###############################################################################

%::
	@mkdir -p Junk
	touch Junk/$(notdir $@)

###############################################################################

-include $(LIBDEP)
-include $(TESTDEP)
