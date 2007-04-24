###########################################################
#
# Example GNUmakefile to build an Application with Marlin
#
# @author Frank Gaede
# @author Jan Engels
# @version $Id: GNUmakefile,v 1.12 2007-04-24 10:46:28 gaede Exp $
#
###########################################################


#----------------------------------------------------------------------------
#        user variables
#----------------------------------------------------------------------------

# name of the program
PROGNAME = MarlinUtil

# source files for the binary
# leave this commented out if your Application does not produce binary files
BINSRCS = testmarlinutil.cc

# marlin working directory
ifndef MARLINWORKDIR
 MARLINWORKDIR=$(MARLIN)
 export MARLINWORKDIR
endif

# additional user include/libraries
ifndef USERINCLUDES
 USERINCLUDES += -I$(MARLIN)/include
 USERINCLUDES += -D USE_CLHEP -I$(CLHEP)/include
 USERINCLUDES += -D USE_GEAR -I$(GEAR)/src/cpp/include
 USERINCLUDES += -I$(GSL_HOME)/include
 USERINCLUDES += -I$(LCIO)/src/cpp/include
endif

ifndef USERLIBS
 USERLIBS += -L$(MARLIN)/lib -lMarlin
 USERLIBS += -L$(CLHEP)/lib -lCLHEP
 USERLIBS += -L$(GEAR)/lib -lgear -lgearxml
 USERLIBS += -L$(GSL_HOME)/lib -lgsl -lgslcblas
 USERLIBS += -L$(LCIO)/lib -llcio
 USERLIBS += -L$(LCIO)/sio/lib -lsio -lz
endif

#----------------------------------------------------------------------------
#        do not modify anything below in order to use top level makefile
#----------------------------------------------------------------------------

# directories
INCDIR = ../include
LIBDIR = ../lib
BINDIR = $(MARLINWORKDIR)/bin
DOCDIR = $(MARLINWORKDIR)/doc/$(PROGNAME)
TMPDIR = $(MARLINWORKDIR)/tmp/$(PROGNAME)
OBJDIR = $(TMPDIR)/obj
DEPDIR = $(TMPDIR)/dep

# library file without path
LIBNAME = lib$(PROGNAME).a

# library file with path
LIB = $(LIBDIR)/$(LIBNAME)

ifdef BINSRCS
# binary files with path
BIN = $(patsubst %.cc,$(BINDIR)/%,$(BINSRCS))

# object files for the binaries without path
BINOBJS = $(patsubst %.cc,%.o,$(BINSRCS))

# object files for the binary with path
BINOBJ = $(patsubst %,$(OBJDIR)/%,$(BINOBJS))
endif

# source files without path
SRCS = $(filter-out $(BINSRCS),$(wildcard *.cc))

# object files without path
OBJS = $(patsubst %.cc,%.o,$(SRCS))

# object files with path
OBJ = $(patsubst %,$(OBJDIR)/%,$(OBJS))

# dependency files without path
DEPS = $(patsubst %.o,%.d,$(OBJS) $(BINOBJS))

# dependency files with path
DEP = $(patsubst %,$(DEPDIR)/%,$(DEPS))


#---- include ANN library: http://www.cs.umd.edu/~mount/ANN/
ANNSRCS = $(wildcard ann/*.cpp)
ANNOBJS = $(patsubst ann/%.cpp,%.o,$(ANNSRCS))
ANNOBJ = $(patsubst %,$(OBJDIR)/%,$(ANNOBJS))
ANNDEPS = $(patsubst %.o,%.d,$(ANNOBJS))
ANNDEP = $(patsubst %,$(DEPDIR)/%,$(ANNDEPS))

OBJ += $(ANNOBJ)
DEP += $(ANNDEP)

#bla: 
#	@echo $(ANNSRCS) 
#	@echo $(ANNOBJS) 
#	@echo $(ANNOBJ) 
#	@echo $(ANNDEPS) 
#	@echo $(ANNDEP) 

# c++ compiler
ifdef MY_CXX
  CXX = $(MY_CXX)
else
  CXX = g++
endif

# compiler flags
ifdef MARLINDEBUG
 CXXFLAGS = -g -Wall -ansi
else
 CXXFLAGS = -O3 -Wall -ansi
endif

# pre-processor flags
CPPFLAGS = -I$(INCDIR)
CPPFLAGS += $(USERINCLUDES)
CPPFLAGS += `. $(MARLINWORKDIR)/bin/marlin_includes.sh`

# libraries
#LIBS = -Wl,-whole-archive -L $(LIBDIR) -l$(PROGNAME) -Wl,-no-whole-archive
LIBS += $(USERLIBS)
LIBS += `. $(MARLINWORKDIR)/bin/marlin_libs.sh`

# linker
LD = $(CXX)

# linker flags
LDFLAGS = -Wl,--demangle

# archiever
AR = ar cru

# check if a directory exists
DIR_EXISTS = test -d

.PHONY: all lib doc clean distclean

all: lib

ifdef BINSRCS
.PHONY: bin
bin: lib $(BINOBJ) $(BIN)
endif

lib: $(LIB)

doc:
	@$(DIR_EXISTS) $(DOCDIR) || mkdir -p $(DOCDIR) ; \
	export DOX_PROGNAME=$(PROGNAME) ; \
	export DOX_DOCDIR=$(DOCDIR) ; \
	cd ../doc ; doxygen ;

clean:
	rm -f $(OBJ) $(BINOBJ) $(BIN) $(LIB)

distclean:
	rm -rf $(BIN) $(LIB) $(TMPDIR) $(DOCDIR)

# tell make to include dependencies
-include $(DEP)

# binary rule
$(BIN) : $(OBJ) $(BINOBJ)
	@$(DIR_EXISTS) $(BINDIR) || mkdir -p $(BINDIR)
	$(LD) $(OBJDIR)/$(notdir $@).o $(LDFLAGS) $(LIBS) -o $@

# objects rule
$(OBJDIR)/%.o : %.cc
	@$(DIR_EXISTS) $(OBJDIR) || mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# objects rule
$(OBJDIR)/%.o : ann/%.cpp
	@$(DIR_EXISTS) $(OBJDIR) || mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# library rule
$(LIB) : $(OBJ) $(ANNOBJ)
	@$(DIR_EXISTS) $(LIBDIR) || mkdir -p $(LIBDIR)
	$(AR) $@ $?

# dependency rule
$(DEPDIR)/%.d : %.cc 
	@$(DIR_EXISTS) $(DEPDIR) || mkdir -p $(DEPDIR)
	@echo "Creating dependencies for $<" ; \
	set -e; $(CXX) -MM $(CPPFLAGS) $< \
	| sed 's|\($*\)\.o[ :]*|$(OBJDIR)/\1.o $(DEPDIR)/$*.d : |g'\
	> $@; [ -s $@ ] || rm -f $@

# end
