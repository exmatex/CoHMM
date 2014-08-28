ARCH := intel64
M_UNAME := $(shell uname -m)
ifeq ($(M_UNAME), i686)
ARCH := ia32
endif
#VPATH = ..
SRCDIR=src
##############CHARM##########
ifeq ($(SET), charm)
#CHARM_ROOT=$(HOME)/charm
ifeq (,$(CHARM_ROOT))
$(info Please establish Charmm++ environment variables before using this Makefile.)
$(info E.g. by running setting CHARM_ROOT pr run 'module load charm++')
$(error CHARM_ROOT is not set)
endif
CXXFLAGS=-DCHARM
CHARMBIN=$(CHARM_ROOT)/bin
CHARMINC=$(CHARM_ROOT)/include
CHARMC=$(CHARMBIN)/charmc
CXX=$(CHARMC)
OBJDIR=charm_obj
BINDIR=charm_bin
#############CNC############
else ifeq ($(SET), cnc) 
CXXFLAGS=-DCNC
CNC=icpc
CXX=$(CNC)
OBJDIR=cnc_obj
BINDIR=cnc_bin
#############OMP############
else ifeq ($(SET), omp)
CXXFLAGS=-DOMP
OMP=g++
CXX=$(OMP)
OBJDIR=omp_obj
BINDIR=omp_bin
else ifeq ($(SET), circle)
CXXFLAGS=-DCIRCLE
CXX=mpic++
LIBCIRCLELIBS=$(shell pkg-config --libs libcircle)
LIBCIRCLE_CFLAGS=$(shell pkg-config --cflags libcircle)
ifeq ($(LIBCIRCLELIBS), )
$(error Set LIBCIRCLELIBS or run 'module load libcircle' first)
endif
OBJDIR=circle_obj
BINDIR=circle_bin
else ifeq ($(SET), )
ifneq "$(MAKECMDGOALS)" "clean"
$(error please SET=cnc, SET=charm, SET=circle or SET=omp)
endif
endif

ifeq ($(CXX), $(CNC))
ifeq (,$(CNCROOT))
$(info Please establish CnC environment variables before using this Makefile.)
$(info E.g. by running cncvars.sh or cncvars.csh)
$(error CNCROOT is not set)
endif
endif

#HIREDIS_INCLUDES=$(HOME)/hiredis
ifeq ($(HIREDIS_INCLUDES), )
$(error Set HIREDIS_INCLUDES or run 'module load hiredis' first)
endif
HIREDISLIB=$(HIREDIS_INCLUDES)/../lib
HIREDISINC=$(HIREDIS_INCLUDES)/hiredis
HIREDIS_CFLAG=-I$(HIREDISINC)
HIREDIS_LDFLAG=-L$(HIREDISLIB) -lhiredis

ifeq ($(MKLROOT), )
$(error Set MKLROOT or run 'module load mkl' first)
endif
MKL=$(MKLROOT)
MKLINC=$(MKL)/include
MKLLIB=$(MKL)/lib/$(ARCH)
MKL_CFLAG= -I$(MKLINC) 

#MKL flags
ifeq ($(CXX), $(CHARMC))
  MKL_LDFLAG= -L$(MKLLIB) -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential
else ifeq ($(SET), circle) 
  MKL_LDFLAG= -L$(MKLLIB) -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential
else ifeq ($(CXX), $(OMP)) 
  MKL_LDFLAG= -L$(MKLLIB) -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential
  OMP_CFLAGS=-fopenmp
  OMP_LDFLAGS=-fopenmp
else ifeq ($(CXX), $(CNC)) 
  MKL_LDFLAG= -L$(MKLLIB) -mkl=sequential
  CNC_LDFLAG=-L$(CNCROOT)/lib/$(ARCH) -lcnc -ltbb -ltbbmalloc
  CNC_CFLAG= -D_DIST_ -I$(CNCROOT)/include -std=c++0x
  #CNC_CFLAG=-I$(CNCROOT)/include -std=c++0x
endif

#BOOST_INCLUDES=/projects/opt/boost/1.55.0/include
#BOOST_INCLUDES=$(HOME)/boost_1_54_0_build/include
ifeq ($(BOOST_INCLUDES), )
$(error Set BOOST_INCLUDES or run 'module load boost' first)
endif
BOOST_CFLAG=-I$(BOOST_INCLUDES)
#We use boost header only so far
#BOOSTLIB=$(HIREDIS_INCLUDES)/../lib
#BOOST_LDFLAG=-L$(BOOSTINC)

COMD=$(PWD)/COMD_lib
COMDINC=$(COMD)/src-lib
SUBDIRS=$(COMDINC)
ifeq ($(wildcard $(COMDINC)/CoMD_lib.h), )
$(error COMDINC=$(COMDINC) seems to point to the wrong directory)
endif
COMDLIB=$(COMD)
COMD_CFLAG=-I$(COMDINC)
COMD_LDFLAG=-L$(COMDLIB) -lCoMD_2D
#COMD_LDFLAG=-L$(COMDLIB) -lcomd

OBJS:=$(addprefix $(OBJDIR)/, 2DKriging.o kriging.o flux.o redisBuckets.o output.o input.o)
CXXFLAGS+=$(HIREDIS_CFLAG) $(MKL_CFLAG) $(COMD_CFLAG) $(BOOST_CFLAG) -g
LDFLAGS=$(HIREDIS_LDFLAG) $(MKL_LDFLAG) $(COMD_LDFLAG) -lm -lrt
ifeq ($(CXX), $(CHARMC))
$(info compiling Charm files)
OBJS+=$(OBJDIR)/main_charm.o
else ifeq ($(CXX), $(CNC)) 
$(info compiling CnC files)
OBJS+=$(OBJDIR)/main_generic.o
LDFLAGS+=$(CNC_LDFLAG) 
CXXFLAGS+=$(CNC_CFLAG)
else ifeq ($(CXX), $(OMP)) 
$(info compiling OpenMP files)
OBJS+=$(OBJDIR)/main_generic.o
CXXFLAGS+=$(OMP_CFLAGS)
LDFLAGS+=$(OMP_LDFLAGS)
else ifeq ($(SET), circle)
OBJS+=$(OBJDIR)/main_generic.o
CXXFLAGS+=$(LIBCIRCLE_CFLAGS)
LDFLAGS+=$(LIBCIRCLELIBS) -lboost_serialization
endif

#target
NAME=$(BINDIR)/2D_Kriging
default: all
all: $(SUBDIRS) $(OBJDIR) $(NAME)

ifeq ($(CXX), $(CHARMC))
$(OBJDIR)/%.d: $(SRCDIR)/%.cpp
	@#1. sed:  put one file per line * * -> *\\\n*
	@#2. sed gcc -MG does not know that missing files will be in $(SRCDIR)
	@# no path -> SRCDIR
	@echo g++ -MM -MG -MT $(OBJDIR)/$*.o $(CXXFLAGS) -I$(SRCDIR) -I$(CHARMINC) $< \> $@
	@g++ -MM -MG -MT $(OBJDIR)/$*.o $(CXXFLAGS) -I$(SRCDIR) -I$(CHARMINC) $< | \
	sed 's/\([^[:space:]]\) \([^\\[:space:]]\)/\1 \\\n \2/g' | \
	sed '/^[^/]*.\(def\|decl\)\./s@[^[:space:]]@$(SRCDIR)/&@' > $@

#Charmm++ ci files
$(SRCDIR)/%.decl.h $(SRCDIR)/%.def.h: $(SRCDIR)/%.ci
	@#charmc writes to pwd only
	cd $(<D) && $(CXX) $(<F)
else
#deps rule for c files
$(OBJDIR)/%.d: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MM -MF $@ -MT $(OBJDIR)/$*.o $<
endif

DEPS=$(OBJS:.o=.d)
ifneq "$(MAKECMDGOALS)" "clean"
-include $(DEPS)
endif

#make subdirs for objects and executable
$(NAME): | $(BINDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJS): | $(OBJDIR)
$(DEPS): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

##--- Executable ---##

$(NAME): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)
ifeq ($(CXX), $(CHARMC))
	mv charmrun $(BINDIR)/
endif

#GNU make implicit rule
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

.PHONY: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) $(MFLAGS) -C $@

subdirclean:
	@for i in $(SUBDIRS); do \
	  echo $(MAKE) $(MFLAGS) -C $$i clean; \
	  $(MAKE) $(MFLAGS) -C $$i clean || exit 1; \
	done

clean: subdirclean
	rm -f $(SRCDIR)/*.decl.h $(SRCDIR)/*.def.h charmrun
	rm -f *.vtk *.dat core.* 
	rm -f $(OBJS) $(DEPS) $(NAME) $(OBJDIR)/main_*.[od]

