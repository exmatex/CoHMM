ARCH := intel64
M_UNAME := $(shell uname -m)
ifeq ($(M_UNAME), i686)
ARCH := ia32
endif
#VPATH = ..
##############CHARM##########
ifeq ($(SET), charm)
OPTS=-DCHARM
CHARMDIR=/home/droehm/charm
CHARMC=$(CHARMDIR)/bin/charmc $(OPTS)
CXX=$(CHARMC)
#############CNC############
else ifeq ($(SET), cnc) 
OPTS=-DCNC
CNC=icpc $(OPTS)
CXX=$(CNC)
#############OMP############
else ifeq ($(SET), omp)
OPTS=-DOMP
OMP=g++ $(OPTS)
CXX=$(OMP)

else ifeq ($(SET), ) 
$(error please SET=cnc, SET=charm or SET=omp)
endif

ifeq ($(CXX), $(CNC))
ifeq (,$(CNCROOT))
$(info Please estblish CnC environment variables before using this Makefile.)
$(info E.g. by running cncvars.sh or cncvars.csh)
$(error CNCROOT is not set)
endif
endif

#Darwin Flags
HIREDIS=/home/droehm/hiredis
HIREDISLIB=$(HIREDIS)
HIREDISINC=$(HIREDIS)
HIREDIS_CFLAG=-I$(HIREDISINC)
HIREDIS_LDFLAG=-L$(HIREDISLIB) -lhiredis

MKL=/projects/opt/intel/compilers/composer_xe_2013.4.183/mkl
MKLINC=$(MKL)/include
MKLLIB=$(MKL)/lib/intel64
MKL_CFLAG= -I$(MKLINC) 

#MKL flags
ifeq ($(CXX), $(CHARMC))
  MKL_LDFLAG= -L$(MKLLIB) -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential
else ifeq ($(CXX), $(OMP)) 
  OPENMP_FLAG=-fopenmp
  MKL_LDFLAG= -L$(MKLLIB) -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential
else ifeq ($(CXX), $(CNC)) 
  MKL_LDFLAG= -L$(MKLLIB) -mkl=sequential
  CNC_LDFLAG=-L$(CNCROOT)/lib/$(ARCH) -lcnc -ltbb -ltbbmalloc
  #CNC_CFLAG=-I$(CNCROOT)/include -D_DIST_ -std=c++0x
  CNC_CFLAG=-I$(CNCROOT)/include -std=c++0x
endif

#BOOST=/projects/opt/boost/1.55.0
BOOST=/home/droehm/boost_1_54_0_build
BOOSTINC=$(BOOST)/include
BOOSTLIB=$(BOOST)/lib
BOOST_CFLAG=-I$(BOOSTINC)
BOOST_LDFLAG=-L$(BOOSTINC)

COMD=/home/droehm/2014/CNC/COMD_lib
COMDINC=$(COMD)/src-lib
COMDLIB=$(COMD)
COMD_CFLAG=-I$(COMDINC)
COMD_LDFLAG=-L$(COMDLIB) -lCoMD_2D
#COMD_LDFLAG=-L$(COMDLIB) -lcomd

ifeq ($(CXX), $(CHARMC))
LDFLAGS= $(HIREDIS_LDFLAG) $(MKL_LDFLAG) $(COMD_LDFLAG) -lm -lrt
CXXFLAGS= $(HIREDIS_CFLAG) $(MKL_CFLAG) $(COMD_CFLAG) $(BOOST_CFLAG) -g $(LDFLAGS) 
else ifeq ($(CXX), $(CNC)) 
LDFLAGS= $(HIREDIS_LDFLAG) $(MKL_LDFLAG) $(COMD_LDFLAG) $(CNC_LDFLAG) -lm -lrt
CXXFLAGS= $(HIREDIS_CFLAG) $(MKL_CFLAG) $(COMD_CFLAG) $(BOOST_CFLAG) $(CNC_CFLAG) -g  
else ifeq ($(CXX), $(OMP)) 
LDFLAGS= $(HIREDIS_LDFLAG) $(MKL_LDFLAG) $(COMD_LDFLAG) -lm -lrt $(OPENMP_FLAG)
CXXFLAGS= $(HIREDIS_CFLAG) $(MKL_CFLAG) $(COMD_CFLAG) $(BOOST_CFLAG) -g 
endif

#target

default: all
all: 2D_Kriging

##--- Executable ---##

ifeq ($(CXX), $(CHARMC))
$(info compiling Charm files)

2D_Kriging: 2DKriging.o main.o kriging.o flux.o redisBuckets.o output.o input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

##--- Main Chare ---##

main.o : main.C main.h main.decl.h main.def.h krigingMod.decl.h input.hpp

main.decl.h main.def.h : main.ci
	$(CXX) main.ci

krigingMod.decl.h krigingMod.def.h : krigingMod.ci
	$(CXX) krigingMod.ci

2DKriging.o: 2DKriging.cpp krigingMod.decl.h krigingMod.def.h main.decl.h main.h

else ifeq ($(CXX), $(CNC)) 
$(info compiling CnC files)

2D_Kriging: 2DKriging.o main_cnc.o kriging.o flux.o redisBuckets.o output.o input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

main_cnc.o : main_cnc.cpp main_cnc.hpp input.hpp

2DKriging.o: 2DKriging.cpp 2DKriging.hpp

else ifeq ($(CXX), $(OMP)) 
$(info compiling OpenMP files)

2D_Kriging: 2DKriging.o main_cnc.o kriging.o flux.o redisBuckets.o output.o input.o
	$(CXX) -o $@ $^ $(LDFLAGS)

main_cnc.o : main_cnc.cpp main_cnc.hpp input.hpp

2DKriging.o: 2DKriging.cpp 2DKriging.hpp

endif

flux.o: flux.cpp

input.o: input.cpp

output.o: output.cpp

kriging.o: kriging.cpp

redisBuckets.o: redisBuckets.cpp

clean:
	rm -f ./*decl.h ./*def.h ./*.o ./*.vtk ./*.dat ./core.* ./2D_Kriging ./charmrun

