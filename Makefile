# (ON/OFF)
DO_OPENMP = OFF
# Note that the (clang) C compiler in OS X
# Xcode 5.x does not yet include OpenMP.
# On OS X, one may instead use gcc provided
# by Homebrew (http://brew.sh/).

# (ON/OFF)
DO_GNUPLOT = OFF
# If enabled, Gnuplot will be used to
# generate postscript plots in 'output/'.

CC = cc 
CFLAGS = -ICoMDLib
LDFLAGS = -LCoMDLib -lCoMD


ifeq ($(DO_GNUPLOT), ON)
CFLAGS += -DGNUPLOT='"/usr/bin/gnuplot -persist"'
endif

ifeq ($(DO_OPENMP), ON)
CFLAGS += -fopenmp
endif

all : cohmm

cohmm : cohmm.c
	$(MAKE) -C CoMDLib
	${CC} cohmm.c $(CFLAGS) $(LDFLAGS) -lm -g -O3 -o cohmm

clean :
	rm -f cohmm
	rm -fr cohmm.dSYM
	$(MAKE) -C CoMDLib clean
