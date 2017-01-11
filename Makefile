# Usage: -e DBG=1 | -e GPROF=1

# Compiler/Linker variables and default values
CC=gcc
CPP=g++
IFLAGS=
#CFLAGS=-Wall -O3 -ffast-math -fexpensive-optimizations -funroll-all-loops -finline-functions -pedantic -std=c99
CFLAGS=-Winline -Wall -O3 -msse4 -msse3 -msse2 -funroll-all-loops -fexpensive-optimizations -finline-functions -pedantic -std=c99 -I. -fopenmp
CXXFLAGS=-msse4 -std=c++11 -O2 -I. -fopenmp
# -funroll-all-loops 
LDFLAGS=-lm
INSTALL_DIR=/usr/local/bin
EXE_NAME=pindel

# Handle specific installation choices, first two redefine default flags, and so must be first
# -e DBG=1	: for debugging version
ifdef DBG
	CFLAGS=-g -msse2 -msse3 -msse4 -I. -fno-inline -std=c99 -fopenmp
	CXXFLAGS=-msse4 -std=c++11 -g -I.
endif
# -e GPROF=1	: for profiling with gprof
ifdef GPROF
	CFLAGS=-g -msse2 -msse3 -msse4 -fno-inline -pg -I. -std=c99
	CXXFLAGS=-msse4 -std=c++11 -g -pg -I.
endif

# Other Directives
# -e VERSION=2.0: for setting version
ifndef VERSION
	VERSION=1.0
endif

VPATH = src/
CFLAGS := $(CFLAGS) $(IFLAGS)

# Local variables
srcs = $(wildcard *.c)
hds = $(wildcard *.h)
objs = $(srcs:.c=.o) nbhd.o
# deps = $(srcs:.c=.d)

$(EXE_NAME) : $(objs)
	$(CPP) -I. -fopenmp -o $(EXE_NAME) $(objs) $(LDFLAGS)

nbhd.o: nbhd.cpp nbhd.h
	$(CPP) -c $(CXXFLAGS) -o nbhd.o nbhd.cpp

.PHONY : clean 

clean:
	rm -f *.o *.d $(EXE_NAME)
