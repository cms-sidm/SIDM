# Makefile for MSTW Fortran and C++ code.
# Comments to Graeme.Watt(at)cern.ch.

# If using g77, uncomment these two lines:
#FC=g77
#LDLIBS=-lg2c

# If using gfortran, uncomment these two lines:
FC=gfortran
LDLIBS=-lgfortran

# C and C++ compilers:
CC=gcc
CXX=g++

# Compiler flags:
FFLAGS=-Wall -W
CXXFLAGS=-Wall -W

EXEC=example example_cpp

all: $(EXEC)

#####################################
# Fortran version:

example : example.o mstwpdf.o alphaS.o
	$(FC) $^ -o $@

example.o : example.f
	$(FC) -c $(FFLAGS) $< -o $@

mstwpdf.o : mstwpdf.f
	$(FC) -c $(FFLAGS) $< -o $@

alphaS.o : alphaS.f
	$(FC) -c $(FFLAGS) $< -o $@

#####################################
# C++ version:

example_cpp : example_cpp.o mstwpdf_cpp.o alphaS.o
	$(CXX) $^ -o $@ $(LDLIBS)

# Or without linking to the Fortran alphaS code.
# (and comment out relevant parts of example.cc.)
#example_cpp : example_cpp.o mstwpdf_cpp.o
#	$(CXX) $^ -o $@

example_cpp.o : example.cc mstwpdf.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

mstwpdf_cpp.o : mstwpdf.cc mstwpdf.h
	$(CXX) -c $(CXXFLAGS) $< -o $@

#####################################
# Mathematica version:

# Build MathLink interface to Fortran alphaS code.
# First try "mcc -o alphaS.exe alphaS.c alphaS.tm alphaS.f -lg2c".
# (Use -lg2c with g77 and -lgfortran with gfortran.)
# If this doesn't work, compile in steps as shown below.
# Set variables (for Linux or Mac OS X) then do "make alphaS.exe".
# See MathLink documentation for more information.

# Linux (example shown for lxplus.cern.ch, modify as appropriate).
# Find next two variables with $InstallationDirectory and $SystemID.
#MDIR=/afs/cern.ch/project/parc/math70
#SYS=Linux-x86-64
#MLDIR=$(MDIR)/SystemFiles/Links/MathLink/DeveloperKit/$(SYS)/CompilerAdditions
#MLLIB=ML64i3 # set to ML32i3 if using a 32-bit system
#EXTRALIBS=-lm -lpthread -lrt

# Mac OS X (example shown for $SystemID=MacOSX-x86, modify as appropriate).
MDIR=/Applications/Mathematica.app
MLDIR=$(MDIR)/SystemFiles/Links/MathLink/DeveloperKit/CompilerAdditions
MLLIB=MLi3
EXTRALIBS=

alphaS.exe : alphaS.c alphaS.tm alphaS.f
	$(MLDIR)/mprep alphaS.tm -o alphaStm.c
	$(CC) -I$(MLDIR) -c alphaStm.c -o alphaStm.o
	$(CC) -I$(MLDIR) -c alphaS.c -o alphaSc.o
	$(FC) -c alphaS.f -o alphaS.o
	$(CXX) -o alphaS.exe alphaS.o alphaSc.o alphaStm.o -L$(MLDIR) -l$(MLLIB) $(LDLIBS) $(EXTRALIBS)

#####################################

clean :
	-rm -f *.o alphaStm.c alphaS.exe $(EXEC)
