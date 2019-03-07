#-*- mode: makefile; mode: font-lock; vc-back-end: Git -*-

SHELL = /bin/sh

# Comms version - set to serial for serial or mpi for message passing
COMMS_ARCH = mpi

# Where you want the binary
prefix     = $(HOME)
bindir     = $(prefix)/bin

# Set version of comms to use accordingly
ifeq ($(COMMS_ARCH),serial)
COMMS_OBJ = comms_serial.o
else
COMMS_OBJ = comms_mpi.o
endif

# Define objects in dependency order
OBJECTS   = constants.o timer.o userparams.o util.o random.o \
	    data_structures.o $(COMMS_OBJ) molint.o io.o init.o \
	    mc_moves.o main.o

F90       = mpif90
LD        = mpif90
FFLAGS    = -O3 -Wall -Wextra 

.PRECIOUS: %.o
.PHONY:  clean

%: %.o
%.o: %.f90
	$(F90) $(FFLAGS) -c -o $@ $<

%.o: %.F90
	$(F90) $(FFLAGS) -c -o $@ $<

%.o: %.f95
	$(F90) $(FFLAGS) -c -o $@ $<
%: %.o
	$(F90) $(FFLAGS) -o $@ $^

mc_water :  $(OBJECTS)

	$(LD) -o $(bindir)/mw_water_ls $(OBJECTS) $(FFLAGS) 

clean : 

	rm -f *.mod *.d *.il *.o work.*
	rm -f $(bindir)/mw_water_ls

