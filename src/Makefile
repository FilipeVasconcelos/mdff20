$(info -----------------------  MDFF20 compiling info ------------------------)
UNAME:=$(shell uname -a)   
HOST:=$(shell hostname)
$(info compiling on $(HOST))
$(info -----------------------------------------------------------------------)

CC=gcc
#CC=mpicc
#CC=clang
$(info CC  $(shell which $(CC) ))

CFLAGS=-O3 -msse4 -mavx2 -Wall
$(info CFLAGS   $(CFLAGS))
LDFLAGS = -lm # bibliotheque math 
$(info LDFLAGS   $(LDFLAGS))
INC=constants.h rand.h read_posff.h \
    cell.h config.h md.h kinetic.h thermo.h field.h
OBJ=mdff20.o io.o mpi_mdff.o constants.o \
    read_posff.o cell.o config.o\
    nmlj.o integration.o field.o md.o rand.o\
    cell.o thermo.o tools.o kinetic.o

BRANCH=$(shell git branch | grep "*" | awk '{print $$2}' )
$(info branch   $(BRANCH))       

# --------------------------------------------------------------------------- 
#   -DHOST      : hostname
#   -DBRANCH    : git branch
#   -DMPI       : 
#   possible debug/timing 
# --------------------------------------------------------------------------- 
ifeq ($(CC),mpicc)
    CPPDEFS=-DMPI -DHOST="$(HOST)" -DBRANCH="$(BRANCH)" 
else
    CPPDEFS=-DHOST="$(HOST)" -DBRANCH="$(BRANCH)"
endif
$(info CPPDEFS  $(CPPDEFS))

.PHONY: clean

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(CPPDEFS) -o $@ -c $<

mdff20.x: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf *.o
$(info -----------------------------------------------------------------------)
