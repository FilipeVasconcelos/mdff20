$(info -----------------------  MDFF20 compiling info ------------------------)
UNAME:=$(shell uname -a)   
HOST:=$(shell hostname)
$(info compiling on $(HOST))
$(info -----------------------------------------------------------------------)
# installation directory
INSTALL_DIR=../bin
$(info installation dir $(shell readlink -f $(INSTALL_DIR)) )
# List of executable to be built within the package
EXE= mdff20.x

#CC=gcc
CC=mpicc
#CC=clang
$(info CC  $(shell which $(CC) ))

#CFLAGS=-O3 -msse4 -mavx2 -Wall -ffast-math -fopenmp
CFLAGS=-O2 -msse4 -mavx2 -ffast-math -fopenmp
$(info CFLAGS   $(CFLAGS))
LDFLAGS = -lm -fopenmp # bibliotheque math 
#LDFLAGS = -lm 
$(info LDFLAGS   $(LDFLAGS))
#INC= constants.h rand.h read_posff.h \
    cell.h config.h md.h kinetic.h thermo.h field.h

OBJ=mdff20.o io.o mpi_mdff.o timing.o constants.o \
    global.o pbc.o verlet.o config.o cell.o \
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
    CPPDEFS=-DMPI -DOMP -DHOST="$(HOST)" -DBRANCH="$(BRANCH)" 
else
    CPPDEFS=-DOMP -DHOST="$(HOST)" -DBRANCH="$(BRANCH)"
endif
#CPPDEFS+=-DDEBUG
$(info CPPDEFS  $(CPPDEFS))

.PHONY: clean

%.o: %.c $(DEPS)
	$(CC) -I./ $(CFLAGS) $(CPPDEFS) -o $@ -c $<

mdff20.x: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

install:
	cp $(EXE) $(INSTALL_DIR)

clean:
	rm -rf *.o


$(info -----------------------------------------------------------------------)
