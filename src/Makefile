#CC=gcc
CC=mpicc
#CC=clang
CFLAGS=-I./
OFLAGS=-O2 -ffast-math -march=native -mfpmath=sse  -mmmx -mavx 
LDFLAGS = -lm # bibliotheque math 
DEPS=constants.h rand.h read_posff.h \
     cell.h config.h md.h kinetic.h thermo.h
OBJ=mdff20.o constants.o read_posff.o cell.o config.o\
    nmlj.o integration.o field.o md.o rand.o cell.o thermo.o tools.o kinetic.o


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OFLAGS)

mdff20: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(OFLAGS) $(LDFLAGS)

clean:
	rm -rf *.o
