CC:=gcc
CFLAGS:=-O2 -Wall -g -pg -ffp-contract=off
shared: CFLAGS += -fPIC

navier: alloc.o boundary.o init.o main.o output.o simulation.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

shared: alloc.o boundary.o init.o shared.o output.o simulation.o
	$(CC) $(CFLAGS) -lm -shared -o libnavier.so $^

.PHONY: all
all:    navier

.PHONY: clean
clean:
	rm -f *.o navier libnavier.so

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

boundary.o       : datadef.h
init.o           : datadef.h
main.o           : alloc.h boundary.h datadef.h init.h simulation.h
simulation.o     : datadef.h init.h
shared.o     : datadef.h boundary.h init.h simulation.h