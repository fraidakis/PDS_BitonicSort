CC = mpicc
CFLAGS = -I/usr/lib/x86_64-linux-gnu/openmpi/include
LDFLAGS = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi

all: start

start: start.o
    $(CC) -o start start.o $(LDFLAGS)

start.o: start.c
    $(CC) $(CFLAGS) -c start.c

clean:
    rm -f start start.o