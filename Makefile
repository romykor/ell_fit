CC 		= g++

OPT3 	= -O3
CFLAGS 	= -Wall

LDFLAGS = -lm -lquadmath

all: elfit4q

elfit: elfit4q.cpp
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OPT3) -fopenmp $+ $(LDFLAGS) -o $@

clean:
	rm -rf elfit4q *.o *.out *.err *.prv *.pcf *.row *.sym

