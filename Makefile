CC 		= g++

OPT3 	= -O3
CFLAGS 	= -Wall

LDFLAGS = -lm -lquadmath

all: elfit

elfit: elfit1.cpp
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OPT3) -fopenmp $+ $(LDFLAGS) -o $@

clean:
	rm -rf elfit *.o *.out *.err *.prv *.pcf *.row *.sym

