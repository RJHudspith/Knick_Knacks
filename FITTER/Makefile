## Makefile for the fit code

INCDIR=$(PWD)/HEADERS/
LDFLAGS=-lm -lgsl -lgslcblas -lfftw3
CFLAGS=-Wall -O3 -ffast-math -I$(INCDIR) -Werror -fopenmp
CC=gcc

## list the directories that are available
LIBDIRS=FITS IO ANALYSIS UTILS ./FITS/ALPHAS/ ./FITS/STATPOT/ ./FITS/CHIRAL_SUSY/ ./FITS/CORRELATORS/ ./

## create the usual SOURCES, DEPENDENCIES and HEADERS
SRC := $(foreach dir,$(LIBDIRS),$(wildcard $(dir)/*.c))
OBJ := $(patsubst %.c,%.o,$(SRC))
DEPS := $(INCDIR)/*.h
TILDED := $(foreach dir,$(LIBDIRS),$(wildcard $(dir)/*~))

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDFLAGS)

GLUFIT: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean

## should really do an expansion here ...
clean:
	rm -v $(OBJ) $(TILDED) $(INCDIR)/*~ GLUFIT ./LOCAL/* *.agr