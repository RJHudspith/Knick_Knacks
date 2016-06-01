CC=clang

INCDIR=$(PWD)/Include/
CFLAGS=-Wall -O3 -I$(INCDIR)
LDFLAGS=-lm
LIBDIRS=./

## create the usual SOURCES, DEPENDENCIES and HEADERS
SRC := $(foreach dir,$(LIBDIRS),$(wildcard $(dir)/*.c))
OBJ := $(patsubst %.c,%.o,$(SRC))
DEPS := $(INCDIR)/*.h
TILDED := $(foreach dir,$(LIBDIRS),$(wildcard $(dir)/*~))

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

INT: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean

clean:
	rm -v $(OBJ) $(TILDED) $(INCDIR)/*~ INT
