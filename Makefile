# CC = gcc-11
CC = gcc
C_VERSION = c99
OPTIMIZE = O3
# DEBUG = -ggdb3
# PROF = -pg

# C_INCLUDE_PATH=$(PATH):/usr/local/include

CC_FLAGS = -Wall -Werror -pedantic -std=$(C_VERSION) -$(OPTIMIZE) $(DEBUG) $(PROF)

BINS_DIR = ..

TARGETS = key_search_optimized

all: $(TARGETS)

key_search_optimized: key_search_optimized.c sbox32.h
	$(CC) -DUSE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS=1 $(CC_FLAGS) -fopenmp -L/usr/local/lib -lm4ri -o kso key_search_optimized.c

clean:
	rm -f *.o
	cd $(BINS_DIR) && rm -f $(TARGETS)
