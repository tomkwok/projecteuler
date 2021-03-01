CC?=gcc
CFLAGS=
# CFLAGS=-DDEBUG
LFLAGS=-lm

projecteuler:
	$(CC) $(CFLAGS) -o ./projecteuler ./projecteuler.c $(LFLAGS)

run: projecteuler
	./projecteuler

clean:
	rm ./projecteuler

all: projecteuler

.PHONY: all run clean
