CC=gcc
CFLAGS=

projecteuler:
	$(CC) $(CFLAGS) -lm -o ./projecteuler ./projecteuler.c

run: projecteuler
	./projecteuler

clean:
	rm ./projecteuler

.PHONY: all run clean
