CC=gcc
CFLAGS=

projecteuler:
	$(CC) $(CFLAGS) -o ./projecteuler ./projecteuler.c

run: projecteuler
	./projecteuler

clean:
	rm ./projecteuler

.PHONY: all run clean
