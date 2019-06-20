# makefile for this mini-project
#

CC = gcc
CFLAGS = -std=c99 -O3 -lm
RAPTOR = encoder.o galois.o bipartite.o decoder.o pivoting.o gaussian.o

raptorTest: $(RAPTOR) test.c
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


.PHONY: clean

clean:
	rm -f *.o raptorTest
