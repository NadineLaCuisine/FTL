CC=/usr/bin/g++
# CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread -pipe
LDFLAGS=-pthread -pipe


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=ftl

all: $(EXEC)

ftl: ftl.o  utils.o mapping.o
	$(CC) -o $@ $^ $(LDFLAGS)

ftl.o: main.cpp  utils.h mapping.h
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

mapping.o: mapping.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
