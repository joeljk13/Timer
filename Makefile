CC = g++
CFLAGS = -save-temps -std=c++11 -O2 -fwhole-program # -ggdb3

all: cpp-test

c-test: test.c
	$(CC) $(CFLAGS) -o $@ $<

cpp-test: test.cpp timer.hpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f c-test cpp-test *.o *.ii *.s
