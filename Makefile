CC = gcc
CXX = g++
CFLAGS = -save-temps -std=gnu99 -O2 # -ggdb3
CXXFLAGS = -save-temps -std=c++11 -O2 -fwhole-program # -ggdb3

all: cpp-test c-test

c-test: test.c timer.h timer.c
	$(CC) $(CFLAGS) -o $@ timer.c $<

cpp-test: test.cpp timer.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f c-test cpp-test *.o *.ii *.i *.s
