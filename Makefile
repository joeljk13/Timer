CC = g++
CFLAGS = -save-temps -std=c++11 -O0 -fwhole-program -ggdb3

c-test: test.c
	$(CC) $(CFLAGS) -o $@ $<

cpp-test: test.cpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm c-test cpp-test *.o *.ii *.s
