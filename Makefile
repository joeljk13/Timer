OPTIMIZATION = -Os -fwhole-program

PARAMS =

LINK =

SOURCES = test.cpp

OBJECTS = $(SOURCES:.cpp=.o)

# All objects will depend on all headers, for simplicity
HEADERS = timer.hpp

EXECNAME = Test

CC = g++

ERR = -fmax-errors=3 -Wpedantic

CPARAMS = $(PARAMS) $(ERR) -std=c++11 -save-temps \
    -fverbose-asm # -masm=intel (Doesn't work with boost)

LPARAMS = $(PARAMS)

timer: $(EXECNAME)
	:

$(EXECNAME): $(OBJECTS)
	$(CC) -o $@ $(OPTIMIZATION) $(OBJECTS) $(LPARAMS) $(LINK)

%.o: %.cpp $(HEADERS)
	$(CC) -c -o $@ $(CPARAMS) $(OPTIMIZATION) $<

clean:
	rm Test *.o *.ii *.s
