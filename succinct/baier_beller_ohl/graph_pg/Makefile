ifndef SDSLLITE
    $(error SDSLLITE is undefined)
endif

CC       = g++
CC_FLAGS = -O9 -std=c++11 -DNDEBUG -ffast-math -funroll-loops -msse4.2 -Wall -I${SDSLLITE}/include -L${SDSLLITE}/lib
CCLIB    = -lsdsl -ldivsufsort -ldivsufsort64
SOURCES  =$(wildcard *.cpp)
EXECS    =$(SOURCES:.cpp=)
TARGETS  =$(SOURCES:.cpp=.x)

all: $(TARGETS)
%.x:	%.cpp
	$(CC) $(CC_FLAGS) $< -o $@ $(CCLIB)

clean:
	rm -f $(TARGETS)
