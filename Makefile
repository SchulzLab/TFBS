OBJECTS=Matrix.o Reader.o Controller.o

PROGNAME=Reader

CXXFLAGS= -std=c++1y -ggdb -Wall -pedantic -pipe

.PHONY: all clean

all: ${PROGNAME}

%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $<

${PROGNAME}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o $@ 

clean:
	rm -f ${OBJECTS} ${PROGNAME}
