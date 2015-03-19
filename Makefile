OBJECTS=Matrix.o Reader.o Controller.o svm.o libsvm_interface.o

PROGNAME=reader

CXXFLAGS= -std=c++0x -ggdb -Wall -pedantic -pipe -fopenmp

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
