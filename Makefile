# whole tool
OBJECTSCOMP=Matrix.o Reader.o Controller.o svm.o libsvm_interface.o tfbs_main.o
PROGNAMECOMP=tfbs

# reader tool
OBJECTSREAD=Matrix.o Reader.o Controller_reader.o reader_main.o
PROGNAMEREAD=reader



CXXFLAGS= -std=c++0x -ggdb -Wall -pedantic -pipe -fopenmp

.PHONY: all clean

all: ${PROGNAMECOMP}

%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $<

${PROGNAMECOMP}: ${OBJECTSCOMP}
	${CXX} ${CXXFLAGS} ${OBJECTSCOMP} -o $@ 

${PROGNAMEREAD}: ${OBJECTSREAD}
	${CXX} ${CXXFLAGS} ${OBJECTSREAD} -o $@ 

clean:
	rm -f ${OBJECTSCOMP} ${PROGNAMECOMP} ${OBJECTSREAD} ${PROGNAMEREAD}
