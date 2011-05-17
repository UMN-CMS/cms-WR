CXXFLAGS := `root-config --cflags` -O2 -g -Wall

LDFLAGS := `root-config --libs`

%.o: %.cc
	${CXX} -c $*.cc -o $@ ${CXXFLAGS} 

%.o: %.C
	${CXX} -c $*.C -o $@ ${CXXFLAGS} 

all : exe

exe: makeLimitFile2010.exe makeLimitFile2011.exe makeLimitFile2Y.exe

makeLimitFile2010.exe: makeLimitFile.o makeLimitFile2010.o
	${CXX} -g -o $@ $^  ${LDFLAGS}

makeLimitFile2011.exe: makeLimitFile.o makeLimitFile2011.o
	${CXX} -g -o $@ $^  ${LDFLAGS}

makeLimitFile2Y.exe: makeLimitFile.o makeLimitFile2Y.o
	${CXX} -g -o $@ $^  ${LDFLAGS}

makeLimitFile.o: systematics.h

