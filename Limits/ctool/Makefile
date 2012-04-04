CXXFLAGS := `root-config --cflags` -O2 -g -Wall

LDFLAGS := `root-config --libs`

%.o: %.cc
	${CXX} -c $*.cc -o $@ ${CXXFLAGS} 

%.o: %.C
	${CXX} -c $*.C -o $@ ${CXXFLAGS} 

all : exe

exe: makeLimitFile.exe makeLimitFileELEC.exe

makeLimitFile.exe: makeLimitFile.o makeLimitFileMain.o systematics.o ratedb.o
	${CXX} -g -o $@ $^  ${LDFLAGS}


makeLimitFileELEC.exe: makeLimitFileELEC.o makeLimitFileMainELEC.o systematics.o
	${CXX} -g -o $@ $^  ${LDFLAGS}

makeLimitFile.o: systematics.h systematics.cc makeLimitFile.hh ratedb.hh

makeLimitFileMain.o: makeLimitFile.hh systematics.h

makeLimitFileELEC.o: systematics.h systematics.cc makeLimitFile.hh

