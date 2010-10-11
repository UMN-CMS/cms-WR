CXXFLAGS := `root-config --cflags` -O2 -g

LDFLAGS := `root-config --libs` -lRooFitCore -lRooStats -lRooFit -lMinuit -lFoam

SOURCES := SimpleCardfile.cc doLimitSetting.C doLimitSettingMain.cc

HEADERS := SimpleCardfile.hh

OBJECTS := $(patsubst %.C,%.o,$(patsubst %.cc,%.o,${SOURCES}))

%.o: %.cc
	${CXX} -c $*.cc -o $@ ${CXXFLAGS} 

%.o: %.C
	${CXX} -c $*.C -o $@ ${CXXFLAGS} -DMAKE_MAIN

all : exe

exe: doLimitSetting.exe

doLimitSetting.exe: ${OBJECTS} 
	${CXX} -g -o $@ $^  ${LDFLAGS}


doLimitSetting.o: doLimitSetting.hh

doLimitSettingMain.o: SimpleCardfile.hh doLimitSetting.hh

SimpleCardfile.o: SimpleCardfile.hh
