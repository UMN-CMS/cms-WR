CXXFLAGS := `root-config --cflags` -O2 -g -Wall

LDFLAGS := `root-config --libs` -lRooFitCore -lRooStats -lRooFit -lMinuit -lFoam

SOURCES := SimpleCardfile.cc doLimitSetting.C doLimitSettingMain.cc doLimitElectron.cc doLimitElectronMain.cc

HEADERS := SimpleCardfile.hh

%.o: %.cc
	${CXX} -c $*.cc -o $@ ${CXXFLAGS} 

%.o: %.C
	${CXX} -c $*.C -o $@ ${CXXFLAGS} -DMAKE_MAIN

all : exe

exe: doLimitSetting.exe doLimitElectron.exe

doLimitSetting.exe: SimpleCardfile.o doLimitSetting.o doLimitSettingMain.o
	${CXX} -g -o $@ $^  ${LDFLAGS}

doLimitElectron.exe: SimpleCardfile.o doLimitElectron.o doLimitElectronMain.o
	${CXX} -g -o $@ $^  ${LDFLAGS}

doLimitSetting.o: doLimitSetting.hh

doLimitSettingMain.o: SimpleCardfile.hh doLimitSetting.hh

doLimitElectron.o: doLimitElectron.hh

doLimitElectronMain.o: SimpleCardfile.hh doLimitElectron.hh

SimpleCardfile.o: SimpleCardfile.hh
