#g++ hybridInstructionalMacro.C -o hybridMacro -I/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/roofit/5.34.22-cms3/include/ -pthread -std=c++11 -Wno-deprecated-declarations -m64 -I/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0_pre9/external/slc6_amd64_gcc491/bin/../../../../../../lcg/root/5.34.22-cms3/include -m64 -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_0_pre9/external/slc6_amd64_gcc491/bin/../../../../../../lcg/root/5.34.22-cms3/lib -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lRooFit -lRooFitCore -lMinuit -lRooStats -L/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/roofit/5.34.22-cms3/lib

#use this makefile to compile hybridInstructionalMacro.C into an executable file named hybridMacro
CXXFLAGS=-pthread -std=c++11 -Wno-deprecated-declarations -m64
INCDIR=/cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/
ROOFIT_INC=roofit/6.02.00-odfocd5/
ROOT_INC=root/6.02.00-odfocd5/
LIBS=-lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -rdynamic -lRooFit -lRooFitCore -lMinuit -lRooStats

hybridMacro: hybridInstructionalMacro.C
	g++ $< -o $@ -I$(INCDIR)$(ROOFIT_INC)include $(CXXFLAGS) -I$(INCDIR)$(ROOT_INC)include -L$(INCDIR)$(ROOT_INC)lib $(LIBS) -L$(INCDIR)$(ROOFIT_INC)lib


