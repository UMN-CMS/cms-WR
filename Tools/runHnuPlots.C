void runHnuPlots(){

	gSystem->Load("libFWCoreFWLite.so");
	AutoLibraryLoader::enable();
	gROOT->ProcessLine(".L hnuPlots.C++");
	executePlot2012();

}
