void runCombinedMiniPlotter(){
	gROOT->ProcessLine(".L combinedMiniPlotter.C+");
	gROOT->ProcessLine("combinedMiniPlotter()");
}
