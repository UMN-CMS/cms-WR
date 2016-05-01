void runCombinedMiniPlotterMuMu(){
	gROOT->ProcessLine(".L combinedMiniPlotterMuMu.C+");
	gROOT->ProcessLine("combinedMiniPlotterMuMu()");
	gROOT->ProcessLine(".q");
}
