void runCombinedMiniPlotterForDYTandPMuMu(){
	gROOT->ProcessLine(".L combinedMiniPlotterForDYTandPMuMu.C+");
	gROOT->ProcessLine("combinedMiniPlotterForDYTandPMuMu()");
}
