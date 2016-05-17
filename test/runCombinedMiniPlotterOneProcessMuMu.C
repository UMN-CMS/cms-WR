void runCombinedMiniPlotterOneProcessMuMu(){
	gROOT->ProcessLine(".L combinedMiniPlotterOneProcessMuMu.C+");
	gROOT->ProcessLine("combinedMiniPlotterOneProcessMuMu()");
}
