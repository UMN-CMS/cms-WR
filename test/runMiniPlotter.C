void runMiniPlotter(){
	gROOT->ProcessLine(".L miniPlotter.C+");
	gROOT->ProcessLine("miniPlotter()");
}
