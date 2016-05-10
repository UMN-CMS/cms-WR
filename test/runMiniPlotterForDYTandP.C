void runMiniPlotterForDYTandP(){
	gROOT->ProcessLine(".L miniPlotterForDYTandP.C+");
	gROOT->ProcessLine("miniPlotterForDYTandP()");
	gROOT->ProcessLine(".q");
}
