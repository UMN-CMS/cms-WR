void runDumpTreePlots(){
	gROOT->ProcessLine(".L dumpTreePlots.C+");
	dumpTreePlots();
}
