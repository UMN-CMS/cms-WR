void runFlavorSideband(){
	gROOT->ProcessLine(".L flavorSideband.C");
	gROOT->ProcessLine("flavorSideband()");
	gROOT->ProcessLine(".q");
}
