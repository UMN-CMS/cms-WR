void runCutEff(){
	gROOT->ProcessLine(".L cutEfficiency.C+");
	cutEfficiency();
}
