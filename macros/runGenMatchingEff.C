void runGenMatchingEff(){
	gROOT->ProcessLine(".L genMatchingEfficiency.C+");
	genMatchingEfficiency();
}
