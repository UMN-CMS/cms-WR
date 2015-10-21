void runTrialRooDatasetMacro(){
	gROOT->ProcessLine(".L trialRooDatasetMacro.C");
	gROOT->ProcessLine("macro()");
}
