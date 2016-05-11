void runCalculateDyScaleFactors(){
	gROOT->ProcessLine(".L calculateDyScaleFactors.C+");
	gROOT->ProcessLine("calculateDyScaleFactors()");
	gROOT->ProcessLine(".q");
}
