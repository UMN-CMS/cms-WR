void runMacroSandBox(){
	gROOT->ProcessLine(".L macroSandBox.C+");
	macroSandBox();
}
