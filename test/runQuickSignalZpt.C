void runQuickSignalZpt(){
	gROOT->ProcessLine(".L quickSignalZpt.C+");
	gROOT->ProcessLine("quickSignalZpt()");
	gROOT->ProcessLine(".q");
}
