#include "commonRooFitBkgndEstimation.C"
#include "../interface/rooFitFxns.h"

void fitEMuTTBarStudyBinningEMuRatios(){
	TChain * ttbaremu = new TChain("Tree_Iter0");
	ttbaremu->Add("/afs/cern.ch/work/s/skalafut/public/WR_starting2015/forPeterRootFiles/selected_tree_TT_flavoursidebandEMu.root");

}///end fitEMuTTBarStudyBinningEMuRatios()

