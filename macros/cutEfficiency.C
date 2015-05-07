#include <TFile.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TBranch.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TObjArray.h>
#include "TCollection.h"
#include <TCut.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//#define DEBUG

void calcEffAndUnc(Float_t Nevts, Float_t kevts, Float_t & eff, Float_t & uncert){
	eff = kevts/Nevts;
	uncert = (1/Nevts)*sqrt(kevts*(1 - (kevts/Nevts) ));
}///end calcEffAndUnc()

/**
 * use this macro to calculate the efficiency of different cuts at reco or gen lvl
 */
void cutEfficiency(){
	///RECO chains and cut efficiences
	TChain * matchedRecoNoCuts = new TChain("matchedRecoAnalyzerOne/matchedRecoObjectsNoCuts","");
	matchedRecoNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");
	TChain * matchedRecoPtEtaCuts = new TChain("matchedRecoAnalyzerTwo/matchedRecoObjectsWithPtEtaCuts","");
	matchedRecoPtEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");
	TChain * matchedRecoPtEtaDileptonMassCuts = new TChain("matchedRecoAnalyzerThree/matchedRecoObjectsWithPtEtaAndDileptonMassCuts","");
	matchedRecoPtEtaDileptonMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");
	TChain * matchedRecoPtEtaDileptonMassDrCuts = new TChain("matchedRecoAnalyzerFour/matchedRecoObjectsWithPtEtaDileptonMassAndDrCuts","");
	matchedRecoPtEtaDileptonMassDrCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");
	TChain * matchedRecoPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");

	Float_t numEvtsWhereLeadingIsHardest = (Float_t) (matchedRecoPtEtaCuts->Draw("dijetMassGen","leadingIsHardest>0"));
	Float_t leadingIsHardestEff=0, dileptonMassEff=0;
	Float_t leadingIsHardestEffUnc=10, dileptonMassEffUnc=10;
	calcEffAndUnc((Float_t) (matchedRecoPtEtaCuts->GetEntriesFast()), numEvtsWhereLeadingIsHardest, leadingIsHardestEff, leadingIsHardestEffUnc);
	calcEffAndUnc(numEvtsWhereLeadingIsHardest,(Float_t) (matchedRecoPtEtaDileptonMassCuts->GetEntries()), dileptonMassEff, dileptonMassEffUnc);

	cout<<"leadingIsHardestEff = \t"<< leadingIsHardestEff <<"\t uncertainty = \t"<< leadingIsHardestEffUnc <<endl;
	cout<<"dileptonMassEff = \t"<< dileptonMassEff <<"\t uncertainty = \t"<< dileptonMassEffUnc <<endl;


	///GEN chains and cut efficiencies
	TChain * matchedNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	matchedNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaCuts = new TChain("matchedGenAnalyzerTwo/matchedGenObjectsWithPtEtaCuts","");
	matchedPtEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaDileptonMassCuts = new TChain("matchedGenAnalyzerThree/matchedGenObjectsWithPtEtaAndDileptonMassCuts","");
	matchedPtEtaDileptonMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaDileptonMassDrCuts = new TChain("matchedGenAnalyzerFour/matchedGenObjectsWithPtEtaDileptonMassAndDrCuts","");
	matchedPtEtaDileptonMassDrCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedGenAnalyzerFive/matchedGenObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	matchedPtEtaDileptonMassDrFourObjMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");

	Float_t numGenEvtsWhereLeadingIsHardest = (Float_t) (matchedPtEtaCuts->Draw("dijetMassGen","leadingIsHardest>0"));
	Float_t leadingIsHardestEffGen=0, dileptonMassEffGen=0;
	Float_t leadingIsHardestEffGenUnc=10, dileptonMassEffGenUnc=10;
	calcEffAndUnc((Float_t) (matchedPtEtaCuts->GetEntriesFast()), numGenEvtsWhereLeadingIsHardest, leadingIsHardestEffGen, leadingIsHardestEffGenUnc);
	calcEffAndUnc(numGenEvtsWhereLeadingIsHardest,(Float_t) (matchedPtEtaDileptonMassCuts->GetEntries()), dileptonMassEffGen, dileptonMassEffGenUnc);

	cout<<"leadingIsHardestEffGen = \t"<< leadingIsHardestEffGen <<"\t uncertainty = \t"<< leadingIsHardestEffGenUnc <<endl;
	cout<<"dileptonMassEffGen = \t"<< dileptonMassEffGen <<"\t uncertainty = \t"<< dileptonMassEffGenUnc <<endl;



}///end cutEffGeniciency()
