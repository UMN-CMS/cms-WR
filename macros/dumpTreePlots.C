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

//using namespace std;

//#define DEBUG

///use this fxn to grab a histogram made in a TTree->Draw() call, draw this histo
///with a thicker line, and save it to a .png file 
void saveSingleHisto(TString canvName,TString histName,TString histTitle,TString outFile){
	Bool_t isPlottingGeV = false;
	
	if(histTitle.Contains("pt") || histTitle.Contains("Mass")) isPlottingGeV = true;

	gStyle->SetOptStat("emriou");
	TCanvas * canv = new TCanvas(canvName,canvName,600,600);
	canv->cd();
	TH1F * hist = (TH1F*) gROOT->FindObject(histName);
	hist->SetTitle(histTitle);
	hist->SetLineColor(1);
	hist->SetLineWidth(3);
	
	char temp[130];
	if(isPlottingGeV) sprintf(temp,"Events / %.2f GeV",hist->GetXaxis()->GetBinWidth(1));
	else sprintf(temp,"Events / %.2f ",hist->GetXaxis()->GetBinWidth(1));
	hist->GetYaxis()->SetTitle(temp);

	hist->Draw();
	canv->SaveAs(outFile,"recreate");
}///end saveSingleHisto()

///use this fxn to take one tuple and dump all branches into plots with unique names
void SaveTreePlots(TChain * chain, TString outputFileName){
	TObjArray * branches = chain->GetListOfBranches();
	TIter brItr(branches);
	for(TBranch * aBranch=(TBranch*) brItr.Next(); aBranch!=NULL; aBranch=(TBranch*) brItr.Next()){
		TString brName = aBranch->GetName();
#ifdef DEBUG
		std::cout<<"brName = \t"<< brName <<std::endl;
#endif
		if(brName.CompareTo("ptGenEle")==0 || brName.CompareTo("etaGenEle")==0 || brName.CompareTo("phiGenEle")==0 || 
				brName.CompareTo("ptGenJet")==0 || brName.CompareTo("etaGenJet")==0 || brName.CompareTo("phiGenJet")==0 ){
			///NOTE the max value of 2 is set knowing that the max number of array elements per entry is 2
			///if the tree structure changes to allow arrays with more than 2 entries, this loop upper limit value
			///will have to change!
			for(unsigned int i=0; i<2; i++){
				TString num = std::to_string(i);
				TString updatedBrName = brName+"["+num+"]"+">>"+brName+num+"Hist";
				chain->Draw(updatedBrName);

				///now updatedBrName has the array name and element of interest, either [0] or [1]
				saveSingleHisto(brName+num,brName+num+"Hist",brName+"["+num+"]",outputFileName+"_"+brName+"["+num+"]"+".png");

			}///end for loop

		}///end filter to do something different if the branch name corresponds to an array branch
		
		else{
			TString updatedBrName = brName+">>"+brName+"Hist";
			TString histName = brName+"Hist";
#ifdef DEBUG
			std::cout<<"updatedBrName = \t"<< updatedBrName <<std::endl;
#endif

			chain->Draw(updatedBrName);
			
			TString outFilePath = outputFileName+"_"+brName+".png";
#ifdef DEBUG
			std::cout<<"outFilePath = \t"<< outFilePath <<std::endl;
#endif
			saveSingleHisto(brName,histName,brName,outFilePath);
			
		}

	}///end loop over branches in TChain pointer named chain

}///end SaveTreePlots()


///use this macro to take a file with multiple directories, one tree per directory, and dump
///all of the plots in all trees into pdf files with unique names
void dumpTreePlots(){
	///chains made without pdgId matching
	/*
	TChain * noCuts = new TChain("genAnalyzerOne/genObjectsNoCuts","");
	noCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	TChain * ptEtaCuts = new TChain("genAnalyzerTwo/genObjectsWithPtEtaCuts","");
	ptEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	TChain * ptEtaDileptonMassCuts = new TChain("genAnalyzerThree/genObjectsWithPtEtaAndDileptonMassCuts","");
	ptEtaDileptonMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	*/

	///chains made with pdgId matching
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


	TString plotDir_noCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/noCuts";
	TString plotDir_withPtEtaCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/withPtEtaCuts";
	TString plotDir_withPtEtaDileptonMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/withPtEtaDileptonMassCuts";

	TString plotDir_matched_noCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_noCuts/noCuts";
	TString plotDir_matched_withPtEtaCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaCuts/withPtEtaCuts";
	TString plotDir_matched_withPtEtaDileptonMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaDileptonMassCuts/withPtEtaDileptonMassCuts";
	TString plotDir_matched_withPtEtaDileptonMassDrCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaDileptonMassDrCuts/withPtEtaDileptonMassDrCuts";
	TString plotDir_matched_withPtEtaDileptonMassDrFourObjMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaDileptonMassDrFourObjMassCuts/withPtEtaDileptonMassDrFourObjMassCuts";




	//SaveTreePlots(noCuts, plotDir_noCuts);
	//SaveTreePlots(ptEtaCuts, plotDir_withPtEtaCuts);
	//SaveTreePlots(ptEtaDileptonMassCuts, plotDir_withPtEtaDileptonMassCuts);
	
	/*SaveTreePlots(matchedNoCuts, plotDir_matched_noCuts);
	SaveTreePlots(matchedPtEtaCuts, plotDir_matched_withPtEtaCuts);
	SaveTreePlots(matchedPtEtaDileptonMassCuts, plotDir_matched_withPtEtaDileptonMassCuts);
	SaveTreePlots(matchedPtEtaDileptonMassDrCuts, plotDir_matched_withPtEtaDileptonMassDrCuts);
	SaveTreePlots(matchedPtEtaDileptonMassDrFourObjMassCuts, plotDir_matched_withPtEtaDileptonMassDrFourObjMassCuts);
	*/


	///chains made with deltaR matching btwn reco and gen
	TChain * matchedGenJetsToGenQuarksNoCuts = new TChain("matchGenJetsToGenQuarksNoCutsNewPath/matchedGenJetsNoCutsTree","");
	matchedGenJetsToGenQuarksNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	
	TChain * matchedRecoJetsToGenJetsNoCuts = new TChain("matchRecoJetsToGenJetsNoCutsNewPath/matchedRecoJetsNoCutsTree","");
	matchedRecoJetsToGenJetsNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	
	TChain * matchedRecoEleToLeadingGenEleNoCuts = new TChain("matchRecoEleToLeadingGenEleNoCutsNewPath/matchedLeadingRecoEleNoCutsTree","");
	matchedRecoEleToLeadingGenEleNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	TChain * matchedRecoEleToSubleadingGenEleNoCuts = new TChain("matchRecoEleToSubleadingGenEleNoCutsNewPath/matchedSubleadingRecoEleNoCutsTree","");
	matchedRecoEleToSubleadingGenEleNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	
	TChain * matchedRecoNoCuts = new TChain("matchedRecoAnalyzerOne/matchedRecoObjectsNoCuts","");
	matchedRecoNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	TChain * matchedRecoPtEtaCuts = new TChain("matchedRecoAnalyzerTwo/matchedRecoObjectsWithPtEtaCuts","");
	matchedRecoPtEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	TChain * matchedRecoPtEtaDileptonMassCuts = new TChain("matchedRecoAnalyzerThree/matchedRecoObjectsWithPtEtaAndDileptonMassCuts","");
	matchedRecoPtEtaDileptonMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	TChain * matchedRecoPtEtaDileptonMassDrCuts = new TChain("matchedRecoAnalyzerFour/matchedRecoObjectsWithPtEtaDileptonMassAndDrCuts","");
	matchedRecoPtEtaDileptonMassDrCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	TChain * matchedRecoPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");


	TString plotDir_reco_matched_noCuts_dR_genJetsToGenQuarks = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_genJetsToGenQuarks/noCuts_genJetsToGenQuarks";
	TString plotDir_reco_matched_noCuts_dR_recoElesToGenLeadingEles = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_recoElesToGenLeadingEles/noCuts_recoElesToGenLeadingEles";
	TString plotDir_reco_matched_noCuts_dR_recoElesToGenSubleadingEles = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_recoElesToGenSubleadingEles/noCuts_recoElesToGenSubleadingEles";
	TString plotDir_reco_matched_noCuts_dR_recoJetsToGenJets = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_recoJetsToGenJets/noCuts_recoJetsToGenJets";
	
	TString plotDir_reco_matched_noCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts/noCuts";
	TString plotDir_reco_matched_withPtEtaCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaCuts/withPtEtaCuts";
	TString plotDir_reco_matched_withPtEtaDileptonMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaDileptonMassCuts/withPtEtaDileptonMassCuts";
	TString plotDir_reco_matched_withPtEtaDileptonMassDrCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaDileptonMassDrCuts/withPtEtaDileptonMassDrCuts";
	TString plotDir_reco_matched_withPtEtaDileptonMassDrFourObjMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaDileptonMassDrFourObjMassCuts/withPtEtaDileptonMassDrFourObjMassCuts";

	///make plots of dR btwn gen jets and gen quarks, reco jets and gen jets, and reco eles and gen eles
	SaveTreePlots(matchedGenJetsToGenQuarksNoCuts,plotDir_reco_matched_noCuts_dR_genJetsToGenQuarks);
	SaveTreePlots(matchedRecoJetsToGenJetsNoCuts,plotDir_reco_matched_noCuts_dR_recoJetsToGenJets);
	SaveTreePlots(matchedRecoEleToLeadingGenEleNoCuts,plotDir_reco_matched_noCuts_dR_recoElesToGenLeadingEles);
	SaveTreePlots(matchedRecoEleToSubleadingGenEleNoCuts,plotDir_reco_matched_noCuts_dR_recoElesToGenSubleadingEles);

	///make plots of reco jets and eles matched to GEN objects after different levels of cuts
	///the matching is done before any cuts are applied at GEN or reco lvl
	SaveTreePlots(matchedRecoNoCuts, plotDir_reco_matched_noCuts);
	SaveTreePlots(matchedRecoPtEtaCuts, plotDir_reco_matched_withPtEtaCuts);
	SaveTreePlots(matchedRecoPtEtaDileptonMassCuts, plotDir_reco_matched_withPtEtaDileptonMassCuts);
	SaveTreePlots(matchedRecoPtEtaDileptonMassDrCuts, plotDir_reco_matched_withPtEtaDileptonMassDrCuts);
	SaveTreePlots(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, plotDir_reco_matched_withPtEtaDileptonMassDrFourObjMassCuts);



}///end dumpTreePlots()

