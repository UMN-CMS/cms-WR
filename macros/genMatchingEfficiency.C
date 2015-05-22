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

/**
 * input an int and double to define the number of bins, an int to define the minimum number of matched lower lvl objs , two TString
 * commands to pass to TTree::Draw(), two TString histo names, one TString canvas name, one TString
 * to define the output file path where the plot image is saved, one TChain * pointer, one TString efficiency
 * histo title, and one TString efficiency histo name
 *
 */
void calcAndPlotEff(Double_t nBins,Int_t nBinsInt, TString cutString, TString uncutHistPlotArgs, TString cutHistPlotArgs, TString uncutHistName, TString cutHistName, TString canvName, TString outputFile, TChain * chain, TString effHistName, TString effHistTitle){
	Float_t binLowEdges[nBinsInt+1];
	Double_t binContents[nBinsInt];
	for(Long64_t i=0; i<(nBinsInt+1); i++){
		binLowEdges[i] = -8;
	}///end initialization of bin low edges for efficiency histo 

	TCanvas * c999 = new TCanvas("t1","t1",500,500);
	c999->cd();
	///find the number of lower level objects which could be matched to a higher level object
	Double_t numPossibleMatches = (Double_t) (chain->Draw(uncutHistPlotArgs,""));
#ifdef DEBUG
	cout<<"numPossibleMatches = \t"<< numPossibleMatches << endl;
#endif
	c999->Update();

	TH1F * hist = (TH1F*) gROOT->FindObject(uncutHistName);
#ifdef DEBUG
	cout<<"cutString = \t"<< cutString << endl;
#endif

	TCanvas * c9999 = new TCanvas("t2","t2",500,500);
	c9999->cd();
	Double_t check = (Double_t) (chain->Draw(cutHistPlotArgs,cutString));
#ifdef DEBUG
	cout<<"this many evts passed the cut: \t" << check <<endl;
#endif
	c9999->Update();

	TH1F * selectedHist = (TH1F*) gROOT->FindObject(cutHistName);

	///now find the bin edge values which divide hist into nBinsInt bins, each bin with roughly the same number of entries
	Double_t cumulativeBinsVal=0;
	Double_t numeratorVal=0;	///< numeratorVal/cumulativeBinsVal = gen matching efficiency
	Double_t resetThreshold = pow(nBins,-1.0);
#ifdef DEBUG
	cout<<"resetThreshold = \t"<< resetThreshold <<endl;
#endif

	if(resetThreshold < 0.001) return;

	Int_t binNumIndex=0;
	for(Int_t bin=1; bin<=hist->GetNbinsX();bin++){
		if((cumulativeBinsVal/numPossibleMatches) >= resetThreshold){
#ifdef DEBUG
			cout<<"bin = \t"<<bin<<endl;
			cout<<"binContents[binNumIndex-1] = \t"<< binContents[binNumIndex-1] <<endl;
			cout<<"binLowEdges[binNumIndex] = \t"<< binLowEdges[binNumIndex] <<endl;
#endif

			binContents[binNumIndex-1] = (numeratorVal/cumulativeBinsVal);
			binLowEdges[binNumIndex] = hist->GetXaxis()->GetBinLowEdge(bin);
			binNumIndex++;
			cumulativeBinsVal=0;
			numeratorVal=0;
		}
		if(binNumIndex==0){
			binLowEdges[binNumIndex] = hist->GetXaxis()->GetBinLowEdge(bin);	///< when binNumIndex=0, bin=1
#ifdef DEBUG
			cout<<"binNumIndex = 0"<<endl;
			cout<<"binLowEdges[binNumIndex] = \t"<< binLowEdges[binNumIndex] <<endl;
#endif
	
			binNumIndex++;
		}
		
		cumulativeBinsVal += hist->GetBinContent(bin);
		numeratorVal += selectedHist->GetBinContent(bin);

		if(bin==hist->GetNbinsX()){
			binContents[binNumIndex-1] = (numeratorVal/cumulativeBinsVal);
			binLowEdges[binNumIndex] = hist->GetXaxis()->GetBinLowEdge(bin);
		}
	
	}///end loop over bins in hist

#ifdef DEBUG
	cout<<"filled binContents array"<<endl;
#endif


#ifdef DEBUG
	for(Int_t j=0;j<(nBinsInt+1);j++){
		cout<<"bin number \t"<<(j+1)<<"\t has lower edge = \t"<< binLowEdges[j] <<endl;
	}
	for(Int_t j=0;j<nBinsInt;j++){
		cout<<"bin num \t"<<j<<"\t has contents = \t"<< binContents[j]<<endl;
	}
#endif

	///now make a new histo with the number of bins = nBinsInt, and the
	///bin contents equal to the gen matching efficiency
	TCanvas * c1 = new TCanvas(canvName,canvName,650,650);
	c1->cd();
	TH1F * effHist = new TH1F(effHistName,effHistTitle,nBinsInt,binLowEdges);
	effHist->SetMarkerColor(kRed);
	effHist->SetMarkerSize(1.5);
	for(Int_t h=1;h<=nBinsInt;h++){
		effHist->SetBinContent(h, binContents[h-1]);
	}///end loop over bins in rebinnedDenomHist
	effHist->Draw("*H");	///< draw histogram bars to indicate bin sizes (variable), and points at center of each bin
	c1->Update();
	c1->SaveAs(outputFile,"recreate");

}///end calcAndPlotEff()

/**
 * use this macro to calculate the efficiency to match a reco object to a gen obj at a fixed value of deltaR
 * as a fxn of pt, eta, and phi of the gen obj
 *
 */
void genMatchingEfficiency(){
	TChain * matchedGenJetsToGenQuarksNoCuts = new TChain("matchGenJetsToGenQuarksNoCutsNewPath/matchedGenJetsNoCutsTree","");
	matchedGenJetsToGenQuarksNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");
	
	TChain * matchedRecoJetsToGenJetsNoCuts = new TChain("matchRecoJetsToGenJetsNoCutsNewPath/matchedRecoJetsNoCutsTree","");
	matchedRecoJetsToGenJetsNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");
	
	TChain * matchedRecoEleToLeadingGenEleNoCuts = new TChain("matchRecoEleToLeadingGenEleNoCutsNewPath/matchedLeadingRecoEleNoCutsTree","");
	matchedRecoEleToLeadingGenEleNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");

	TChain * matchedRecoEleToSubleadingGenEleNoCuts = new TChain("matchRecoEleToSubleadingGenEleNoCutsNewPath/matchedSubleadingRecoEleNoCutsTree","");
	matchedRecoEleToSubleadingGenEleNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel.root");

	//calcAndPlotEff(Double_t nBins, Int_t nBinsInt, TString cutString, TString uncutHistPlotArgs, TString cutHistPlotArgs, TString uncutHistName, TString cutHistName, TString canvName, TString outputFile, TChain * chain, TString effHistName, TString effHistTitle)
	
	///matching efficiency vs gen pT
	/*
	calcAndPlotEff(9,9,"nWithMatch>0","ptLowerLevel>>ptHistoOne(400,0.,1300.)","ptLowerLevel>>selectedPtHistoOne(400,0.,1300.)","ptHistoOne","selectedPtHistoOne","a1","subleadingRecoToGenEleMatchingEff_vs_pt.png",matchedRecoEleToSubleadingGenEleNoCuts,"matchEffSubleadingEleVsPt","reco to gen subleading electron matching efficiency vs gen pT");
	calcAndPlotEff(9,9,"nWithMatch>0","ptLowerLevel>>ptHistoTwo(400,0.,1300.)","ptLowerLevel>>selectedPtHistoTwo(400,0.,1300.)","ptHistoTwo","selectedPtHistoTwo","a2","leadingRecoToGenEleMatchingEff_vs_pt.png",matchedRecoEleToLeadingGenEleNoCuts,"matchEffLeadingEleVsPt","reco to gen leading electron matching efficiency vs gen pT");
	calcAndPlotEff(9,9,"nWithMatch>1","ptLowerLevel>>ptHistoThree(400,0.,1300.)","ptLowerLevel>>selectedPtHistoThree(400,0.,1300.)","ptHistoThree","selectedPtHistoThree","a3","recoJetsToGenJetsMatchingEff_vs_pt.png",matchedRecoJetsToGenJetsNoCuts,"matchEffRecoJetsVsPt","reco to gen jet matching efficiency vs gen pT");
	calcAndPlotEff(9,9,"nWithMatch>1","ptLowerLevel>>ptHistoFour(400,0.,1300.)","ptLowerLevel>>selectedPtHistoFour(400,0.,1300.)","ptHistoFour","selectedPtHistoFour","a4","genJetsToGenQuarksMatchingEff_vs_pt.png",matchedGenJetsToGenQuarksNoCuts,"matchEffGenJetsVsPt","gen jet to quark matching efficiency vs gen quark pT");

	///matching efficiency vs gen eta
	calcAndPlotEff(8,8,"nWithMatch>0","etaLowerLevel>>etaHistoOne(375,-3.,3.05)","etaLowerLevel>>selectedEtaHistoOne(375,-3.,3.05)","etaHistoOne","selectedEtaHistoOne","b1","subleadingRecoToGenEleMatchingEff_vs_eta.png",matchedRecoEleToSubleadingGenEleNoCuts,"matchEffSubleadingEleVsEta","reco to gen subleading electron matching efficiency vs gen #eta");
	calcAndPlotEff(8,8,"nWithMatch>0","etaLowerLevel>>etaHistoTwo(375,-3.,3.05)","etaLowerLevel>>selectedEtaHistoTwo(375,-3.,3.05)","etaHistoTwo","selectedEtaHistoTwo","b2","leadingRecoToGenEleMatchingEff_vs_eta.png",matchedRecoEleToLeadingGenEleNoCuts,"matchEffLeadingEleVsEta","reco to gen leading electron matching efficiency vs gen #eta");
	calcAndPlotEff(8,8,"nWithMatch>1","etaLowerLevel>>etaHistoThree(375,-3.,3.05)","etaLowerLevel>>selectedEtaHistoThree(375,-3.,3.05)","etaHistoThree","selectedEtaHistoThree","b3","recoJetsToGenJetsMatchingEff_vs_eta.png",matchedRecoJetsToGenJetsNoCuts,"matchEffRecoJetsVsEta","reco to gen jet matching efficiency vs gen #eta");
	calcAndPlotEff(8,8,"nWithMatch>1","etaLowerLevel>>etaHistoFour(375,-3.,3.05)","etaLowerLevel>>selectedEtaHistoFour(375,-3.,3.05)","etaHistoFour","selectedEtaHistoFour","b4","genJetsToGenQuarksMatchingEff_vs_eta.png",matchedGenJetsToGenQuarksNoCuts,"matchEffGenJetsVsEta","gen jet to quark matching efficiency vs gen quark #eta");
	*/

	///matching efficiency vs gen phi 
	calcAndPlotEff(8,8,"nWithMatch>0","phiLowerLevel>>phiHistoOne(375,-3.3,3.35)","phiLowerLevel>>selectedPhiHistoOne(375,-3.3,3.35)","phiHistoOne","selectedPhiHistoOne","d1","subleadingRecoToGenEleMatchingEff_vs_phi.png",matchedRecoEleToSubleadingGenEleNoCuts,"matchEffSubleadingEleVsPhi","reco to gen subleading electron matching efficiency vs gen #phi");
	calcAndPlotEff(8,8,"nWithMatch>0","phiLowerLevel>>phiHistoTwo(375,-3.3,3.35)","phiLowerLevel>>selectedPhiHistoTwo(375,-3.3,3.35)","phiHistoTwo","selectedPhiHistoTwo","d2","leadingRecoToGenEleMatchingEff_vs_phi.png",matchedRecoEleToLeadingGenEleNoCuts,"matchEffLeadingEleVsPhi","reco to gen leading electron matching efficiency vs gen #phi");
	calcAndPlotEff(8,8,"nWithMatch>1","phiLowerLevel>>phiHistoThree(375,-3.3,3.35)","phiLowerLevel>>selectedPhiHistoThree(375,-3.3,3.35)","phiHistoThree","selectedPhiHistoThree","d3","recoJetsToGenJetsMatchingEff_vs_phi.png",matchedRecoJetsToGenJetsNoCuts,"matchEffRecoJetsVsPhi","reco to gen jet matching efficiency vs gen #phi");
	calcAndPlotEff(8,8,"nWithMatch>1","phiLowerLevel>>phiHistoFour(375,-3.3,3.35)","phiLowerLevel>>selectedPhiHistoFour(375,-3.3,3.35)","phiHistoFour","selectedPhiHistoFour","d4","genJetsToGenQuarksMatchingEff_vs_phi.png",matchedGenJetsToGenQuarksNoCuts,"matchEffGenJetsVsPhi","gen jet to quark matching efficiency vs gen quark #phi");
	

}///end genMatchingEfficiency()


