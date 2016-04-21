#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TObject.h"
#include <vector>
#include <iostream>
#include <string>
#include <cstdio>

//#define DEBUG

void quickHistoDivide(){
	TFile * histosOne = new TFile("jorgeHistosFromHisMinitrees.root","READ");
	TFile * histosTwo = new TFile("seanHistosUsingJorgeMinitreesWithRochCorrs.root","READ");
#ifdef DEBUG
	std::cout<<"initialized pointers to input files"<<std::endl;
#endif

	TString powhegHistName = "h_Mll_DY_pow", madgraphHistName = "h_Mll_DY_mad";
	histosOne->cd();
#ifdef DEBUG
	std::cout<<"cd to histosOne input file"<<std::endl;
#endif
	TH1 * powhegHistoOne = (TH1*) histosOne->Get(powhegHistName);
	TH1 * madgraphHistoOne = (TH1*) histosOne->Get(madgraphHistName);
#ifdef DEBUG
	std::cout<<"retrieved two histos from histosOne input file"<<std::endl;
#endif
	
	histosTwo->cd();
#ifdef DEBUG
	std::cout<<"cd to histosTwo input file"<<std::endl;
#endif
	TH1 * powhegHistoTwo = (TH1*) histosTwo->Get(powhegHistName);
	TH1 * madgraphHistoTwo = (TH1*) histosTwo->Get(madgraphHistName);
#ifdef DEBUG
	std::cout<<"retrieved two histos from histosTwo input file"<<std::endl;
#endif
	
	//divide the two powheg histos (and madgraph histos) and plot the division
	TCanvas * c1 = new TCanvas("c1","c1",800,800);
	c1->cd();
	powhegHistoOne->SetTitle("M_MuMu for PowhegM50to120 JorgeMinitrees  JorgeAnalysis divided by SeanAnalysis");
	powhegHistoOne->SetMarkerColor(kBlack);
	powhegHistoOne->SetMarkerStyle(21);
	powhegHistoOne->SetMarkerSize(1);
	powhegHistoOne->Divide(powhegHistoTwo);
#ifdef DEBUG
	std::cout<<"called divide on powheg histos"<<std::endl;
#endif
	powhegHistoOne->Draw("HISTP");
	c1->SaveAs("dimuonMassRatioPowhegM50to120_JorgeMinitrees_JorgeAnalysis_dividedBy_SeanAnalysis.png","RECREATE");
	c1->SaveAs("dimuonMassRatioPowhegM50to120_JorgeMinitrees_JorgeAnalysis_dividedBy_SeanAnalysis.pdf","RECREATE");

	TCanvas * c2 = new TCanvas("c2","c2",800,800);
	c2->cd();
	madgraphHistoOne->SetTitle("M_MuMu for MadIncl JorgeMinitrees  JorgeAnalysis divided by SeanAnalysis");\
	madgraphHistoOne->SetMarkerColor(kBlack);
	madgraphHistoOne->SetMarkerStyle(21);
	madgraphHistoOne->SetMarkerSize(1);
	madgraphHistoOne->Divide(madgraphHistoTwo);
	madgraphHistoOne->Draw("HISTP");
	c2->SaveAs("dimuonMassRatioMadIncl_JorgeMinitrees_JorgeAnalysis_dividedBy_SeanAnalysis.png","RECREATE");
	c2->SaveAs("dimuonMassRatioMadIncl_JorgeMinitrees_JorgeAnalysis_dividedBy_SeanAnalysis.pdf","RECREATE");

}///end quickHistoDivide()

