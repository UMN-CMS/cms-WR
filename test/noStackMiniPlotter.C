#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <string>
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>

/**
 * this macro is derived from quickCalcDyScaleFactors.C, and is a generic plotting macro used to overlay up 
 * to three continuous curves and one set of points on one canvas.  The muon and electron channels are processed 
 * simultaneously.  This macro should be run by itself using minitrees which have been processed by analysis.cpp. 
 * The resulting plots do not have a ratio plot.
 *
 */


#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

TString dir = "../analysisCppOutputRootFilesNewHltSf/";


void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1D*> *hs, Float_t normRescale, bool applyCorr);
void drawPlots(TH1D* hs_DYPowheg, TH1D* hs_DYMadIncl, TH1D* hs_DYAmcIncl, TH1D* hs_data, TString xtitle, TString fname, std::string channel);
void noStackMiniPlotter()
{
	
	TString treeName = "Tree_Iter0";
	TChain * chain_DYPowhegEE = new TChain(treeName,"DYPowhegMassBinnedEE");
	TChain * chain_DYMadInclEE = new TChain(treeName,"DYMadgraphInclusiveEE");
	TChain * chain_DYAmcInclEE = new TChain(treeName,"DYAMCInclusiveEE");
	TChain * chain_dataEE = new TChain(treeName,"DataEE");
	TChain * chain_DYPowhegMuMu = new TChain(treeName,"DYPowhegMassBinnedMuMu");
	TChain * chain_DYMadInclMuMu = new TChain(treeName,"DYMadgraphInclusiveMuMu");
	TChain * chain_DYAmcInclMuMu = new TChain(treeName,"DYAMCInclusiveMuMu");
	TChain * chain_dataMuMu = new TChain(treeName,"DataMuMu");

	chain_DYPowhegEE->Add(dir+"selected_tree_DYPOWHEG_dytagandprobeEE"+mcFileTag+".root");
	chain_DYMadInclEE->Add(dir+"selected_tree_DYMadInclAndHT_dytagandprobeEE"+mcFileTag+".root");
	chain_DYAmcInclEE->Add(dir+"selected_tree_DYAMC_dytagandprobeEE"+mcFileTag+".root");
	chain_dataEE->Add(dir+"selected_tree_data_dytagandprobeEE"+dataFileTag+".root");
	
	chain_DYAmcInclMuMu->Add(dir+"selected_tree_DYAMC_dytagandprobeMuMu"+mcFileTag+".root");
	chain_DYPowhegMuMu->Add(dir+"selected_tree_DYPOWHEG_dytagandprobeMuMu"+mcFileTag+".root");
	chain_DYMadInclMuMu->Add(dir+"selected_tree_DYMadInclAndHT_dytagandprobeMuMu"+mcFileTag+".root");
	
	//temporary change to save disk space
	//chain_DYAmcInclMuMu->Add(dir+"selected_tree_DYMadInclAndHT_dytagandprobeMuMu"+mcFileTag+".root");
	//chain_DYPowhegMuMu->Add(dir+"selected_tree_DYMadInclAndHT_dytagandprobeMuMu"+mcFileTag+".root");
	//chain_DYMadInclMuMu->Add(dir+"selected_tree_DYMadInclAndHT_dytagandprobeMuMu"+mcFileTag+".root");
	chain_dataMuMu->Add(dir+"selected_tree_data_dytagandprobeMuMu"+dataFileTag+".root");

	Selector myEvent_DYPowhegEE;
	Selector myEvent_DYMadInclEE;
	Selector myEvent_DYAmcInclEE;
	Selector myEvent_dataEE;
	Selector myEvent_DYPowhegMuMu;
	Selector myEvent_DYMadInclMuMu;
	Selector myEvent_DYAmcInclMuMu;
	Selector myEvent_dataMuMu;

	myEvent_DYPowhegEE.SetBranchAddresses(chain_DYPowhegEE);
	myEvent_DYMadInclEE.SetBranchAddresses(chain_DYMadInclEE);
	myEvent_DYAmcInclEE.SetBranchAddresses(chain_DYAmcInclEE);
	myEvent_dataEE.SetBranchAddresses(chain_dataEE);
	myEvent_DYPowhegMuMu.SetBranchAddresses(chain_DYPowhegMuMu);
	myEvent_DYMadInclMuMu.SetBranchAddresses(chain_DYMadInclMuMu);
	myEvent_DYAmcInclMuMu.SetBranchAddresses(chain_DYAmcInclMuMu);
	myEvent_dataMuMu.SetBranchAddresses(chain_dataMuMu);


	std::vector<TH1D*> hs_DYPowhegEE;
	MakeHistos(chain_DYPowhegEE, &myEvent_DYPowhegEE, &hs_DYPowhegEE, 1, false);
	std::vector<TH1D*> hs_DYMadInclEE;
	MakeHistos(chain_DYMadInclEE, &myEvent_DYMadInclEE, &hs_DYMadInclEE, 1, false);
	std::vector<TH1D*> hs_DYAmcInclEE;
	MakeHistos(chain_DYAmcInclEE, &myEvent_DYAmcInclEE, &hs_DYAmcInclEE, 1, false);
	std::vector<TH1D*> hs_DYPowhegMuMu;
	MakeHistos(chain_DYPowhegMuMu, &myEvent_DYPowhegMuMu, &hs_DYPowhegMuMu, 1.0, false);
	std::vector<TH1D*> hs_DYMadInclMuMu;
	MakeHistos(chain_DYMadInclMuMu, &myEvent_DYMadInclMuMu, &hs_DYMadInclMuMu, 1.0, false);
	std::vector<TH1D*> hs_DYAmcInclMuMu;
	MakeHistos(chain_DYAmcInclMuMu, &myEvent_DYAmcInclMuMu, &hs_DYAmcInclMuMu, 1.0, false);


	std::vector<TH1D*> hs_dataEE;
	MakeHistos(chain_dataEE, &myEvent_dataEE, &hs_dataEE, 1.0, false);
	std::vector<TH1D*> hs_dataMuMu;
	MakeHistos(chain_dataMuMu, &myEvent_dataMuMu, &hs_dataMuMu, 1.0, false);
	
	//now the vectors of TH1D pointers are filled with pointers to histos with nonzero entries
	unsigned int nPlots = hs_DYPowhegEE.size();

	TString xtitles[] = {"leading lepton p_{T}", "subleading lepton p_{T}", "leading jet p_{T}", "subleading jet p_{T}", "leading lepton #eta", "subleading lepton #eta", "leading jet #eta", "subleading jet #eta", "leading lepton #phi", "subleading lepton #phi", "leading jet #phi", "subleading jet #phi", "dilepton mass", "nPV", "Z #phi", "Z rapidity", "Z p_{T}", "dilepton mass both leptons in barrel", "dilepton mass both leptons in endcap", "dilepton mass one lepton in barrel other lepton in endcap","nPV","Z P_{T} (GeV)","Z P_{T} (GeV)","M_{LLJJ} (GeV)"};

	TString fnames[] = {"l1_pt", "l2_pt", "j1_pt", "j2_pt", "l1_eta", "l2_eta", "j1_eta", "j2_eta", "l1_phi", "l2_phi", "j1_phi", "j2_phi", "Mll", "nPV", "Z_phi", "Z_rapidity", "Z_pt", "Mll_bothEB","Mll_bothEE","Mll_EBEE","nPU","Z_pt_coarseBinning","Z_pt_variableBinning", "Mlljj"};

	int i = 0;
	for(unsigned int i = 0; i < nPlots; i++) {
		std::string s = std::to_string(i);
		drawPlots(hs_DYPowhegEE[i], hs_DYMadInclEE[i], hs_DYAmcInclEE[i], hs_dataEE[i], xtitles[i], fnames[i], "EE");
		drawPlots(hs_DYPowhegMuMu[i], hs_DYMadInclMuMu[i], hs_DYAmcInclMuMu[i], hs_dataMuMu[i], xtitles[i], fnames[i], "MuMu");

	}

}//end noStackMiniPlotter()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1D*> *hs, Float_t normRescale, bool applyCorr)
{
	TH1D *h_lepton_pt0 = new TH1D("h_lepton_pt0", "", 100, -20, 600);
	TH1D *h_lepton_eta0 = new TH1D("h_lepton_eta0", "", 50, -3.2, 3.2);
	TH1D *h_lepton_phi0 = new TH1D("h_lepton_phi0", "", 50, -3.2, 3.2);
	TH1D *h_lepton_pt1 = new TH1D("h_lepton_pt1", "", 120, -15, 400);
	TH1D *h_lepton_eta1 = new TH1D("h_lepton_eta1", "", 50, -3, 3);
	TH1D *h_lepton_phi1 = new TH1D("h_lepton_phi1", "", 50, -3.15, 3.15);

	TH1D *h_jet_pt0 = new TH1D("h_jet_pt0", "", 50, -5, 300);
	TH1D *h_jet_eta0 = new TH1D("h_jet_eta0", "", 50, -3.2, 3.2);
	TH1D *h_jet_phi0 = new TH1D("h_jet_phi0", "", 50, -3.2, 3.2);
	TH1D *h_jet_pt1 = new TH1D("h_jet_pt1", "", 50, -5, 300);
	TH1D *h_jet_eta1 = new TH1D("h_jet_eta1", "", 50, -3.2, 3.2);
	TH1D *h_jet_phi1 = new TH1D("h_jet_phi1", "", 50, -3.2, 3.2);

	TH1D *h_dilepton_mass = new TH1D("h_dilepton_mass", "", 80, 70., 110.);	//coarser binning
	//TH1D *h_dilepton_mass = new TH1D("h_dilepton_mass", "", 280, 70., 110.);	//default
	//TH1D *h_dilepton_mass = new TH1D("h_dilepton_mass", "", 560, 70., 110.);	//finer binning
	TH1D *h_nPV = new TH1D("pileup", "nPV distribution", 50., 0, 50.);
	TH1D *h_nPU = new TH1D("pileup", "true nPU distribution for MC", 50., 0, 50.);

	TH1D *h_Z_phi = new TH1D("h_Z_phi", "", 50, -3.2, 3.2);
	TH1D *h_Z_rapidity = new TH1D("h_Z_rapidity", "", 70, -5., 8.);
	TH1D *h_Z_pt = new TH1D("h_Z_pt", "", 35, -10., 300.);
	TH1D *h_Z_ptCoarseBinning = new TH1D("h_Z_ptCoarseBinning", "", 10, -10., 1000.);

	//Z pT with variable bin widths
	Float_t bins[] = {0,50,100,150,200,300,400,500,600,2000};
	Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
	TH1D *h_Z_ptVariableBinning = new TH1D("h_Z_ptVariableBinning", "", binnum, bins);

	TH1D *h_dilepton_massBothEB = new TH1D("h_dilepton_massBothEB", "", 280, 55., 125.);
	TH1D *h_dilepton_massBothEE = new TH1D("h_dilepton_massBothEE", "", 280, 55., 125.);
	TH1D *h_dilepton_massEBEE = new TH1D("h_dilepton_massEBEE", "", 280, 55., 125.);
	
	Float_t binsWR[] = { 100, 200, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1900, 2500};	//wider binsWR work better at high WR mass, include overflow evts in last bin shown on plot
	Int_t  binnumWR = sizeof(binsWR)/sizeof(Float_t) - 1;
	TH1D * h_WR_mass = new TH1D("h_WR_mass","",binnumWR, binsWR);
	
	Long64_t nEntries = chain->GetEntries();

	cout << nEntries << endl;

	TString chTitle(chain->GetTitle());
	//std::cout<<"chTitle=\t"<< chTitle << std::endl;
	Double_t hltCorrFactor=1.0;
	
	for(int ev = 0; ev < nEntries; ++ev) {
		chain->GetEntry(ev);
		if(myEvent->dilepton_mass > upperMllCut || myEvent->dilepton_mass < lowerMllCut) continue;
		if(myEvent->lead_lepton_pt < leadLeptonPtCut || myEvent->sublead_lepton_pt < subleadLeptonPtCut || std::fabs(myEvent->sublead_lepton_eta) > leptonEtaCut || std::fabs(myEvent->lead_lepton_eta) > leptonEtaCut) continue;
		if(myEvent->lead_jet_pt < leadJetPtCut || myEvent->sublead_jet_pt < subLeadJetPtCut) continue;

		TLorentzVector leadLeptonFourMom, subleadLeptonFourMom, zFourMom;
		leadLeptonFourMom.SetPtEtaPhiE(myEvent->lead_lepton_pt, myEvent->lead_lepton_eta, myEvent->lead_lepton_phi, myEvent->lead_lepton_pt);
		subleadLeptonFourMom.SetPtEtaPhiE(myEvent->sublead_lepton_pt, myEvent->sublead_lepton_eta, myEvent->sublead_lepton_phi, myEvent->sublead_lepton_pt);
		zFourMom = leadLeptonFourMom + subleadLeptonFourMom;

		h_Z_pt->Fill(zFourMom.Pt(), (myEvent->weight)*hltCorrFactor);
		h_Z_phi->Fill(zFourMom.Phi(), (myEvent->weight)*hltCorrFactor);
		h_Z_rapidity->Fill(zFourMom.Rapidity(), (myEvent->weight)*hltCorrFactor);
		h_Z_ptCoarseBinning->Fill(zFourMom.Pt(), (myEvent->weight)*hltCorrFactor);
		h_Z_ptVariableBinning->Fill(zFourMom.Pt(), (myEvent->weight)*hltCorrFactor);

		h_lepton_pt0->Fill(myEvent->lead_lepton_pt, (myEvent->weight)*hltCorrFactor);
		h_lepton_pt1->Fill(myEvent->sublead_lepton_pt, (myEvent->weight)*hltCorrFactor);
		h_lepton_eta0->Fill(myEvent->lead_lepton_eta, (myEvent->weight)*hltCorrFactor);
		h_lepton_eta1->Fill(myEvent->sublead_lepton_eta, (myEvent->weight)*hltCorrFactor);
		h_lepton_phi0->Fill(myEvent->lead_lepton_phi, (myEvent->weight)*hltCorrFactor);
		h_lepton_phi1->Fill(myEvent->sublead_lepton_phi, (myEvent->weight)*hltCorrFactor);

		h_jet_pt0->Fill(myEvent->lead_jet_pt, (myEvent->weight)*hltCorrFactor);
		h_jet_pt1->Fill(myEvent->sublead_jet_pt, (myEvent->weight)*hltCorrFactor);
		h_jet_eta0->Fill(myEvent->lead_jet_eta, (myEvent->weight)*hltCorrFactor);
		h_jet_eta1->Fill(myEvent->sublead_jet_eta, (myEvent->weight)*hltCorrFactor);
		h_jet_phi0->Fill(myEvent->lead_jet_phi, (myEvent->weight)*hltCorrFactor);
		h_jet_phi1->Fill(myEvent->sublead_jet_phi, (myEvent->weight)*hltCorrFactor);

		h_dilepton_mass->Fill(myEvent->dilepton_mass, (myEvent->weight)*hltCorrFactor);
		h_nPV->Fill(myEvent->nPV, (myEvent->weight)*hltCorrFactor);
		h_nPU->Fill(myEvent->nPU, (myEvent->weight)*hltCorrFactor);

		//these three ifs are mutually exclusive
		if(std::fabs(myEvent->lead_lepton_eta) < 1.44 && std::fabs(myEvent->sublead_lepton_eta) < 1.44) h_dilepton_massBothEB->Fill(myEvent->dilepton_mass, (myEvent->weight)*hltCorrFactor);
		if(std::fabs(myEvent->lead_lepton_eta) > 1.56 && std::fabs(myEvent->sublead_lepton_eta) > 1.56) h_dilepton_massBothEE->Fill(myEvent->dilepton_mass, (myEvent->weight)*hltCorrFactor);
		if( (std::fabs(myEvent->lead_lepton_eta) > 1.56 && std::fabs(myEvent->sublead_lepton_eta) < 1.44) || (std::fabs(myEvent->lead_lepton_eta) < 1.44 && std::fabs(myEvent->sublead_lepton_eta) > 1.56) ) h_dilepton_massEBEE->Fill(myEvent->dilepton_mass, (myEvent->weight)*hltCorrFactor);

		TLorentzVector leadJetFourMom, subleadJetFourMom, wrCandFourMom;
		leadJetFourMom.SetPtEtaPhiE(myEvent->lead_jet_pt, myEvent->lead_jet_eta, myEvent->lead_jet_phi, myEvent->lead_jet_pt);
		subleadJetFourMom.SetPtEtaPhiE(myEvent->sublead_jet_pt, myEvent->sublead_jet_eta, myEvent->sublead_jet_phi, myEvent->sublead_jet_pt);
		wrCandFourMom = leadJetFourMom + subleadJetFourMom + leadLeptonFourMom + subleadLeptonFourMom;
		h_WR_mass->Fill(wrCandFourMom.M(), (myEvent->weight)*hltCorrFactor);
	
	}

	h_Z_ptCoarseBinning->Scale(normRescale);
	h_Z_ptVariableBinning->Scale(normRescale);
	h_Z_pt->Scale(normRescale);
	h_Z_phi->Scale(normRescale);
	h_Z_rapidity->Scale(normRescale);

	h_lepton_pt0->Scale(normRescale);
	h_lepton_pt1->Scale(normRescale);
	h_lepton_eta0->Scale(normRescale);
	h_lepton_eta1->Scale(normRescale);
	h_lepton_phi0->Scale(normRescale);
	h_lepton_phi1->Scale(normRescale);

	h_jet_pt0->Scale(normRescale);
	h_jet_pt1->Scale(normRescale);
	h_jet_eta0->Scale(normRescale);
	h_jet_eta1->Scale(normRescale);
	h_jet_phi0->Scale(normRescale);
	h_jet_phi1->Scale(normRescale);

	h_dilepton_mass->Scale(normRescale);
	h_nPV->Scale(normRescale);
	h_nPU->Scale(normRescale);

	h_dilepton_massBothEB->Scale(normRescale);
	h_dilepton_massBothEE->Scale(normRescale);
	h_dilepton_massEBEE->Scale(normRescale);

	h_WR_mass->Scale(normRescale);

	///this order of push_back calls should not be changed
	hs->push_back(h_lepton_pt0);
	hs->push_back(h_lepton_pt1);
	hs->push_back(h_jet_pt0);
	hs->push_back(h_jet_pt1);
	hs->push_back(h_lepton_eta0);
	hs->push_back(h_lepton_eta1);
	hs->push_back(h_jet_eta0);
	hs->push_back(h_jet_eta1);
	hs->push_back(h_lepton_phi0);
	hs->push_back(h_lepton_phi1);
	hs->push_back(h_jet_phi0);
	hs->push_back(h_jet_phi1);
	hs->push_back(h_dilepton_mass);
	hs->push_back(h_nPV);
	hs->push_back(h_Z_phi);
	hs->push_back(h_Z_rapidity);
	hs->push_back(h_Z_pt);
	hs->push_back(h_dilepton_massBothEB);
	hs->push_back(h_dilepton_massBothEE);
	hs->push_back(h_dilepton_massEBEE);
	hs->push_back(h_nPU);
	hs->push_back(h_Z_ptCoarseBinning);
	
	//divide bin contents of h_Z_ptVariableBinning by width of each bin
	//get num bins in histo
	Int_t nBins = h_Z_ptVariableBinning->GetNbinsX();
	for(Int_t j=1; j<=nBins; j++){
		//include the overflows in the very last bin shown on the plot
		if(j==nBins){
			Double_t origBinContents = h_Z_ptVariableBinning->GetBinContent(j);
			Double_t overflowContents = h_Z_ptVariableBinning->GetBinContent(j+1);
			std::cout<<"overflow contents =\t"<< overflowContents << std::endl;	//sanity check
			h_Z_ptVariableBinning->SetBinContent(j, origBinContents+overflowContents);
		}//end work to include overflows in last bin shown on plot
		//in each bin, divide the bin contents by the bin width
		Double_t oldBinContents = h_Z_ptVariableBinning->GetBinContent(j);
		Double_t oldBinErrors = h_Z_ptVariableBinning->GetBinError(j);
		Double_t binWidth = h_Z_ptVariableBinning->GetBinWidth(j);
		h_Z_ptVariableBinning->SetBinContent(j, oldBinContents/binWidth);
		h_Z_ptVariableBinning->SetBinError(j, oldBinErrors/binWidth);
	}//end loop over bins in histo

	//add overflow events from h_WR_mass into last bin shown on plot, and divide bin contents and bin errors by bin widths
	Int_t nBinsWR = h_WR_mass->GetNbinsX();
	for(Int_t j=1; j<=nBinsWR; j++){
		//include the overflows in the very last bin shown on the plot
		if(j==nBinsWR){
			Double_t origBinContents = h_WR_mass->GetBinContent(j);
			Double_t overflowContents = h_WR_mass->GetBinContent(j+1);
			//std::cout<<"overflow contents =\t"<< overflowContents << std::endl;	//sanity check
			h_WR_mass->SetBinContent(j, origBinContents+overflowContents);
		}//end including overflows in last bin shown on plot

		//in each bin, divide the bin contents by the bin width
		Double_t oldBinContents = h_WR_mass->GetBinContent(j);
		Double_t oldBinErrors = h_WR_mass->GetBinError(j);
		Double_t binWidth = h_WR_mass->GetBinWidth(j);
		h_WR_mass->SetBinContent(j, oldBinContents/binWidth);
		h_WR_mass->SetBinError(j, oldBinErrors/binWidth);
	}//end loop over bins in h_WR_mass

	hs->push_back(h_Z_ptVariableBinning);
	hs->push_back(h_WR_mass);
	
}

void drawPlots(TH1D* hs_DYPowheg, TH1D* hs_DYMadIncl, TH1D* hs_DYAmcIncl, TH1D* hs_data, TString xtitle, TString fname, std::string channel)
{


	//gStyle->SetOptStat("eou");
	gStyle->SetOptStat("");
	TLegend *leg = new TLegend( 0.60, 0.60, 0.90, 0.90 ) ;
	leg->AddEntry( hs_DYPowheg, "DY Powheg" ) ;
	leg->AddEntry( hs_DYMadIncl, "ABC" ) ;
	leg->AddEntry( hs_DYAmcIncl, "DEF" ) ;
	leg->AddEntry( hs_data, "Data");
	leg->SetFillColor( kWhite ) ;

	TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
	mycanvas->cd();
	hs_DYPowheg->SetLineColor(kRed);
	hs_DYPowheg->SetLineWidth(3);
	hs_DYMadIncl->SetLineColor(kBlack);
	hs_DYMadIncl->SetLineWidth(3);
	hs_DYAmcIncl->SetLineColor(kBlue);
	hs_DYAmcIncl->SetLineWidth(3);
	hs_data->SetMarkerStyle(20);
	hs_data->SetMarkerSize(1);
	hs_data->SetMarkerColor(kBlack);

	/*for ratio plot*/
	Double_t eps = 0.001;
	TPad* p1 = new TPad("p1", "p1", 0, 0.25, 1, 1, 0);
	p1->Draw();
	TPad* p2 = new TPad("p2", "p2", 0, 0.1, 1, 0.25 + eps, 0);
	p2->Draw();
	p1->SetBottomMargin(0);
	p2->SetTopMargin(0);
	p1->cd();
	/**/
	hs_data->SetStats(1);
	hs_DYPowheg->SetStats(1);
	TH1D *ratio_Powheg = (TH1D*)hs_data->Clone();
	TH1D *ratio_Mad = (TH1D*)hs_data->Clone();
	TH1D *ratio_Amc = (TH1D*)hs_data->Clone();
	//TString plotTitle = "CMS Private   #surds = 13 TeV #int lumi = 2.6 fb^{-1}";
	TString plotTitle = "CMS Private          2.6 fb^{-1} (13 TeV)";
	hs_DYPowheg->SetTitle(plotTitle);
	hs_data->SetTitle(plotTitle);
	TString ytitle = "Events";
	//TString ytitle = "Events/(";
	//ytitle += (hs_data->GetXaxis()->GetNbins());
	//ytitle += (hs_data->GetXaxis()->GetBinWidth(2));
	//ytitle += " GeV)";
	hs_DYAmcIncl->GetYaxis()->SetTitle(ytitle.Data());
	hs_DYPowheg->GetYaxis()->SetTitle(ytitle.Data());
	hs_DYMadIncl->GetYaxis()->SetTitle(ytitle.Data());
	hs_data->GetYaxis()->SetTitle(ytitle.Data());
	hs_data->GetYaxis()->SetTitleOffset(1.5);
	if(fname.EqualTo("Mlljj")) hs_DYMadIncl->GetYaxis()->SetTitle("Events/GeV"), hs_data->GetYaxis()->SetTitle("Events/GeV");
	hs_DYAmcIncl->GetYaxis()->SetTitleOffset(1.5);
	hs_DYMadIncl->GetYaxis()->SetTitleOffset(1.5);
	hs_DYPowheg->GetYaxis()->SetTitleOffset(1.5);
	//hs_data->GetYaxis()->SetTitleSize(0.042);
	//hs_DYAmcIncl->GetYaxis()->SetTitleSize(0.042);
	
	hs_data->GetXaxis()->SetTitleSize(0.04);
	//hs_data->SetLabelSize(0.04, "x");
	//hs_data->SetLabelSize(0.025, "y");
	hs_data->GetXaxis()->SetTitle(xtitle.Data());
	hs_DYAmcIncl->GetXaxis()->SetTitleSize(0.04);
	hs_DYAmcIncl->SetLabelSize(0.04, "x");
	hs_DYAmcIncl->GetXaxis()->SetTitle(xtitle.Data());
	hs_DYPowheg->GetXaxis()->SetTitleSize(0.04);
	hs_DYPowheg->SetLabelSize(0.04, "x");
	hs_DYPowheg->SetLabelSize(0.025, "y");
	hs_DYPowheg->GetXaxis()->SetTitle(xtitle.Data());
	hs_DYMadIncl->GetXaxis()->SetTitleSize(0.04);
	hs_DYMadIncl->SetLabelSize(0.04, "x");
	hs_DYMadIncl->GetXaxis()->SetTitle(xtitle.Data());
	if(fname.EqualTo("Z_pt") ) hs_data->GetXaxis()->SetTitle("Z P_{T} [GeV]"), hs_DYAmcIncl->GetXaxis()->SetTitle("Z P_{T} [GeV]");
	
	hs_data->Draw("ep");
	hs_DYMadIncl->Draw("histo same");
#ifdef showAMC	
	hs_DYAmcIncl->Draw("histo same");
#endif
	hs_data->Draw("epsame");
	hs_DYAmcIncl->GetXaxis()->SetTitle(xtitle.Data());
	ratio_Amc->GetXaxis()->SetTitle(xtitle.Data());
	if(fname.EqualTo("Mll")) ratio_Amc->GetXaxis()->SetTitle("M_{LL} [GeV]"), ratio_Mad->GetXaxis()->SetTitle("M_{LL} [GeV]");
	if(fname.EqualTo("Z_pt")) ratio_Amc->GetXaxis()->SetTitle("Z P_{T} [GeV]"), ratio_Mad->GetXaxis()->SetTitle("Z P_{T} [GeV]");
	if(fname.EqualTo("nPV")) ratio_Amc->GetXaxis()->SetTitle("nPV");
	if(fname.EqualTo("nPU")) ratio_Amc->GetXaxis()->SetTitle("nPU");
	if(fname.EqualTo("Mlljj")) ratio_Amc->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), ratio_Mad->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
	
	Float_t labelSize = 0.25;
	ratio_Amc->GetXaxis()->SetTickSize(0.40);
	ratio_Amc->GetXaxis()->SetTitleSize(labelSize);
	ratio_Amc->SetLabelSize(labelSize - 0.07, "x");
	
	ratio_Mad->GetXaxis()->SetTickSize(0.40);
	ratio_Mad->GetXaxis()->SetTitleSize(labelSize);
	ratio_Mad->SetLabelSize(labelSize - 0.07, "x");
	
	ratio_Powheg->GetXaxis()->SetTitle(xtitle.Data());
	ratio_Powheg->GetXaxis()->SetTickSize(0.40);
	ratio_Powheg->GetXaxis()->SetTitleSize(labelSize);
	ratio_Powheg->SetLabelSize(labelSize - 0.07, "x");
	leg->Draw();
	//mycanvas->cd();
	p2->cd();	//for ratio plot
	ratio_Powheg->Sumw2();
	ratio_Powheg->SetStats(0);
	ratio_Mad->Sumw2();
	ratio_Mad->SetStats(0);
	ratio_Amc->Sumw2();
	ratio_Amc->SetStats(0);

	ratio_Powheg->Divide(hs_DYPowheg);
	ratio_Powheg->SetMarkerStyle(20);
	ratio_Powheg->SetMarkerColor(kRed);
	ratio_Powheg->SetLabelSize(labelSize - 0.07, "y");
	ratio_Powheg->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Powheg->GetYaxis()->SetNdivisions(505);

	ratio_Mad->Divide(hs_DYMadIncl);
	ratio_Mad->SetMarkerStyle(21);
	ratio_Mad->SetMarkerColor(kBlack);
	ratio_Mad->SetLabelSize(labelSize - 0.07, "y");
	ratio_Mad->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Mad->GetYaxis()->SetNdivisions(505);

	ratio_Amc->Divide(hs_DYAmcIncl);
	ratio_Amc->SetMarkerStyle(22);
	ratio_Amc->SetMarkerColor(kBlue);
	ratio_Amc->SetLabelSize(labelSize - 0.07, "y");
	ratio_Amc->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Amc->GetYaxis()->SetNdivisions(505);

	ratio_Mad->Draw("p");
#ifdef showAMC	
	ratio_Amc->Draw("p");
#endif
	//ratio_Powheg->Draw("p");	//comment if only plotting AMC
	float xmax = ratio_Amc->GetXaxis()->GetXmax();
	float xmin = ratio_Amc->GetXaxis()->GetXmin();
	TF1 *f1 = new TF1("f1", "1", xmin, xmax);
	//ratio_Powheg->Draw("p");
	ratio_Mad->Draw("psame");

#ifdef showAMC	
	ratio_Amc->Draw("psame");
#endif
	f1->Draw("same");
	mycanvas->cd();
	/**/

	//no Ratio
	//TString cuts = "_minLeadLeptPt_" + to_string(minLeadLeptonPt) +"_minSubleadLeptPt_" + to_string(minSubleadLeptonPt) + "_minLeadJetPt_" + to_string(leadJetPtCut) + "_minSubleadJetPt_" + to_string(subLeadJetPtCut)+"_noRatioPlot_" + channel;

	//only MadHT ratio
	TString cuts = "_minLeadLeptPt_" + to_string(minLeadLeptonPt) +"_minSubleadLeptPt_" + to_string(minSubleadLeptonPt) + "_minLeadJetPt_" + to_string(leadJetPtCut) + "_minSubleadJetPt_" + to_string(subLeadJetPtCut)+"_onlyMadRatioPlot_" + channel;

#ifdef showAMC
	//only MadHT and AMC ratios
	cuts = "_minLeadLeptPt_" + to_string(minLeadLeptonPt) +"_minSubleadLeptPt_" + to_string(minSubleadLeptonPt) + "_minLeadJetPt_" + to_string(leadJetPtCut) + "_minSubleadJetPt_" + to_string(subLeadJetPtCut) + "_onlyAMCandMadHtRatioPlots_" + channel;
#endif

	if(requireSeparatedLeptonsAndJets) cuts += "_withLeptonJetDrCuts";
	if(!requireSeparatedLeptonsAndJets) cuts += "_withoutLeptonJetDrCuts";
	if(useMllReweighted) cuts += "_mcIsMllReweighted";


	TString fn = fname + cuts;

	/**/
	if(fname.EqualTo("Mll") == true){
		//if(fname.EqualTo("nPV") == true || fname.EqualTo("nPU") == true) fn = fname + "_DYTagAndProbeNoJetCuts_noRatio_" + channel;
		/*
		mycanvas->Print((fn + "_highestZoomRatio.pdf").Data());
		mycanvas->Print((fn + "_highestZoomRatio.png").Data());

		///do the zoomed plots only for Z_pt and M_ll distributions

		if(fname.EqualTo("Mll") == true){
			//reset the Y axis scale on the ratio plot
			ratio_Powheg->GetYaxis()->SetRangeUser(0.9, 1.1);
			ratio_Mad->GetYaxis()->SetRangeUser(0.9, 1.1);
			ratio_Amc->GetYaxis()->SetRangeUser(0.9, 1.1);

			mycanvas->Update();
			mycanvas->Print((fn + "_highZoomRatio.pdf").Data());
			mycanvas->Print((fn + "_highZoomRatio.png").Data());

			//reset the Y axis scale on the ratio plot
			ratio_Powheg->GetYaxis()->SetRangeUser(0.85, 1.15);
			ratio_Mad->GetYaxis()->SetRangeUser(0.85, 1.15);
			ratio_Amc->GetYaxis()->SetRangeUser(0.85, 1.15);

			mycanvas->Update();
			mycanvas->Print((fn + "_mediumZoomRatio.pdf").Data());
			mycanvas->Print((fn + "_mediumZoomRatio.png").Data());

			//reset the Y axis scale on the ratio plot
			ratio_Powheg->GetYaxis()->SetRangeUser(0.8, 1.2);
			ratio_Mad->GetYaxis()->SetRangeUser(0.8, 1.2);
			ratio_Amc->GetYaxis()->SetRangeUser(0.8, 1.2);

			mycanvas->Update();
			mycanvas->Print((fn + "_lowZoomRatio.pdf").Data());
			mycanvas->Print((fn + "_lowZoomRatio.png").Data());

			//reset the Y axis scale on the ratio plot
			ratio_Powheg->GetYaxis()->SetRangeUser(0.7, 1.3);
			ratio_Mad->GetYaxis()->SetRangeUser(0.7, 1.3);
			ratio_Amc->GetYaxis()->SetRangeUser(0.7, 1.3);

			mycanvas->Update();
			mycanvas->Print((fn + "_lowestZoomRatio.pdf").Data());
			mycanvas->Print((fn + "_lowestZoomRatio.png").Data());
		}///end zoom ratio plots for Mll and Zpt distributions
		*/

		//reset the Y axis scale on the ratio plot
		ratio_Powheg->GetYaxis()->SetRangeUser(0.61, 1.39);
		ratio_Mad->GetYaxis()->SetRangeUser(0.61, 1.39);
		ratio_Amc->GetYaxis()->SetRangeUser(0.61, 1.39);

		if(fname.EqualTo("Z_pt_variableBinning") == true) fn += "_overflowsAddedToLastBin";
		mycanvas->Update();
		mycanvas->Print((fn + "_noZoomRatio.pdf").Data());
		mycanvas->Print((fn + "_noZoomRatio.png").Data());
		mycanvas->Print((fn + "_noZoomRatio.C").Data());

		p1->SetLogy();	//for plotting ratio plot
		//mycanvas->SetLogy();	//when not plotting ratio plot
		mycanvas->Print((fn + "_log.pdf").Data());
		mycanvas->Print((fn + "_log.png").Data());
		mycanvas->Print((fn + "_log.C").Data());

	}
	/**/
	
	/*
	mycanvas->Update();

	if(fname.EqualTo("Z_pt") == true || fname.EqualTo("nPV") == true){
		mycanvas->Print((fn + "_noZoomRatio.pdf").Data());
		mycanvas->Print((fn + "_noZoomRatio.png").Data());

		mycanvas->SetLogy();	//when not plotting ratio plot
		//p1->SetLogy();	//for plotting ratio plot
		mycanvas->Print((fn + "_log.pdf").Data());
		mycanvas->Print((fn + "_log.png").Data());
	}
	*/
	
	mycanvas->Close();

}
