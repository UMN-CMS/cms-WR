#include "TH1F.h"
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
// #include "ExoAnalysis/cmsWR/src/Selector.cc"
// #include "ExoAnalysis/cmsWR/src/miniTreeEvent.cc"
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif
//#define DBG
TString dir = "../analysisCppOutputRootFiles/", mcFileTag = "_withMllWeight";

void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t normRescale);
void drawPlots(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, TString xtitle, TString fname, std::string channel);
void quickSignalCompareDYMC()
{

	TString treeName = "Tree_Iter0";
	//TChain * chain_DYPowhegInclEE = new TChain(treeName,"DYPowhegInclusiveEE");
	TChain * chain_DYPowhegEE = new TChain(treeName,"DYPowhegM50toInfEE");
	TChain * chain_DYMadInclEE = new TChain(treeName,"DYMadgraphHTEE");
	TChain * chain_DYAmcInclEE = new TChain(treeName,"DYAMCInclusiveEE");
	TChain * chain_dataEE = new TChain(treeName,"DataEE");
	TChain * chain_DYPowhegMuMu = new TChain(treeName,"DYPowhegM50toInfMuMu");
	TChain * chain_DYMadInclMuMu = new TChain(treeName,"DYMadgraphHTMuMu");
	TChain * chain_DYAmcInclMuMu = new TChain(treeName,"DYAMCInclusiveMuMu");
	TChain * chain_dataMuMu = new TChain(treeName,"DataMuMu");

	//chain_DYPowhegInclEE->Add(dir+"selected_tree_DYPOWINCL_dytagandprobeEE"+mcFileTag+".root");
	chain_DYPowhegEE->Add(dir+"selected_tree_DYPOWHEG_signal_eeEE"+mcFileTag+".root");
	chain_DYMadInclEE->Add(dir+"selected_tree_DYMADHT_signal_eeEE"+mcFileTag+".root");
	chain_DYAmcInclEE->Add(dir+"selected_tree_DYAMC_signal_eeEE"+mcFileTag+".root");
	chain_dataEE->Add(dir+"selected_tree_DYAMC_signal_eeEE"+mcFileTag+".root");
	
	chain_DYPowhegMuMu->Add(dir+"selected_tree_DYPOWHEG_signal_mumuMuMu"+mcFileTag+".root");
	chain_DYMadInclMuMu->Add(dir+"selected_tree_DYMADHT_signal_mumuMuMu"+mcFileTag+".root");
	chain_DYAmcInclMuMu->Add(dir+"selected_tree_DYAMC_signal_mumuMuMu"+mcFileTag+".root");
	chain_dataMuMu->Add(dir+"selected_tree_DYAMC_signal_mumuMuMu"+mcFileTag+".root");


	//Selector myEvent_DYPowhegInclEE;
	Selector myEvent_DYPowhegEE;
	Selector myEvent_DYMadInclEE;
	Selector myEvent_DYAmcInclEE;
	Selector myEvent_dataEE;
	Selector myEvent_DYPowhegMuMu;
	Selector myEvent_DYMadInclMuMu;
	Selector myEvent_DYAmcInclMuMu;
	Selector myEvent_dataMuMu;

	//myEvent_DYPowhegInclEE.SetBranchAddresses(chain_DYPowhegInclEE);
	myEvent_DYPowhegEE.SetBranchAddresses(chain_DYPowhegEE);
	myEvent_DYMadInclEE.SetBranchAddresses(chain_DYMadInclEE);
	myEvent_DYAmcInclEE.SetBranchAddresses(chain_DYAmcInclEE);
	myEvent_dataEE.SetBranchAddresses(chain_dataEE);
	myEvent_DYPowhegMuMu.SetBranchAddresses(chain_DYPowhegMuMu);
	myEvent_DYMadInclMuMu.SetBranchAddresses(chain_DYMadInclMuMu);
	myEvent_DYAmcInclMuMu.SetBranchAddresses(chain_DYAmcInclMuMu);
	myEvent_dataMuMu.SetBranchAddresses(chain_dataMuMu);


	std::vector<TH1F*> hs_DYPowhegEE;
	MakeHistos(chain_DYPowhegEE, &myEvent_DYPowhegEE, &hs_DYPowhegEE,1.0);
	std::vector<TH1F*> hs_DYMadInclEE;
	MakeHistos(chain_DYMadInclEE, &myEvent_DYMadInclEE, &hs_DYMadInclEE, 1.0);
	std::vector<TH1F*> hs_DYAmcInclEE;
	MakeHistos(chain_DYAmcInclEE, &myEvent_DYAmcInclEE, &hs_DYAmcInclEE, 1.0);
	std::vector<TH1F*> hs_DYPowhegMuMu;
	MakeHistos(chain_DYPowhegMuMu, &myEvent_DYPowhegMuMu, &hs_DYPowhegMuMu, 1.0);
	std::vector<TH1F*> hs_DYMadInclMuMu;
	MakeHistos(chain_DYMadInclMuMu, &myEvent_DYMadInclMuMu, &hs_DYMadInclMuMu, 1.0);
	std::vector<TH1F*> hs_DYAmcInclMuMu;
	MakeHistos(chain_DYAmcInclMuMu, &myEvent_DYAmcInclMuMu, &hs_DYAmcInclMuMu, 1.0);

	std::vector<TH1F*> hs_dataEE;
	MakeHistos(chain_dataEE, &myEvent_dataEE, &hs_dataEE, 1.0);
	std::vector<TH1F*> hs_dataMuMu;
	MakeHistos(chain_dataMuMu, &myEvent_dataMuMu, &hs_dataMuMu, 1.0);


	//now the vectors of TH1F pointers are filled with pointers to histos with nonzero entries
	unsigned int nPlots = hs_DYPowhegEE.size();

	TString xtitles[] = {"leading lepton p_{T} [GeV]", "subleading lepton p_{T} [GeV]", "leading jet p_{T} [GeV]", "subleading jet p_{T} [GeV]", "leading lepton #eta", "subleading lepton #eta", "leading jet #eta", "subleading jet #eta", "leading lepton #phi", "subleading lepton #phi", "leading jet #phi", "subleading jet #phi", "dilepton mass [GeV]", "Z p_{T} [GeV]","M_{lljj} [GeV]"};

	TString fnames[] = {"l1_pt", "l2_pt", "j1_pt", "j2_pt", "l1_eta", "l2_eta", "j1_eta", "j2_eta", "l1_phi", "l2_phi", "j1_phi", "j2_phi", "Mll", "Z_pt","Mlljj"};

	int i = 0;
	for(unsigned int i = 0; i < nPlots; i++) {
		std::string s = std::to_string(i);
		drawPlots(hs_DYPowhegEE[i], hs_DYMadInclEE[i], hs_DYAmcInclEE[i], hs_dataEE[i], xtitles[i], fnames[i], "EE");
		drawPlots(hs_DYPowhegMuMu[i], hs_DYMadInclMuMu[i], hs_DYAmcInclMuMu[i], hs_dataMuMu[i], xtitles[i], fnames[i], "MuMu");

	}

}//end quickSignalCompareDYMC()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t normRescale)
{
	TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0", "", 20, 0, 1000);
	TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0", "", 10, -3.2, 3.2);
	TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0", "", 10, -3.2, 3.2);
	TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1", "", 10, 0, 600);
	TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1", "", 15, -3.2, 3.2);
	TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1", "", 15, -3.2, 3.2);

	TH1F *h_jet_pt0 = new TH1F("h_jet_pt0", "", 20, 0, 1000);
	TH1F *h_jet_eta0 = new TH1F("h_jet_eta0", "", 10, -3.2, 3.2);
	TH1F *h_jet_phi0 = new TH1F("h_jet_phi0", "", 50, -3.2, 3.2);
	TH1F *h_jet_pt1 = new TH1F("h_jet_pt1", "", 10, 0, 600);
	TH1F *h_jet_eta1 = new TH1F("h_jet_eta1", "", 15, -3.2, 3.2);
	TH1F *h_jet_phi1 = new TH1F("h_jet_phi1", "", 15, -3.2, 3.2);

	TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass", "", 20, 100, 1500);
	TH1F *h_nPV = new TH1F("h_nPV", "", 100, 0, 100);

	TH1F *h_Z_phi = new TH1F("h_Z_phi", "", 15, -3.2, 3.2);
	TH1F *h_Z_rapidity = new TH1F("h_Z_rapidity", "", 15, -5., 8.);
	TH1F *h_Z_pt = new TH1F("h_Z_pt", "", 10, 0., 800.);

  	TH1F *h_WR_mass = new TH1F("h_WR_mass","",20,400,3000);	//fixed bin width

	Long64_t nEntries = chain->GetEntries();
	cout << nEntries << endl;
	
	for(int ev = 0; ev < nEntries; ++ev) {
		chain->GetEntry(ev);
		if(myEvent->WR_mass < 600) continue;
		if(myEvent->dilepton_mass < 200) continue;

		TLorentzVector leadLeptonFourMom, subleadLeptonFourMom, zFourMom;
		leadLeptonFourMom.SetPtEtaPhiE(myEvent->lead_lepton_pt, myEvent->lead_lepton_eta, myEvent->lead_lepton_phi, myEvent->lead_lepton_pt);
		subleadLeptonFourMom.SetPtEtaPhiE(myEvent->sublead_lepton_pt, myEvent->sublead_lepton_eta, myEvent->sublead_lepton_phi, myEvent->sublead_lepton_pt);
		zFourMom = leadLeptonFourMom + subleadLeptonFourMom;

		h_Z_pt->Fill(zFourMom.Pt(), (myEvent->weight));
		h_Z_phi->Fill(zFourMom.Phi(), (myEvent->weight));
		h_Z_rapidity->Fill(zFourMom.Rapidity(), (myEvent->weight));

		h_lepton_pt0->Fill(myEvent->lead_lepton_pt, (myEvent->weight));
		h_lepton_pt1->Fill(myEvent->sublead_lepton_pt, (myEvent->weight));
		h_lepton_eta0->Fill(myEvent->lead_lepton_eta, (myEvent->weight));
		h_lepton_eta1->Fill(myEvent->sublead_lepton_eta, (myEvent->weight));
		h_lepton_phi0->Fill(myEvent->lead_lepton_phi, (myEvent->weight));
		h_lepton_phi1->Fill(myEvent->sublead_lepton_phi, (myEvent->weight));

		h_jet_pt0->Fill(myEvent->lead_jet_pt, (myEvent->weight));
		h_jet_pt1->Fill(myEvent->sublead_jet_pt, (myEvent->weight));
		h_jet_eta0->Fill(myEvent->lead_jet_eta, (myEvent->weight));
		h_jet_eta1->Fill(myEvent->sublead_jet_eta, (myEvent->weight));
		h_jet_phi0->Fill(myEvent->lead_jet_phi, (myEvent->weight));
		h_jet_phi1->Fill(myEvent->sublead_jet_phi, (myEvent->weight));

		h_dilepton_mass->Fill(myEvent->dilepton_mass, (myEvent->weight));
		h_nPV->Fill(myEvent->nPV, (myEvent->weight));
		h_WR_mass->Fill(myEvent->WR_mass, (myEvent->weight));
	
	}

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
	//hs->push_back(h_nPV);
	//hs->push_back(h_Z_phi);
	//hs->push_back(h_Z_rapidity);
	hs->push_back(h_Z_pt);
	hs->push_back(h_WR_mass);

}

void drawPlots(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, TString xtitle, TString fname, std::string channel)
{

	//data should be set to DYAMC
	gStyle->SetOptStat("");
	TLegend *leg = new TLegend( 0.60, 0.60, 0.9, 0.90 ) ;
	leg->AddEntry( hs_DYPowheg, "DY Powheg" ) ;
	leg->AddEntry( hs_DYMadIncl, "DY MAD HT binned" ) ;
	leg->AddEntry( hs_DYAmcIncl, "DY AMC" ) ;
	//leg->AddEntry( histos[2][0], "10 x WR 2600" ) ;
	//leg->AddEntry( hs_data, "Data");
	leg->SetFillColor( kWhite ) ;

	TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
	hs_DYPowheg->SetLineColor(kRed);
	hs_DYPowheg->SetLineWidth(3);
	hs_DYMadIncl->SetLineColor(kBlack);
	hs_DYMadIncl->SetLineWidth(3);
	hs_DYAmcIncl->SetLineColor(kBlue);
	hs_DYAmcIncl->SetLineWidth(3);
	hs_data->SetMarkerStyle(20);
	hs_data->SetMarkerSize(1);
	hs_data->SetMarkerColor(kBlack);

	Double_t eps = 0.001;
	TPad* p1 = new TPad("p1", "p1", 0, 0.25, 1, 1, 0);
	p1->Draw();
	TPad* p2 = new TPad("p2", "p2", 0, 0.1, 1, 0.25 + eps, 0);
	p2->Draw();
	p1->SetBottomMargin(0);
	p2->SetTopMargin(0);
	p1->cd();
	hs_data->SetStats(1);
	hs_DYPowheg->SetStats(1);
	TH1F *ratio_Powheg = (TH1F*)hs_data->Clone();
	TH1F *ratio_Mad = (TH1F*)hs_data->Clone();
	TH1F *ratio_Amc = (TH1F*)hs_data->Clone();
	TString plotTitle = "CMS Private   #surds = 13 TeV #int lumi = 2.6 fb^{-1}";
	hs_DYPowheg->SetTitle(plotTitle);
	hs_data->SetTitle(plotTitle);
	//hs_data->Draw("ep");
	hs_DYPowheg->Draw("histo");
	hs_DYMadIncl->Draw("histo same");
	hs_DYAmcIncl->Draw("histo same");
	//hs_data->Draw("epsame");
	TString ytitle = "Events";
	//TString ytitle = "Events/(";
	//ytitle += (hs_data->GetXaxis()->GetNbins());
	//ytitle += ")";
	hs_DYAmcIncl->GetYaxis()->SetTitle(ytitle.Data());
	/*if(fname.EqualTo("Mll") == true || fname.EqualTo("Z_pt") == true){
		//show the data over MC ratio on each Mll and Z_pt plot
		xtitle += " dataOvrAMC = ";
		xtitle += (to_string(dataOvrAmc)).c_str();
		xtitle += " dataOvrPow = ";
		xtitle += (to_string(dataOvrPowhegMassBinned)).c_str();
		xtitle += " dataOvrMad = ";
		xtitle += (to_string(dataOvrMad)).c_str();
	}*/
	hs_DYAmcIncl->GetXaxis()->SetTitle(xtitle.Data());

  	if(fname.EqualTo("Mlljj")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), hs_DYPowheg->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), hs_DYPowheg->GetYaxis()->SetTitle("Events"), hs_data->GetYaxis()->SetTitle("Events"), hs_DYAmcIncl->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), hs_DYAmcIncl->GetYaxis()->SetTitle("Events"), hs_DYMadIncl->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), hs_DYMadIncl->GetYaxis()->SetTitle("Events");

	ratio_Powheg->GetXaxis()->SetTitle(xtitle.Data());
	ratio_Powheg->GetXaxis()->SetTickSize(0.40);
	ratio_Powheg->GetXaxis()->SetTitleSize(0.25);
	ratio_Powheg->SetLabelSize(0.18, "x");
	leg->Draw();
	mycanvas->cd();
	p2->cd();	///<change to ratio TPad
	ratio_Powheg->Sumw2();
	ratio_Powheg->SetStats(0);
	ratio_Mad->Sumw2();
	ratio_Mad->SetStats(0);
	ratio_Amc->Sumw2();
	ratio_Amc->SetStats(0);

	ratio_Powheg->Divide(hs_DYPowheg);
	ratio_Powheg->SetMarkerStyle(20);
	ratio_Powheg->SetMarkerColor(kRed);
	ratio_Powheg->SetLabelSize(0.18, "y");
	ratio_Powheg->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Powheg->GetYaxis()->SetNdivisions(505);

	ratio_Mad->Divide(hs_DYMadIncl);
	ratio_Mad->SetMarkerStyle(21);
	ratio_Mad->SetMarkerColor(kBlack);
	ratio_Mad->SetLabelSize(0.18, "y");
	ratio_Mad->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Mad->GetYaxis()->SetNdivisions(505);

	ratio_Amc->Divide(hs_DYAmcIncl);
	ratio_Amc->SetMarkerStyle(22);
	ratio_Amc->SetMarkerColor(kBlue);
	ratio_Amc->SetLabelSize(0.18, "y");
	ratio_Amc->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Amc->GetYaxis()->SetNdivisions(505);

	ratio_Mad->Draw("p");
	ratio_Amc->Draw("p");
	ratio_Powheg->Draw("p");
	float xmax = ratio_Powheg->GetXaxis()->GetXmax();
	float xmin = ratio_Powheg->GetXaxis()->GetXmin();
	TF1 *f1 = new TF1("f1", "1", xmin, xmax);
	ratio_Powheg->Draw("p");
	ratio_Mad->Draw("psame");
	ratio_Amc->Draw("psame");
	f1->Draw("same");
	mycanvas->cd();

	TString fn = fname + "_" + channel.c_str() + "_signalRegionWithMllWeight_ratioWRTDYAMC";

	mycanvas->Print((fn + "_highestZoomRatio.pdf").Data());
	mycanvas->Print((fn + "_highestZoomRatio.png").Data());

	/*
	if(fname.EqualTo("Mll") == true || fname.EqualTo("Z_pt") == true || fname.EqualTo("Mlljj") == true){
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
	ratio_Powheg->GetYaxis()->SetRangeUser(0.6, 1.4);
	ratio_Mad->GetYaxis()->SetRangeUser(0.6, 1.4);
	ratio_Amc->GetYaxis()->SetRangeUser(0.6, 1.4);

	mycanvas->Update();
	mycanvas->Print((fn + "_noZoomRatio.pdf").Data());
	mycanvas->Print((fn + "_noZoomRatio.png").Data());
	
	p1->SetLogy();
	mycanvas->Print((fn + "_log.pdf").Data());
	mycanvas->Print((fn + "_log.png").Data());
 
	mycanvas->Close();

}
