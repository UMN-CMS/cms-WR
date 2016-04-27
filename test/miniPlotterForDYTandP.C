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

bool useJorgeMinitrees = false;
bool isReweighted = false;
Selector::tag_t channel = Selector::MuMu;
//Selector::tag_t channel = Selector::EE;

void writeIntegralsToTxtFile(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, Float_t minMll, Float_t maxMll, Float_t minSubleadLeptonPt, Float_t minLeadLeptonPt, Float_t maxLeptonEta, Float_t minLeadJetPt, Int_t minNJets, Float_t minSubleadJetPt);
void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t leadJetPtCut, Float_t leadLeptonPtCut, Float_t subleadLeptonPtCut, Float_t upperMllCut, Float_t lowerMllCut, Float_t leptonEtaCut, Int_t minNJets, Float_t subleadJetPtCut, Float_t normRescale);
void drawPlots(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, TString xtitle, TString fname, Float_t minMll, Float_t maxMll, Float_t minSubleadLeptonPt, Float_t minLeadLeptonPt, Float_t maxLeptonEta, Float_t minLeadJetPt, Int_t minNJets, Float_t minSubleadJetPt);
void miniPlotterForDYTandP()
{

	TString treeName = "treeDyCheck";
	TChain * chain_DYPowheg = new TChain(treeName,"DYPowhegM50to120");
	TChain * chain_DYMadIncl = new TChain(treeName,"DYMadgraphInclusive");
	TChain * chain_DYAmcIncl = new TChain(treeName,"DYAMCInclusive");
	TChain * chain_data = new TChain(treeName,"Data");

	switch (channel) {
	case Selector::EE:
		chain_DYPowheg->Add("../rootFiles/selected_tree_DYPOWHEG_dytagandprobeEE_v15noRoch.root");
		chain_DYMadIncl->Add("../rootFiles/selected_tree_DYMAD_dytagandprobeEE_v15noRoch.root"); // 0 - Electrons
		chain_DYAmcIncl->Add("../rootFiles/selected_tree_DYAMC_dytagandprobeEE_v15noRoch.root");
		chain_data->Add("../rootFiles/selected_tree_data_dytagandprobeEE_v14.root");
		break;
	case Selector::MuMu:
		if(!useJorgeMinitrees) chain_DYPowheg->Add("../rootFiles/selected_tree_DYPOWHEG_dytagandprobeMuMu_v15NoToysOrXsxnWeightInAnalysisCpp.root");
		if(!useJorgeMinitrees) chain_DYMadIncl->Add("../rootFiles/selected_tree_DYMAD_dytagandprobeMuMu_v15NoToysOrXsxnWeightInAnalysisCpp.root"); // 1 - Muons
		if(!useJorgeMinitrees) chain_DYAmcIncl->Add("../rootFiles/selected_tree_DYAMC_dytagandprobeMuMu_v15NoToysOrXsxnWeightInAnalysisCpp.root");
		if(!useJorgeMinitrees) chain_data->Add("../rootFiles/selected_tree_data_dytagandprobeMuMu_v15NoToysInAnalysisCpp.root");
		
		if(useJorgeMinitrees) chain_DYPowheg->Add("../rootFiles/selected_tree_DYPOWHEG_dytagandprobeMuMu_jorgeMinitreeWithoutToysOrXsxnWeightInAnalysisCpp.root");
		if(useJorgeMinitrees) chain_DYMadIncl->Add("../rootFiles/selected_tree_DYMAD_dytagandprobeMuMu_jorgeMinitreeWithoutToysOrXsxnWeightInAnalysisCpp.root"); // 1 - Muons
		if(useJorgeMinitrees) chain_DYAmcIncl->Add("../rootFiles/selected_tree_DYAMC_dytagandprobeMuMu_jorgeMinitreeWithoutToysOrXsxnWeightInAnalysisCpp.root");
		if(useJorgeMinitrees) chain_data->Add("../rootFiles/selected_tree_data_dytagandprobeMuMu_jorgeMinitreeWithoutToysInAnalysisCpp.root");
		
		break;
	default:
		std::cout << "Unknown tag" << std::endl;
	}

	Selector myEvent_DYPowheg;
	Selector myEvent_DYMadIncl;
	Selector myEvent_DYAmcIncl;
	Selector myEvent_data;

	myEvent_DYPowheg.SetBranchAddresses(chain_DYPowheg);
	myEvent_DYMadIncl.SetBranchAddresses(chain_DYMadIncl);
	myEvent_DYAmcIncl.SetBranchAddresses(chain_DYAmcIncl);
	myEvent_data.SetBranchAddresses(chain_data);

	Float_t intLumi = 2640.523267;
	Float_t dyPowZMassBin = (1975./2836876.);
	Float_t dyMad = (5991./9042031.);
	Float_t dyAmc = (5915./19548618.);
	if(channel == Selector::EE) dyPowZMassBin = (1997./49921782.);
	Float_t minLeadJetPt = -10, minLeadLeptonPt = 33, minSubleadLeptonPt = 20, maxMll = 123, minMll = 57, maxLeptonEta = 2.4, minSubleadJetPt = -10;
	Int_t minNumJets = -1;
	std::vector<TH1F*> hs_DYPowheg;
	MakeHistos(chain_DYPowheg, &myEvent_DYPowheg, &hs_DYPowheg, minLeadJetPt, minLeadLeptonPt, minSubleadLeptonPt, maxMll, minMll, maxLeptonEta, minNumJets, minSubleadJetPt, intLumi*dyPowZMassBin);
	std::vector<TH1F*> hs_DYMadIncl;
	MakeHistos(chain_DYMadIncl, &myEvent_DYMadIncl, &hs_DYMadIncl, minLeadJetPt, minLeadLeptonPt, minSubleadLeptonPt, maxMll, minMll, maxLeptonEta, minNumJets, minSubleadJetPt, intLumi*dyMad);
	std::vector<TH1F*> hs_DYAmcIncl;
	MakeHistos(chain_DYAmcIncl, &myEvent_DYAmcIncl, &hs_DYAmcIncl, minLeadJetPt, minLeadLeptonPt, minSubleadLeptonPt, maxMll, minMll, maxLeptonEta, minNumJets, minSubleadJetPt, intLumi*dyAmc);

	std::vector<TH1F*> hs_data;
	MakeHistos(chain_data, &myEvent_data, &hs_data, minLeadJetPt, minLeadLeptonPt, minSubleadLeptonPt, maxMll, minMll, maxLeptonEta, minNumJets, minSubleadJetPt, 1.0);

	unsigned int nPlots = hs_DYPowheg.size();

	// hs_data[13]->SetLineColor(kRed);
	// hs_data[13]->Draw();
	// hs_DYMadIncl[13]->Draw("same");

	TString xtitles[] = {"leading lepton p_{T}", "subleading lepton p_{T}", "leading jet p_{T}", "subleading jet p_{T}", "leading lepton #eta", "subleading lepton #eta", "leading jet #eta", "subleading jet #eta", "leading lepton #phi", "subleading lepton #phi", "leading jet #phi", "subleading jet #phi", "dilepton mass", "nPV", "Z #phi", "Z rapidity", "Z p_{T}","nPV no PU weight","nJets"};

	TString fnames[] = {"l1_pt", "l2_pt", "j1_pt", "j2_pt", "l1_eta", "l2_eta", "j1_eta", "j2_eta", "l1_phi", "l2_phi", "j1_phi", "j2_phi", "Mll", "nPV", "Z_phi", "Z_rapidity", "Z_pt","nPV_noPUwgt","nJets"};


	int i = 0;
	for(unsigned int i = 0; i < nPlots; i++) {
		std::string s = std::to_string(i);
		drawPlots(hs_DYPowheg[i], hs_DYMadIncl[i], hs_DYAmcIncl[i], hs_data[i], xtitles[i], fnames[i], minMll, maxMll, minSubleadLeptonPt, minLeadLeptonPt, maxLeptonEta, minLeadJetPt, minNumJets, minSubleadJetPt);
	}

}//end miniPlotterForDYTandP()
void writeIntegralsToTxtFile(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, Float_t minMll, Float_t maxMll, Float_t minSubleadLeptonPt, Float_t minLeadLeptonPt, Float_t maxLeptonEta, Float_t minLeadJetPt, Int_t minNJets, Float_t minSubleadJetPt)
{

	//define the integral ranges around the Z peak
	Double_t zCentr = 91.1876;
	std::vector<Double_t> zMassRanges;
	zMassRanges.push_back(5);
	zMassRanges.push_back(10);
	zMassRanges.push_back(15);
	zMassRanges.push_back(20);
	zMassRanges.push_back(25);
	zMassRanges.push_back(30);
	unsigned int num = zMassRanges.size();

	//write the integral of each dilepton_mass histo over a certain range into a txt file
	std::string zPeakEvtCountFile = "ZPeakRegionIntegralsFromSean.txt";
	ofstream writeToZpeakFile(zPeakEvtCountFile.c_str(), ofstream::app);
	writeToZpeakFile << "\t" << std::endl;
	writeToZpeakFile << "\t" << std::endl;
	std::string weightString = "Mll normalization IS NOT applied";
	if(isReweighted) weightString = "Mll normalization IS applied";
	if(useJorgeMinitrees) weightString += "  used Jorges minitrees with Roch corrections for MC and data";
	if(!useJorgeMinitrees) weightString += "  used Shervins V15 minitrees for MC and data";
	writeToZpeakFile << "MC processed WITHOUT systematics re-evaluation or Xsxn normalization in analysis.cpp" <<std::endl;
	writeToZpeakFile << weightString << std::endl;
	writeToZpeakFile << "min dilepton mass =\t"<< minMll <<"\tmax dilepton mass =\t"<< maxMll <<"\tmin lead jet pt =\t"<< minLeadJetPt <<"\tmin sublead jet pt =\t"<< minSubleadJetPt << "\tmin num jets =\t"<< minNJets << std::endl;
	writeToZpeakFile << "min sublead lepton pt =\t"<< minSubleadLeptonPt <<"\tmin lead lepton pt =\t"<< minLeadLeptonPt << "\tmax lepton |eta| =\t"<< maxLeptonEta << std::endl;
	std::string channelHeader = (channel == Selector::EE) ? "ELECTRON" : "MUON";
	writeToZpeakFile << channelHeader << std::endl;
	writeToZpeakFile << "\t"<< "DATA\t"<<"AMCNLO\t"<<"MADGRAPH\t"<<"POWHEGM50to120"<< std::endl;

	for(unsigned int i=0; i<num; i++){
		Int_t lowBin = hs_data->GetXaxis()->FindBin(zCentr - zMassRanges[i]);
		Int_t highBin = hs_data->GetXaxis()->FindBin(zCentr + zMassRanges[i]);
	
		/*	

#ifdef DBG
		std::cout<<"scanning through histo named\t"<< hs_data->GetName() << std::endl;
		std::cout<<"looking for values =\t" << zCentr-zMassRanges[i] << "\tand =\t"<< zCentr+zMassRanges[i] << std::endl;
		std::cout<<"max bin=\t"<< hs_data->GetNbinsX() << "\thas central value =\t"<< hs_data->GetXaxis()->GetBinCenter(hs_data->GetNbinsX()) << std::endl;
		std::cout<<"lowBin center = "<< hs_data->GetXaxis()->GetBinCenter(lowBin) <<"\thighBin center = "<< hs_data->GetXaxis()->GetBinCenter(highBin) << std::endl;
#endif
*/

		writeToZpeakFile <<"Z peak +/- "<<zMassRanges[i]<<" GeV\t"<< hs_data->Integral(lowBin, highBin) <<"\t"<< hs_DYAmcIncl->Integral(lowBin, highBin) <<"\t"<< hs_DYMadIncl->Integral(lowBin, highBin) <<"\t"<< hs_DYPowheg->Integral(lowBin, highBin) << std::endl;

	}//end loop over Z mass ranges
	
	writeToZpeakFile.close();

}//end writeIntegralsToTxtFile()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t leadJetPtCut, Float_t leadLeptonPtCut, Float_t subleadLeptonPtCut, Float_t upperMllCut, Float_t lowerMllCut, Float_t leptonEtaCut, Int_t minNJets, Float_t subleadJetPtCut, Float_t normRescale)
{

	Float_t reweight = 1.0;	///<reweight the MC by this scale factor to account for discrepancies in ZToll dilepton mass distribution
	if(isReweighted){
		TString chTitle = chain->GetTitle();
		if( !(chTitle.Contains("Data")) ){
			reweight = (channel == Selector::EE) ? (1) : (1);
		}
		std::cout<<"reweight =\t"<< reweight << std::endl;
	}//end if(isReweighted)
	TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0", "", 100, -20, 600);
	TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0", "", 50, -3.2, 3.2);
	TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0", "", 50, -3.2, 3.2);
	TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1", "", 120, -15, 400);
	TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1", "", 50, -3, 3);
	TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1", "", 50, -3.15, 3.15);

	TH1F *h_jet_pt0 = new TH1F("h_jet_pt0", "", 50, -5, 300);
	TH1F *h_jet_eta0 = new TH1F("h_jet_eta0", "", 50, -3.2, 3.2);
	TH1F *h_jet_phi0 = new TH1F("h_jet_phi0", "", 50, -3.2, 3.2);
	TH1F *h_jet_pt1 = new TH1F("h_jet_pt1", "", 50, -5, 300);
	TH1F *h_jet_eta1 = new TH1F("h_jet_eta1", "", 50, -3.2, 3.2);
	TH1F *h_jet_phi1 = new TH1F("h_jet_phi1", "", 50, -3.2, 3.2);

	//TH1F *h_WR_mass = new TH1F("h_WR_mass", "", 50, 0, 2000);
	float dilepton_max = 125.;
	if(channel == Selector::EMu)
		dilepton_max = 1000;
	TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass", "", 280, 55, dilepton_max);
	TH1F *h_nPV = new TH1F("h_nPV", "", 100, 0, 100);
	TH1F *h_nPV_noPUWeight = new TH1F("h_nPV_noPUWeight", "", 100, 0, 100);

	TH1F *h_Z_phi = new TH1F("h_Z_phi", "", 50, -3.2, 3.2);
	TH1F *h_Z_rapidity = new TH1F("h_Z_rapidity", "", 70, -5., 8.);
	TH1F *h_Z_pt = new TH1F("h_Z_pt", "", 80, -10., 300.);

	TH1F *h_nJets = new TH1F("h_nJets", "", 26, -1, 25);
	
	Long64_t nEntries = chain->GetEntries();

	cout << nEntries << endl;
	
	//TString chTitle = chain->GetTitle();
	
	for(int ev = 0; ev < nEntries; ++ev) {
		chain->GetEntry(ev);
		if(myEvent->dilepton_mass > upperMllCut || myEvent->dilepton_mass < lowerMllCut) continue;
		if(myEvent->njets < minNJets) continue;
		if(myEvent->lead_lepton_pt < leadLeptonPtCut || myEvent->sublead_lepton_pt < subleadLeptonPtCut || std::fabs(myEvent->sublead_lepton_eta) > leptonEtaCut || std::fabs(myEvent->lead_lepton_eta) > leptonEtaCut) continue;
		if(myEvent->lead_jet_pt < leadJetPtCut || myEvent->sublead_jet_pt < subleadJetPtCut) continue;

		//if( !(chTitle.Contains("Data")) ){
		//	(myEvent->weight) *= (2640.523267/2640.);
		//}

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

		//h_WR_mass->Fill(myEvent->WR_mass, (myEvent->weight));
		h_dilepton_mass->Fill(myEvent->dilepton_mass, (myEvent->weight));
		h_nPV->Fill(myEvent->nPV, (myEvent->weight));
		h_nPV_noPUWeight->Fill(myEvent->nPV, (myEvent->weight) / (myEvent->pu_weight));
		h_nJets->Fill(myEvent->njets, (myEvent->weight));
	
	}

	h_Z_pt->Scale(reweight*normRescale);
	h_Z_phi->Scale(reweight*normRescale);
	h_Z_rapidity->Scale(reweight*normRescale);

	h_lepton_pt0->Scale(reweight*normRescale);
	h_lepton_pt1->Scale(reweight*normRescale);
	h_lepton_eta0->Scale(reweight*normRescale);
	h_lepton_eta1->Scale(reweight*normRescale);
	h_lepton_phi0->Scale(reweight*normRescale);
	h_lepton_phi1->Scale(reweight*normRescale);

	h_jet_pt0->Scale(reweight*normRescale);
	h_jet_pt1->Scale(reweight*normRescale);
	h_jet_eta0->Scale(reweight*normRescale);
	h_jet_eta1->Scale(reweight*normRescale);
	h_jet_phi0->Scale(reweight*normRescale);
	h_jet_phi1->Scale(reweight*normRescale);

	//h_WR_mass->Scale(reweight*normRescale);
	h_dilepton_mass->Scale(reweight*normRescale);
	h_nPV->Scale(reweight*normRescale);
	h_nPV_noPUWeight->Scale(reweight*normRescale);
	h_nJets->Scale(reweight*normRescale);


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
	//hs->push_back(h_WR_mass);
	hs->push_back(h_dilepton_mass);
	hs->push_back(h_nPV);
	hs->push_back(h_Z_phi);
	hs->push_back(h_Z_rapidity);
	hs->push_back(h_Z_pt);
	hs->push_back(h_nPV_noPUWeight);
	hs->push_back(h_nJets);

}

void drawPlots(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, TString xtitle, TString fname, Float_t minMll, Float_t maxMll, Float_t minSubleadLeptonPt, Float_t minLeadLeptonPt, Float_t maxLeptonEta, Float_t minLeadJetPt, Int_t minNJets, Float_t minSubleadJetPt)
{

	if(fname.EqualTo("Mll") == true) writeIntegralsToTxtFile(hs_DYPowheg,hs_DYMadIncl,hs_DYAmcIncl,hs_data,minMll,maxMll,minSubleadLeptonPt,minLeadLeptonPt,maxLeptonEta, minLeadJetPt, minNJets, minSubleadJetPt);

	gStyle->SetOptStat("eou");
	TLegend *leg = new TLegend( 0.80, 0.50, 0.98, 0.70 ) ;
	leg->AddEntry( hs_DYPowheg, "DY Powheg" ) ;
	leg->AddEntry( hs_DYMadIncl, "DY MAD Incl" ) ;
	//leg->AddEntry( hs_DYAmcIncl, "DY AMC Incl" ) ;
	//leg->AddEntry( histos[2][0], "10 x WR 2600" ) ;
	leg->AddEntry( hs_data, "Data");
	leg->SetFillColor( kWhite ) ;

	TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
	hs_DYPowheg->SetLineColor(kRed);
	hs_DYPowheg->SetLineWidth(3);
	hs_DYMadIncl->SetLineColor(kBlack);
	hs_DYMadIncl->SetLineWidth(3);
	//hs_DYAmcIncl->SetLineColor(kBlue);
	//hs_DYAmcIncl->SetLineWidth(3);
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
	hs_DYPowheg->SetTitle("CMS Preliminary");
	hs_data->SetTitle("CMS Preliminary");
	//th->Draw("histo");
	//hs_data->Draw("epsame");
	hs_data->Draw("ep");
	hs_DYPowheg->Draw("histo same");
	hs_DYMadIncl->Draw("histo same");
	//hs_DYAmcIncl->Draw("histo same");
	hs_data->Draw("epsame");
	TString ytitle = "Events/(";
	ytitle += (hs_data->GetXaxis()->GetNbins());
	ytitle += ")";
	//hs_DYAmcIncl->GetYaxis()->SetTitle(ytitle.Data());
	//hs_DYAmcIncl->GetXaxis()->SetTitle(xtitle.Data());

	ratio_Powheg->GetXaxis()->SetTitle(xtitle.Data());
	ratio_Powheg->GetXaxis()->SetTickSize(0.40);
	ratio_Powheg->GetXaxis()->SetTitleSize(0.2);
	ratio_Powheg->SetLabelSize(0.1, "x");
	leg->Draw();
	mycanvas->cd();
	p2->cd();	///<change to ratio TPad
	ratio_Powheg->Sumw2();
	ratio_Powheg->SetStats(0);
	ratio_Mad->Sumw2();
	ratio_Mad->SetStats(0);
	//ratio_Amc->Sumw2();
	//ratio_Amc->SetStats(0);

	ratio_Powheg->Divide(hs_DYPowheg);
	ratio_Powheg->SetMarkerStyle(20);
	ratio_Powheg->SetMarkerColor(kRed);
	ratio_Powheg->SetLabelSize(0.1, "y");
	ratio_Powheg->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_Powheg->GetYaxis()->SetNdivisions(505);

	ratio_Mad->Divide(hs_DYMadIncl);
	ratio_Mad->SetMarkerStyle(21);
	ratio_Mad->SetMarkerColor(kBlack);
	ratio_Mad->SetLabelSize(0.1, "y");
	ratio_Mad->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_Mad->GetYaxis()->SetNdivisions(505);

	//ratio_Amc->Divide(hs_DYAmcIncl);
	//ratio_Amc->SetMarkerStyle(22);
	//ratio_Amc->SetMarkerColor(kBlue);
	//ratio_Amc->SetLabelSize(0.1, "y");
	//ratio_Amc->GetYaxis()->SetRangeUser(0.5, 1.5);
	//ratio_Amc->GetYaxis()->SetNdivisions(505);

	ratio_Mad->Draw("p");
	//ratio_Amc->Draw("p");
	ratio_Powheg->Draw("p");
	float xmax = ratio_Powheg->GetXaxis()->GetXmax();
	float xmin = ratio_Powheg->GetXaxis()->GetXmin();
	TF1 *f1 = new TF1("f1", "1", xmin, xmax);
	ratio_Powheg->Draw("p");
	ratio_Mad->Draw("psame");
	//ratio_Amc->Draw("psame");
	f1->Draw("same");
	mycanvas->cd();

	TString fn = "";
	TString cuts = "_minLeadLeptPt_" + to_string(minLeadLeptonPt) +"_minSubleadLeptPt_" + to_string(minSubleadLeptonPt) + "_minLeadJetPt_" + to_string(minLeadJetPt) + "_minSubleadJetPt_" + to_string(minSubleadJetPt) + "_minNJets_" + to_string(minNJets);
	if(isReweighted) cuts += "_isReweighted";
	else cuts += "_isNotReweighted";

	if(channel == Selector::EE)
		fn = fname + cuts + "_eeChannel";
	if(channel == Selector::MuMu)
		fn = fname + cuts + "_mumuChannel";

	mycanvas->Print((fn + ".pdf").Data());
	mycanvas->Print((fn + ".png").Data());
	p1->SetLogy();
	mycanvas->Print((fn + "_log.pdf").Data());
	mycanvas->Print((fn + "_log.png").Data());

	mycanvas->Close();
	if(fname.EqualTo("Mll") == true){
		TFile f("seanHistos.root","RECREATE");
		hs_data->SetName("h_Mll_DATA"), hs_DYAmcIncl->SetName("h_Mll_DY"), hs_DYPowheg->SetName("h_Mll_DY_pow"), hs_DYMadIncl->SetName("h_Mll_DY_mad");
		hs_data->Write();
		hs_DYAmcIncl->Write();
		hs_DYPowheg->Write();
		hs_DYMadIncl->Write();
	}
}
