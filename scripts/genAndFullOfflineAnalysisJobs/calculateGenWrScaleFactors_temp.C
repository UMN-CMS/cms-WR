#include <TFile.h>
#include <TList.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TProfile.h>
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
#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

using namespace std;

void calculateGenWrScaleFactorsNULL(){
	
	///make a 2D plot with MNu as the vertical axis, MWR as the horizontal axis, and the fraction of GEN WR->eejj evts which
	///pass all offline cuts at each point in the (MWR, MNu) space
	///also print a table of these values

	///all input .root files should be in the same directory, and have file names which differ only in the WR and Nu mass values
	string dir= "PTHTOTREES/";
	string fileBegin = "analyzed_genWrToEEJJFullOfflineAnalysis_WR_";
	string fileEnd = "_QNUM.root";
	string fileMiddle = "_NU_";
	gStyle->SetTitleOffset(1.6,"Y");
	gStyle->SetOptStat("");

	string cutEfficiencyVsMassFile = "offlineEfficienciesVsMasses.txt";
	ofstream writeToEfficiencyFile(cutEfficiencyVsMassFile.c_str(),ofstream::app);
	Float_t passingPercentage=-1;
	TH2F * twoDimAcceptanceHist = new TH2F("twoDimAccHist","Percentage of events passing all GEN cuts",58,700,3600,72,0,3600);

	int maxWrMass = MXWRMSS, increment = INCRMT, minWrMass = MNWRMSS, minNuMass = MNNUMSS;
	//starting value for wrMass equals the minimum wr mass
	for(int wrMass=minWrMass; wrMass<=maxWrMass ; wrMass+=increment){
		///loop over WR mass values
		
		for(int nuMass=minNuMass; nuMass<wrMass ; nuMass+=increment){
			///loop over Nu mass values

			///define input root file name
			string pfn = dir+fileBegin+to_string(wrMass)+fileMiddle+to_string(nuMass)+fileEnd;
			
			///define one TChain to count the number of events generated, and record the desired WR and Nu masses
			///define another TChain to count the number of evts passing all offline cuts at GEN lvl (no HLT or ID)
			TChain * genInfo = new TChain("genNuAnalyzerOne/genNu");
			genInfo->Add(pfn.c_str());
			TChain * afterOfflineCuts = new TChain("genMatchedParticleAnalyzerThree/genLeptonsAndJetsWithAllCuts");
			afterOfflineCuts->Add(pfn.c_str());

			///calculate percentage of evts which pass GEN cuts, and store this percentage along with the nu and wr mass values in a txt file
			passingPercentage = (100)*((Float_t) afterOfflineCuts->GetEntries()/genInfo->GetEntries());
			writeToEfficiencyFile << passingPercentage <<"\tpercent of events with WR mass=\t"<< wrMass <<"\tand Nu mass=\t"<< nuMass <<"\tpass all cuts"<< endl;
			
			///fill a bin in the 2D histo with a weight equal to passingPercentage
			twoDimAcceptanceHist->Fill((Double_t) wrMass, (Double_t) nuMass,passingPercentage);
			genInfo->Delete();
			afterOfflineCuts->Delete();

		}///end loop over Nu mass values

	}///end loop over WR mass values	
	
	writeToEfficiencyFile.close();


	///set axis labels and draw the histo, save an image of it (both pdf and png formats), and close the txt file
	twoDimAcceptanceHist->GetXaxis()->SetTitle("WR Mass [GeV]");
	twoDimAcceptanceHist->GetYaxis()->SetTitle("Nu Mass [GeV]");
	TCanvas * r1 = new TCanvas("r1","r1",900,700);
	r1->cd();
	twoDimAcceptanceHist->Draw("COLZ");
	r1->SaveAs("twoDimGenWrAcceptances_afterAllCuts.png","recreate");
	r1->SaveAs("twoDimGenWrAcceptances_afterAllCuts.pdf","recreate");
	
}///end macroSandBox()

