/**
 * this is the main program for the analysis
 * please document it here
 */

#include "TChain.h"
#include "TFile.h"
// #include "TLorentzVector.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom3.h"
// #include <iostream>
// // #include "ExoAnalysis/cmsWR/interface/FitRooDataSet.h"

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"
#include "ExoAnalysis/cmsWR/interface/Selector.h"
// #include "../interface/FitRooDataSet.h"
// #include "../src/Selector.cc"
// #include "../src/miniTreeEvent.cc"

#include <vector>
#include <string>
#include <fstream>
// #include "../interface/miniTreeEvent.h"
#include "FitRooDataSet.h"
#include "rooFitFxns.h"
#include "ToyThrower.h"
#include "ExoAnalysis/cmsWR/interface/configReader.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

int main(void)
{

	float integratedLumi = 2.4; ///\todo should not be hard-coded!!!
	using namespace RooFit;
	std::cout << "******************************* Analysis ******************************" << std::endl;
	std::cout << "[WARNING] no weights associated to jets yet" << std::endl;

	configReader myReader("configs/2015-v1.conf");
	myReader.getNorm1fb("TTJets_DiLept_v2");

// KEY: TDirectoryFile   miniTree_flavoursideband;1      miniTree_flavoursideband
// KEY: TDirectoryFile   miniTree_lowdileptonsideband;1  miniTree_lowdileptonsideband
// KEY: TDirectoryFile   miniTree_signal_ee;1    miniTree_signal_ee
// KEY: TDirectoryFile   miniTree_signal_mumu;1  miniTree_signal_mumu
// KEY: TDirectoryFile   miniTree_dytagandprobe;1        miniTree_dytagandprobe
// KEY: TDirectoryFile   zToEEAnalyzer;1 zToEEAnalyzer
// KEY: TDirectoryFile   zToMuMuAnalyzer;1       zToMuMuAnalyzer

	std::vector<std::string> TTchainNames;
	TTchainNames.push_back("TTJets_DiLept_v1");
	TTchainNames.push_back("TTJets_DiLept_v2");
	TChain *chain = (myReader.getMiniTreeChain(TTchainNames, "miniTree_dytagandprobe"));
	std::cout << chain->GetEntries() << std::endl;
	// if you want to check if the config file is read correctly:
#ifdef DEBUG
	std::cout << myReader << std::endl;
#endif

	TString mode = "test";
	// Plotting trees
	TFile f("selected_tree_" + mode + ".root", "recreate");
	TTree * t1 = new TTree("t1", "");

	miniTreeEvent myEvent;

	myEvent.SetBranchAddresses(chain);
	Selector selEvent;
	selEvent.SetBranches(t1);

	Int_t nToys = 1;
	Int_t nEntries = chain->GetEntries();

	int isData = 1; // Fill in with 1 or 0 based on information from the trees
	TRandom3 Rand;
	const int Total_Number_of_Systematics_Smear = 1;// electron scale(MC)
	const int Total_Number_of_Systematics_Up_Down = 4;// muon id, muon iso, electron scale(data) and jet energy scale
	float Random_Numbers_for_Systematics_Smear[Total_Number_of_Systematics_Smear] = {0.};
	float Random_Numbers_for_Systematics_Up_Down[Total_Number_of_Systematics_Up_Down] = {0.};


	std::vector<std::string> List_Systematics;
	List_Systematics.push_back("smear");
	std::string word;
	std::ifstream Syst_File;

	std::cout << "[INFO] Reading systematics to be evaluated" << std::endl;
	std::string systFileName = "Systematics_To_Be_Evaluated.txt";
	Syst_File.open(systFileName);
	if(Syst_File.is_open() == false) {
		std::cerr << "[ERROR] File " << systFileName << " not opened. Check if the file exists" << std::endl;
		return 1;
	}

	while(Syst_File.peek() != EOF && Syst_File.good()) {
		if(Syst_File.peek() == 10) { // 10 = \n
			Syst_File.get();
			continue;
		}

		if(Syst_File.peek() == 35) { // 35 = #
			Syst_File.ignore(1000, 10); // ignore the rest of the line until \n
			continue;
		}

		Syst_File >> word;
		std::cout << word << std::endl;
		List_Systematics.push_back(word);
	}
	std::cout << "[INFO] nSyst = " << List_Systematics.size() << std::endl;
	if(List_Systematics.size() == 0) {
		std::cerr << "[ERROR] No systematics defined!" << std::endl;
		return 1;

	}
	std::cout << "[INFO] Running nToys = " << nToys << std::endl;
	for(int i = 0; i < nToys; i++) {
		Rand.SetSeed(i + 1);
		RooRealVar massWR("fourObjectMass", "fourObjectMass", 600, 6500);
		RooRealVar evtWeight("evtWeight", "evtWeight", -2, 2);
		RooArgSet vars(massWR, evtWeight);
		RooDataSet * tempDataSet = new RooDataSet("temp", "temp", vars);

		for(int ev = 0; ev < nEntries; ev++) {
			chain->GetEntry(ev);

			Selector tmp_selEvent(myEvent);
			selEvent = tmp_selEvent;

#ifdef DEBUG
			cout << "RUN=" << myEvent.run << endl;
			cout << "Mu" << endl;
			for(auto m : * (myEvent.muons_p4))
				cout << m.Pt() << " " << m.Eta() << endl;
			cout << "Jet" << endl;
			for(auto m : * (myEvent.jets_p4))
				cout << m.Pt() << " " << m.Eta() << endl;
#endif

			for(int Rand_Smear_Iter = 0; Rand_Smear_Iter < Total_Number_of_Systematics_Smear; Rand_Smear_Iter++)
				Random_Numbers_for_Systematics_Smear[Rand_Smear_Iter] = Rand.Gaus(0.0, 1.);

			//ToyThrower(myEvent, Random_Numbers_for_Systematics_Smear, Random_Numbers_for_Systematics_Up_Down, i + 1, List_Systematics, isData);



			// Select events with one good WR candidate
			// Tags:
			// 0 -- EEJJ Channel
			// 1 -- MuMuJJ Channel
			// 2 -- EMuJJ Channel



			if(selEvent.isPassing(Selector::EMu) && selEvent.dilepton_mass > 200) {
				float weight = selEvent.weight * myReader.getNorm1fb(selEvent.datasetName) * integratedLumi; // the weight is the event weight * single object weights

				massWR.setVal(selEvent.WR_mass);
				evtWeight.setVal(weight);
				tempDataSet->add(vars);

				t1->Fill();

			}

		}
		assert(tempDataSet->sumEntries() > 0);

//		t1->Write();

		tempDataSet->Print();

		RooFitResult * tempFitRslt = NULL;
		fitRooDataSet(tempFitRslt, tempDataSet, expPdfRooAbsPdf);

		std::cout << "Res=" << std::endl;
		expPdfRooAbsPdf->Print();

	}

	return 0;

}
