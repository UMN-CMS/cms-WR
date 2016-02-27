/**
 * this is the main program for the analysis
 * please document it here
 */

#include "TChain.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom3.h"
// #include <iostream>
// // #include "ExoAnalysis/cmsWR/interface/FitRooDataSet.h"

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"
#include "ExoAnalysis/cmsWR/interface/Selector.h"

#include <vector>
#include <string>
#include <fstream>
#include "FitRooDataSet.h"
#include "rooFitFxns.h"
#include "ToyThrower.h"
#include "ExoAnalysis/cmsWR/interface/configReader.h"


//#define DEBUG

int main(void)
{
	char name[100];// to name the tree per iteration of the Toy
	float integratedLumi = 2.52e3; ///\todo should not be hard-coded!!!
	using namespace RooFit;
	std::cout << "******************************* Analysis ******************************" << std::endl;
	std::cout << "[WARNING] no weights associated to jets yet" << std::endl;

	configReader myReader("configs/2015-v1.conf");


#ifdef DEBUG
	std::cout << myReader << std::endl;
#endif

	Selector::tag_t channel = Selector::EMu;
	std::vector<TString> modes;

	modes.push_back("ttbar");
	modes.push_back("DYMuMu");
	modes.push_back("DYEE");
	modes.push_back("WJets");
	modes.push_back("WZ");
	modes.push_back("ZZ");
	modes.push_back("data_EMu");
	modes.push_back("data_EE");
	modes.push_back("data_MuMu");
	channel = Selector::EE;

	TString tree_channel = "";
	// Select the channel to be studied //
	if(channel == Selector::EE)
		tree_channel = "_signal_ee";
	else if(channel == Selector::MuMu)
		tree_channel = "_signal_mumu";
	else if(channel == Selector: EMu)
		tree_channel = "_flavoursideband";
	else {
		tree_channel = "";
		std::cerr << "[ERROR] No channel defined" << std::endl;
		return 1;
	}

	TString treeName = "miniTree" + tree_channel;

	for(auto m : modes) {
		int isData = 0; // Fill in with 1 or 0 based on information from the trees
		std::vector<std::string> TTchainNames;
		TString mode = m;
		bool run_toys = true;

		// Select the dataset to run over //
		if(mode.EqualTo("ttbar")) {
			//TTchainNames.push_back("TTJets_DiLept_v1");
			TTchainNames.push_back("TTJets_DiLept_v2");
		} else if(mode.EqualTo("DYMuMu")) {
			TTchainNames.push_back("DYToMuMu_powheg_50to120");
			TTchainNames.push_back("DYToMuMu_powheg_120to200");
			TTchainNames.push_back("DYToMuMu_powheg_200to400");
			TTchainNames.push_back("DYToMuMu_powheg_400to800");
			TTchainNames.push_back("DYToMuMu_powheg_800to1400");
			TTchainNames.push_back("DYToMuMu_powheg_1400to2300");
			TTchainNames.push_back("DYToMuMu_powheg_2300to3500");
			TTchainNames.push_back("DYToMuMu_powheg_3500to4500");
			TTchainNames.push_back("DYToMuMu_powheg_4500to6000");
			TTchainNames.push_back("DYToMuMu_powheg_6000toInf");
			if(channel == Selector::EMu)
				run_toys = false;
		} else if(mode.EqualTo("DYEE")) {
			TTchainNames.push_back("DYToEE_powheg_50to120");
			TTchainNames.push_back("DYToEE_powheg_120to200");
			TTchainNames.push_back("DYToEE_powheg_200to400");
			TTchainNames.push_back("DYToEE_powheg_400to800");
			TTchainNames.push_back("DYToEE_powheg_800to1400");
			TTchainNames.push_back("DYToEE_powheg_1400to2300");
			TTchainNames.push_back("DYToEE_powheg_2300to3500");
			TTchainNames.push_back("DYToEE_powheg_3500to4500");
			TTchainNames.push_back("DYToEE_powheg_4500to6000");
			TTchainNames.push_back("DYToEE_powheg_6000toInf");
			if(channel == Selector::EMu)
				run_toys = false;
		} else if(mode.EqualTo("WJets")) {
			TTchainNames.push_back("WJetsLNu");
			run_toys = false;
		} else if(mode.EqualTo("WZ")) {
			TTchainNames.push_back("WZ");
			run_toys = false;
		} else if(mode.EqualTo("ZZ")) {
			TTchainNames.push_back("ZZ");
			run_toys = false;
		} else if(mode.EqualTo("data_EMu")) {
			isData = 1;
			TTchainNames.push_back("MuEG_RunC");
			TTchainNames.push_back("MuEG_RunD_v3");
			TTchainNames.push_back("MuEG_RunD_v4");
		} else if(mode.EqualTo("data_EE")) {
			isData = 1;
			TTchainNames.push_back("DoubleEG_RunC");
			TTchainNames.push_back("DoubleEG_RunD_v3");
			TTchainNames.push_back("DoubleEG_RunD_v4");
		} else if(mode.EqualTo("data_MuMu")) {
			isData = 1;
			TTchainNames.push_back("SingleMu_RunC");
			TTchainNames.push_back("SingleMu_RunD_v3");
			TTchainNames.push_back("SingleMu_RunD_v4");
		}
		// Select the channel to be studied //
		if(channel == 0)
			tree_channel = "_lowdileptonsideband";//"_signal_ee";
		else if(channel == 1)
			tree_channel = "_lowdileptonsideband";//"_signal_mumu";
		else if(channel == 2)
			tree_channel = "_flavoursideband";
		else
			tree_channel = "";

		TChain *c = (myReader.getMiniTreeChain(TTchainNames, ("miniTree" + tree_channel).Data()));
		std::cout << c->GetEntries() << std::endl;

		// if you want to check if the config file is read correctly:
#ifdef DEBUG
		std::cout << myReader.getNorm1fb("TTJets_DiLept_v1") << std::endl;
#endif

		// Plotting trees
		TFile f("selected_tree_" + mode + tree_channel + std::to_string(channel) + ".root", "recreate");
		// store the fitted results for every toy in a tree
		TTree * tf1 = new TTree("tf1", "");
		Float_t normalization;
		std::vector<Float_t> fit_parameters, fit_parameter_errors;
		tf1->Branch("Normalization", &normalization);
		tf1->Branch("FitParameters", &fit_parameters);
		tf1->Branch("FitParameterErrors", &fit_parameter_errors);

		miniTreeEvent myEvent;

		myEvent.SetBranchAddresses(c);
		Selector selEvent;

		const Int_t nToys = 1;
		TTree * t1[nToys];
		Int_t nEntries = c->GetEntries();

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
			for(int Rand_Up_Down_Iter = 0; Rand_Up_Down_Iter < Total_Number_of_Systematics_Up_Down; Rand_Up_Down_Iter++)
				Random_Numbers_for_Systematics_Up_Down[Rand_Up_Down_Iter] = Rand.Gaus(0.0, 1.);
			RooRealVar massWR("fourObjectMass", "fourObjectMass", 600, 6500);
			RooRealVar evtWeight("evtWeight", "evtWeight", -2, 2);
			RooArgSet vars(massWR, evtWeight);
			RooDataSet * tempDataSet = new RooDataSet("temp", "temp", vars);
			sprintf(name, "Tree_Iter%i", i);
			t1[i] = new TTree(name, "");
			selEvent.SetBranches(t1[i]);

			for(int ev = 0; ev < nEntries; ev++) {

				if(ev % 50000 == 1) std::cout << std::endl << 100 * ev / nEntries << " % ..." << std::endl;
				c->GetEntry(ev);

#ifdef DEBUG
				std::cout << "RUN=" << myEvent.run << std::endl;
				std::cout << "Mu" << std::endl;
				for(auto m : * (myEvent.muons_p4))
					std::cout << m.Pt() << " " << m.Eta() << std::endl;
				std::cout << "Jet" << std::endl;
				for(auto m : * (myEvent.jets_p4))
					std::cout << m.Pt() << " " << m.Eta() << std::endl;
#endif

				for(int Rand_Smear_Iter = 0; Rand_Smear_Iter < Total_Number_of_Systematics_Smear; Rand_Smear_Iter++)
					Random_Numbers_for_Systematics_Smear[Rand_Smear_Iter] = Rand.Gaus(0.0, 1.);

				ToyThrower( &myEvent, Random_Numbers_for_Systematics_Smear, Random_Numbers_for_Systematics_Up_Down, i + 1, List_Systematics, isData);

				Selector tmp_selEvent(myEvent);
				selEvent = tmp_selEvent;

				// Select events with one good WR candidate
				// Tags:
				// 0 -- EEJJ Channel
				// 1 -- MuMuJJ Channel
				// 2 -- EMuJJ Channel

#ifdef DEBUG
				if(selEvent.isPassing(channel)) {
					std::cout << std::endl << selEvent.dilepton_mass << std::endl;
					std::cout << "Mu " << ev  << std::endl;
					for(auto m : * (myEvent.muons_p4))
						std::cout << m.Pt() << " " << m.Eta() << std::endl;
					std::cout << "Ele" << std::endl;
					for(auto m : * (myEvent.electrons_p4))
						std::cout << m.Pt() << " " << m.Eta() << std::endl;
					std::cout << "Jet" << std::endl;
					for(auto m : * (myEvent.jets_p4))
						std::cout << m.Pt() << " " << m.Eta() << std::endl;
				}
#endif

				if(selEvent.isPassing(channel)) {

					bool dilepton_cut = (channel == Selector::EMu) ? (selEvent.dilepton_mass > 200) : (selEvent.dilepton_mass < 200);

					if(!dilepton_cut)
						continue;

					//std::cout<<selEvent.weight << " " <<myReader.getNorm1fb(selEvent.datasetName) << std::endl;

					if(isData == 0)
						selEvent.weight = myEvent.weight * myReader.getNorm1fb(selEvent.datasetName) * integratedLumi; // the weight is the event weight * single object weights
					else
						selEvent.weight = 1.0;
					massWR.setVal(selEvent.WR_mass);
					evtWeight.setVal(selEvent.weight);
					tempDataSet->add(vars);
					selEvent.nPV = myEvent.nPV;


					t1[i]->Fill();

				}

			}

			if(i == 0) {
				t1[i]->Write();
				tempDataSet->Write();
			}


			tempDataSet->Print();

			if(mode.EqualTo("ttbar")) {

				assert(tempDataSet->sumEntries() > 0);
				expPower.setVal(-0.004);

				RooFitResult * tempFitRslt = NULL;
				// fit dataset to given PDF
				fitRooDataSet(tempFitRslt, tempDataSet, expPdfRooAbsPdf);

				// std::cout << "Res=" << std::endl;
				// expPdfRooAbsPdf->Print();

				// dataset normalization is the number of entries in the dataset
				normalization = tempDataSet->sumEntries();
				// set of variables in the PDF
				RooArgSet *vset = expPdfRooAbsPdf->getVariables();
				// these variables are stored in vectors
				// clear these vector for each iteration
				fit_parameters.clear();
				fit_parameter_errors.clear();

				// loop over RooRealVars in the set
				TIterator * iter = vset->createIterator();
				TObject * var = iter->Next();
				RooRealVar *var_pdf;
				while (var) {
					// ignore the M_WR variable
					if(strcmp(var->GetName(), "fourObjectMass") != 0) {
						var_pdf = (RooRealVar*)vset->find(var->GetName());
						// store the value of the fitted parameters and the corresponding errors
						fit_parameters.push_back(var_pdf->getVal());
						fit_parameter_errors.push_back(var_pdf->getError());
					}
					var = iter->Next();
				}
				// fill the tree with the normalization, parameters, and errors
				tf1->Fill();
			}
			if(!run_toys)
				break;
		}
		// only write the fitted branch for the modes that make sense (ttbar, DY, and data)
		if(mode.EqualTo("ttbar")) // add the other modes later
			tf1->Write();
	}
	return 0;

}
