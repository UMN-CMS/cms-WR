/**
 * this is the main program for the analysis
 * please document it here
 */

#include "TH1F.h"
#include "TChain.h"
#include "TFile.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom3.h"
// #include <iostream>
// // #include "ExoAnalysis/cmsWR/interface/FitRooDataSet.h"

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"
#include "ExoAnalysis/cmsWR/interface/AnalysisResult.h"
#include "ExoAnalysis/cmsWR/interface/Selector.h"

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include "FitRooDataSet.h"
#include "rooFitFxns.h"
#include "ToyThrower.h"
#include "analysisTools.h"
#include "configReader.h"

#include <unordered_set>

#include "Calibration/ZFitter/interface/EnergyScaleCorrection_class.hh"

#include <TStopwatch.h>
#define _ENDSTRING std::string::npos
//#define DEBUG
//#define DEBUGG
//process only 50000 events when DEBUGG is defined

/**
TT
DY TANDP POWHEG AMC MAD POWINCL
W
WZ
ZZ
data TANDP
*/

/** \class chainNames
	\brief this class helps in finding the right tree name based on the sample, sideband and channel one wants to analyze
*/

class chainNames
{

public:
	chainNames(): ///< default constructor
		all_modes(  // list of all possible modes
	{"TT", "W", "WZ", "ZZ", "data", "DYPOWHEG", "DYMADHT", "DYAMC", "DYMAD", "DYPOWINCL", "signal"
	}
	)
	{
	};

	bool isData(std::string mode)
	{
		if(mode == "data") return true;
		return false;
	}

	std::vector<std::string> getChainNames(std::string mode, Selector::tag_t channel, bool isTagAndProbe)
	{

		std::vector<std::string> TTchainNames;
		if(checkValidMode(mode) == false) {
			std::cerr << "[ERROR]" << std::endl;
			return TTchainNames;
		}
		if(mode == "TT") {
			//TTchainNames.push_back("TTJets_DiLept_v1");
			TTchainNames.push_back("TTJets_DiLept_v2");
		} else if(mode.find("DY") != _ENDSTRING) {
			//if(mode.Contains("TANDP") ) tree_channel = "_dytagandprobe";
			std::string tagName = "";
			if(channel == Selector::EE) tagName = "EE";
			if(channel == Selector::MuMu) tagName = "MuMu";
			if(channel == Selector::EMu) { ///\todo to be fixes, it should be possible
				std::cout << "ERROR looking for DY in EMu channel" << std::endl;
				return TTchainNames;
			}
			if(mode.find("POWHEG") != _ENDSTRING) {
				TTchainNames.push_back("DYTo" + tagName + "_powheg_50to120");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_120to200");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_200to400");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_400to800");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_800to1400");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_1400to2300");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_2300to3500");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_3500to4500");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_4500to6000");
				TTchainNames.push_back("DYTo" + tagName + "_powheg_6000toInf");
			} else if(mode.find("AMC") != _ENDSTRING) {
				//amc at nlo inclusive sample gen dilepton mass greater than 50 GeV
				TTchainNames.push_back("DYJets_amctnlo");
			} else if(mode.find("MAD") != _ENDSTRING) {
				//madgraph inclusive sample gen dilepton mass greater than 50 GeV
				TTchainNames.push_back("DYJets_madgraph");
			} else if(mode.find("POWINCL") != _ENDSTRING && channel == Selector::EE) {
				TTchainNames.push_back("DYToEE_powheg");
			} else if(mode.find("MADHT") != _ENDSTRING) {
				TTchainNames.push_back("DYJets_madgraph_ht100to200");
				TTchainNames.push_back("DYJets_madgraph_ht200to400");
				TTchainNames.push_back("DYJets_madgraph_ht400to600");
				TTchainNames.push_back("DYJets_madgraph_ht600toInf");
			}
		} else if(mode == "W") {
			TTchainNames.push_back("WJetsLNu");
		} else if(mode == "WZ") {
			TTchainNames.push_back("WZ");
		} else if(mode == "ZZ") {
			TTchainNames.push_back("ZZ");
		} else if(mode == "data") {
			std::string dataTag = "";
			if(channel == Selector::EMu)  dataTag = "MuEG";
			if(channel == Selector::EE)   dataTag = "DoubleEG";
			if(channel == Selector::MuMu) dataTag = "SingleMu";
			TTchainNames.push_back(dataTag + "_RunC");
			TTchainNames.push_back(dataTag + "_RunD_v3");
			TTchainNames.push_back(dataTag + "_RunD_v4");
		}
		if(mode.find("WRto") != _ENDSTRING) {
			TTchainNames.push_back(mode);
		}
		return TTchainNames;
	};


	std::string getTreeName(Selector::tag_t channel, bool isTagAndProbe, bool isLowDiLepton)
	{
		std::string tree_channel = "";

		// Select the channel to be studied //
		if(isLowDiLepton && channel != Selector::EMu)
			tree_channel = "_lowdileptonsideband";
		else if(isTagAndProbe)
			tree_channel = "_dytagandprobe";
		else if(channel == Selector::EE)
			tree_channel = "_signal_ee";
		else if(channel == Selector::MuMu)
			tree_channel = "_signal_mumu";
		else if(channel == Selector::EMu)
			tree_channel = "_flavoursideband";
		else {
			tree_channel = "";
			std::cerr << "[ERROR] No channel defined" << std::endl;
		}

		return  tree_channel;
	};

	bool checkValidMode(std::string mode)
	{
		if(mode.find("WRto") != _ENDSTRING) {
			return true;
		}

		if(all_modes.count(mode) == 0) {
			std::cerr << "[ERROR] Mode " << mode << " not part of the standard modes:" << std::endl;
			for(auto allowed_mode : all_modes) std::cerr << "        " << allowed_mode << std::endl;
			return false;
		}
		return true;
	};

private:
	std::unordered_set<std::string> all_modes;

};




int main(int ac, char* av[])
{
	namespace po = boost::program_options;

	std::vector<std::string> modes;
	chainNames chainNames_;

	std::string channel_str, outDir, outFileTag;
	float integratedLumi;
	Int_t nToys;
	bool debug;
	bool isTagAndProbe, isLowDiLepton, saveToys, ignoreDyScaleFactors;
	int nStatToys;
	int signalN;
	int seed;
	// Declare the supported options.
	po::options_description required("Mandatory command line options");
	required.add_options()
	("mode,m", po::value<std::vector<std::string> >(&modes)->required(), "Set mode to use:\n")
	("channel,c", po::value<std::string>(&channel_str)->required(), "Set Channel (EE, MuMu, EMu)")
	;

	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("lumi,l", po::value<float>(&integratedLumi)->default_value(2640.523267), "Integrated luminosity")
	("toys,t", po::value<int>(&nToys)->default_value(0), "Number of Toys")
	("seed,s", po::value<int>(&seed)->default_value(123456), "Starting seed")
	("saveToys", po::bool_switch(&saveToys)->default_value(false), "Save t1 tree vector for every toy iteration")
	("outputDir,d", po::value<std::string>(&outDir)->default_value(""), "output dir for file with plotting trees")
	("outputFileTag,f", po::value<std::string>(&outFileTag)->default_value(""), "tag name added to output file with plotting trees")
	("ignoreDyScaleFactors", po::value<bool>(&ignoreDyScaleFactors)->default_value(true), "Ignore DyScaleFactors defined in configs directory")
	("verbose,v", po::bool_switch(&debug)->default_value(false), "Turn on debug statements")
	("isTagAndProbe", po::bool_switch(&isTagAndProbe)->default_value(false), "use the tag&probe tree variants")
	("isLowDiLepton", po::bool_switch(&isLowDiLepton)->default_value(false), "low di-lepton sideband")
	("nStatToys", po::value<int>(&nStatToys)->default_value(0), "throw N toys for stat uncertainty.")
	("signalN", po::value<int>(&signalN)->default_value(0), "throw N toys for stat uncertainty.")
	;

	po::variables_map vm;
	po::options_description all("all");
	all.add(desc).add(required);
	po::store(po::parse_command_line(ac, av, all), vm);


	if (vm.count("help")) {
		std::cout << required << "\n";
		std::cout << desc << "\n";
		return 1;
	}

	try {
		po::notify(vm);
	} catch(const po::required_option & e) {
		std::cerr << "[ERROR] "  << e.what() << std::endl;
		std::cerr << desc << std::endl;
		return 1;
	}


	EnergyScaleCorrection_class eSmearer("Calibration/ZFitter/data/scales_smearings/76X_16DecRereco_2015_Etunc");

	//------------------------------ check if modes given in the command line are allowed
	for(auto s : modes ) {
		if(chainNames_.checkValidMode(s) == false) return 1;
	}


	//------------------------------ translate the channel option into the selector type
	Selector::tag_t channel;
	if(channel_str == "EE")
		channel = Selector::EE;
	else if(channel_str == "MuMu")
		channel = Selector::MuMu;
	else if(channel_str == "EMu")
		channel = Selector::EMu;
	else {
		std::cerr << "[ERROR] Channel " << channel_str << " not recognized" << std::endl;
		std::cerr << desc << std::endl;
		return 1;
	}

	configReader myReader("configs/2015-v1.conf");
	if(debug) std::cout << myReader << std::endl;

	std::cout << "[INFO] Selected modes: \n";
	unsigned int msize = modes.size();
	modes.erase( std::remove( modes.begin(), modes.end(), "signal" ), modes.end() );
	if(modes.size() != msize) {
		int di = 0;
		for(std::string datasetName : myReader.getDatasetNames()) {
			if(datasetName.find("WRto" + channel_str + "JJ_") != _ENDSTRING) {
				di++;
				if(!signalN || signalN == di)
					modes.push_back(datasetName);
			}
		}
	}
	for(auto s : modes) {
		std::cout << "       - " << s << "\n";
	}
	std::cout << std::endl;

	char name[100];// to name the tree per iteration of the Toy
	using namespace RooFit;
	std::cout << "******************************* Analysis ******************************" << std::endl;
	std::cout << "[WARNING] no weights associated to jets yet" << std::endl;

	myReader.setupDyMllScaleFactor("configs/dyScaleFactors.txt");


	std::map<int, std::pair<int, int> > mass_cut = getMassCutMap();
	std::vector<int> mass_vec = getMassVec();

	std::string treeName = "miniTree" + chainNames_.getTreeName(channel, isTagAndProbe, isLowDiLepton);
	unsigned long long zMass60to120EvtCount = 0;	///<count the number of evts from each dataset with 60 < dilepton_mass < 120 which pass loose selector cuts
	unsigned long long zMass65to115EvtCount = 0;
	unsigned long long zMass70to110EvtCount = 0;
	unsigned long long zMass75to105EvtCount = 0;
	unsigned long long zMass80to100EvtCount = 0;
	unsigned long long zMass85to95EvtCount = 0;

	for(auto mode : modes) {
		bool isData = chainNames_.isData(mode);
		if(isData) eSmearer.doScale = true;
		else eSmearer.doSmearings = true;

		TChain *c = myReader.getMiniTreeChain(chainNames_.getChainNames(mode, channel, isTagAndProbe), treeName);
#ifdef DEBUG
		c->Print();
#endif
		std::cout << "[INFO] Entries: " <<  c->GetEntries() << std::endl;
		if(c->GetEntries() == 0) {
			std::cerr << "[ERROR] No entries in chain... something went wrong" << std::endl;
			return 1;
		}

		// if you want to check if the config file is read correctly:
		if(debug) std::cout << myReader.getNorm1fb("TTJets_DiLept_v1") << std::endl;

		// Plotting trees
		std::string chnlName = channel_str, usingWeights = "";
		if(!ignoreDyScaleFactors && mode.find("DY") != _ENDSTRING) usingWeights = "_withMllWeight";
		TFile f((outDir + "selected_tree_" + mode + chainNames_.getTreeName(channel, isTagAndProbe, isLowDiLepton) + chnlName + usingWeights + outFileTag + ".root").c_str(), "recreate");
		f.WriteObject(&mass_vec, "signal_mass");
		// store the fitted results for every toy in a tree
		AnalysisResult result;
		result.nmasses = mass_vec.size();

		TTree * syst_tree = new TTree("syst_tree", "");
		TTree * central_value_tree = new TTree("central_value_tree", "");

		//First loop will always be central value
		result.SetBranches(syst_tree);
		result.SetBranches(central_value_tree);

		miniTreeEvent myEvent;

		myEvent.SetBranchAddresses(c);
		Selector selEvent;

		std::vector<TTree *> t1(nToys + 1, NULL);
		TTree * tDyCheck = new TTree("treeDyCheck", "");
		ULong64_t nEntries = c->GetEntries();
#ifdef DEBUGG
		nEntries = 500;
#endif

		TRandom3 Rand;
		const int Total_Number_of_Systematics_Smear = 1;// electron scale(MC)
		const int Total_Number_of_Systematics_Up_Down = 4;// muon id, muon iso, electron scale(data) and jet energy scale
		float Random_Numbers_for_Systematics_Smear[Total_Number_of_Systematics_Smear] = {0.};
		float Random_Numbers_for_Systematics_Up_Down[Total_Number_of_Systematics_Up_Down] = {0.};


		std::vector<std::string> List_Systematics;
		//List_Systematics.push_back("smear");
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

		TStopwatch ts;

		std::vector< miniTreeEvent> myEventVector;
		ts.Start();
		std::cout << "Loading events (nEvents = " << nEntries << "): [ 0%]" << std::flush;
		unsigned long long int nEntries_100 = nEntries / 100;
		for(unsigned long long int ev = 0; ev < nEntries; ++ev) {
			if(nEntries > 100 && ev % nEntries_100 == 1) {
				std::cout << "\b\b\b\b\b[" << std::setw (2) <<  (int)(ev / nEntries_100) << "%]" << std::flush;
				std::cout << myEventVector.size() << std::endl;
			}
			c->GetEntry(ev);
			Selector sel(myEvent);
			if((isTagAndProbe == true && (myEvent.electrons_p4->size() > 1 || myEvent.muons_p4->size() > 1) ) || sel.isPassingPreselect()) {
				unsigned int nEle = myEvent.electrons_p4->size();
				for(unsigned int iEle = 0; iEle < nEle; ++iEle) {
					TLorentzVector& p4 = (*myEvent.electrons_p4)[iEle];
					if(isData) { //only scales are corrected
						(*myEvent.electron_scale)[iEle] = eSmearer.ScaleCorrection(myEvent.run, fabs(p4.Eta()) < 1.479, (*myEvent.electron_r9)[iEle], p4.Eta(), p4.Et());
						(*myEvent.electron_smearing)[iEle] = 0.;
					} else { // only the smearings are corrected
						(*myEvent.electron_scale)[iEle] = 1.;
						(*myEvent.electron_smearing)[iEle] = eSmearer.getSmearingSigma(myEvent.run, fabs(p4.Eta()) < 1.479, (*myEvent.electron_r9)[iEle], p4.Eta(), p4.Et(), EnergyScaleCorrection_class::kRho, 0);
					}
				}
				myEventVector.push_back(myEvent);
			}
		}
		nEntries = myEventVector.size();
		nEntries_100 = nEntries / 100;
		ts.Stop();
		ts.Print();

		std::cout << "[INFO] Running nToys = " << nToys << std::endl;
		bool loop_one = true;
		int seed_i = seed + 1;

		// electron systematics!
		bool Flag_Smear_Electron_Scale = false;
		for(unsigned int iii = 0; iii < List_Systematics.size(); iii++) {
			if(List_Systematics[iii] == "Smear_Electron_Scale") Flag_Smear_Electron_Scale = true;
		}


		for(int i = 0; i < nToys + 1; ++i, ++seed_i) {

			Rand.SetSeed(seed_i);
			//for central values, we take the central value of Mu ID/ISO efficiencies and dont smear for JES systematics
			// Roch and Electron scales are smeared with a pre-defined seed(1), to give conistent results.
			if(loop_one) {
				Random_Numbers_for_Systematics_Up_Down[0] = 0.;//Mu Eff ID
				Random_Numbers_for_Systematics_Up_Down[1] = 0.;//Mu Eff ISO
				Random_Numbers_for_Systematics_Up_Down[2] = 0.; //Rand.Gaus(0.0, 1.);// Electron Scale(Data)
				Random_Numbers_for_Systematics_Up_Down[3] = 0.;//JES
			} else {
				for(int Rand_Up_Down_Iter = 0; Rand_Up_Down_Iter < Total_Number_of_Systematics_Up_Down; Rand_Up_Down_Iter++)
					Random_Numbers_for_Systematics_Up_Down[Rand_Up_Down_Iter] = Rand.Gaus(0.0, 1.);
			}
			RooDataSet * tempDataSet = new RooDataSet("temp", "temp", Fits::vars);
			sprintf(name, "Tree_Iter%i", i);
			t1[i] = new TTree(name, "");
			selEvent.SetBranches(t1[i]);
			if(loop_one) selEvent.SetBranches(tDyCheck);

			ts.Stop();
			ts.Print();
			ts.Start();

			std::cout << "Processing events (nEvents = " << nEntries << "): [ 0%]" << std::flush;

			unsigned long long int ev = 0;

			//------------------------------ scale random numbers: one set of numbers per toy, common to all events and electrons
			// it's the Random_Numbers_for_Systematics_Up_Down[2]

			for(auto myEvent : myEventVector) {


				if(nEntries > 100 && ev % nEntries_100 == 1) {
//					std::cout << "Processing events (nEvents = " << nEntries << "): [ 0%]" << std::flush;
					std::cout << "\b\b\b\b\b[" << std::setw (2) <<  (int)(ev / nEntries_100) << "%]" << std::flush;
				}
				//std::cout << "\b\b\b\b" << (int)( ev/nEntries_100) << " %" << std::endl;

//#ifdef DEBUG
				if (debug) {
					std::cout << "RUN=" << myEvent.run << std::endl;
					std::cout << "Mu" << std::endl;
					for(auto m : * (myEvent.muons_p4))
						std::cout << m.Pt() << " " << m.Eta() << std::endl;
					std::cout << "Jet" << std::endl;
					for(auto m : * (myEvent.jets_p4))
						std::cout << m.Pt() << " " << m.Eta() << std::endl;
				}
//#endif
				for(int Rand_Smear_Iter = 0; Rand_Smear_Iter < Total_Number_of_Systematics_Smear; Rand_Smear_Iter++)
					Random_Numbers_for_Systematics_Smear[Rand_Smear_Iter] = Rand.Gaus(0.0, 1.);
				ToyThrower( &myEvent, Random_Numbers_for_Systematics_Smear, Random_Numbers_for_Systematics_Up_Down, seed_i, List_Systematics, isData);


				unsigned int nEle = myEvent.electrons_p4->size();
				for(unsigned int iEle = 0; iEle < nEle; ++iEle) {
					TLorentzVector& p4 = (*myEvent.electrons_p4)[iEle];

					if(isData) { //only scales are corrected
#ifdef DEBUG
						std::cout << p4.Et() << "-->";
#endif
						p4 *= (*myEvent.electron_scale)[iEle];
						// since the Random_Numbers_for_Systematics is the same for all the categories, then the uncertainties are considered fully correlated
						if(loop_one!=false && Flag_Smear_Electron_Scale) p4 *= 1 + Random_Numbers_for_Systematics_Up_Down[2] * eSmearer.ScaleCorrectionUncertainty(myEvent.run, fabs(p4.Eta()) < 1.479, 0., p4.Eta(), p4.Et());
#ifdef DEBUG
						std::cout << p4.Et() << "\t" << p4.Eta() << "\t" << (*myEvent.electron_r9)[iEle] << "\t" << (*myEvent.electron_scale)[iEle] << "\t" << (eSmearer.ScaleCorrectionUncertainty(myEvent.run, fabs(p4.Eta()) < 1.479, 0., p4.Eta(), p4.Et())) << std::endl;
#endif
					} else { // only the smearings are corrected
						double r = (loop_one==false && Flag_Smear_Electron_Scale) ? Rand.Gaus(0., 1.) : 0;
						(*myEvent.electron_smearing)[iEle] = eSmearer.getSmearingSigma(myEvent.run, fabs(p4.Eta()) < 1.479, 0., p4.Eta(), p4.Et(), EnergyScaleCorrection_class::kRho, r);
#ifdef DEBUG
						std::cout << p4.Et() << "-->";
#endif
						p4 *= Rand.Gaus(1., (*myEvent.electron_smearing)[iEle]);
#ifdef DEBUG
						std::cout << p4.Et() << "\t" << p4.Eta() << "\t" << (*myEvent.electron_r9)[iEle] << "\t" <<	(*myEvent.electron_smearing)[iEle] << "\t" <<  eSmearer.getSmearingSigma(myEvent.run, fabs(p4.Eta()) < 1.479, 0., p4.Eta(), p4.Et(), EnergyScaleCorrection_class::kRho, 0) << "\t" << eSmearer.getSmearingSigma(myEvent.run, fabs(p4.Eta()) < 1.479, 0., p4.Eta(), p4.Et(), EnergyScaleCorrection_class::kRho, r) << std::endl;
#endif

					}
				}

				Selector tmp_selEvent(myEvent);
				selEvent = tmp_selEvent;
				// Select events with one good WR candidate
				// Tags:
				// 0 -- EEJJ Channel
				// 1 -- MuMuJJ Channel
				// 2 -- EMuJJ Channel

				// if(debug && selEvent.isPassing(channel)) {
				// 	std::cout << std::endl << selEvent.dilepton_mass << std::endl;
				// 	std::cout << "Mu " << ev  << std::endl;
				// 	for(auto m : * (myEvent.muons_p4))
				// 		std::cout << m.Pt() << " " << m.Eta() << std::endl;
				// 	std::cout << "Ele" << std::endl;
				// 	for(auto m : * (myEvent.electrons_p4))
				// 		std::cout << m.Pt() << " " << m.Eta() << std::endl;
				// 	std::cout << "Jet" << std::endl;
				// 	for(auto m : * (myEvent.jets_p4))
				// 		std::cout << m.Pt() << " " << m.Eta() << std::endl;
				// }

				if(loop_one && selEvent.isPassingLooseCuts(channel)) {
					if(isData == false) {
						selEvent.weight *= myReader.getNorm1fb(selEvent.datasetName) * integratedLumi; // the weight is the event weight * single object weights

						//multiply by an additional weight when processing DY samples
						if(mode.find("DY") != _ENDSTRING && !ignoreDyScaleFactors) {
							selEvent.weight *= myReader.getDyMllScaleFactor(channel_str, mode);
						}
					} else {
						selEvent.weight = 1;
#ifdef DEBUGG
						std::cout << "selEvent.weight=\t" << selEvent.weight << std::endl;
#endif
						assert(selEvent.weight == 1);
					}
#ifdef DEBUGG
					std::cout << "selEvent.weight=\t" << selEvent.weight << std::endl;
					std::cout << "integratedLumi=\t" << integratedLumi << std::endl;
					std::cout << "myReader.getNorm1fb(selEvent.datasetName)=\t" << myReader.getNorm1fb(selEvent.datasetName) << std::endl;
					std::cout << "myEvent.weight=\t" << myEvent.weight << std::endl;
#endif

					if(selEvent.dilepton_mass > 61.2 && selEvent.dilepton_mass < 121.2 && loop_one) ++zMass60to120EvtCount;
					if(selEvent.dilepton_mass > 66.2 && selEvent.dilepton_mass < 116.2 && loop_one) ++zMass65to115EvtCount;
					if(selEvent.dilepton_mass > 71.2 && selEvent.dilepton_mass < 111.2 && loop_one) ++zMass70to110EvtCount;
					if(selEvent.dilepton_mass > 76.2 && selEvent.dilepton_mass < 106.2 && loop_one) ++zMass75to105EvtCount;
					if(selEvent.dilepton_mass > 81.2 && selEvent.dilepton_mass < 101.2 && loop_one) ++zMass80to100EvtCount;
					if(selEvent.dilepton_mass > 86.2 && selEvent.dilepton_mass < 96.2 && loop_one)  ++zMass85to95EvtCount ;
					if(loop_one) tDyCheck->Fill();
				}

				if(selEvent.isPassing(channel)) {

					if (channel == Selector::EMu && selEvent.dilepton_mass < 200) continue;


					if(isData == false) {
						selEvent.weight *= myReader.getNorm1fb(selEvent.datasetName) * integratedLumi; // the weight is the event weight * single object weights

						//multiply by an additional weight when processing DY samples
						if(mode.find("DY") != _ENDSTRING && !ignoreDyScaleFactors) {
							selEvent.weight *= myReader.getDyMllScaleFactor(channel_str, mode);
						}
					} else {
						selEvent.weight = 1;
						assert(selEvent.weight == 1);
					}

					Fits::massWR.setVal(selEvent.WR_mass);
					Fits::evtWeight.setVal(selEvent.weight);
					tempDataSet->add(Fits::vars);

					t1[i]->Fill();

				}
				++ev;
			}//end loop over all input evts, and adding events to the RooDataSet pointer named tempDataSet
			ts.Stop();
			ts.Print();
			ts.Start();

			if(loop_one) std::cout << zMass60to120EvtCount << "\tevents from the dataset named\t" << selEvent.datasetName << "\tpass isPassingLooseCuts and have 61.2 < dilepton_mass < 121.2" << std::endl;
			if(loop_one) std::cout << zMass65to115EvtCount << "\tevents from the dataset named\t" << selEvent.datasetName << "\tpass isPassingLooseCuts and have 66.2 < dilepton_mass < 116.2" << std::endl;
			if(loop_one) std::cout << zMass70to110EvtCount << "\tevents from the dataset named\t" << selEvent.datasetName << "\tpass isPassingLooseCuts and have 71.2 < dilepton_mass < 111.2" << std::endl;
			if(loop_one) std::cout << zMass75to105EvtCount << "\tevents from the dataset named\t" << selEvent.datasetName << "\tpass isPassingLooseCuts and have 76.2 < dilepton_mass < 106.2" << std::endl;
			if(loop_one) std::cout << zMass80to100EvtCount << "\tevents from the dataset named\t" << selEvent.datasetName << "\tpass isPassingLooseCuts and have 81.2 < dilepton_mass < 101.2" << std::endl;
			if(loop_one) std::cout << zMass85to95EvtCount << "\tevents from the dataset named\t" << selEvent.datasetName << "\tpass isPassingLooseCuts and have 86.2 < dilepton_mass < 96.2" << std::endl;

			///make a permanent RooDataSet which has the same information as tempDataSet, but with events which are weighted according to the var named evtWeight
			RooDataSet * permanentWeightedDataSet = new RooDataSet("permanentWeightedDataSet", "permanentWeightedDataSet", tempDataSet, Fits::vars, "", Fits::evtWeight.GetName());
			// Count number of events in each mass range to store in tree.
			TH1F * hWR_mass = new TH1F("hWR_mass", "hWR_mass", 140, 0, 7000);
			t1[i]->Draw("WR_mass>>hWR_mass", "weight", "goff");
			double error = 0;
			for(size_t mass_i = 0; mass_i < mass_vec.size(); mass_i++) {
				auto range = mass_cut[mass_vec.at(mass_i)];
				result.events_in_range[mass_i] = hWR_mass->IntegralAndError(hWR_mass->FindBin(range.first), hWR_mass->FindBin(range.second), error);
				result.error_in_range[mass_i] = float(error);
			}
			delete hWR_mass;

			if(saveToys) t1[i]->Write();
			if(loop_one) {
				if(!saveToys) t1[i]->Write();
				permanentWeightedDataSet->Write();
				tDyCheck->Write();
				delete tDyCheck;
			}
			delete t1[i];

			permanentWeightedDataSet->Print();

			if((isTagAndProbe == false) && (mode == "TT" || mode.find("DY") != _ENDSTRING || (mode == "data" && channel == Selector::EMu) ) ) {
				assert(permanentWeightedDataSet->sumEntries() > 0);
				Fits::expPower.setVal(-0.004);
				RooFitResult * tempFitRslt = NULL;
				// fit dataset to given PDF
				fitRooDataSet(tempFitRslt, permanentWeightedDataSet, Fits::expPdfRooAbsPdf);	///< expPdfRooAbsPdf defined over massWR 600 to 6500

				// dataset normalization is the number of entries in the dataset with fourObjectMass (name of Fits::massWR) above 600
				result.normalization = permanentWeightedDataSet->sumEntries("fourObjectMass > 600");
				// set of variables in the PDF
				RooArgSet *vset = Fits::expPdfRooAbsPdf->getVariables();
				// loop over RooRealVars in the set
				TIterator * iter = vset->createIterator();
				TObject * var = iter->Next();
				RooRealVar *var_pdf;
				result.nparam = 0;
				while (var) {
					// ignore the M_WR variable
					if(strcmp(var->GetName(), "fourObjectMass") != 0) {
						var_pdf = (RooRealVar*)vset->find(var->GetName());
						// store the value of the fitted parameters and the corresponding errors
						result.fit_parameters[result.nparam] = var_pdf->getVal();
						result.fit_parameter_errors[result.nparam++] = var_pdf->getError();
					}
					var = iter->Next();
				}

				for(size_t mass_i = 0; mass_i < mass_vec.size(); mass_i++) {
					auto range = mass_cut[mass_vec.at(mass_i)];
					double integral =  NormalizedIntegral(Fits::expPdfRooAbsPdf, Fits::massWR, range.first, range.second);
					result.fit_integral_in_range[mass_i] = integral;
				}

				if(nStatToys && loop_one) {
					std::cout << "Doing " << nStatToys << " stat Toys" << std::endl;
					AnalysisResult stat_result;
					stat_result.normalization = result.normalization;
					stat_result.nmasses = mass_vec.size();
					TTree * stat_tree = new TTree("stat_tree", "");
					stat_result.SetBranches(stat_tree);

					for(int stat_i = 0; stat_i < nStatToys; stat_i++) {
						// loop over RooRealVars in the set, vary them and store
						TIterator * iter = vset->createIterator();
						TObject * var = iter->Next();
						RooRealVar *var_pdf;
						stat_result.nparam = 0;
						while (var) {
							// ignore the M_WR variable
							if(strcmp(var->GetName(), "fourObjectMass") != 0) {
								var_pdf = (RooRealVar*)vset->find(var->GetName());
								var_pdf->setVal(Rand.Gaus(var_pdf->getVal(), var_pdf->getError()));
								stat_result.fit_parameters[stat_result.nparam++] = var_pdf->getVal();
							}
							var = iter->Next();
						}
						//Calculate integrals for each fit
						for(size_t mass_i = 0; mass_i < mass_vec.size(); mass_i++) {
							auto range = mass_cut[mass_vec.at(mass_i)];
							double integral =  NormalizedIntegral(Fits::expPdfRooAbsPdf, Fits::massWR, range.first, range.second);
							stat_result.fit_integral_in_range[mass_i] = integral;
						}

						stat_tree->Fill();

					}
					stat_tree->Write();
					delete stat_tree;

				}
			}

			// fill the tree with the normalization, parameters, and errors
			if(loop_one) central_value_tree->Fill();
			else syst_tree->Fill();

			loop_one = false;
		}
		syst_tree->Write();
		central_value_tree->Write();
	}
	return 0;

}
