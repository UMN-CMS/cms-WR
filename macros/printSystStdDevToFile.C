///ROOT includes
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TBranch.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPad.h>

///C++ includes
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

/**
 * read the syst_tree TTree from three sets of root files with user defined path names.  Each set of root files corresponds to signal and bkgnd MC which have been processed
 * by analysis.cpp with one systematic enabled (for which toys are thrown by analysis.cpp and ToyThrower).  One set has Smear_Jet_Scale enabled, another set has
 * Smear_Electron_Scale enabled, and the other set has both Smear_Muon_Scale and Smear_Muon_ID_Iso enabled.
 *
 * Each set will have four files: one low mass WR signal, one medium mass WR signal, one high mass WR signal, and the DY root file.  The ttbar MC and EMu data root files are
 * not needed.
 *
 * This macro looks at all sets of root files, and prints the syst uncertainty due to the quantity which was varied.  For each WR signal file, only one lepton syst uncertainty is
 * printed, corresponding to the optimized mass window of that particular WR signal.  For each DY file, three lepton uncertainties will be printed - one in a low mass window, one
 * in a medium mass window, and one in a high mass window.  For each DY file, three jet uncertainties will be printed.
 * No plots are made by this macro, just a txt file with numbers.
 *
 * analysisTools.h defined in the interface directory is not used here because it does not provide a quick way to match both the lepton channel and WR mass used to define the mass window.
 * Instead of using analysisTools.h, the mass_cuts.txt file is read directly, and the WR mass is specified by the user.
 *
 */

///given the location of the mass_cuts.txt file, and the WR mass and lepton channel as a string, this function returns the correct index
///of NEventsInRange to use for plotting.  When NEventsInRange is plotted with this index, the RMS of the distribution
///is equal to the systematic uncertainty.
int getIndexForSystematicUncertainty(string wrMass, string leptonChannel, string pathToMassCutsFile){
	int indexToReturn = 0;
	ifstream massWindowReader(pathToMassCutsFile.c_str());
	while(massWindowReader.peek() != EOF && massWindowReader.good() ){
		if(massWindowReader.peek() == 35){
			massWindowReader.ignore(1000, '\n');
			continue;
		}//end discard of lines which start with a number sign

		string channelFromFile, wrMassForMassWindow, massWindowLowerBound, massWindowUpperBound;
		massWindowReader >> channelFromFile >> wrMassForMassWindow >> massWindowLowerBound >> massWindowUpperBound;
		if(channelFromFile == leptonChannel && wrMassForMassWindow == wrMass) break;	///<leave the loop once a match is found

		indexToReturn++;
	}//end loop over lines in mass cuts file
	massWindowReader.close();

	return indexToReturn;
}//end getIndexForSystematicUncertainty()

void printSystStdDevToFile(){

	///user defined path to txt file which lists WR mass window ranges
	string pathToMassWindowsFile = "../configs/mass_cuts.txt";

	///user defined low, medium, and high WR mass points
	///make sure each mass point is listed in the mass cuts file
	int lowMass = 800, mediumMass = 2200, highMass = 4000;
	int wrMassPoints[] = {800,2200,4000};
	vector<int> wrMassVect(wrMassPoints, wrMassPoints + sizeof(wrMassPoints)/sizeof(int) );

	///user defined paths to root file dirs     the combination absPath + relPath must be an existing directory
	///and other strings related to TTree file access
	string absPathToMainRootFileDir = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/";
	string relDirPaths[] = {"800toysOnlyJetScaleSyst/", "800toysOnlySmearEleScaleSyst/", "410toysOnlyMuScaleIdSyst/"};
	vector<string> relDirPathsVector(relDirPaths, relDirPaths + sizeof(relDirPaths)/sizeof(string) );
	string uncertTagNames[] = {"jet energy scale smearing","electron energy scale smearing","muon energy scale smearing and ID Iso smearing"};
	string treeName = "syst_tree", eeChannelInFileNames = "eeEE", mumuChannelInFileNames = "mumuMuMu";
	string dyFileName = "selected_tree_DYAMC_signal_";
	//string wrMuMuLowMassFileName = "selected_tree_WRtoMuMuJJ_" + to_string(lowMass) + "_" + to_string(lowMass/2) + "_signal_" + mumuChannelInFileNames;

	///user defined path to output txt file created by this macro, which lists systematic uncertainties for different processes
	///in both lepton channels, and several mass windows
	string pathToSystUncOutputFile = "systematicUncertaintiesFromAnalysisCpp.txt";
	ofstream writeToSystFile(pathToSystUncOutputFile.c_str(),ofstream::trunc);
	
	///now that all the user defined vars have been declared, calculate the systematic uncertainty for the specified mass windows for WR and DY in both lepton channels
	int numWrMasses = wrMassVect.size(), numSystematics = relDirPathsVector.size();
	for(int m=0; m<numWrMasses; m++){
		for(int s=0; s<numSystematics; s++){
			//for each WR mass and source of uncertainty, determine the systematic uncertainty and write it to a txt file
			TChain * dyMuMu = new TChain(treeName.c_str());
			dyMuMu->Add( (absPathToMainRootFileDir +  ).c_str() );

		}//end loop over different systematics to evaluate

	}//end loop over different WR mass points
	



}//end printSystStdDevToFile()
