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
#include <array>
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

//#define DEBUG

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
int getIndexForSystematicUncertainty(string wrMass, string inputLeptonChannel, string pathToMassCutsFile){
	int indexToReturn = 0;
	ifstream massWindowReader(pathToMassCutsFile.c_str());
	while(massWindowReader.peek() != EOF && massWindowReader.good() ){
		if(massWindowReader.peek() == 35){
			massWindowReader.ignore(1000, '\n');
			continue;
		}//end discard of lines which start with a number sign

		string channelFromFile, wrMassForMassWindow, massWindowLowerBound, massWindowUpperBound;
		massWindowReader >> channelFromFile >> wrMassForMassWindow >> massWindowLowerBound >> massWindowUpperBound;
		string copyInputWrMass = wrMass;
		
		//WR mass read in from mass cuts file has a leading zero on 800 GeV mass, 0800
		if(wrMass == "800") copyInputWrMass = "0800";
		if(inputLeptonChannel.find(channelFromFile) != string::npos && wrMassForMassWindow == copyInputWrMass) break;	///<leave the loop once a match is found

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
	//int wrMassPoints[] = {800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3600,3800,4000,4200,4400,4600,4800,5000,5200,5600,5800,6000};
	//int wrMassPoints[] = {1000,1600,2200,2800,3600};
	int wrMassPoints[] = {1600,2200,2800,3600};
	vector<int> wrMassVect(wrMassPoints, wrMassPoints + sizeof(wrMassPoints)/sizeof(int) );

	///user defined paths to root file dirs     the combination absPath + relPath must be an existing directory
	///and other strings related to TTree file access
	string absPathToMainRootFileDir = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/";
	
	//relDirPathsVect and uncertTagNamesVect must have the same size
	string relDirPaths[] = {"3200toysOnlyJetScaleSyst/", "3200toysSmearEleScaleSyst/", "3200toysOnlyMuScaleIdSyst/","3200toysAllSyst/"};
	string uncertTagNames[] = {"jet energy scale","electron energy scale","muon energy scale and ID Iso","all sources"};
	//string relDirPaths[] = {"3200toysAllSyst/"};
	//string uncertTagNames[] = {"all sources"};
	vector<string> relDirPathsVect(relDirPaths, relDirPaths + sizeof(relDirPaths)/sizeof(string) );
	vector<string> uncertTagNamesVect(uncertTagNames, uncertTagNames + sizeof(uncertTagNames)/sizeof(string) );
	string treeName = "syst_tree", eeChannelInFileNames = "eeEE", mumuChannelInFileNames = "mumuMuMu";
	string dyFileName = "selected_tree_DYAMC_signal_";
	string emuDataFileName = "selected_tree_data_flavoursidebandEMu";

	///user defined path to output txt file created by this macro, which lists systematic uncertainties for different processes
	///in both lepton channels, and several mass windows
	string pathToSystUncOutputFile = "systUncertaintiesFromAnalysisCpp.txt";
	ofstream writeToSystFile(pathToSystUncOutputFile.c_str(),ofstream::trunc);
	
	writeToSystFile<<"#WR mass\tprocess channel\tsyst source\tmean\tsyst uncert\tstat uncert\tsyst plus stat uncert percent\tsyst sqd plus stat sqd over nEvts\tmean sqd over syst sqd plus stat sqd"<< endl;
	
	///now that all the user defined vars have been declared, calculate the systematic uncertainty for the specified mass windows for WR and DY in both lepton channels
	int numWrMasses = wrMassVect.size(), numSystematics = relDirPathsVect.size();
	if( numSystematics != uncertTagNamesVect.size() ){
		cout<<"relDirPathsVect and uncertTagNamesVect do not have the same number of elements"<<endl;
		exit(1);
	}
	for(int m=0; m<numWrMasses; m++){
		for(int s=0; s<numSystematics; s++){
			//for each WR mass and source of uncertainty, determine the systematic uncertainty and write it to a txt file
			map<string,TChain*> chainPointers;
			if(uncertTagNamesVect[s].find("jet") != string::npos){
				TChain * dyMuMu = new TChain(treeName.c_str(),"DYMuMu");
				dyMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + mumuChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttMuMu = new TChain(treeName.c_str(),"TTMuMu");
				ttMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +".root" ).c_str() );
				TChain * wrMuMu = new TChain(treeName.c_str(),"WRMuMu");
				wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
				TChain * dyEE = new TChain(treeName.c_str(),"DYEE");
				dyEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + eeChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttEE = new TChain(treeName.c_str(),"TTEE");
				ttEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"Copy.root" ).c_str() );
				TChain * wrEE = new TChain(treeName.c_str(),"WREE");
				wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
				//chainPointers["MuMuDY"] = dyMuMu;
				//chainPointers["MuMuWR"] = wrMuMu;
				//chainPointers["EEDY"] = dyEE;
				//chainPointers["EEWR"] = wrEE;
				chainPointers["EETT"] = ttEE;
				chainPointers["MuMuTT"] = ttMuMu;
	
			}

			if(uncertTagNamesVect[s].find("electron") != string::npos){
				TChain * dyEE = new TChain(treeName.c_str(),"DYEE");
				dyEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + eeChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttEE = new TChain(treeName.c_str(),"TTEE");
				ttEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"Copy.root" ).c_str() );
				TChain * wrEE = new TChain(treeName.c_str(),"WREE");
				wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
				//chainPointers["EEDY"] = dyEE;
				//chainPointers["EEWR"] = wrEE;
				chainPointers["EETT"] = ttEE;
			}

			if(uncertTagNamesVect[s].find("muon") != string::npos){
				TChain * dyMuMu = new TChain(treeName.c_str(),"DYMuMu");
				dyMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + mumuChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttMuMu = new TChain(treeName.c_str(),"TTMuMu");
				ttMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +".root" ).c_str() );
				TChain * wrMuMu = new TChain(treeName.c_str(),"WRMuMu");
				wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
				//chainPointers["MuMuDY"] = dyMuMu;
				//chainPointers["MuMuWR"] = wrMuMu;
				chainPointers["MuMuTT"] = ttMuMu;
			}

			else{
				TChain * dyMuMu = new TChain(treeName.c_str(),"DYMuMu");
				dyMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + mumuChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttMuMu = new TChain(treeName.c_str(),"TTMuMu");
				ttMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +".root" ).c_str() );
	
				TChain * wrMuMu = new TChain(treeName.c_str(),"WRMuMu");
				wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
				TChain * dyEE = new TChain(treeName.c_str(),"DYEE");
				dyEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + eeChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttEE = new TChain(treeName.c_str(),"TTEE");
				ttEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"Copy.root" ).c_str() );
	
				TChain * wrEE = new TChain(treeName.c_str(),"WREE");
				wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
				//chainPointers["MuMuDY"] = dyMuMu;
				//chainPointers["MuMuWR"] = wrMuMu;
				//chainPointers["EEDY"] = dyEE;
				//chainPointers["EEWR"] = wrEE;
				chainPointers["EETT"] = ttEE;
				chainPointers["MuMuTT"] = ttMuMu;
			}

#ifdef DEBUG
			cout<<"#WR mass\tprocess channel\tsyst source\tmean\tsyst uncert\t stat uncert\tsyst and stat uncert percent\tsyst sqd plus stat sqd over nEvts\tmean sqd over syst sqd plus stat sqd"<< endl;
#endif

			for(map<string,TChain*>::const_iterator chMapIt=chainPointers.begin(); chMapIt!=chainPointers.end(); chMapIt++){
				//get the integer corresponding to the correct index position in NEventsInRange
				string nEventsIndex = to_string( getIndexForSystematicUncertainty( to_string(wrMassVect[m]), chMapIt->first, pathToMassWindowsFile) );
				string drawArg = "NEventsInRange["+nEventsIndex+"]>>";
				string histName = (chMapIt->first) + "tempHist";
				(chMapIt->second)->Draw( (drawArg + histName + "()").c_str() );
				TH1F* tempHist = (TH1F*) gROOT->FindObject(histName.c_str());
				Double_t stdDev = tempHist->GetRMS(), stdDevErr = tempHist->GetRMSError(), mean = tempHist->GetMean();
				string statUncDrawArg = "ErrorEventsInRange["+nEventsIndex+"]>>";
				string statUncHistName = (chMapIt->first) + "StatTempHist";
				(chMapIt->second)->Draw( (statUncDrawArg + statUncHistName + "()").c_str() );
				TH1F* statUncTempHist = (TH1F*) gROOT->FindObject(statUncHistName.c_str());
				Double_t statUnc = statUncTempHist->GetMean();
#ifdef DEBUG
				cout<< wrMassVect[m] << "\t"<< chMapIt->first << "\t" << uncertTagNamesVect[s] << "\t"<< mean << "\t" << stdDev << "\t"<< statUnc <<"\t"<< 100*sqrt(stdDev*stdDev + statUnc*statUnc)/mean << endl;
				cout<<"index used in NEventsInRange and ErrorEventsInRange=\t"<< nEventsIndex <<endl;
				cout<<"\t"<<endl;
#endif
				if((chMapIt->first).find("TT") != string::npos){
					//rescale mean, stdDev and statUnc by the emu data scale factor
					Double_t rescale = ((chMapIt->first).find("EE") != string::npos) ? 0.414 : 0.657;
					mean *= rescale;
					stdDev *= rescale;
					statUnc *= rescale;
				}
				writeToSystFile << wrMassVect[m] << "\t"<< chMapIt->first << "\t" << uncertTagNamesVect[s] << "\t" << mean << "\t" << stdDev << "\t"<< statUnc <<"\t"<< 100*sqrt(stdDev*stdDev + statUnc*statUnc)/mean << "\t" << (stdDev*stdDev + statUnc*statUnc)/mean << "\t" << mean*mean/(stdDev*stdDev + statUnc*statUnc) << endl;
			}//end loop over TChains

			//loop over chainPointers elements, and clear allocated memory
			for(map<string,TChain*>::const_iterator mapIt=chainPointers.begin(); mapIt!=chainPointers.end(); mapIt++){ (mapIt->second)->Delete();}
			chainPointers.clear();

		}//end loop over different systematics to evaluate

	}//end loop over different WR mass points
	
	writeToSystFile.close();
	

}//end printSystStdDevToFile()
