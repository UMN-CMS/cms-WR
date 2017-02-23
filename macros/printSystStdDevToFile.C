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

//#define useDYMAD
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

///this method is designed to read evts from one index of NEventsInRange from two TChains for many toys, divide
///the evts for each toy, then make a distribution of the ratio vs the number of toys (NOT vs LLJJ mass)
///divide entries from chainOne by entries from chainTwo
void calcAndDrawRatio(TChain * chainOne, TChain * chainTwo, string userIndex, string tagOne, TString outFileName, TString plotTitle, TString xAxisTitle, TString yAxisTitle, Float_t minX, Float_t maxX){

	string aliasTwo = "chTwo";
	///declare chainTwo as a friend of chainOne
	chainOne->AddFriend(chainTwo, aliasTwo.c_str());

	///first draw weighted evts vs toys from both chains
	string drawArg = "(NEventsInRange["+userIndex+"]/"+aliasTwo+".NEventsInRange["+userIndex+"])>>";
	//cout<<"drawArg= "<< drawArg << endl;
	string histName = tagOne + "tempHist";
	chainOne->Draw( (drawArg + histName + "(1000,0.,1.)").c_str() );
	//chainOne->Draw( ("(NEventsInRange[0]/chTwo.NEventsInRange[0])>>" + histName + "(500,0.,1.)").c_str() );
	TH1F* tempHist = (TH1F*) gROOT->FindObject(histName.c_str());

	//now draw ratio of the two histos
	gStyle->SetOptStat("emrou");
	TCanvas * cRatio = new TCanvas("cRatio","",800,800);
	cRatio->cd();
	tempHist->SetLineColor(kBlack);
	tempHist->SetLineWidth(3);
	tempHist->SetTitle(plotTitle);
	tempHist->GetXaxis()->SetRangeUser(minX, maxX);
	tempHist->GetYaxis()->SetTitle(yAxisTitle);
	tempHist->GetYaxis()->SetTitleOffset(1.3);
	tempHist->GetXaxis()->SetTitle(xAxisTitle);
	tempHist->Draw();
	cRatio->Print((outFileName+".pdf").Data());
	cRatio->Print((outFileName+".png").Data());
	cRatio->Print((outFileName+".C").Data());
	cRatio->Close();

}//end calcAndDrawRatio()


void getExpEvtsUnwgtEvtsAndUncsSqdUserIndex(Double_t & expEvts, Double_t & unwgtEvts, Double_t & statUncSqd, Double_t & systUncSqd, string processAndChannel, TChain * chain, string userIndex){
	
	string nEventsIndex = userIndex;
	
	//syst uncert
	string drawArg = "NEventsInRange["+nEventsIndex+"]>>";
	string histName = processAndChannel + "tempHist";
	chain->Draw( (drawArg + histName + "()").c_str() );
	TH1F* tempHist = (TH1F*) gROOT->FindObject(histName.c_str());
	expEvts = tempHist->GetMean();
	Double_t systUnc = tempHist->GetRMS(), stdDevErr = tempHist->GetRMSError();
	systUncSqd = systUnc*systUnc;

	//stat uncert
	string statUncDrawArg = "ErrorEventsInRange["+nEventsIndex+"]>>";
	string statUncHistName = processAndChannel + "StatTempHist";
	chain->Draw( (statUncDrawArg + statUncHistName + "()").c_str() );
	TH1F* statUncTempHist = (TH1F*) gROOT->FindObject(statUncHistName.c_str());
	Double_t statUnc = statUncTempHist->GetMean();
	statUncSqd = statUnc*statUnc;

	//num unweighted evts
	string unweightEvtsDrawArg = "UnweightedNEventsInRange["+nEventsIndex+"]>>";
	string unweightEvtsHistName = processAndChannel + "UnweightEvtTempHist";
	chain->Draw( (unweightEvtsDrawArg + unweightEvtsHistName + "()").c_str() );
	TH1F* unweightEvtsTempHist = (TH1F*) gROOT->FindObject(unweightEvtsHistName.c_str());
	unwgtEvts = unweightEvtsTempHist->GetMean();

}//end getExpEvtsUnwgtEvtsAndUncsSqdUserIndex()


void getExpEvtsUnwgtEvtsAndUncsSqd(Double_t & expEvts, Double_t & unwgtEvts, Double_t & statUncSqd, Double_t & systUncSqd, int wrMassPoint, string processAndChannel, string pathToMassWindowsDef, TChain * chain){
	
	string nEventsIndex = to_string( getIndexForSystematicUncertainty( to_string(wrMassPoint), processAndChannel, pathToMassWindowsDef) );
	
	//syst uncert
	string drawArg = "NEventsInRange["+nEventsIndex+"]>>";
	string histName = processAndChannel + "tempHist";
	chain->Draw( (drawArg + histName + "()").c_str() );
	TH1F* tempHist = (TH1F*) gROOT->FindObject(histName.c_str());
	expEvts = tempHist->GetMean();
	Double_t systUnc = tempHist->GetRMS(), stdDevErr = tempHist->GetRMSError();
	systUncSqd = systUnc*systUnc;

	//stat uncert
	string statUncDrawArg = "ErrorEventsInRange["+nEventsIndex+"]>>";
	string statUncHistName = processAndChannel + "StatTempHist";
	chain->Draw( (statUncDrawArg + statUncHistName + "()").c_str() );
	TH1F* statUncTempHist = (TH1F*) gROOT->FindObject(statUncHistName.c_str());
	Double_t statUnc = statUncTempHist->GetMean();
	statUncSqd = statUnc*statUnc;

	//num unweighted evts
	string unweightEvtsDrawArg = "UnweightedNEventsInRange["+nEventsIndex+"]>>";
	string unweightEvtsHistName = processAndChannel + "UnweightEvtTempHist";
	chain->Draw( (unweightEvtsDrawArg + unweightEvtsHistName + "()").c_str() );
	TH1F* unweightEvtsTempHist = (TH1F*) gROOT->FindObject(unweightEvtsHistName.c_str());
	unwgtEvts = unweightEvtsTempHist->GetMean();

}//end getExpEvtsUnwgtEvtsAndUncs()

void printSystStdDevToFile(){

	///user defined path to txt file which lists WR mass window ranges
	string pathToMassWindowsFile = "../configs/mass_cuts.txt";

	///user defined low, medium, and high WR mass points
	///make sure each mass point is listed in the mass cuts file
	//int wrMassPoints[] = {800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3600,3800,4000,4200,4400,4600,4800,5000,5200,5600,5800,6000};
	int wrMassPoints[] = {1600,2200,2800};
	vector<int> wrMassVect(wrMassPoints, wrMassPoints + sizeof(wrMassPoints)/sizeof(int) );

	///user defined paths to root file dirs     the combination absPath + relPath must be an existing directory
	///and other strings related to TTree file access
	string absPathToMainRootFileDir = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/";
	
	//relDirPathsVect and uncertTagNamesVect must have the same size
	string relDirPaths[] = {"3200toysAllSystSmoothedWindowsFebrTwentyOne/", "3200toysJetSystNewMassWindows/","3200toysAllLeptSystNewMassWindows/"};
	string uncertTagNames[] = {"all sources","jet syst","mixed lepton"};
	//string relDirPaths[] = {"3200toysAllSystSmoothedWindowsFebrTwentyOne/"};
	//string uncertTagNames[] = {"all sources"};
	vector<string> relDirPathsVect(relDirPaths, relDirPaths + sizeof(relDirPaths)/sizeof(string) );
	vector<string> uncertTagNamesVect(uncertTagNames, uncertTagNames + sizeof(uncertTagNames)/sizeof(string) );
	string treeName = "syst_tree", eeChannelInFileNames = "eeEE", mumuChannelInFileNames = "mumuMuMu";
	string dyFileName = "selected_tree_DYMadInclAndHT_signal_", dyFileNameTwo = "";
	string dyEETagName="EEDYMadInclAndHT", dyEETwoTagName="";
	string dyMuMuTagName="MuMuDYMadInclAndHT", dyMuMuTwoTagName="";
#ifdef useDYMAD
	dyFileName = "selected_tree_DYMadInclAndHT_signal_";
	dyFileNameTwo = "selected_tree_DYMAD_signal_";
	dyEETagName="EEDYHT", dyEETwoTagName="EEDYIncl";
	dyMuMuTagName="MuMuDYHT", dyMuMuTwoTagName="MuMuDYIncl";
#endif
	
	string emuDataFileName = "selected_tree_data_flavoursidebandEMu";

	///user defined path to output txt file created by this macro, which lists systematic uncertainties for different processes
	///in both lepton channels, and several mass windows
	string pathToSystUncOutputFile = "systUncertaintiesFromAnalysisCpp.txt";
	ofstream writeToSystFile(pathToSystUncOutputFile.c_str(),ofstream::trunc);
	
	//writeToSystFile<<"#WR mass\tprocess channel\tsyst source\tmean\tsyst sqd plus stat sqd over nEvts\tmean sqd over syst sqd plus stat sqd"<< endl;
	writeToSystFile<<"#WR mass\tprocess channel\tsyst source\tmean\tstat uncert\tsyst uncert\tunweighted evts"<< endl;
	
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
			if(uncertTagNamesVect[s].find("mixed") != string::npos){
				TChain * dyEE = new TChain(treeName.c_str(),"DYEE");
				dyEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + eeChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttEE = new TChain(treeName.c_str(),"TTEE");
				ttEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"EE.root" ).c_str() );
				TChain * wrEE = new TChain(treeName.c_str(),"WREE");
				wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
	
				TChain * ttMuMu = new TChain(treeName.c_str(),"TTMuMu");
				ttMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"MuMu.root" ).c_str() );
				TChain * dyMuMu = new TChain(treeName.c_str(),"DYMuMu");
				dyMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + mumuChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * wrMuMu = new TChain(treeName.c_str(),"WRMuMu");
				wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
	
				chainPointers["EEDY"] = dyEE;
				chainPointers["EEWR"] = wrEE;
				chainPointers["EETop"] = ttEE;
				chainPointers["MuMuTop"] = ttMuMu;
				chainPointers["MuMuDY"] = dyMuMu;
				chainPointers["MuMuWR"] = wrMuMu;
			}
		
			else if(uncertTagNamesVect[s].find("jet") != string::npos){
				TChain * dyMuMu = new TChain(treeName.c_str(),"DYMuMu");
				dyMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + mumuChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttMuMu = new TChain(treeName.c_str(),"TTMuMu");
				ttMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"MuMu.root" ).c_str() );
				TChain * wrMuMu = new TChain(treeName.c_str(),"WRMuMu");
				wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
				TChain * dyEE = new TChain(treeName.c_str(),"DYEE");
				dyEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + eeChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttEE = new TChain(treeName.c_str(),"TTEE");
				ttEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"EE.root" ).c_str() );
				TChain * wrEE = new TChain(treeName.c_str(),"WREE");
				wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
				chainPointers["MuMuDY"] = dyMuMu;
				chainPointers["MuMuWR"] = wrMuMu;
				chainPointers["EEDY"] = dyEE;
				chainPointers["EEWR"] = wrEE;
				chainPointers["EETop"] = ttEE;
				chainPointers["MuMuTop"] = ttMuMu;
			}

			else if(uncertTagNamesVect[s].find("electron") != string::npos){
				TChain * dyEE = new TChain(treeName.c_str(),"DYEE");
				dyEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + eeChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttEE = new TChain(treeName.c_str(),"TTEE");
				ttEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"EE.root" ).c_str() );
				TChain * wrEE = new TChain(treeName.c_str(),"WREE");
				wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
				chainPointers["EEDY"] = dyEE;
				chainPointers["EEWR"] = wrEE;
				chainPointers["EETop"] = ttEE;
			}

			else if(uncertTagNamesVect[s].find("muon") != string::npos){
				TChain * dyMuMu = new TChain(treeName.c_str(),"DYMuMu");
				dyMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + mumuChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttMuMu = new TChain(treeName.c_str(),"TTMuMu");
				ttMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"MuMu.root" ).c_str() );
				TChain * wrMuMu = new TChain(treeName.c_str(),"WRMuMu");
				wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
				chainPointers["MuMuDY"] = dyMuMu;
				chainPointers["MuMuWR"] = wrMuMu;
				chainPointers["MuMuTop"] = ttMuMu;
			}

			else if(uncertTagNamesVect[s].find("all") != string::npos){ //all systematics enabled
				TChain * dyMuMu = new TChain(treeName.c_str(),"DYMuMu");
				dyMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + mumuChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttMuMu = new TChain(treeName.c_str(),"TTMuMu");
				ttMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"MuMu.root" ).c_str() );
				TChain * wrMuMu = new TChain(treeName.c_str(),"WRMuMu");
				wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
				
				TChain * dyEE = new TChain(treeName.c_str(),"DYEE");
				dyEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + dyFileName + eeChannelInFileNames +"_withMllWeight.root" ).c_str() );
				TChain * ttEE = new TChain(treeName.c_str(),"TTEE");
				ttEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + emuDataFileName +"EE.root" ).c_str() );
				TChain * wrEE = new TChain(treeName.c_str(),"WREE");
				wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_" + to_string(wrMassVect[m]/2) + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
	
				if(wrMassVect[m] == 1800 ){
					wrMuMu->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoMuMuJJ_" + to_string(wrMassVect[m]) + "_1400"  + "_signal_" + mumuChannelInFileNames +".root" ).c_str() );
					wrEE->Add( (absPathToMainRootFileDir + relDirPathsVect[s] + "selected_tree_WRtoEEJJ_" + to_string(wrMassVect[m]) + "_1400" + "_signal_" + eeChannelInFileNames +".root" ).c_str() );
	
				}
				
				chainPointers["MuMuDY"] = dyMuMu;
				chainPointers["MuMuWR"] = wrMuMu;
				chainPointers["EEDY"] = dyEE;
				chainPointers["EEWR"] = wrEE;
				chainPointers["EETop"] = ttEE;
				chainPointers["MuMuTop"] = ttMuMu;
			}

#ifdef DEBUG
			cout<<"#WR mass\tprocess channel\tsyst source\tmean\tstat uncert\t syst uncert\tsyst and stat uncert percent\tsyst sqd plus stat sqd over nEvts\tmean sqd over syst sqd plus stat sqd"<< endl;
#endif

			for(map<string,TChain*>::const_iterator chMapIt=chainPointers.begin(); chMapIt!=chainPointers.end(); chMapIt++){
				//get the integer corresponding to the correct index position in NEventsInRange
				//expected num evts and syst uncert
#ifdef useDYMAD
				if((chMapIt->first).find("DY") != string::npos) continue;

#endif
				Double_t mean, statDevSqd, systDevSqd, evtsNoWgts;

				getExpEvtsUnwgtEvtsAndUncsSqd(mean, evtsNoWgts, statDevSqd, systDevSqd, wrMassVect[m], chMapIt->first, pathToMassWindowsFile, chMapIt->second);


				if((chMapIt->first).find("Top") != string::npos){
					//rescale mean, unweighted evts and uncertainties by the emu data scale factor
					Double_t rescale = ((chMapIt->first).find("EE") != string::npos) ? 0.4194 : 0.6563;
					mean *= rescale;
					systDevSqd *= (rescale*rescale);
					statDevSqd *= (rescale*rescale);
					evtsNoWgts *= rescale;
				}
				//writeToSystFile << wrMassVect[m] << "\t"<< chMapIt->first << "\t" << uncertTagNamesVect[s] << "\t" << mean << "\t" << (statDevSqd + systDevSqd)/mean << "\t" << mean*mean/(statDevSqd + systDevSqd) << endl;
				writeToSystFile << wrMassVect[m] << "\t"<< chMapIt->first << "\t" << uncertTagNamesVect[s] << "\t" << mean << "\t" << sqrt(statDevSqd) << "\t"<< sqrt(systDevSqd) <<"\t"<< evtsNoWgts << endl;
	

			}//end loop over TChains at one mass point and both lepton channels
			writeToSystFile << "\t" << endl;
			writeToSystFile << "\t" << endl;


#ifdef useDYMAD

			//sum expected evts and unweighted evts, add stat and syst uncertainties in quadrature
			//do electron channel first
			Double_t meanEEDYIncl, statUncSqdEEDYIncl, systUncSqdEEDYIncl, unweightEvtsEEDYIncl;
			Double_t meanEEDYHT, statUncSqdEEDYHT, systUncSqdEEDYHT, unweightEvtsEEDYHT;

			getExpEvtsUnwgtEvtsAndUncsSqd(meanEEDYHT, unweightEvtsEEDYHT, statUncSqdEEDYHT, systUncSqdEEDYHT, wrMassVect[m], dyEETagName, pathToMassWindowsFile, chainPointers[dyEETagName]);
			getExpEvtsUnwgtEvtsAndUncsSqd(meanEEDYIncl, unweightEvtsEEDYIncl, statUncSqdEEDYIncl, systUncSqdEEDYIncl, wrMassVect[m], dyEETwoTagName, pathToMassWindowsFile, chainPointers[dyEETwoTagName]);

			writeToSystFile << wrMassVect[m] << "\t" << "DYTotEE" << "\t" << uncertTagNamesVect[s] << "\t" << (meanEEDYHT+meanEEDYIncl) << "\t" << sqrt(statUncSqdEEDYHT+statUncSqdEEDYIncl) << "\t"<< sqrt(systUncSqdEEDYHT+systUncSqdEEDYIncl) <<"\t"<< (unweightEvtsEEDYHT+unweightEvtsEEDYIncl) << endl;
			//writeToSystFile << wrMassVect[m] << "\t" << "DYTotEE" << "\t" << uncertTagNamesVect[s] << "\t" << (meanEEDYHT+meanEEDYIncl) << "\t" << (statUncSqdEEDYHT+statUncSqdEEDYIncl+systUncSqdEEDYHT+systUncSqdEEDYIncl)/(meanEEDYHT+meanEEDYIncl) <<"\t"<< (meanEEDYHT+meanEEDYIncl)*(meanEEDYHT+meanEEDYIncl)/(statUncSqdEEDYHT+statUncSqdEEDYIncl+systUncSqdEEDYHT+systUncSqdEEDYIncl) << endl;


			///////////////////////////////////////////
			//now do muon channel
			Double_t meanMuMuDYIncl, statUncSqdMuMuDYIncl, systUncSqdMuMuDYIncl, unweightEvtsMuMuDYIncl;
			Double_t meanMuMuDYHT, statUncSqdMuMuDYHT, systUncSqdMuMuDYHT, unweightEvtsMuMuDYHT;

			getExpEvtsUnwgtEvtsAndUncsSqd(meanMuMuDYHT, unweightEvtsMuMuDYHT, statUncSqdMuMuDYHT, systUncSqdMuMuDYHT, wrMassVect[m], dyMuMuTagName, pathToMassWindowsFile, chainPointers[dyMuMuTagName]);
			getExpEvtsUnwgtEvtsAndUncsSqd(meanMuMuDYIncl, unweightEvtsMuMuDYIncl, statUncSqdMuMuDYIncl, systUncSqdMuMuDYIncl, wrMassVect[m], dyMuMuTwoTagName, pathToMassWindowsFile, chainPointers[dyMuMuTwoTagName]);

			writeToSystFile << wrMassVect[m] << "\t" << "DYTotMuMu" << "\t" << uncertTagNamesVect[s] << "\t" << (meanMuMuDYHT+meanMuMuDYIncl) << "\t" << sqrt(statUncSqdMuMuDYHT+statUncSqdMuMuDYIncl) << "\t"<< sqrt(systUncSqdMuMuDYHT+systUncSqdMuMuDYIncl) <<"\t"<< (unweightEvtsMuMuDYHT+unweightEvtsMuMuDYIncl) << endl;
			//writeToSystFile << wrMassVect[m] << "\t" << "DYTotMuMu" << "\t" << uncertTagNamesVect[s] << "\t" << (meanMuMuDYHT+meanMuMuDYIncl) << "\t" << (statUncSqdMuMuDYHT+statUncSqdMuMuDYIncl+systUncSqdMuMuDYHT+systUncSqdMuMuDYIncl)/(meanMuMuDYHT+meanMuMuDYIncl) <<"\t"<< (meanMuMuDYHT+meanMuMuDYIncl)*(meanMuMuDYHT+meanMuMuDYIncl)/(statUncSqdMuMuDYHT+statUncSqdMuMuDYIncl+systUncSqdMuMuDYHT+systUncSqdMuMuDYIncl) << endl;

#endif

			//loop over chainPointers elements, and clear allocated memory
			for(map<string,TChain*>::const_iterator mapIt=chainPointers.begin(); mapIt!=chainPointers.end(); mapIt++){ (mapIt->second)->Delete();}
			chainPointers.clear();

		}//end loop over different systematics to evaluate

	}//end loop over different WR mass points
	
	writeToSystFile.close();
	

}//end printSystStdDevToFile()
