//use this macro to make plots of and related to HEEP ID efficiency, RECO efficiency, and HLT efficiency for electrons in data and MC

#include <TFile.h>
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
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TObjArray.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include "macroSandBox.C"

using namespace std;

void calculateEffAndSFwithTwoChains(TChain * numerCh, TChain * denomCh, TString branchToScan, TString baseSelection, TString tighterSelection, Float_t & numerChEff, Float_t & numerChEffUnc, Float_t & denomChEff, Float_t & denomChEffUnc, Float_t & SF, Float_t & SFUnc, string title, string xLabel, string outFileTag){
	calculateEfficiencyWithOneChain(numerCh,branchToScan,baseSelection,tighterSelection,numerChEffUnc,numerChEff);
	calculateEfficiencyWithOneChain(denomCh,branchToScan,baseSelection,tighterSelection,denomChEffUnc,denomChEff);
	Float_t tempNumerEff = numerChEff, tempNumerEffUnc = numerChEffUnc, tempDenomEff = denomChEff, tempDenomEffUnc = denomChEffUnc;
	
	SF = (tempNumerEff/tempDenomEff);
	
	///SF uncertainty formula derived through error propagation on a ratio of two quantities with uncertainties
	SFUnc = sqrt( (1/pow(tempDenomEff,2))*(pow(tempNumerEffUnc,2)) + (pow(tempNumerEff,2)/pow(tempDenomEff,4))*(pow(tempDenomEffUnc,2)) );

	///plot the distribution used for scanning from each chain with base cuts and tighter cuts on one histo
	///two curves per plot, two plots per call to calculateEffAndSFwithTwoChains()
	map<string,TChain*> numerPlottingChainMap, denomPlottingChainMap;
	string baseCut(baseSelection), tightCut(tighterSelection), branch(branchToScan);
	numerPlottingChainMap[branch+">>dataBeforeCut_dileptonMass(50,70.,110.),("+baseCut+")"] = numerCh;
	numerPlottingChainMap[branch+">>dataPassingCut_dileptonMass(50,70.,110.),("+tightCut+")"] = numerCh;
	
	makeAndSaveMultipleCurveOverlayHisto(numerPlottingChainMap,"c1",0.7,0.6,0.98,0.9,false,title,xLabel,outFileTag+"_Data",false);

	denomPlottingChainMap[branch+">>DYMCBeforeCut_dileptonMass(50,70.,110.),("+baseCut+")"] = denomCh;
	denomPlottingChainMap[branch+">>DYMCPassingCut_dileptonMass(50,70.,110.),("+tightCut+")"] = denomCh;
	
	makeAndSaveMultipleCurveOverlayHisto(denomPlottingChainMap,"c2",0.7,0.6,0.98,0.9,false,title,xLabel,outFileTag+"_DYMC",false);
	
}//end calculateSFwithTwoChains()

void egmTnPEfficiency(){
	TString pathToFiles = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/egmTnPRootFiles/";
	TString inclusiveReco = "GsfElectronToSC/fitter_tree", inclusiveHeep = "GsfElectronToRECO/fitter_tree", inclusiveHlt = "GsfElectronToTrigger/fitter_tree";
	TChain * inclusiveRecoEffData = new TChain(inclusiveReco);
	TChain * inclusiveHeepEffData = new TChain(inclusiveHeep);
	TChain * inclusiveHltEtEffData = new TChain(inclusiveHlt);
	TChain * inclusiveHltIdEffData = new TChain(inclusiveHlt);
	
	TChain * inclusiveRecoEffDYMC = new TChain(inclusiveReco);
	TChain * inclusiveHeepEffDYMC = new TChain(inclusiveHeep);
	TChain * inclusiveHltEtEffDYMC = new TChain(inclusiveHlt);
	TChain * inclusiveHltIdEffDYMC = new TChain(inclusiveHlt);
	
	inclusiveRecoEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfRecoUsingSC4path_tagPassesHighETProbePassesLowET_thirtyInputFiles.root");
	inclusiveHeepEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfHeepUsingSC4path_minProbePtTwentyFive_thirtyInputFiles.root");
	inclusiveHltEtEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfHighETCutUsingSC4path_minProbePtTwentyFive_thirtyInputFiles.root");
	inclusiveHltIdEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfWP60CutsUsingSC4path_minProbePtTwentyFive_thirtyInputFiles.root");
	inclusiveRecoEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfRecoUsingSC4path_tagPassesHighETProbePassesLowET_thirtyInputFiles.root");
	inclusiveHeepEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfHeepUsingSC4path_minProbePtTwentyFive_thirtyInputFiles.root");
	inclusiveHltEtEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfHighETCutUsingSC4path_minProbePtTwentyFive_thirtyInputFiles.root");
	inclusiveHltIdEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfWP60CutsUsingSC4path_minProbePtTwentyFive_thirtyInputFiles.root");

	///calculate inclusive efficiencies and SFs which do not depend on ET
	string inclusiveEfficiencyFile = "inclusiveEleRecoHeepHltEfficiency.txt";
	ofstream writeToInclusiveEfficiencyFile(inclusiveEfficiencyFile.c_str(),ofstream::trunc);
	Float_t dataInclRecoEff = 0, dataInclHeepEff=0, dataInclHltEtEff=0, dataInclHltIdEff=0;
	Float_t dataInclRecoEffUnc = -1, dataInclHeepEffUnc=-1, dataInclHltEtEffUnc=-1, dataInclHltIdEffUnc=-1;
	Float_t dymcInclRecoEff = 0, dymcInclHeepEff=0, dymcInclHltEtEff=0, dymcInclHltIdEff=0;
	Float_t dymcInclRecoEffUnc = -1, dymcInclHeepEffUnc=-1, dymcInclHltEtEffUnc=-1, dymcInclHltIdEffUnc=-1;
	Float_t datadymcRecoSf = 0, datadymcHeepSf=0, datadymcHltEtSf=0, datadymcHltIdSf=0;
	Float_t datadymcRecoSfUnc = -1, datadymcHeepSfUnc=-1, datadymcHltEtSfUnc=-1, datadymcHltIdSfUnc=-1;
	
	calculateEffAndSFwithTwoChains(inclusiveRecoEffData, inclusiveRecoEffDYMC, "mass", "mass>80 && mass<100 && probe_sc_et>35", "mass>80 && mass<100 && probe_sc_et>35 && passingRECO>0", dataInclRecoEff, dataInclRecoEffUnc, dymcInclRecoEff, dymcInclRecoEffUnc, datadymcRecoSf, datadymcRecoSfUnc, "Ele RECO efficiency","M_{LL} (GeV)", "_recoEfficiency");
	calculateEffAndSFwithTwoChains(inclusiveHeepEffData, inclusiveHeepEffDYMC, "mass", "mass>80 && mass<100 && probe_Ele_et>35", "mass>80 && mass<100 && probe_Ele_et>35 && passingVeto>0", dataInclHeepEff, dataInclHeepEffUnc, dymcInclHeepEff, dymcInclHeepEffUnc, datadymcHeepSf, datadymcHeepSfUnc, "Ele HEEP ID efficiency","M_{LL} (GeV)", "_heepEfficiency");
	calculateEffAndSFwithTwoChains(inclusiveHltEtEffData, inclusiveHltEtEffDYMC, "mass", "mass>80 && mass<100 && probe_sc_et>35", "mass>80 && mass<100 && probe_sc_et>35 && passingHLT>0", dataInclHltEtEff, dataInclHltEtEffUnc, dymcInclHltEtEff, dymcInclHltEtEffUnc, datadymcHltEtSf, datadymcHltEtSfUnc, "Ele HLT E_{T} efficiency","M_{LL} (GeV)", "_hltEtEfficiency");
	calculateEffAndSFwithTwoChains(inclusiveHltIdEffData, inclusiveHltIdEffDYMC, "mass", "mass>80 && mass<100 && probe_sc_et>35", "mass>80 && mass<100 && probe_sc_et>35 && passingHLT>0", dataInclHltIdEff, dataInclHltIdEffUnc, dymcInclHltIdEff, dymcInclHltIdEffUnc, datadymcHltIdSf, datadymcHltIdSfUnc, "Ele HLT ID efficiency","M_{LL} (GeV)", "_hltIdEfficiency");
	
	writeToInclusiveEfficiencyFile << "#PER ELE EFFICIENCY IN REAL DATA" << std::endl;
	writeToInclusiveEfficiencyFile << "Reco\t" << dataInclRecoEff<< "\t+/-\t"<< dataInclRecoEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "Heep\t" << dataInclHeepEff<< "\t+/-\t"<< dataInclHeepEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HltEt\t" << dataInclHltEtEff<< "\t+/-\t"<< dataInclHltEtEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HltId\t" << dataInclHltIdEff<< "\t+/-\t"<< dataInclHltIdEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "#PER ELE EFFICIENCY IN DYMC" << std::endl;
	writeToInclusiveEfficiencyFile << "Reco\t" << dymcInclRecoEff<< "\t+/-\t"<< dymcInclRecoEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "Heep\t" << dymcInclHeepEff<< "\t+/-\t"<< dymcInclHeepEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HltEt\t" << dymcInclHltEtEff<< "\t+/-\t"<< dymcInclHltEtEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HltId\t" << dymcInclHltIdEff<< "\t+/-\t"<< dymcInclHltIdEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "#PER ELE DATA/MC EFFICIENCY RATIO" << std::endl;
	writeToInclusiveEfficiencyFile << "Reco\t" << datadymcRecoSf << "\t+/-\t"<< datadymcRecoSfUnc << std::endl;
	writeToInclusiveEfficiencyFile << "Heep\t" << datadymcHeepSf << "\t+/-\t"<< datadymcHeepSfUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HltEt\t" << datadymcHltEtSf << "\t+/-\t"<< datadymcHltEtSfUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HltId\t" << datadymcHltIdSf << "\t+/-\t"<< datadymcHltIdSfUnc << std::endl;
	writeToInclusiveEfficiencyFile.close();
	
	//calculateEffAndSFwithTwoChains(TChain * numerCh, TChain * denomCh, TString branchToScan, TString baseSelection, TString tighterSelection, Float_t & numerChEff, Float_t & numerChEffUnc, Float_t & denomChEff, Float_t & denomChEffUnc, Float_t & SF, Float_t & SFUnc, string title, string xLabel, string outFileTag){
	










	delete inclusiveRecoEffData;
	delete inclusiveHeepEffData;
	delete inclusiveHltEtEffData;
	delete inclusiveHltIdEffData;
	delete inclusiveRecoEffDYMC;
	delete inclusiveHeepEffDYMC;
	delete inclusiveHltEtEffDYMC;
	delete inclusiveHltIdEffDYMC;

}///end egmTnPEfficiency()

