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
	
	makeAndSaveMultipleCurveOverlayHisto(numerPlottingChainMap,"c1",0.7,0.6,0.98,0.9,false,title,xLabel,outFileTag+"_Data",false, -10);

	denomPlottingChainMap[branch+">>DYMCBeforeCut_dileptonMass(50,70.,110.),("+baseCut+")"] = denomCh;
	denomPlottingChainMap[branch+">>DYMCPassingCut_dileptonMass(50,70.,110.),("+tightCut+")"] = denomCh;
	
	makeAndSaveMultipleCurveOverlayHisto(denomPlottingChainMap,"c2",0.7,0.6,0.98,0.9,false,title,xLabel,outFileTag+"_DYMC",false, -10);
	
}//end calculateSFwithTwoChains()

void egmTnPEfficiency(){
	TString pathToFiles = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/egmTnPRootFiles/";
	TString inclusiveReco = "GsfElectronToSC/fitter_tree", inclusiveHeep = "GsfElectronToRECO/fitter_tree", inclusiveHlt = "GsfElectronToTrigger/fitter_tree";
	TChain * inclusiveRecoEffData = new TChain(inclusiveReco);
	TChain * inclusiveHeepEffData = new TChain(inclusiveHeep);
	TChain * inclusiveHltEtEffData = new TChain(inclusiveHlt);
	TChain * inclusiveHltIdEffData = new TChain(inclusiveHlt);
	TChain * inclusiveDblEle33HltIdEffData = new TChain(inclusiveHlt);

	TChain * inclusiveRecoEffDYMC = new TChain(inclusiveReco);
	TChain * inclusiveHeepEffDYMC = new TChain(inclusiveHeep);
	TChain * inclusiveHltEtEffDYMC = new TChain(inclusiveHlt);
	TChain * inclusiveHltIdEffDYMC = new TChain(inclusiveHlt);
	TChain * inclusiveDblEle33HltIdEffDYMC = new TChain(inclusiveHlt);


	inclusiveRecoEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfRecoUsingSC4path_tagPassesHighETProbePassesLowET_allWrSkims.root");
	inclusiveHeepEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfHeepUsingSC4path_tagPassesHighETProbePassesLowET_allWrSkims.root");
	inclusiveHltEtEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfHighETUsingSC4path_tagPassesHighETProbePassesLowET_allWrSkims.root");
	inclusiveHltIdEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfWP60IdUsingSC4path_tagPassesTightLegProbePassesHighET_allWrSkims.root");
	inclusiveDblEle33HltIdEffData->Add(pathToFiles+"TnPTree_data_efficiencyOfUnseededDoubleEle33CaloTrkIdCuts_tagEt50PassesL1SeededDblEle33LegProbeEt40PassesEtUnseeded_allWrSkims.root");

	//DYAMC Incl files
	inclusiveRecoEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfRecoUsingSC4path_tagPassesHighETProbePassesLowET_allWrSkims.root");
	inclusiveHeepEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfHeepUsingSC4path_tagPassesHighETProbePassesLowET_allWrSkims.root");
	inclusiveHltEtEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfHighETUsingSC4path_tagPassesHighETProbePassesLowET_allWrSkims.root");
	inclusiveHltIdEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfWP60IdUsingSC4path_tagPassesTightLegProbePassesHighET_allWrSkims.root");
	inclusiveDblEle33HltIdEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfUnseededDoubleEle33CaloTrkIdCuts_tagEt50PassesL1SeededDblEle33LegProbeEt40PassesEtUnseeded_allWrSkims.root");

	//DYMadIncl files
	//inclusiveRecoEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfRecoUsingSC4path_tagPassesHighETProbePassesLowET_allDyMadInclSkims.root");
	//inclusiveHeepEffDYMC->Add(pathToFiles+"TnPTree_mc_efficiencyOfHeepUsingSC4path_tagPassesHighETProbePassesLowET_allDyMadInclSkims.root");
	
	///calculate inclusive efficiencies and SFs which do not depend on ET
	string inclusiveEfficiencyFile = "inclusiveEleRecoHeepHltEfficiency.txt";
	ofstream writeToInclusiveEfficiencyFile(inclusiveEfficiencyFile.c_str(),ofstream::trunc);
	Float_t dataInclRecoEff = 0, dataInclHeepEff=0, dataInclHltEtEff=0, dataInclHltIdEff=0;
	Float_t dataInclRecoEffUnc = -1, dataInclHeepEffUnc=-1, dataInclHltEtEffUnc=-1, dataInclHltIdEffUnc=-1;
	Float_t dymcInclRecoEff = 0, dymcInclHeepEff=0, dymcInclHltEtEff=0, dymcInclHltIdEff=0;
	Float_t dymcInclRecoEffUnc = -1, dymcInclHeepEffUnc=-1, dymcInclHltEtEffUnc=-1, dymcInclHltIdEffUnc=-1;
	Float_t datadymcRecoSf = 0, datadymcHeepSf=0, datadymcHltEtSf=0, datadymcHltIdSf=0;
	Float_t datadymcRecoSfUnc = -1, datadymcHeepSfUnc=-1, datadymcHltEtSfUnc=-1, datadymcHltIdSfUnc=-1;
	Float_t dataInclDblEle33HltIdEff = 0, dymcInclDblEle33HltIdEff = 0, datadymcDblEle33HltIdSf = 0, datadymcDblEle33HltIdSfUnc = -1, dataInclDblEle33HltIdEffUnc = -1, dymcInclDblEle33HltIdEffUnc = -1;
	
	calculateEffAndSFwithTwoChains(inclusiveRecoEffData, inclusiveRecoEffDYMC, "mass", "mass>80 && mass<100 && probe_sc_et>35", "mass>80 && mass<100 && probe_sc_et>35 && passingRECO>0", dataInclRecoEff, dataInclRecoEffUnc, dymcInclRecoEff, dymcInclRecoEffUnc, datadymcRecoSf, datadymcRecoSfUnc, "Ele RECO efficiency","M_{LL} (GeV)", "_recoEfficiency");
	calculateEffAndSFwithTwoChains(inclusiveHeepEffData, inclusiveHeepEffDYMC, "mass", "mass>80 && mass<100 && probe_Ele_et>35", "mass>80 && mass<100 && probe_Ele_et>35 && passingVeto>0", dataInclHeepEff, dataInclHeepEffUnc, dymcInclHeepEff, dymcInclHeepEffUnc, datadymcHeepSf, datadymcHeepSfUnc, "Ele HEEP ID efficiency","M_{LL} (GeV)", "_heepEfficiency");
	calculateEffAndSFwithTwoChains(inclusiveHltEtEffData, inclusiveHltEtEffDYMC, "mass", "mass>80 && mass<100 && probe_sc_et>35", "mass>80 && mass<100 && probe_sc_et>35 && passingHLT>0", dataInclHltEtEff, dataInclHltEtEffUnc, dymcInclHltEtEff, dymcInclHltEtEffUnc, datadymcHltEtSf, datadymcHltEtSfUnc, "Ele HLT E_{T} efficiency","M_{LL} (GeV)", "_hltEtEfficiency");
	calculateEffAndSFwithTwoChains(inclusiveHltIdEffData, inclusiveHltIdEffDYMC, "mass", "mass>80 && mass<100 && probe_sc_et>35", "mass>80 && mass<100 && probe_sc_et>35 && passingHLT>0", dataInclHltIdEff, dataInclHltIdEffUnc, dymcInclHltIdEff, dymcInclHltIdEffUnc, datadymcHltIdSf, datadymcHltIdSfUnc, "Ele HLT ID efficiency","M_{LL} (GeV)", "_hltIdEfficiency");
	calculateEffAndSFwithTwoChains(inclusiveDblEle33HltIdEffData, inclusiveDblEle33HltIdEffDYMC, "mass", "mass>80 && mass<100 && probe_sc_et>50", "mass>80 && mass<100 && probe_sc_et>50 && passingHLT>0", dataInclDblEle33HltIdEff, dataInclDblEle33HltIdEffUnc, dymcInclDblEle33HltIdEff, dymcInclDblEle33HltIdEffUnc, datadymcDblEle33HltIdSf, datadymcDblEle33HltIdSfUnc, "HLT DoubleEle33 ID efficiency","M_{LL} (GeV)", "_dblEle33HltIdEfficiency");

	writeToInclusiveEfficiencyFile << "#Reco and Heep and HLT Ele*WP60_SC4_Mass55 high ET cut eff calculated in data (MC) evts which fired HLT_Ele30(25)WP60_SC4_Mass55, one ele has pt>35 and is dR<0.2 matched to object passing high ET HLT cut, a second ele with pt>35 is dR<0.2 matched to HLT object passing SC4 cut, both eles have |eta|<2.4, and their dilepton mass is btwn 80 and 100" << std::endl;
	writeToInclusiveEfficiencyFile << "#HLT Ele*WP60_SC4_Mass55 HltId eff calculated in data (MC) evts which fired HLT_Ele30(25)WP60_SC4_Mass55, one ele has pt>35 and is dR<0.2 matched to object passing entire Ele*WP60 leg, a second ele with pt>35 is dR<0.2 matched to HLT object passing SC4 cut, both eles have |eta|<2.4, and their dilepton mass is btwn 80 and 100" << std::endl;
	writeToInclusiveEfficiencyFile << "#HLT DoubleEle33_CaloIdL_GsfTrkIdVL HltId eff calculated in data and MC evts which fired HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*, one ele has pt>50 and is dR<0.2 matched to object passing entire L1Seeded leg of DoubleEle33, a second ele with pt>50 is dR<0.2 matched to HLT object passing ET33 cut on unseeded leg, both eles have |eta|<2.4, and their dilepton mass is btwn 80 and 100" << std::endl;
	writeToInclusiveEfficiencyFile << "#the uncertainties quoted here are purely statistical" << std::endl;
	writeToInclusiveEfficiencyFile << "#" << std::endl;
	writeToInclusiveEfficiencyFile << "#PER ELE EFFICIENCY IN REAL DATA" << std::endl;
	writeToInclusiveEfficiencyFile << "Reco\t" << dataInclRecoEff<< "\t+/-\t"<< dataInclRecoEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "Heep\t" << dataInclHeepEff<< "\t+/-\t"<< dataInclHeepEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT Ele30WP60_SC4_Mass55 HltEt\t" << dataInclHltEtEff<< "\t+/-\t"<< dataInclHltEtEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT Ele30WP60_SC4_Mass55 HltId\t" << dataInclHltIdEff<< "\t+/-\t"<< dataInclHltIdEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT DoubleEle33_CaloIdL_GsfTrkIdVL HltId\t" << dataInclDblEle33HltIdEff<< "\t+/-\t"<< dataInclDblEle33HltIdEffUnc << std::endl;
	
	writeToInclusiveEfficiencyFile << "#PER ELE EFFICIENCY IN DYMC" << std::endl;
	writeToInclusiveEfficiencyFile << "Reco\t" << dymcInclRecoEff<< "\t+/-\t"<< dymcInclRecoEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "Heep\t" << dymcInclHeepEff<< "\t+/-\t"<< dymcInclHeepEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT Ele25WP60_SC4_Mass55 HltEt\t" << dymcInclHltEtEff<< "\t+/-\t"<< dymcInclHltEtEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT Ele25WP60_SC4_Mass55 HltId\t" << dymcInclHltIdEff<< "\t+/-\t"<< dymcInclHltIdEffUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT DoubleEle33_CaloIdL_GsfTrkIdVL HltId\t" << dymcInclDblEle33HltIdEff<< "\t+/-\t"<< dymcInclDblEle33HltIdEffUnc << std::endl;

	writeToInclusiveEfficiencyFile << "#PER ELE DATA/MC EFFICIENCY RATIO" << std::endl;
	writeToInclusiveEfficiencyFile << "Reco\t" << datadymcRecoSf << "\t+/-\t"<< datadymcRecoSfUnc << std::endl;
	writeToInclusiveEfficiencyFile << "Heep\t" << datadymcHeepSf << "\t+/-\t"<< datadymcHeepSfUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT Ele*WP60_SC4_Mass55 HltEt\t" << datadymcHltEtSf << "\t+/-\t"<< datadymcHltEtSfUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT Ele*WP60_SC4_Mass55 HltId\t" << datadymcHltIdSf << "\t+/-\t"<< datadymcHltIdSfUnc << std::endl;
	writeToInclusiveEfficiencyFile << "HLT DoubleEle33_CaloIdL_GsfTrkIdVL HltId\t" << datadymcDblEle33HltIdSf << "\t+/-\t"<< datadymcDblEle33HltIdSfUnc << std::endl;
	writeToInclusiveEfficiencyFile.close();
	

	delete inclusiveRecoEffData;
	delete inclusiveHeepEffData;
	delete inclusiveHltEtEffData;
	delete inclusiveHltIdEffData;
	delete inclusiveDblEle33HltIdEffData;
	delete inclusiveRecoEffDYMC;
	delete inclusiveHeepEffDYMC;
	delete inclusiveHltEtEffDYMC;
	delete inclusiveHltIdEffDYMC;
	delete inclusiveDblEle33HltIdEffDYMC;

}///end egmTnPEfficiency()

