#include <TFile.h>
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
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

///use this macro to take a file with multiple directories, one tree per directory, and dump
///all of the plots in all trees into pdf files with unique names
void SaveTreePlots(TChain * chain, string outputFileName){
	TObjArray * branches = chain->GetListOfBranches();
	TIter brItr(branches);
	for(TBranch * aBranch=(TBranch*) brItr.Next(); aBranch!=NULL; aBranch=(TBranch*) brItr.Next()){
		string brName = aBranch->GetName();
		if(brName.compare("ptGenEle")==0 || brName.compare("etaGenEle")==0 || brName.compare("phiGenEle")==0 || 
				brName.compare("ptGenJet")==0 || brName.compare("etaGenJet")==0 || brName.compare("phiGenJet")==0 ){
			for(unsigned int i=0; i<2; i++){
				string updatedBrName = "";
				updatedBrName += brName;
				updatedBrName += "[";
				updatedBrName += to_string(i);
				updatedBrName += "]";

				///now updatedBrName has the array name and element of interest, either [0] or [1]

				TCanvas * c1 = new TCanvas(updatedBrName.c_str(),updatedBrName.c_str(),800,800);
				c1->cd();
				chain->Draw(updatedBrName.c_str());
				string combName = chain->GetName();
				combName += "_";
				combName += updatedBrName;
				c1->SaveAs((outputFileName + "_" + combName+".png").c_str(), "recreate");

			}///end for loop

		}///end filter to do something different if the branch name corresponds to an array branch
		else{
			TCanvas * c1 = new TCanvas(brName.c_str(),brName.c_str(),800,800);
			c1->cd();
			chain->Draw(brName.c_str());
			string combName = chain->GetName();
			combName += "_";
			combName += brName;
			c1->SaveAs((outputFileName + "_" + combName+".png").c_str(), "recreate");

		}
	}///end loop over branches in tempChain to initialize keys and values in cutBranchesMap

}///end SaveTreePlots()

void dumpTreePlots(){
	TChain * noCuts = new TChain("genAnalyzerOne/genObjectsNoCuts","");
	noCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	TChain * ptEtaCuts = new TChain("genAnalyzerTwo/genObjectsWithPtEtaCuts","");
	ptEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	TChain * ptEtaDileptonMassCuts = new TChain("genAnalyzerThree/genObjectsWithPtEtaAndDileptonMassCuts","");
	ptEtaDileptonMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");

	string plotDir_noCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/noCuts";
	string plotDir_withPtEtaCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/withPtEtaCuts";
	string plotDir_withPtEtaDileptonMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/withPtEtaDileptonMassCuts";

	//SaveTreePlots(noCuts, plotDir_noCuts);
	//SaveTreePlots(ptEtaCuts, plotDir_withPtEtaCuts);
	SaveTreePlots(ptEtaDileptonMassCuts, plotDir_withPtEtaDileptonMassCuts);


}///end dumpTreePlots()
