#include "TLorentzVector.h"
#include "TMath.h"
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
#include <TEventList.h>
#include <TEntryList.h>
#include <TEntryListArray.h>
#include "dumpTreePlots.C"
#include "printSystStdDevToFile.C"

using namespace std;

//#define sigRgnWrAndBkgnds
//#define multiStepCutEffsWRandBkgnds
//#define makeShiftedMCPileupFiles
//#define nMinusOneCutEffsGenAndReco
//#define studyGenWrKinematicsVsWrAndNuMasses
//#define genAndRecoWrPlotsMinimalCuts
//#define privateGenEff
//#define twoDimPlotGenWrAcceptance
//#define recoAndGenHLTEfficiency
//#define genPlotsUsingWRDecayProducts
//#define compareCentrallyProducedToPrivateWrSignal
//#define lowMassSkimmedBkgndOnRealData
//#define lowMassFlavorSidebandBkgndOnData
//#define checkWellSeparatedGenPtBins
//#define PtRatioProfiles
//#define RecoGenOverlays
//#define StudyEffectOfMassPairs
//#define bkgndOverlaidOnMatchedSignal
//#define sOverBsensitivity
#define showMassWindows
//#define printNewDySyst
//#define DYHTPlot

//#define DEBUG
//#define DEBUGEVTWEIGHTMTHD
//#define DEBUGVECTOR


//given a TChain and a branch name, calculate and return the average value of all entries in the branch
Double_t calculateBranchMean(TChain * chain, std::string branchName, Double_t weight){
	chain->Draw((branchName+">>tempHist()").c_str(), (to_string(weight)+"*(1>0)").c_str() );
	TH1D * tempHist = (TH1D*) gROOT->FindObject("tempHist");
	Double_t mean = tempHist->GetMean(1);	//1 for x axis, 2 for y axis, 3 for z axis
#ifdef DEBUG
	cout<< "mean = "<< mean << endl;
#endif
	return mean;
}//end calculateBranchMean()


//given a TChain and a branch name, calculate and return the std dev of all entries in the branch
Double_t calculateBranchStdDev(TChain * chain, std::string branchName, Double_t weight){
	chain->Draw((branchName+">>tempHist()").c_str(), (to_string(weight)+"*(1>0)").c_str() );
	TH1D * tempHist = (TH1D*) gROOT->FindObject("tempHist");
	Double_t stdDev = tempHist->GetStdDev(1);	//1 for x axis, 2 for y axis, 3 for z axis
#ifdef DEBUG
	cout<< "stdDev = "<< stdDev << endl;
#endif
	return stdDev;
}//end calculateBranchStdDev()



//calculate four object mass given pt, eta, phi of four objects
//cannot get this function to work within a TTree Draw call
Float_t lljj(Float_t ptOne, Float_t etaOne, Float_t phiOne, Float_t ptTwo, Float_t etaTwo, Float_t phiTwo, Float_t ptThree, Float_t etaThree, Float_t phiThree, Float_t ptFour, Float_t etaFour, Float_t phiFour){
#ifdef DEBUG
	cout<<" "<<endl;
	cout<<"enterred calcFourObjMass fxn"<<endl;
#endif
	TLorentzVector One, Two, Three, Four;
#ifdef DEBUG
	cout<<"declared four TLorentzVectors"<<endl;
#endif
	One.SetPtEtaPhiM(ptOne, etaOne, phiOne, 0);
	Two.SetPtEtaPhiM(ptTwo, etaTwo, phiTwo, 0);
	Three.SetPtEtaPhiM(ptThree, etaThree, phiThree, 0);
	Four.SetPtEtaPhiM(ptFour, etaFour, phiFour, 0);
#ifdef DEBUG
	cout<<"initialized pt, eta, phi, and mass vals for four TLorentzVectors"<<endl;
	cout<<"One.Pt()=\t"<< One.Pt() <<endl;
	cout<<"Two.Pt()=\t"<< Two.Pt() <<endl;
	cout<<"Three.Pt()=\t"<< Three.Pt() <<endl;
	cout<<"Four.Pt()=\t"<< Four.Pt() <<endl;
#endif
	Float_t mass = (One+Two+Three+Four).M();
#ifdef DEBUG
	cout<<"four object mass=\t"<< mass <<endl;
#endif
	return mass;

}//end lljj()

///for calculating deltaR given eta and phi of two objects
float deltaR(float etaOne, float phiOne, float etaTwo, float phiTwo){
	float pi=3.14159;
	float deta = (etaOne - etaTwo);
#ifdef DEBUG
	cout<<"deta =\t"<< deta <<endl;
#endif
	float dphi = fabs(phiOne-phiTwo);
	if(dphi > pi) dphi -= 2*pi;
#ifdef DEBUG
	cout<<"dphi =\t"<< dphi <<endl;
	cout<<"about to leave deltaR function"<<endl;
#endif
	return sqrt(deta*deta + dphi*dphi);

}//end deltaR()

//calculates the dilepton mass given the pt, eta, and phi of two objects
Float_t calcDileptonMass(Float_t ptOne, Float_t etaOne, Float_t phiOne, Float_t ptTwo, Float_t etaTwo, Float_t phiTwo){
	Float_t mll = 0.;
	Float_t mllSqd = 2*ptOne*ptTwo*(TMath::CosH(etaOne - etaTwo) - TMath::Cos(phiOne - phiTwo));
	if(mllSqd > 0.) mll = TMath::Sqrt(mllSqd);
	return mll;
}//end calcDileptonMass()

//returns 0 if cand is has same eta and pt as either One or Two, returns 1 if no match is found
Int_t matching(Float_t candEta, Float_t candPt, Float_t etaOne, Float_t ptOne, Float_t etaTwo, Float_t ptTwo){
	if(candPt == ptOne && candEta == etaOne) return 0.;
	if(candPt == ptTwo && candEta == etaTwo) return 0.;
	return 1.;

}//end matching()

///use this fxn to calculate the efficiency of a selection applied to a TChain, and return the efficiency value
void calculateEfficiencyWithOneChain(TChain * chain,TString branchToScan,TString baseSelection,TString tighterSelection,Float_t & effUnc,Float_t & eff){
	chain->SetScanField(0);
	Long64_t evtsBeforeCut = chain->Scan(branchToScan,baseSelection);
	Long64_t evtsAfterCut = chain->Scan(branchToScan,tighterSelection);
	eff = ((Float_t) evtsAfterCut/evtsBeforeCut);
	effUnc = eff*sqrt(((Float_t) 1/evtsBeforeCut) + ((Float_t) 1/evtsAfterCut));

}///end calculateEfficiencyWithOneChain()

///use this fxn to plot and save one histogram using one branch from one TTree (or TChain)
void makeAndSaveSingleHistoFromTree(TChain * chain,string canvName,string treeDrawArgs,string histName,string histTitle,string xAxisTitle,string outputFileName){
	gStyle->SetOptStat("");
	TCanvas * cc = new TCanvas(canvName.c_str(),canvName.c_str(),750,700);
	cc->cd();
	chain->Draw(treeDrawArgs.c_str());
	TH1F * tempHist = (TH1F*) gROOT->FindObject(histName.c_str());
	tempHist->SetTitle(histTitle.c_str());
	tempHist->GetXaxis()->SetTitle(xAxisTitle.c_str());
	tempHist->SetLineWidth(3);
	tempHist->SetLineColor(1);

	//adjust the horizontal axis range
	//Double_t mean = tempHist->GetMean(1);
	//tempHist->GetXaxis()->SetRangeUser((mean)/2,1.3*(mean));
	resetXaxisLimits(tempHist);	///<defined in dumpTreePlots.C
	tempHist->Draw("");
	cc->SaveAs(outputFileName.c_str(),"recreate");
	tempHist->Delete();

}///end makeAndSaveSingleHistoFromTree()


///use this fxn to plot and save one histogram using one branch from one TTree (or TChain) with cuts
void makeAndSaveSingleHistoFromTreeWithCuts(TChain * chain,string canvName,string cuts,string treeDrawArgs,string histName,string histTitle,string xAxisTitle,string outputFileName,Float_t efficiencyThreshold,bool saveCutEffToFile,string cutEffsFilePath,string branchName,bool multiCutEff, bool autoResetLimits){
#ifdef DEBUG
	cout<<"  "<<endl;
	cout<<"  "<<endl;
	cout<<"  "<<endl;
	cout<<"in makeAndSaveSingleHistoFromTreeWithCuts fxn"<<endl;
#endif
	gStyle->SetOptStat("");
	TCanvas * cc = new TCanvas(canvName.c_str(),canvName.c_str(),750,700);
	cc->cd();
#ifdef DEBUG
	cout<<"about to draw tree"<<endl;
	cout<<"treeDrawArgs=\t"<< treeDrawArgs <<endl;
	cout<<"cuts=\t"<< cuts <<endl;
#endif
	chain->Draw(treeDrawArgs.c_str(),cuts.c_str());
#ifdef DEBUG
	cout<<"successfully drew tree, about to look for histo with FindObject"<<endl;
#endif
	TH1F * tempHist = (TH1F*) gROOT->FindObject(histName.c_str());
	tempHist->SetTitle(histTitle.c_str());
	tempHist->GetXaxis()->SetTitle(xAxisTitle.c_str());
	tempHist->SetLineWidth(3);
	tempHist->SetLineColor(1);


#ifdef DEBUG
	cout<<"about to reset x axis limits"<<endl;
	cout<<"autoResetLimits = "<< autoResetLimits <<endl;
	cout<<"hist Min =\t"<< tempHist->GetXaxis()->GetBinCenter(1) <<endl;
	cout<<"hist Max =\t"<< tempHist->GetXaxis()->GetBinCenter(tempHist->GetNbinsX()) <<endl;
#endif
	//adjust the horizontal axis range
	if(autoResetLimits) resetXaxisLimits(tempHist);	///< defined in dumpTreePlots.C
	if(saveCutEffToFile && !multiCutEff){
		Float_t cutEffDenom = (Float_t) chain->GetEntries(cuts.c_str());
		Float_t cutEffNumer = (Float_t) chain->GetEntries( (cuts + " && " + branchName + " > " + to_string(efficiencyThreshold) ).c_str() );
		ofstream writeToCutEffFile(cutEffsFilePath.c_str(),ofstream::app);
		writeToCutEffFile << chain->GetTitle() << "  &  " << histName << "  &  "<< (cutEffNumer/cutEffDenom) << "  DBLSLSH" << endl;
		cout<<"efficiency of lower bound cut with threshold=\t"<< efficiencyThreshold <<"\t=\t"<< (cutEffNumer/cutEffDenom) <<endl;
	}

	if(saveCutEffToFile && multiCutEff){
		//if multiCutEff is true, then the field named branchName will have several branch names and the value of each cut"
		Float_t cutEffDenom = (Float_t) chain->GetEntries(cuts.c_str());
		Float_t cutEffNumer = (Float_t) chain->GetEntries( (cuts + " && " + branchName).c_str() );
		ofstream writeToCutEffFile(cutEffsFilePath.c_str(),ofstream::app);
		writeToCutEffFile << chain->GetTitle() << "  &  " << branchName << "  &  "<< (cutEffNumer/cutEffDenom) << "  DBLSLSH" << endl;
		cout<<"efficiency of the cuts \t"<< branchName <<"\t=\t"<< (cutEffNumer/cutEffDenom) <<endl;
	}


	if(xAxisTitle.find_first_of(" not ") != string::npos && tempHist->GetNbinsX() < 5) tempHist->Draw("HISTTEXT90");	///<for histos showing how often a leading or subleading GEN particle is not matched to the expected GEN mother
	else tempHist->Draw("");
#ifdef DEBUG
	cout<<"drew histo"<<endl;
#endif
	cc->SaveAs((outputFileName+".png").c_str(),"recreate");
	cc->SaveAs((outputFileName+".pdf").c_str(),"recreate");
	tempHist->Delete();
	cc->Close();

}///end makeAndSaveSingleHistoFromTreeWithCuts()


///use this fxn to plot and save one histogram using one branch from one TTree (or TChain) with a fitted curve
void makeAndSaveSingleHistoFromTreeWithFit(TChain * chain,string canvName,string cuts,string treeDrawArgs,string histName,string histTitle,string xAxisTitle,string outputFileName, TF1 * fitFxn){
	gStyle->SetOptStat("");
	gStyle->SetOptFit(1112);
	TCanvas * cc = new TCanvas(canvName.c_str(),canvName.c_str(),750,700);
	cc->cd();
	chain->Draw(treeDrawArgs.c_str(),cuts.c_str());
	TH1F * tempHist = (TH1F*) gROOT->FindObject(histName.c_str());
	tempHist->Fit(fitFxn,"QR");
	//tempHist->Fit(fitFxn,"Q");
	tempHist->SetTitle(histTitle.c_str());
	tempHist->GetXaxis()->SetTitle(xAxisTitle.c_str());
	tempHist->SetLineWidth(3);
	tempHist->SetLineColor(1);

	//adjust the horizontal axis range
	resetXaxisLimits(tempHist);	///< defined in dumpTreePlots.C
	tempHist->Scale(1/tempHist->Integral());
	Float_t histMax = tempHist->GetBinContent(tempHist->GetMaximumBin());
	Float_t fitMax = fitFxn->GetMaximum(400,7000);
	
	//rescale the hist area so that the hist max equals the fit max divided by 7
	tempHist->Scale((fitMax/10)/histMax);
	
	//rescale the vertical axis
	histMax = tempHist->GetBinContent(tempHist->GetMaximumBin());
	//Float_t newMax = ((fitMax > histMax) ? 1.3*fitMax : 1.3*histMax);
	Float_t newMax = (1.3*histMax);
	
	tempHist->SetMaximum(newMax);
#ifdef DEBUG
	cout<<"histMax =\t"<< histMax <<endl;
	cout<<"fitMax =\t"<< fitMax <<endl;
	cout<<"location of fitMax=\t"<< fitFxn->GetMaximumX(400,7000) <<endl;
	cout<<"newMax =\t"<< newMax <<endl;
	cout<<"tempHist integral=\t"<< tempHist->Integral() <<endl;
	cout<<"fit integral=\t"<< fitFxn->Integral(0,8000) <<endl;
	cout<<"fit param 0=\t"<< fitFxn->GetParameter(0) <<endl;
	cout<<"fit param 1=\t"<< fitFxn->GetParameter(1) <<endl;
#endif
	tempHist->Draw("");
	fitFxn->SetLineColor(kRed);
	fitFxn->Draw("LSAME");
	cc->SaveAs((outputFileName+".png").c_str(),"recreate");
	cc->SaveAs((outputFileName+".pdf").c_str(),"recreate");
	tempHist->Delete();
	cc->Close();

}///end makeAndSaveSingleHistoFromTreeWithFit()


///use this fxn to plot and save one histogram using two branches from one TTree (or TChain) and a friend TTree or TChain
void makeAndSaveSingleHistoFromTreeWithFriend(TChain * chain,TChain * friendChain,string friendAlias,string canvName,string treeDrawArgs,string histName,string histTitle,string xAxisTitle,string outputFileName){
	gStyle->SetOptStat("");
	TCanvas * cc = new TCanvas(canvName.c_str(),canvName.c_str(),750,700);
	cc->cd();
	chain->AddFriend(friendChain,friendAlias.c_str());
	chain->Draw(treeDrawArgs.c_str());
	TH1 * tempHist = (TH1*) gROOT->FindObject(histName.c_str());
	tempHist->SetTitle(histTitle.c_str());
	tempHist->GetXaxis()->SetTitle(xAxisTitle.c_str());
	tempHist->SetLineWidth(3);
	tempHist->SetLineColor(1);

	//adjust the horizontal axis range
	Double_t mean = tempHist->GetMean(1);
	tempHist->GetXaxis()->SetRangeUser((mean)/3,1.6*(mean));
	tempHist->Draw("");
	cc->SaveAs(outputFileName.c_str(),"recreate");
	tempHist->Delete();

}///end makeAndSaveSingleHistoFromTreeWithFriend()



/**
 *
 * use this fxn to add std::pair<Double_t, string> containers to a vector such that the 0th element of the vector
 * has the lowest Double_t value, and the last element of the vector has the highest Double_t value
 *
 */
void addPairToOrderedVector(vector<pair<Double_t,string> > & orderedVector, pair<Double_t,string> & elementToAdd){
#ifdef DEBUGVECTOR
	cout<<"in addPairToOrderedVector"<<endl;
#endif

	if(orderedVector.size() == 0) orderedVector.push_back(elementToAdd);
	else{
		for(vector<pair<Double_t,string> >::iterator it=orderedVector.begin(); it!=orderedVector.end() ; it++){
#ifdef DEBUGVECTOR
			cout<<"elementToAdd.first =\t"<< elementToAdd.first << endl;
			cout<<"elementToAdd.second =\t"<< elementToAdd.second << endl;
			cout<<"(*it).first =\t"<< (*it).first << endl;
			cout<<"(*it).second =\t"<< (*it).second << endl;
#endif
			if(elementToAdd.first > (*it).first) continue;
			///if the continue does not kick in, then elementToAdd.first is <= the Double_t value stored in the element currently pointed to by it
			///add elementToAdd to orderedVector using insert
			orderedVector.insert(it,elementToAdd);
			break;
		}///end for loop over vector elements
	}///end else
#ifdef DEBUGVECTOR
	cout<<"added element to orderedVector"<<endl;
#endif

}///end addPairToOrderedVector()


/**
 * use this fxn to compute the evt weight for one evt in a TChain
 * return the weight as a Float_t value
 */
Float_t getEvtWeight(Float_t mcEvWgtSign, Int_t numVertices, map<string,vector<Float_t> > mcPuWeights, string mcName){
#ifdef DEBUGEVTWEIGHTMTHD
	cout<<"in getEvtWeight() fxn"<<endl;
	cout<<"num vertices = "<< numVertices <<endl;
#endif
	Float_t weight = 0;
	for(map<string,vector<Float_t> >::const_iterator puWgtIt=mcPuWeights.begin(); puWgtIt!=mcPuWeights.end(); puWgtIt++){
		///loop over the unique keys in mcPuWeights and find the one which matches the MC process name (ttBar, WZ, etc)
		if( (puWgtIt->second).size() <= numVertices ) return 0;
		if( (puWgtIt->first).find(mcName)!= string::npos ) weight = mcEvWgtSign*( (puWgtIt->second).at(numVertices) );
	}

	return weight;
}///end getEvtWeight()


/**
 * use this fxn to compute the PU weight based on the number of vertices in real data compared to MC
 * this fxn returns a map<string, vector<Float_t> > object.  The vector stores the weights (the index of
 * the vector corresponds to the number of vertices in the event), and the string key indicates the
 * bkgnd source.  Examples of the string key are ttBar and dyPlusJets. 
 */
map<string,vector<Float_t> > computePileupWeights(map<string,TChain*> bkgndChainMap, TChain * dataChain){
#ifdef DEBUG
	cout<<"in computePileupWeights fxn"<<endl;
#endif

	map<string,vector<Float_t> > wgtsMap;
	Long64_t dataEvts = dataChain->GetEntries();
	string val="50";

	///loop over chains in bkgndChainMap
	for(map<string,TChain*>::const_iterator bkgIt=bkgndChainMap.begin(); bkgIt!=bkgndChainMap.end(); bkgIt++){
		TH1F *h_data = new TH1F("h_data","",stoi(val),0,stof(val) );
		TH1F *hBkgnd = new TH1F("hBkgnd","",stoi(val),0,stof(val) );

		Int_t nvertices;
		Int_t nvertices_data;
		(bkgIt->second)->SetBranchAddress("nVertices",&nvertices);
		dataChain->SetBranchAddress("nVertices",&nvertices_data);

		Long64_t bkgndEvts = (bkgIt->second)->GetEntries();
		for (Long64_t ev = 0; ev < bkgndEvts; ev++) {
			(bkgIt->second)->GetEntry(ev);
			hBkgnd->Fill(nvertices);
		}
		for (Long64_t ev = 0; ev < dataEvts; ev++) {
			dataChain->GetEntry(ev);
			h_data->Fill(nvertices_data);
		}

		h_data->Scale(1/h_data->Integral());
		hBkgnd->Scale(1/hBkgnd->Integral());

		vector<Float_t> PUW(stoi(val)+1);

		for(unsigned int i = 0;i<PUW.size() ;i++){
			if(hBkgnd->GetBinContent(i) != 0) PUW[i] = h_data->GetBinContent(i)/hBkgnd->GetBinContent(i);
			else PUW[i] = 1.0;
		}

#ifdef DEBUG
		cout<<"filled a vector named PUW with weights for \t"<<bkgIt->first<<"\t bkgnd MC evts"<<endl;
		for(Int_t i=0; i<20; i+=5){
			cout<<"MC events with "<< i <<" vertices are weighted with this value: \t"<< PUW[i] <<endl;
		}///end loop over a subset of elements in weights vector
#endif
		///add the PU weights to the map
		wgtsMap[bkgIt->first] = PUW;

		hBkgnd->Delete();
		h_data->Delete();
	}///end loop over chains in bkgndChainMap

	return wgtsMap;
}///end computePileupWeights()

/**
 * use this fxn to overlay kinematic distributions from MC events using a THStack object
 * similar to overlayPointsOnStackedHistos()
 */
void makeStackedHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,map<string,Float_t> crossSxnsMap,Float_t intLumi,map<string,Float_t> nEvtsMap,TString treeCuts,Bool_t doLogYaxis){
	TCanvas * canv = new TCanvas(canvName,canvName,700,700);
	canv->cd();
	if(doLogYaxis) canv->SetLogy(1);
	TLegend * leg = new TLegend(legXmin,legYmin,legXmax,legYmax);
	map<string,TH1F*> stackedHistoMap;	///< links TString keys to TH1F histos which will ultimately be stacked
	
	///loop over elements in inputChainMap
	for(map<string,TChain*>::const_iterator chMapIt=inputChainMap.begin(); chMapIt!=inputChainMap.end(); chMapIt++){
		size_t openParenth = (chMapIt->first).find_first_of('('), lastChevron = (chMapIt->first).find_last_of('>');
		string uncutHistoName(chMapIt->first);
		///now initialize a new string, get rid of the content in uncutHistoName before '>>' and after '(',
		///and store the substring in the new string object
		string oneHistoName( uncutHistoName.substr(lastChevron+1,openParenth-lastChevron-1) );
		string afterLastChevron( uncutHistoName.substr(lastChevron+1) );
		(chMapIt->second)->Draw((chMapIt->first).c_str(), treeCuts);
		stackedHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
	}///end loop over elements in inputChainMap

	///now add all TH1F objects in stackedHistoMap into one THStack object
	///1 = black for colors
	int colors[] = {2,4,5,8,12,30,40,45};
	vector<int> colorVect(colors,colors + sizeof(colors)/sizeof(int) );
	Int_t i=0;
	THStack * histoStack = new THStack("","");
	vector<pair<Double_t,string> > vectOfBkgnds;	///< use this to link the # of rescaled bkgnd evts to each bkgnd source
	pair<Double_t,string> newPair;
	if(stackedHistoMap.size() > colorVect.size() ) cout<<"not enough unique colors in MultipleCurveOverlayHisto fxn!"<<endl;
#ifdef DEBUG
	std::cout<<"made a THStack object with null name and title"<<std::endl;
#endif

	Double_t bkgndIntegral;	///< use this to track the integral of all bkgnd histos in the THStack object
	for(map<string,TH1F*>::const_iterator histIt=stackedHistoMap.begin(); histIt!=stackedHistoMap.end(); histIt++){
		
		///set the line and fill color, and normalize the MC histos
		(histIt->second)->SetLineColor(colorVect[i]);
		(histIt->second)->SetFillColor(colorVect[i]);
#ifdef DEBUG
		std::cout<<"set line and fill color for the histo \t"<< histIt->first <<std::endl;
		std::cout<<"before rescaling there are "<< (histIt->second)->Integral() <<" entries in the histo mentioned immediately above"<<std::endl;
#endif

		for(map<string,Float_t>::const_iterator mcIt=nEvtsMap.begin(); mcIt!=nEvtsMap.end(); mcIt++){
			if((histIt->first).find(mcIt->first) != string::npos) (histIt->second)->Scale(crossSxnsMap[mcIt->first]*intLumi/nEvtsMap[mcIt->first]);
		}///end loop which normalizes the MC histos to integrated lumi of real data

		newPair = make_pair((histIt->second)->Integral(), histIt->first);
		addPairToOrderedVector(vectOfBkgnds, newPair);
		bkgndIntegral += (histIt->second)->Integral();
#ifdef DEBUG
		std::cout<<"rescaled the histo \t"<< histIt->first <<std::endl;
		std::cout<<"after rescaling there are "<< (histIt->second)->Integral() <<" entries in the histo mentioned immediately above"<<std::endl;
#endif
	
		//histoStack->Add((histIt->second));	///< add each histo in stackedHistoMap to the THStack object

#ifdef DEBUG
		std::cout<<"added \t"<< histIt->first <<"\t to the stacked histo object"<<std::endl;
#endif
		
		if(histIt==stackedHistoMap.begin()){
			histoStack->SetTitle( (histIt->second)->GetTitle() );
			string xLabel="";
			if((histIt->first).find("Mass") !=string::npos || (histIt->first).find("ptEle") !=string::npos || (histIt->first).find("ptJet") !=string::npos){
				xLabel += "GeV";
			}///end if(histo is plotting a variable with dimension of energy)
#ifdef DEBUG
			std::cout<<"xLabel = \t"<< xLabel <<std::endl;
#endif
			(histIt->second)->GetXaxis()->SetTitle(xLabel.c_str());
		}///end filter to set histogram X axis label
		
		size_t lastChevron = (histIt->first).find_last_of('>');
		size_t underscorePos = (histIt->first).find_first_of("_",lastChevron);
#ifdef DEBUG
		std::cout<<"lastChevron located at \t"<< lastChevron<<std::endl;
		std::cout<<"underscorePos located at \t"<< underscorePos<<std::endl;
#endif
		string legEntryName = (histIt->first).substr(lastChevron+1,underscorePos-lastChevron-1);
#ifdef DEBUG
		cout<<"histo name = \t"<< histIt->first <<endl;
		cout<<"legEntry has name = \t"<< legEntryName << endl;
#endif
		leg->AddEntry(histIt->second,legEntryName.c_str());
		i++;
	}///end loop to set line colors of histos, and add entries to TLegend

	for(unsigned int i=0; i<vectOfBkgnds.size(); i++){
		///calling second on the ith element of vectOfBkgnds returns a string which should be a key in stackedHistoMap
		histoStack->Add( stackedHistoMap[vectOfBkgnds[i].second] );
	}///end loop over pairs in vectOfBkgnds

	if(doLogYaxis){
		Double_t originalMax = histoStack->GetMaximum();
		histoStack->SetMaximum(25*originalMax);
		histoStack->SetMinimum(0.1);
	}
	histoStack->Draw("hist");

#ifdef DEBUG
	std::cout<<"stacked histo has been drawn"<<std::endl;
#endif
	
	string outputFile;
	map<string,TH1F*>::const_iterator hIt = stackedHistoMap.begin();
	size_t firstChevron = (hIt->first).find_first_of('>');
	outputFile = (hIt->first).substr(0,firstChevron);
	outputFile += ".png";
	
	leg->Draw();
	canv->SaveAs(outputFile.c_str(),"recreate");

}///end makeStackedHisto()

/**
 * use this fxn to overlay kinematic distributions from real data (points) onto MC (filled histos) using a THStack object and an overlaid TH1F object
 * the keys of crossSxnsMap and nEvtsMap are names of physics processes, like ttBar and dyPlusJets
 * the keys of inputChainMap contains the histogram plotting argument to use with TChain->Draw("plottingArgs")
 * the areas of the MC histograms are normalized to the integrated luminosity of real data
 * this is done by calling Scale(crossSxn * integratedLumi / numEvts) on each of the histos from MC data
 * if there are N unique keys in inputChainMap, then there are N-1 unique keys in the maps called crossSxnsMap and nEvtsMap
 *
 * the TChain to real data should contain the phrase "Data" in the inputChainMap key
 * the TChain to each MC sample should contain a key which is used in crossSxnsMap and nEvtsMap
 *
 */
void overlayPointsOnStackedHistos(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,map<string,Float_t> crossSxnsMap,Float_t intLumi,map<string,Float_t> nEvtsMap,TString treeCuts,Bool_t doPuReweighting,map<string,vector<Float_t> > puWeights,Bool_t doLogYaxis){
	TCanvas * canv = new TCanvas(canvName,canvName,750,700);
	canv->cd();
	if(doLogYaxis) canv->SetLogy(1);
	TLegend * leg = new TLegend(legXmin,legYmin,legXmax,legYmax);
	map<string,TH1F*> stackedHistoMap;	///< links TString keys to TH1F histos which will ultimately be stacked
	map<string,TH1F*> overlaidHistoMap; ///< links one TString key to the TH1F histo coming from real data
	for(map<string,TChain*>::const_iterator chMapIt=inputChainMap.begin(); chMapIt!=inputChainMap.end(); chMapIt++){
		size_t openParenth = (chMapIt->first).find_first_of('('), lastChevron = (chMapIt->first).find_last_of('>');
		string uncutHistoName(chMapIt->first);
		///now initialize a new string, get rid of the content in uncutHistoName before '>>' and after '(',
		///and store the substring in the new string object
		string oneHistoName( uncutHistoName.substr(lastChevron+1,openParenth-lastChevron-1) );
		string afterLastChevron( uncutHistoName.substr(lastChevron+1) );
		if((chMapIt->first).find("Data") == string::npos ){
			if(doPuReweighting){
				///to apply the GEN evt and PU weights:
				//1. identify the quantity which should be plotted
				//2. make a histogram with the appropriate limits and bins
				//3. link a local variable to the appropriate branch in the TChain
				//4. link two more local variables to the evWeightSign (Float_t) and nVertices (Int_t) branches
				//5. loop over all evts in the TChain, and fill the histogram with the appropriate weights (PU and gen evt weights)
				//string branchNames[] = {"ptEle[0]","ptEle[1]","dileptonMass","nVertices"};
				//string histoEndings[] = {"_leadLeptonPt(40,0.,200.)","_subleadLeptonPt(20,0.,100.)","_dileptonMass(50,0.,250.)","_nVertices(45,0.,45.)"};
				size_t firstChevron = (chMapIt->first).find_first_of('>'), closeParenthesis = afterLastChevron.find_first_of(')');
				size_t firstComma = afterLastChevron.find_first_of(','), secondComma = afterLastChevron.find_last_of(',');
				size_t firstUnderscore = afterLastChevron.find_first_of('_');
				size_t updatedOpenParenthesis = afterLastChevron.find_first_of('(');
				string targetBrName( (chMapIt->first).substr(0,firstChevron) );
				string nBinsString( (afterLastChevron).substr(updatedOpenParenthesis+1,firstComma - updatedOpenParenthesis-1) );
				string minString( (afterLastChevron).substr(firstComma+1,secondComma - firstComma-1) );
				string maxString( (afterLastChevron).substr(secondComma+1,closeParenthesis - secondComma-1) );
				string bkgndString( (afterLastChevron).substr(0,firstUnderscore) );
				Int_t vertices;
				Float_t evWgtSign;
				//Float_t mLL;
				//Float_t pts[2];
				(chMapIt->second)->SetBranchAddress("nVertices",&vertices);
				(chMapIt->second)->SetBranchAddress("evWeightSign",&evWgtSign);
				//(chMapIt->second)->SetBranchAddress("dileptonMass",&mLL);
				//(chMapIt->second)->SetBranchAddress("ptEle",pts);
				TH1F * hTemp = new TH1F(oneHistoName.c_str(),targetBrName.c_str(), stoi(nBinsString), stof(minString), stof(maxString) );
				
#ifdef DEBUG
				cout<<"two local vars have been linked to nVertices and evWeightSign branches"<<endl;
				cout<<"targetBrName = \t"<< targetBrName <<endl;
				cout<<"nBinsString = \t"<< nBinsString <<endl;
				cout<<"minString = \t"<< minString <<endl;
				cout<<"maxString = \t"<< maxString <<endl;
				cout<<"bkgndString = \t"<< bkgndString <<endl;
				cout<<"declared a histogram with correct name, number of bins, and min and max x axis values"<<endl;
#endif

				Long64_t totalEntries = (chMapIt->second)->GetEntries();
				Float_t wgt;
				Bool_t filledHisto = false;
				if(targetBrName.find("nVert") != string::npos || targetBrName.find("nJets") != string::npos || targetBrName.find("nLeptons") != string::npos || targetBrName.find("runNum") != string::npos || targetBrName.find("evtNum") != string::npos || targetBrName.find("leadingIs") != string::npos ){
					///use an Int_t variable for plotting nVertices, nJets, nLeptons, runNumber, evtNumber, and leadingIsHardest
					Int_t desiredInt;
					Bool_t useDesiredInt = false;
#ifdef DEBUG
					cout<<"in the Int_t filter"<<endl;
#endif
					if(targetBrName.find("nVertices") == string::npos){
						(chMapIt->second)->SetBranchAddress(targetBrName.c_str(), &desiredInt);
						useDesiredInt = true;
					}

					///now loop over all entries in the TChain, and fill the histogram with appropriate weights
					for(Long64_t ev=0; ev<totalEntries; ev++){
						(chMapIt->second)->GetEntry(ev);
						///now nVertices and evWeightSign in this particular event ev can be accessed
						//if(!(mLL<=120. && mLL>=60. && pts[0]>40 && pts[1]>40 ) ) continue;	///<hard coded cut to check Z peak region
						wgt = getEvtWeight(evWgtSign, vertices, puWeights, bkgndString);
#ifdef DEBUG
						if(ev == 10000) cout<<"on event # "<< ev <<endl;
						if(ev == 100000) cout<<"on event # "<< ev <<endl;
						if(ev == 300000) cout<<"on event # "<< ev <<endl;
#endif
						if(useDesiredInt) hTemp->Fill(desiredInt, wgt);
						else hTemp->Fill(vertices, wgt);
					}///end loop over entries in TChain
					filledHisto = true;

				}///end single Int_t branch filter
				if(targetBrName.find('[') == string::npos && filledHisto == false){
					///the branch of interest stores a single Float_t value
					Float_t desiredFloat;
					(chMapIt->second)->SetBranchAddress(targetBrName.c_str(), &desiredFloat);
					for(Long64_t ev=0; ev<totalEntries; ev++){
						(chMapIt->second)->GetEntry(ev);
						//if(!(mLL<=120. && mLL>=60. && pts[0]>40 && pts[1]>40 ) ) continue;	///<hard coded cut to check Z peak region
#ifdef DEBUG
						if(ev == 10000) cout<<"on event # "<< ev <<endl;
						if(ev == 100000) cout<<"on event # "<< ev <<endl;
						if(ev == 300000) cout<<"on event # "<< ev <<endl;
#endif
						wgt = getEvtWeight(evWgtSign, vertices, puWeights, bkgndString);
						hTemp->Fill(desiredFloat, wgt);
					}
					filledHisto = true;
	
				}///end single Float_t branch filter
				if(filledHisto == false){
					///the branch of interest stores a 2 element array of Float_t values
					///get the true branch name, and the element of the array which should be plotted (first or second)
					Float_t desiredArray[2];
					size_t openBracket = targetBrName.find_first_of('['), closeBracket = targetBrName.find_first_of(']');
					string trueBrName( targetBrName.substr(0,openBracket) );
					string arrElementString( targetBrName.substr(openBracket+1,closeBracket - openBracket-1) );
#ifdef DEBUG
					cout<<"trueBrName = \t"<< trueBrName <<endl;
					cout<<"arrElementString = \t"<< arrElementString <<endl;
#endif
					(chMapIt->second)->SetBranchAddress(trueBrName.c_str(), desiredArray);
					for(Long64_t ev=0; ev<totalEntries; ev++){
						(chMapIt->second)->GetEntry(ev);
						//if(!(mLL<=120. && mLL>=60. && pts[0]>40 && pts[1]>40 ) ) continue;	///<hard coded cut to check Z peak region
#ifdef DEBUG
						if(ev == 10000) cout<<"on event # "<< ev <<endl;
						if(ev == 100000) cout<<"on event # "<< ev <<endl;
						if(ev == 300000) cout<<"on event # "<< ev <<endl;
#endif
						wgt = getEvtWeight(evWgtSign, vertices, puWeights, bkgndString);
						hTemp->Fill(desiredArray[stoi(arrElementString)], wgt);
					}

				}///end Float_t[] array branch filter
			
			}///end if(doPuReweighting == true)
			else if(doPuReweighting == false) (chMapIt->second)->Draw((chMapIt->first).c_str(), treeCuts);


			stackedHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
		}
		else if((chMapIt->first).find("Data") != string::npos ){
			(chMapIt->second)->Draw((chMapIt->first).c_str(), treeCuts);
			overlaidHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
#ifdef DEBUG
			std::cout<<"there are \t"<< (overlaidHistoMap[chMapIt->first])->Integral() <<"\t entries in the real data histo"<<std::endl;
#endif


		}

#ifdef DEBUG
		std::cout<<"input chain map key = \t"<< chMapIt->first <<std::endl;
#endif
	
		/*
		if((chMapIt->first).find("Data") == string::npos ){
			stackedHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
		}///end search for histos made from MC datasets
		
		else if((chMapIt->first).find("Data") != string::npos){
			overlaidHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
#ifdef DEBUG
			std::cout<<"there are \t"<< (overlaidHistoMap[chMapIt->first])->Integral() <<"\t entries in the real data histo"<<std::endl;
#endif


		}
		*/
	}///end loop over elements in inputChainMap

#ifdef DEBUG
	std::cout<<"there are \t"<< stackedHistoMap.size() <<"\t elements in stackedHistoMap"<<std::endl;
	std::cout<<"there are \t"<< overlaidHistoMap.size() <<"\t elements in overlaidHistoMap"<<std::endl;
#endif

	///now add all TH1F objects in stackedHistoMap into one THStack object
	///1 = black for colors
	int colors[] = {2,4,5,8,12,30,40,45};
	vector<int> colorVect(colors,colors + sizeof(colors)/sizeof(int) );
	Int_t i=0;
	THStack * histoStack = new THStack("","");
	vector<pair<Double_t,string> > vectOfBkgnds;	///< use this to link the # of rescaled bkgnd evts to each bkgnd source
	pair<Double_t,string> newPair;
	if(stackedHistoMap.size() > colorVect.size() ) cout<<"not enough unique colors in MultipleCurveOverlayHisto fxn!"<<endl;
#ifdef DEBUG
	std::cout<<"made a THStack object with null name and title"<<std::endl;
#endif

	Double_t bkgndIntegral;	///< use this to track the integral of all bkgnd histos in the THStack object
	for(map<string,TH1F*>::const_iterator histIt=stackedHistoMap.begin(); histIt!=stackedHistoMap.end(); histIt++){
		
		///set the line and fill color, and normalize the MC histos
		(histIt->second)->SetLineColor(colorVect[i]);
		(histIt->second)->SetFillColor(colorVect[i]);
#ifdef DEBUG
		std::cout<<"set line and fill color for the histo \t"<< histIt->first <<std::endl;
		std::cout<<"before rescaling there are "<< (histIt->second)->Integral() <<" entries in the histo mentioned immediately above"<<std::endl;
#endif

		for(map<string,Float_t>::const_iterator mcIt=nEvtsMap.begin(); mcIt!=nEvtsMap.end(); mcIt++){
			if((histIt->first).find(mcIt->first) != string::npos) (histIt->second)->Scale(crossSxnsMap[mcIt->first]*intLumi/nEvtsMap[mcIt->first]);
		}///end loop which normalizes the MC histos to integrated lumi of real data

		newPair = make_pair((histIt->second)->Integral(), histIt->first);
		addPairToOrderedVector(vectOfBkgnds, newPair);
		bkgndIntegral += (histIt->second)->Integral();
#ifdef DEBUG
		std::cout<<"rescaled the histo \t"<< histIt->first <<std::endl;
		std::cout<<"after rescaling there are "<< (histIt->second)->Integral() <<" entries in the histo mentioned immediately above"<<std::endl;
#endif
	
		//histoStack->Add((histIt->second));	///< add each histo in stackedHistoMap to the THStack object

#ifdef DEBUG
		std::cout<<"added \t"<< histIt->first <<"\t to the stacked histo object"<<std::endl;
#endif
		
		if(histIt==stackedHistoMap.begin()){
			histoStack->SetTitle( (histIt->second)->GetTitle() );
			string xLabel="";
			if((histIt->first).find("Mass") !=string::npos || (histIt->first).find("ptEle") !=string::npos || (histIt->first).find("ptJet") !=string::npos){
				xLabel += "GeV";
			}///end if(histo is plotting a variable with dimension of energy)
#ifdef DEBUG
			std::cout<<"xLabel = \t"<< xLabel <<std::endl;
#endif
			(histIt->second)->GetXaxis()->SetTitle(xLabel.c_str());
		}///end filter to set histogram X axis label
		
		size_t lastChevron = (histIt->first).find_last_of('>');
		size_t underscorePos = (histIt->first).find_first_of("_",lastChevron);
#ifdef DEBUG
		std::cout<<"lastChevron located at \t"<< lastChevron<<std::endl;
		std::cout<<"underscorePos located at \t"<< underscorePos<<std::endl;
#endif
		string legEntryName = (histIt->first).substr(lastChevron+1,underscorePos-lastChevron-1);
#ifdef DEBUG
		cout<<"histo name = \t"<< histIt->first <<endl;
		cout<<"legEntry has name = \t"<< legEntryName << endl;
#endif
		leg->AddEntry(histIt->second,legEntryName.c_str());
		i++;
	}///end loop to set line colors of histos, and add entries to TLegend

	for(unsigned int i=0; i<vectOfBkgnds.size(); i++){
		///calling second on the ith element of vectOfBkgnds returns a string which should be a key in stackedHistoMap
		histoStack->Add( stackedHistoMap[vectOfBkgnds[i].second] );
	}///end loop over pairs in vectOfBkgnds

	if(doLogYaxis){
		Double_t originalMax = histoStack->GetMaximum();
		histoStack->SetMaximum(40*originalMax);
		histoStack->SetMinimum(0.1);
	}
	histoStack->Draw("hist");

#ifdef DEBUG
		std::cout<<"stacked histo has been drawn"<<std::endl;
#endif

	Double_t realDataIntegral;
	string outputFile;
	for(map<string,TH1F*>::const_iterator hIt=overlaidHistoMap.begin(); hIt!=overlaidHistoMap.end(); hIt++){
		size_t lastChevron = (hIt->first).find_last_of('>');
		size_t underscorePos = (hIt->first).find_first_of("_",lastChevron);
		string legEntryName = (hIt->first).substr(lastChevron+1,underscorePos-lastChevron-1);
		leg->AddEntry(hIt->second,legEntryName.c_str(),"ep");
		realDataIntegral += (hIt->second)->Integral();
		(hIt->second)->SetMarkerColor(1);
		(hIt->second)->SetMarkerStyle(21);
		(hIt->second)->SetMarkerSize(1.1);
		(hIt->second)->Draw("same ep");	///< draw the real data histogram as points
		size_t firstChevron = (hIt->first).find_first_of('>');
		outputFile = (hIt->first).substr(0,firstChevron);
		outputFile += ".png";
#ifdef DEBUG
		cout<<"outputFile = \t"<< outputFile <<endl;
#endif

		/*
		if(hIt==overlaidHistoMap.begin()){
			(hIt->second)->Draw();
			size_t firstChevron = (hIt->first).find_first_of('>');
			outputFile = (hIt->first).substr(0,firstChevron);
			outputFile += ".png";
#ifdef DEBUG
			cout<<"outputFile = \t"<< outputFile <<endl;
#endif
		}
		else (hIt->second)->Draw("same");
		*/
	}///end loop which draws histograms
	leg->Draw();
	canv->SaveAs(outputFile.c_str(),"recreate");
	cout<<"the integral of the real data histo equals "<< realDataIntegral <<endl;
	cout<<"the integral of all bkgnd histos equals    "<< bkgndIntegral <<endl;
	cout<<"data/MC = "<< realDataIntegral/bkgndIntegral <<endl;

}///end overlayPointsOnStackedHistos


/**
 * use this fxn to overlay multiple curves onto one TCanvas 
 * the TString key in inputChainMap contains the histogram plotting argument to use with TChain->Draw("plottingArgs")
 * the histogram name does not need to be passed into the fxn as an argument, it will be pulled from the map key
 * if doNormalizationByArea is true, then normalize the plotted histos so that the area under each curve = 1
 * title will be shown at the top of the plot
 * xLabel will be shown just below the horizontal axis
 *
 *
 */
void makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea,string title,string xLabel,string outputFileNameModifier,Bool_t specialGrouping,Float_t cutVal){
	gStyle->SetOptStat("");
	TCanvas * canv = new TCanvas(canvName,canvName,750,700);
	canv->cd();
	TLegend * leg = new TLegend(legXmin,legYmin,legXmax,legYmax);
	map<string,TH1F*> overlayHistoMap;	///< links string keys to TH1F histos which will ultimately be overlaid
	for(map<string,TChain*>::const_iterator chMapIt=inputChainMap.begin(); chMapIt!=inputChainMap.end(); chMapIt++){
		size_t openParenth = (chMapIt->first).find_first_of('('), lastChevron = ((chMapIt->first).find_first_of('>')+1);
		string uncutHistoName(chMapIt->first);
		///now initialize a new string, get rid of the content in uncutHistoName before '>>' and after '(',
		///and store the substring in the new string object
		string oneHistoName( uncutHistoName.substr(lastChevron+1,openParenth-lastChevron-1) );
		if( (chMapIt->first).find_first_of('(') == (chMapIt->first).find_last_of('(') ) (chMapIt->second)->Draw((chMapIt->first).c_str());
		else{	///a cut should be applied in tree Draw
			size_t lastOpenParenth = (chMapIt->first).find_last_of('(');
			size_t lastClosedParenth = (chMapIt->first).find_last_of(')');
			string cut( uncutHistoName.substr(lastOpenParenth+1, lastClosedParenth-lastOpenParenth-1) );
			string drawArg( uncutHistoName.substr(0,lastOpenParenth-1) );
			//cout<<"drawArg=\t"<< drawArg << endl;
			//cout<<"cut=\t"<< cut << endl;
			//cout<<"oneHistoName=\t"<< oneHistoName <<endl;
			(chMapIt->second)->Draw(drawArg.c_str(),cut.c_str());

		}
		
		///save pointers to all histograms into overlayHistoMap
		overlayHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
	}///end loop over elements in inputChainMap
	//canv->Update();

	//cout<<"left first loop over elements in input chain"<<endl;
	///now overlay all TH1F objects in overlayHistoMap onto one TCanvas
	int colors[] = {1,2,4,8,12,25,30,40,45};
	vector<int> colorVect(colors,colors + sizeof(colors)/sizeof(int) );
	//short colorsSpecial[] = {kBlack,kBlack,kBlack,kRed,kRed,kRed,kBlue,kBlue,kBlue};
	short colorsSpecial[] = {kRed,kBlack,kBlue};	//standard
	//short colorsSpecial[] = {kBlack,kBlue,kRed};	//temp for fixes
	vector<short> colorSpecialVect(colorsSpecial,colorsSpecial + sizeof(colorsSpecial)/sizeof(short) );
	//int lineStyleSpecial[] = {2,1,6,1,6,2,1,2,6};	///1 is solid, 2 is small dash, 3 is dotted, 7 is medium dash, 9 is large dash, 10 is large dash dot
	int lineStyleSpecial[] = {1,1,1};	///1 is solid, 2 is small dash, 3 is dotted, 7 is medium dash, 9 is large dash, 10 is large dash dot
	vector<int> lineStyleSpecialVect(lineStyleSpecial,lineStyleSpecial + sizeof(lineStyleSpecial)/sizeof(int));
	Int_t i=0;
	typedef map<string,TH1F*>::const_iterator cMapIt;
	Double_t oldMax = 0;	///<height of tallest peak across all histos in overlayHistoMap
	if(overlayHistoMap.size() > colorVect.size() ) cout<<"not enough unique colors in MultipleCurveOverlayHisto fxn!"<<endl;
	for(cMapIt histIt=overlayHistoMap.begin(); histIt!=overlayHistoMap.end(); histIt++){
		(histIt->second)->SetLineWidth(3);
		if(!specialGrouping) (histIt->second)->SetLineColor(colorVect[i]);
		else{
			cout<<"using colorSpecialVect"<<endl;
			cout<<"NOTE colorSpecialVect and lineStyleSpecialVect are current setup for 3 groups, each with 3 elements"<<endl;
			(histIt->second)->SetLineColor(colorSpecialVect[i]);
			(histIt->second)->SetLineStyle(lineStyleSpecialVect[i]);
		}
		if(doNormalizationByArea){
			Double_t oldIntegral = (histIt->second)->Integral();
			(histIt->second)->Scale(1.0/oldIntegral);
		}///end filter to normalize histo area to 1

		///add the x axis label and title to all histos
		(histIt->second)->GetXaxis()->SetTitle(xLabel.c_str());
		(histIt->second)->SetTitle(title.c_str());
		Double_t tempMax = (histIt->second)->GetBinContent((histIt->second)->GetMaximumBin());
		if(oldMax < tempMax) oldMax = tempMax;
	
		size_t lastChevron = (histIt->first).find_first_of('>')+1;
		size_t underscorePos = (histIt->first).find_first_of("_",lastChevron);
		string legEntryName = (histIt->first).substr(lastChevron+1,underscorePos-lastChevron-1);
#ifdef DEBUG
		cout<<"histo name = \t"<< histIt->first <<endl;
		cout<<"legEntry has name = \t"<< legEntryName << endl;
#endif
		leg->AddEntry(histIt->second,legEntryName.c_str());
		i++;
	}///end loop to set line colors of histos, and add entries to TLegend

	///set the max Y scale on all histos in overlayHistoMap to 1.5 times oldMax
	for(cMapIt mIt=overlayHistoMap.begin(); mIt!=overlayHistoMap.end(); mIt++){
		(mIt->second)->SetMaximum(1.5*oldMax);
		if(doNormalizationByArea){
			(mIt->second)->GetYaxis()->SetTitle("Arbitrary Units");
			(mIt->second)->GetYaxis()->SetTitleOffset(1.4);
		}
	}///end loop over elements in overlayHistoMap
	
	string outputFile, outputFilePdf;
	for(cMapIt hIt=overlayHistoMap.begin(); hIt!=overlayHistoMap.end(); hIt++){
		if(hIt==overlayHistoMap.begin()){
			(hIt->second)->Draw();
			size_t firstChevron = (hIt->first).find_first_of('>');
			outputFile = (hIt->first).substr(0,firstChevron);
			size_t badFwdSlash = outputFile.find_first_of('/');
			if(badFwdSlash != string::npos) outputFile.replace(badFwdSlash,1,"Over");
			outputFile += outputFileNameModifier;
			outputFilePdf += outputFile;
			outputFile += ".png";
			outputFilePdf += ".pdf";
#ifdef DEBUG
			cout<<"outputFile = \t"<< outputFile <<endl;
#endif
		}
		else (hIt->second)->Draw("same");
	}///end loop which draws histograms
	leg->SetTextSize(0.027);
	leg->Draw();
	canv->Update();
	canv->SaveAs(outputFile.c_str(),"recreate");
	canv->SaveAs(outputFilePdf.c_str(),"recreate");
	//canv->Clear();
	canv->Close();
	delete canv;
	delete leg;
	for(cMapIt hIt=overlayHistoMap.begin(); hIt!=overlayHistoMap.end(); hIt++){
		delete (hIt->second);
	}///end loop which draws histograms
	overlayHistoMap.clear();
	
}///end makeAndSaveMultipleCurveOverlayHisto()


/*use this fxn to calculate (reco pT)/(matched gen pT) as a function of gen pT, eta, or phi
 * and plot this curve on a histogram with a red dot at the center of each bin
 * the histo is saved to a png file
 * arrSize is the size of the array which is stored in the tuple 
 * arrIndex is the specific index of the array which should be plotted in the TProfile
 * 
 * the TProfile shows (entries in recoBranchName/entries in matchedGenBranchName) as a fxn of entries in otherMatchedGenBranchName
 * if doTGraph is true, then turn the TProfile into a TGraphErrors and save it
 * 
 * if checkSeparatedBins = true, then this fxn should make and save an image of two histos overlaid onto each other
 * each histo shows all values of reco pT/gen pT for one particular gen pT bin 
 */
void calcAndPlotPtRatioProfile(TChain * chain, Int_t arrSize, Int_t arrIndex, TString profileName, TString profileTitle, Int_t nProfileBins, Float_t minX, Float_t maxX, Float_t minY, Float_t maxY, TString matchedGenBranchName, TString recoBranchName, TString otherMatchedGenBranchName,TString canvName, TString drawOptions, TString xAxisLabel, TString outputFile,bool doTGraph,bool checkSeparatedBins){
	if(checkSeparatedBins && doTGraph){
		cout<<"check input bools to calcAndPlotPtRatioProfile fxn"<<endl;
		return;
	}
	TCanvas * g0 = new TCanvas(canvName,canvName,600,600);
	g0->cd();
	gStyle->SetOptStat("emriou");

	///declare two histos in case checkSeparatedBins is true
	TString lowPtHistName = "";
	lowPtHistName += "lowPt";
	string temp(profileTitle);
	string lowPtHistTitle = temp.substr(0,temp.find("vs")-1);
	//TString lowPtHistTitle = "";
	//lowPtHistTitle += profileTitle;
	TString highPtHistName = "";
	highPtHistName += "highPt";
	Int_t nHistoBins = 15;

	TH1F * lowPtRatioHist = new TH1F(lowPtHistName,lowPtHistTitle.c_str(),nHistoBins,0.4,1.3);
	TH1F * highPtRatioHist = new TH1F(highPtHistName,"",nHistoBins,0.4,1.3);

	///declare two arrays which will store info from the TChain
	Float_t recoArr[arrSize], matchedGenArr[arrSize], otherMatchedGenArr[arrSize];
	
	///initialize the TProfile object, and use the "s" option to set the error bar size equal to the RMS of each bin
	TProfile * profilePtr = new TProfile(profileName, profileTitle, nProfileBins, minX, maxX, minY, maxY,"s");
	
	///link the two arrays to the TChain branches which should be read out
	chain->SetBranchAddress(recoBranchName, &recoArr);
	chain->SetBranchAddress(matchedGenBranchName, &matchedGenArr);
	///do not use otherMatchedGenArr if otherMatchedGenBranchName is identical to matchedGenBranchName 
	if(otherMatchedGenBranchName.CompareTo(matchedGenBranchName)!=0) chain->SetBranchAddress(otherMatchedGenBranchName, &otherMatchedGenArr);
	Long64_t nEntries = chain->GetEntries();
	for(Long64_t i=0;i<nEntries; i++){
		chain->GetEntry(i);
		if(otherMatchedGenBranchName.CompareTo(matchedGenBranchName)!=0){
			profilePtr->Fill(otherMatchedGenArr[arrIndex],recoArr[arrIndex]/matchedGenArr[arrIndex]);
		}
		else{
			///this if statement decides what is executed when the x axis is gen pT
			profilePtr->Fill(matchedGenArr[arrIndex],recoArr[arrIndex]/matchedGenArr[arrIndex]);
			if(matchedGenArr[arrIndex] > 200 && matchedGenArr[arrIndex] < 250) lowPtRatioHist->Fill(recoArr[arrIndex]/matchedGenArr[arrIndex]);
			if(matchedGenArr[arrIndex] > 500 && matchedGenArr[arrIndex] < 550) highPtRatioHist->Fill(recoArr[arrIndex]/matchedGenArr[arrIndex]);
		}
	}///end loop over entries in TChain

	///if doTGraph is true, then read the TProfile points and error bars and plot these points on a TGraphErrors plot
	if(doTGraph){
		Double_t xAxisVals[nProfileBins], xAxisErrorVals[nProfileBins], yAxisVals[nProfileBins], yAxisErrorVals[nProfileBins];
		///loop over all of the bins in profilePtr and save the x coordinate, y coordinate, and y axis error bar size
		///from each bin to the appropriate array
		for(Int_t i=1; i<=nProfileBins; i++){
			xAxisVals[i-1]=profilePtr->GetXaxis()->GetBinCenter(i);
			xAxisErrorVals[i-1]=0;
			yAxisVals[i-1]=profilePtr->GetBinContent(i);
			yAxisErrorVals[i-1]=profilePtr->GetBinError(i);
		}//end loop over bins in TProfile which was filled earlier

		///now that the arrays are filled the TGraphErrors object can be made and drawn
		TGraphErrors * ptRatioGrph = new TGraphErrors(nProfileBins, xAxisVals, yAxisVals, xAxisErrorVals, yAxisErrorVals);
		ptRatioGrph->SetMinimum(0);
		ptRatioGrph->SetMaximum(1.1);
		ptRatioGrph->SetMarkerColor(kRed);
		ptRatioGrph->SetMarkerStyle(21);
		ptRatioGrph->GetXaxis()->SetTitle(xAxisLabel);
		ptRatioGrph->SetTitle(profileTitle);
		ptRatioGrph->Draw("AP");
	}//end doTGraph == true

	///if checkSeparatedBins is true, then update the fill and line properties of the two TH1F histos made earlier
	///and overlay the histos 
	if(checkSeparatedBins){
		gStyle->SetOptStat("");

		//rescale the histos so they have the same integral
		lowPtRatioHist->Scale(1/( lowPtRatioHist->Integral() ));
		highPtRatioHist->Scale(1/( highPtRatioHist->Integral() ));
		
		//rescale the y axis so that the max bin doesn't go above the top boundary of the plot 
		Double_t yMax;
		if(lowPtRatioHist->GetBinContent(lowPtRatioHist->GetMaximumBin()) < highPtRatioHist->GetBinContent(highPtRatioHist->GetMaximumBin()) ) yMax = (1.1)*(highPtRatioHist->GetBinContent(highPtRatioHist->GetMaximumBin()));
		else yMax = (1.1)*(lowPtRatioHist->GetBinContent(lowPtRatioHist->GetMaximumBin()));
		lowPtRatioHist->SetMaximum(yMax);
		lowPtRatioHist->GetXaxis()->SetTitle("reco P_{T}/gen P_{T}");

		//rescale the x axis to cluster the bins around the maximum bin
		Double_t centroid = lowPtRatioHist->GetXaxis()->GetBinCenter(lowPtRatioHist->GetMaximumBin() );
		Double_t rms = lowPtRatioHist->GetRMS(1);
		lowPtRatioHist->GetXaxis()->Set(nHistoBins, centroid - 3*rms, centroid + 1.5*rms);
		highPtRatioHist->GetXaxis()->Set(nHistoBins, centroid - 3*rms, centroid + 1.5*rms);
		
		TLegend * leg = new TLegend(0.15,0.7,0.5,0.9);
		lowPtRatioHist->SetFillColor(kRed);
		lowPtRatioHist->SetLineColor(kRed);
		highPtRatioHist->SetLineColor(kBlack);
		highPtRatioHist->SetLineWidth(3);
		leg->AddEntry(lowPtRatioHist,"200 < gen pT < 250");
		leg->AddEntry(highPtRatioHist,"500 < gen pT < 550");
		lowPtRatioHist->Draw();
		highPtRatioHist->Draw("same");
		leg->Draw();
	}//end checkSeparatedBins == true
	
	else{
		///now draw the TProfile and save an image of the plot
		profilePtr->SetMinimum(0);
		profilePtr->SetMarkerColor(kRed);
		profilePtr->GetXaxis()->SetTitle(xAxisLabel);
		profilePtr->Draw(drawOptions);
	}//end doTGraph == false and checkSeparatedBins == false

	g0->SaveAs(outputFile,"recreate");

}///end calcAndPlotPtRatioProfile()

/**use this fxn to plot (reco pT/gen pT) as a fxn of gen pT, eta, or phi
 * in a variable number of bins (specified by an input integer)
 */ 
/*
void fitAndPlotPtRatioStats(map<string,TChain *> inputChainMap, TString fitFunc, T){

}///end fitAndPlotPtRatioStats()
*/

///use this macro to develop and run new macros which don't have a central theme, other than being useful to the WR analysis
///the first use of this macro was to plot reco pT/gen pT for reco jets and leptons matched to GEN counterparts

void macroSandBox(){
	TChain * matchedRecoPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_RECO_WR_signal_and_bkgnd_files/analysis_recoElectronChannel_single_stage_matching_for_jets_centrally_produced_signal_MWr_2600_MNu_1300.root");
	TChain * copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_RECO_WR_signal_and_bkgnd_files/analysis_recoElectronChannel_single_stage_matching_for_jets_centrally_produced_signal_MWr_2600_MNu_1300.root");

	string dyPlusJetsKey = "DYJets";
	string ttBarKey = "TTBar";
	string wJetsKey = "WJets";
	string wzKey = "WZ";
	string zzKey = "ZZ";
	
	///cross sxn values (picobarns) and dataset sizes (=total number of evts) for bkgnd processes
	map<string,Float_t> xSxnsFiftyNs;
	xSxnsFiftyNs[ttBarKey]=831.76;
	xSxnsFiftyNs[dyPlusJetsKey]=6025.2;
	xSxnsFiftyNs[wJetsKey]=61500;
	xSxnsFiftyNs[wzKey]=66.1;
	xSxnsFiftyNs[zzKey]=15.4;
	map<string,Float_t> numEvtsFiftyNs;
	numEvtsFiftyNs[ttBarKey]=10000;
	numEvtsFiftyNs[dyPlusJetsKey]=10000;	///< madgraph DY+JetsToLL M-50 sample
	numEvtsFiftyNs[wJetsKey]=10000;
	numEvtsFiftyNs[wzKey]=10000;
	numEvtsFiftyNs[zzKey]=10000;
	map<string,Float_t> numEvtsTwentyFiveNs;
	numEvtsTwentyFiveNs[ttBarKey]=19899500;
	numEvtsTwentyFiveNs[dyPlusJetsKey]=9052671;	///< madgraph DY+JetsToLL M-50 sample
	numEvtsTwentyFiveNs[wJetsKey]=72121586;
	numEvtsTwentyFiveNs[wzKey]=991232;
	numEvtsTwentyFiveNs[zzKey]=996168;


#ifdef lowMassFlavorSidebandBkgndOnData
	string modAndTreeName = "recoAnalyzerTwo/recoObjectsWithPtEtaCuts";

	TChain * dyPlusJetsEMuJJLowMassSkim = new TChain(modAndTreeName.c_str());
	dyPlusJetsEMuJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/analyzed_50ns_skims_low_dilepton_mass_emujj/analyzed_DYJets_Madgraph_50ns_skim_low_dilepton_mass_region_emujj.root");
	//dyPlusJetsEMuJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_DYJets_50ns_skim_low_dilepton_mass_region_emujj.root");
	
	TChain * ttBarEMuJJLowMassSkim = new TChain(modAndTreeName.c_str());
	ttBarEMuJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/analyzed_50ns_skims_low_dilepton_mass_emujj/analyzed_TTJets_50ns_skim_low_dilepton_mass_region_emujj.root");
	TChain * wJetsEMuJJLowMassSkim = new TChain(modAndTreeName.c_str());
	wJetsEMuJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/analyzed_50ns_skims_low_dilepton_mass_emujj/analyzed_WJets_50ns_skim_low_dilepton_mass_region_emujj.root");
	TChain * wzEMuJJLowMassSkim = new TChain(modAndTreeName.c_str());
	wzEMuJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/analyzed_50ns_skims_low_dilepton_mass_emujj/analyzed_WZ_50ns_skim_low_dilepton_mass_region_emujj.root");
	TChain * zzEMuJJLowMassSkim = new TChain(modAndTreeName.c_str());
	zzEMuJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/analyzed_50ns_skims_low_dilepton_mass_emujj/analyzed_ZZ_50ns_skim_low_dilepton_mass_region_emujj.root");
	
	TChain * muonEGEMuJJLowMassSkim = new TChain(modAndTreeName.c_str());
	muonEGEMuJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/analyzed_50ns_skims_low_dilepton_mass_emujj/analyzed_MuonEG_50ns_skim_low_dilepton_mass_region_emujj.root");

#ifdef DEBUG
	cout<<"declared TChain to bkgnd and data"<<endl;
#endif

	///compute the PU weights for different bkgnd processes
	//computePileupWeights(map<string,TChain*> bkgndChainMap, TChain * dataChain)
	map<string,TChain*> bkgndMap;
	bkgndMap[ttBarKey]=ttBarEMuJJLowMassSkim;
	bkgndMap[dyPlusJetsKey]=dyPlusJetsEMuJJLowMassSkim;
	bkgndMap[wJetsKey]=wJetsEMuJJLowMassSkim;
	bkgndMap[wzKey]=wzEMuJJLowMassSkim;
	bkgndMap[zzKey]=zzEMuJJLowMassSkim;
	
	map<string,vector<Float_t> > pileupWeights = computePileupWeights(bkgndMap, muonEGEMuJJLowMassSkim);


	Float_t integratedLumi = 41.8;	///< in inverse picobarns


	//string branchNames[] = {"ptEle[0]","ptEle[1]","etaEle[0]","etaEle[1]","ptJet[0]","ptJet[1]","etaJet[0]","etaJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet","leadLeptonThreeObjMass","subleadingLeptonThreeObjMass","dijetMass","nLeptons","nJets","nVertices"};
	string link=">>";
	//string histoEndings[] = {"_leadLeptonPt(25,0.,250.)","_subleadLeptonPt(14,0.,140.)","_leadLeptonEta(30,-3.0,3.0)","_subleadLeptonEta(30,-3.0,3.0)","_leadJetPt(25,0.,250.)","_subleadJetPt(14,0.,140.)","_leadJetEta(30,-3.0,3.0)","_subleadJetEta(30,-3.0,3.0)","_dileptonMass(20,0.,200.)","_fourObjectMass(30,0.,600.)","_dR_leadingLeptonLeadingJet(25,0.,5.)","_dR_leadingLeptonSubleadingJet(25,0.,5.)","_dR_subleadingLeptonLeadingJet(25,0.,5.)","_dR_subleadingLeptonSubleadingJet(25,0.,5.)","_dR_leadingLeptonSubleadingLepton(25,0.,5.)","_dR_leadingJetSubleadingJet(25,0.,5.)","_leadLeptonThreeObjMass(25,0.,500.)","_subleadingLeptonThreeObjMass(25,0.,500.)","_dijetMass(40,0.,400)","_nLeptons(4,0.,4.)","_nJets(6,0.,6.)","_nVertices(35,0.,35.)"};
	
	//string branchNames[] = {"ptEle[0]","ptEle[1]","ptJet[0]","ptJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet"};
	//string histoEndings[] = {"_leadLeptonPt(20,30.,200.)","_subleadLeptonPt(20,30.,110.)","_leadJetPt(25,30.,200.)","_subleadJetPt(20,30.,110.)","_dileptonMass(25,0.,200.)","_fourObjectMass(20,0.,600.)","_dR_leadingLeptonLeadingJet(30,0.,5.)","_dR_leadingLeptonSubleadingJet(30,0.,5.)","_dR_subleadingLeptonLeadingJet(30,0.,5.)","_dR_subleadingLeptonSubleadingJet(30,0.,5.)","_dR_leadingLeptonSubleadingLepton(30,0.,5.)","_dR_leadingJetSubleadingJet(30,0.,5.)"};
	
	//string branchNames[] = {"leadLeptonThreeObjMass","subleadingLeptonThreeObjMass"};
	//string histoEndings[] = {"_leadLeptonThreeObjMass(25,50.,600.)","_subleadingLeptonThreeObjMass(25,50.,600.)"};
	
	//string branchNames[] = {"etaEle[0]","etaEle[1]","etaJet[0]","etaJet[1]"};
	//string histoEndings[] = {"_leadLeptonEta(15,-3.0,3.0)","_subleadLeptonEta(15,-3.0,3.0)","_leadJetEta(15,-3.0,3.0)","_subleadJetEta(15,-3.0,3.0)"};
	
	//string branchNames[] = {"nJets","nLeptons","nVertices"};
	//string histoEndings[] = {"_nJets(12,0.,12.)","_nLeptons(6,0.,6.)","_nVertices(25,0.,50.)"};
	
	//string branchNames[] = {"nLeptons","etaEle[0]","etaEle[1]"};
	//string histoEndings[] = {"_nLeptons(6,0.,6.)","_leadLeptonEta(30,-2.6,2.6)","_subleadLeptonEta(30,-2.6,2.6)"};
	

	//string branchNames[] = {"etaEle[0]","etaEle[1]","ptEle[0]","ptEle[1]","dileptonMass","phiEle[0]","phiEle[1]"};
	//string histoEndings[] = {"_leadLeptonEta(50,-2.5,2.5)","_subleadLeptonEta(50,-2.5,2.5)","_leadLeptonPt(40,0.,200.)","_subleadLeptonPt(20,0.,100.)","_dileptonMass(200,50.,250.)","_leadLeptonPhi(64,-3.2,3.2)","_subleadLeptonPhi(64,-3.2,3.2)"};
	
	string branchNames[] = {"dileptonMass"};
	string histoEndings[] = {"_dileptonMass(20,20.,200.)"};
	

	TString evWeightCut = "(evWeightSign < 0 ? -1. : 1.)";

	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {zzKey,wzKey,wJetsKey,ttBarKey,dyPlusJetsKey,"Data"};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[5]+histoEndings[i]] = muonEGEMuJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[4]+histoEndings[i]] = dyPlusJetsEMuJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[3]+histoEndings[i]] = ttBarEMuJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = wJetsEMuJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = wzEMuJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = zzEMuJJLowMassSkim;
		string cName = "o"+to_string(i);
		overlayPointsOnStackedHistos(placeHolderMap,cName.c_str(),0.5,0.75,0.9,0.89,xSxnsFiftyNs,integratedLumi,numEvtsFiftyNs,evWeightCut,false,pileupWeights,true);
		placeHolderMap.clear();
	}///end loop over branchNames



#endif
///end lowMassFlavorSidebandBkgndOnData



#ifdef lowMassSkimmedBkgndOnRealData
	//overlayPointsOnStackedHistos(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,map<string,Float_t> crossSxnsMap,Float_t intLumi,map<string,Float_t> nEvtsMap, TString treeCuts, Bool_t doPuReweighting, map<string,vector<Float_t> > puWeights,Bool_t doLogYaxis)
	//recoAnalyzerOne/zEleEleNoCuts or recoAnalyzerTwo/recoObjectsWithPtEtaCuts
	string moduleAndTreeName = "recoAnalyzerTwo/recoObjectsWithPtEtaCuts";
	//string moduleAndTreeName = "recoAnalyzerOne/zEleEleNoCuts";

	//skim_check_Zee_peak_twoHEEP_noHLT_Aug25 or skim_low_mass_region_eejj
	TChain * dyPlusJetsEEJJLowMassSkim = new TChain(moduleAndTreeName.c_str());
	dyPlusJetsEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/analyzed_DYJets_Madgraph_50ns_skim_low_mass_region_eejj.root");
	//dyPlusJetsEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_DYJets_50ns_skim_low_mass_region_eejj.root");
	
	TChain * ttBarEEJJLowMassSkim = new TChain(moduleAndTreeName.c_str());
	ttBarEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_TTJets_50ns_skim_low_mass_region_eejj.root");
	TChain * wJetsEEJJLowMassSkim = new TChain(moduleAndTreeName.c_str());
	wJetsEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_WJets_50ns_skim_low_mass_region_eejj.root");
	TChain * wzEEJJLowMassSkim = new TChain(moduleAndTreeName.c_str());
	wzEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/WZ_TuneCUETP8M1_13TeV-pythia8/analyzed_WZ_50ns_skim_low_mass_region_eejj.root");
	TChain * zzEEJJLowMassSkim = new TChain(moduleAndTreeName.c_str());
	zzEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/ZZ_TuneCUETP8M1_13TeV-pythia8/analyzed_ZZ_50ns_skim_low_mass_region_eejj.root");
	
	TChain * doubleEGEEJJLowMassSkim = new TChain(moduleAndTreeName.c_str());
	doubleEGEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/DoubleEG/analyzed_DoubleEG_50ns_skim_low_mass_region_eejj.root");

#ifdef DEBUG
	cout<<"declared TChain to bkgnd and data"<<endl;
#endif

	///compute the PU weights for different bkgnd processes
	//computePileupWeights(map<string,TChain*> bkgndChainMap, TChain * dataChain)
	map<string,TChain*> bkgndMap;
	bkgndMap[ttBarKey]=ttBarEEJJLowMassSkim;
	bkgndMap[dyPlusJetsKey]=dyPlusJetsEEJJLowMassSkim;
	bkgndMap[wJetsKey]=wJetsEEJJLowMassSkim;
	bkgndMap[wzKey]=wzEEJJLowMassSkim;
	bkgndMap[zzKey]=zzEEJJLowMassSkim;
	
	map<string,vector<Float_t> > pileupWeights = computePileupWeights(bkgndMap, doubleEGEEJJLowMassSkim);


	Float_t integratedLumi = 41.8;	///< in picobarns


	//string branchNames[] = {"ptEle[0]","ptEle[1]","etaEle[0]","etaEle[1]","ptJet[0]","ptJet[1]","etaJet[0]","etaJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet","leadLeptonThreeObjMass","subleadingLeptonThreeObjMass","dijetMass","nLeptons","nJets","nVertices"};
	string link=">>";
	//string histoEndings[] = {"_leadLeptonPt(25,0.,250.)","_subleadLeptonPt(14,0.,140.)","_leadLeptonEta(30,-3.0,3.0)","_subleadLeptonEta(30,-3.0,3.0)","_leadJetPt(25,0.,250.)","_subleadJetPt(14,0.,140.)","_leadJetEta(30,-3.0,3.0)","_subleadJetEta(30,-3.0,3.0)","_dileptonMass(20,0.,200.)","_fourObjectMass(30,0.,600.)","_dR_leadingLeptonLeadingJet(25,0.,5.)","_dR_leadingLeptonSubleadingJet(25,0.,5.)","_dR_subleadingLeptonLeadingJet(25,0.,5.)","_dR_subleadingLeptonSubleadingJet(25,0.,5.)","_dR_leadingLeptonSubleadingLepton(25,0.,5.)","_dR_leadingJetSubleadingJet(25,0.,5.)","_leadLeptonThreeObjMass(25,0.,500.)","_subleadingLeptonThreeObjMass(25,0.,500.)","_dijetMass(40,0.,400)","_nLeptons(4,0.,4.)","_nJets(6,0.,6.)","_nVertices(35,0.,35.)"};
	
	//string branchNames[] = {"ptEle[0]","ptEle[1]","ptJet[0]","ptJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet"};
	//string histoEndings[] = {"_leadLeptonPt(20,30.,200.)","_subleadLeptonPt(20,30.,110.)","_leadJetPt(25,30.,200.)","_subleadJetPt(20,30.,110.)","_dileptonMass(25,0.,200.)","_fourObjectMass(20,0.,600.)","_dR_leadingLeptonLeadingJet(30,0.,5.)","_dR_leadingLeptonSubleadingJet(30,0.,5.)","_dR_subleadingLeptonLeadingJet(30,0.,5.)","_dR_subleadingLeptonSubleadingJet(30,0.,5.)","_dR_leadingLeptonSubleadingLepton(30,0.,5.)","_dR_leadingJetSubleadingJet(30,0.,5.)"};
	
	//string branchNames[] = {"leadLeptonThreeObjMass","subleadingLeptonThreeObjMass"};
	//string histoEndings[] = {"_leadLeptonThreeObjMass(25,50.,600.)","_subleadingLeptonThreeObjMass(25,50.,600.)"};
	
	//string branchNames[] = {"etaEle[0]","etaEle[1]","etaJet[0]","etaJet[1]"};
	//string histoEndings[] = {"_leadLeptonEta(15,-3.0,3.0)","_subleadLeptonEta(15,-3.0,3.0)","_leadJetEta(15,-3.0,3.0)","_subleadJetEta(15,-3.0,3.0)"};
	
	//string branchNames[] = {"nJets","nLeptons","nVertices"};
	//string histoEndings[] = {"_nJets(12,0.,12.)","_nLeptons(6,0.,6.)","_nVertices(25,0.,50.)"};
	
	//string branchNames[] = {"nLeptons","etaEle[0]","etaEle[1]"};
	//string histoEndings[] = {"_nLeptons(6,0.,6.)","_leadLeptonEta(30,-2.6,2.6)","_subleadLeptonEta(30,-2.6,2.6)"};
	

	//string branchNames[] = {"etaEle[0]","etaEle[1]","ptEle[0]","ptEle[1]","dileptonMass","phiEle[0]","phiEle[1]"};
	//string histoEndings[] = {"_leadLeptonEta(50,-2.5,2.5)","_subleadLeptonEta(50,-2.5,2.5)","_leadLeptonPt(40,0.,200.)","_subleadLeptonPt(20,0.,100.)","_dileptonMass(200,50.,250.)","_leadLeptonPhi(64,-3.2,3.2)","_subleadLeptonPhi(64,-3.2,3.2)"};
	
	string branchNames[] = {"dileptonMass"};
	string histoEndings[] = {"_dileptonMass(200,50.,250.)"};
	

	TString evWeightCut = "(evWeightSign < 0 ? -1. : 1.)";

	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	//string histoBeginnings[] = {"Data",ttBarKey,dyPlusJetsKey,wJetsKey,wzKey,zzKey};
	string histoBeginnings[] = {zzKey,wzKey,wJetsKey,ttBarKey,dyPlusJetsKey,"Data"};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[5]+histoEndings[i]] = doubleEGEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[4]+histoEndings[i]] = dyPlusJetsEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[3]+histoEndings[i]] = ttBarEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = wJetsEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = wzEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = zzEEJJLowMassSkim;
		string cName = "o"+to_string(i);
		overlayPointsOnStackedHistos(placeHolderMap,cName.c_str(),0.5,0.75,0.9,0.89,xSxnsFiftyNs,integratedLumi,numEvtsFiftyNs,evWeightCut,false,pileupWeights,true);
		placeHolderMap.clear();
	}///end loop over branchNames


	
#endif

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///take all values of reco pT/matched gen pT from two well separated gen pT bins, and plot all of these
	///values on two histograms, and overlay the two histos
	///draw the first histo with a fill color, and the second with no fill color and different color, thicker line 
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef checkWellSeparatedGenPtBins
	///leading lepton
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingLeptonPtRatioVsGenPt","leading electron P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1100.,0.,1.1,"ptMatchingEle","ptEle","ptMatchingEle","leadingLeptonGenPt","*H","P_{T,GEN} (GeV)","leadingLeptonPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",false,true);

	///subsubleading lepton
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingLeptonPtRatioVsGenPt","subleading electron P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1000.,0.,1.1,"ptMatchingEle","ptEle","ptMatchingEle","subleadingLeptonGenPt","*H","P_{T,GEN} (GeV)","subleadingLeptonPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",false,true);

	///leading jet
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingJetPtRatioVsGenPt","leading jet P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1000.,0.,1.1,"ptMatchingJet","ptJet","ptMatchingJet","leadingJetGenPt","*H","P_{T,GEN} (GeV)","leadingJetPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",false,true);

	///subleading jet
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingJetPtRatioVsGenPt","subleading jet P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,600.,0.,1.1,"ptMatchingJet","ptJet","ptMatchingJet","subleadingJetGenPt","*H","P_{T,GEN} (GeV)","subleadingJetPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",false,true);

#endif
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///plot reco pT/matched gen pT for both leptons and jets as a function of gen pT, eta, and phi
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef PtRatioProfiles

	/*
	for(Int_t i=0; i<2; i++){

	}///end loop over bins of gen eta
	matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Draw("etaMatchingEle[0]>>wrDauGenEleEta","abs(etaMatchingEle[0])<0.5");
	matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Fit("gaus","ptEle[0]/ptMatchingEle[0]>>wrDauElePtRatioVsGenEta(20,0.5,1.4)","abs(etaMatchingEle[0])<0.5");
	TH1F * hist = (TH1F*) gROOT->FindObject("wrDauElePtRatioVsGenEta");
	TH1F * genHist = (TH1F*) gROOT->FindObject("wrDauGenEleEta");
	TF1 * fit = hist->GetFunction("gaus");
	Double_t mean[10], rms[10], peak[10], xAxis[10];
	Int_t n=1;
	mean[0] = fit->GetParameter(1);
	rms[0] = fit->GetParameter(2);
	peak[0] = hist->GetXaxis()->GetBinCenter( hist->GetMaximumBin() );
	xAxis[0] = genHist->GetMean(1);	//gets the mean value of the x axis
	TCanvas * c0 = new TCanvas("c0","c0",600,600);
	TCanvas * c1 = new TCanvas("c1","c1",600,600);
	TCanvas * c2 = new TCanvas("c2","c2",600,600);
	c0->cd();
	TGraph * meanGr = new TGraph(n,xAxis,mean);
	meanGr->SetTitle("mean reco P_{T}/gen P_{T} vs gen #eta for WR daugther electron");
	meanGr->GetXaxis()->SetTitle("gen #eta");
	meanGr->GetYaxis()->SetTitle("mean reco P_{T}/gen P_{T}");
	meanGr->Draw("AC*");

	c1->cd();
	TGraph * rmsGr = new TGraph(n,xAxis,rms);
	rmsGr->GetXaxis()->SetTitle("gen #eta");
	rmsGr->GetYaxis()->SetTitle("rms reco P_{T}/gen P_{T}");
	
	rmsGr->Draw("AC*");

	c2->cd();
	TGraph * peakGr = new TGraph(n,xAxis,peak);
	peakGr->GetXaxis()->SetTitle("gen #eta");
	peakGr->GetYaxis()->SetTitle("peak reco P_{T}/gen P_{T}");
	peakGr->Draw("AC*");
	*/
	
	///leading lepton
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingLeptonPtRatioVsGenPt","leading electron P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1100.,0.,1.1,"ptMatchingEle","ptEle","ptMatchingEle","leadingLeptonGenPt","*H","P_{T,GEN} (GeV)","leadingLeptonPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingLeptonPtRatioVsGenEta","leading electron P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.0,3.0,0.,1.1,"ptMatchingEle","ptEle","etaMatchingEle","leadingLeptonGenEta","*H","#eta_{GEN}","leadingLeptonPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingLeptonPtRatioVsGenPhi","leading electron P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingEle","ptEle","phiMatchingEle","leadingLeptonGenPhi","*H","#phi_{GEN}","leadingLeptonPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	///subsubleading lepton
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingLeptonPtRatioVsGenPt","subleading electron P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1000.,0.,1.1,"ptMatchingEle","ptEle","ptMatchingEle","subleadingLeptonGenPt","*H","P_{T,GEN} (GeV)","subleadingLeptonPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingLeptonPtRatioVsGenEta","subleading electron P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.2,3.2,0.,1.1,"ptMatchingEle","ptEle","etaMatchingEle","subleadingLeptonGenEta","*H","#eta_{GEN}","subleadingLeptonPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingLeptonPtRatioVsGenPhi","subleading electron P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingEle","ptEle","phiMatchingEle","subleadingLeptonGenPhi","*H","#phi_{GEN}","subleadingLeptonPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png",true,false);


	///leading jet
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingJetPtRatioVsGenPt","leading jet P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1000.,0.,1.1,"ptMatchingJet","ptJet","ptMatchingJet","leadingJetGenPt","*H","P_{T,GEN} (GeV)","leadingJetPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingJetPtRatioVsGenEta","leading jet P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.0,3.0,0.,1.1,"ptMatchingJet","ptJet","etaMatchingJet","leadingJetGenEta","*H","#eta_{GEN}","leadingJetPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingJetPtRatioVsGenPhi","leading jet P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingJet","ptJet","phiMatchingJet","leadingJetGenPhi","*H","#phi_{GEN}","leadingJetPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png",true,false);


	///subleading jet
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingJetPtRatioVsGenPt","subleading jet P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,600.,0.,1.1,"ptMatchingJet","ptJet","ptMatchingJet","subleadingJetGenPt","*H","P_{T,GEN} (GeV)","subleadingJetPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingJetPtRatioVsGenEta","subleading jet P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.2,3.2,0.,1.1,"ptMatchingJet","ptJet","etaMatchingJet","subleadingJetGenEta","*H","#eta_{GEN}","subleadingJetPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png",true,false);

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingJetPtRatioVsGenPhi","subleading jet P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingJet","ptJet","phiMatchingJet","subleadingJetGenPhi","*H","#phi_{GEN}","subleadingJetPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png",true,false);

#endif



	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///plot pT, eta, and phi for matched reco objects (leptons and jets) on top of the gen objects which they are matched to
	///the gen matching is done using GEN objects (leptons, jets) before any GEN cuts

#ifdef RecoGenOverlays
	//makeAndSaveOverlayHisto(TChain * sigChain,TChain * bkgndChain,TString sigHistPlotArg,TString bkgndHistPlotArg,TString sigHistName,TString bkgndHistName,TString histTitle,TString xAxisTitle,TString canvName,TCut sigFilt,TCut bkgndFilt,TString outputFile,Bool_t isPlottingEnergy,Bool_t isPlottingInverseEnergy,Bool_t normalizeHistos,TString sigLegendLabel,TString bkgndLegendLabel,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax)

	///single lepton pT and eta
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"ptEle[0]>>leadingLeptonPtHist(50,0.,1200.)","ptMatchingEle[0]>>leadingMatchedGenLeptonPtHist(50,0.,1200.)","leadingLeptonPtHist","leadingMatchedGenLeptonPtHist","Leading lepton P_{T} RECO and GEN","P_{T}","leadingLeptonPt","","","leadingLeptonPt_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.3,0.6,0.5,0.8);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"ptEle[1]>>subleadingLeptonPtHist(50,0.,1000.)","ptMatchingEle[1]>>subleadingMatchedGenLeptonPtHist(50,0.,1000.)","subleadingLeptonPtHist","subleadingMatchedGenLeptonPtHist","Subleading lepton P_{T} RECO and GEN","P_{T}","subleadingLeptonPt","","","subleadingLeptonPt_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.65,0.53,0.85,0.68);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"etaEle[0]>>leadingLeptonEtaHist(50,-3.0,3.0)","etaMatchingEle[0]>>leadingMatchedGenLeptonEtaHist(50,-3.0,3.0)","leadingLeptonEtaHist","leadingMatchedGenLeptonEtaHist","Leading lepton #eta RECO and GEN","#eta","leadingLeptonEta","","","leadingLeptonEta_matched_reco_and_gen_all_reco_cuts_applied.png",false,false,false,"RECO","GEN",0.13,0.6,0.3,0.8);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"etaEle[1]>>subleadingLeptonEtaHist(50,-3.0,3.0)","etaMatchingEle[1]>>subleadingMatchedGenLeptonEtaHist(50,-3.0,3.0)","subleadingLeptonEtaHist","subleadingMatchedGenLeptonEtaHist","Subleading lepton #eta RECO and GEN","#eta","subleadingLeptonEta","","","subleadingLeptonEta_matched_reco_and_gen_all_reco_cuts_applied.png",false,false,false,"RECO","GEN",0.13,0.6,0.3,0.8);


	///single jet pT and eta
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"ptJet[0]>>leadingJetPtHist(50,0.,1200.)","ptMatchingJet[0]>>leadingMatchedGenJetPtHist(50,0.,1200.)","leadingJetPtHist","leadingMatchedGenJetPtHist","Leading jet P_{T} RECO and GEN","P_{T}","leadingJetPt","","","leadingJetPt_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.75,0.53,0.95,0.68);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"ptJet[1]>>subleadingJetPtHist(50,0.,750.)","ptMatchingJet[1]>>subleadingMatchedGenJetPtHist(50,0.,750.)","subleadingJetPtHist","subleadingMatchedGenJetPtHist","Subleading jet P_{T} RECO and GEN","P_{T}","subleadingJetPt","","","subleadingJetPt_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.75,0.53,0.95,0.68);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"etaJet[0]>>leadingJetEtaHist(50,-3.0,3.0)","etaMatchingJet[0]>>leadingMatchedGenJetEtaHist(50,-3.0,3.0)","leadingJetEtaHist","leadingMatchedGenJetEtaHist","Leading jet #eta RECO and GEN","#eta","leadingJetEta","","","leadingJetEta_matched_reco_and_gen_all_reco_cuts_applied.png",false,false,false,"RECO","GEN",0.13,0.6,0.3,0.8);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"etaJet[1]>>subleadingJetEtaHist(50,-3.0,3.0)","etaMatchingJet[1]>>subleadingMatchedGenJetEtaHist(50,-3.0,3.0)","subleadingJetEtaHist","subleadingMatchedGenJetEtaHist","Subleading jet #eta RECO and GEN","#eta","subleadingJetEta","","","subleadingJetEta_matched_reco_and_gen_all_reco_cuts_applied.png",false,false,false,"RECO","GEN",0.13,0.6,0.3,0.8);

	///dilepton and dijet mass, both heavy Nu masses (with leading or subleading lepton), and WR mass
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"dileptonMass>>dileptonMassHist(50,0.,2500.)","dileptonMassMatching>>dileptonMassMatchingHist(50,0.,2500.)","dileptonMassHist","dileptonMassMatchingHist","M_{LL} for RECO and GEN","M_{LL}","dileptonMass","","","dileptonMass_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.75,0.53,0.95,0.68);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"dijetMass>>dijetMassHist(50,0.,1400.)","dijetMassMatching>>dijetMassMatchingHist(50,0.,1400.)","dijetMassHist","dijetMassMatchingHist","M_{JJ} for RECO and GEN","M_{JJ}","dijetMass","","","dijetMass_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.75,0.53,0.95,0.68);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"leadLeptonThreeObjMass>>leadLeptonThreeObjMassHist(50,0.,2800.)","leadLeptonThreeObjMassMatching>>leadLeptonThreeObjMassMatchingHist(50,0.,2800.)","leadLeptonThreeObjMassHist","leadLeptonThreeObjMassMatchingHist","M_{leading L JJ} for RECO and GEN","M_{leading L JJ}","leadLeptonThreeObjMass","","","leadLeptonThreeObjMass_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.15,0.53,0.3,0.68);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"subleadingLeptonThreeObjMass>>subleadingLeptonThreeObjMassHist(50,0.,1500.)","subleadingLeptonThreeObjMassMatching>>subleadingLeptonThreeObjMassMatchingHist(50,0.,1500.)","subleadingLeptonThreeObjMassHist","subleadingLeptonThreeObjMassMatchingHist","M_{subleading L JJ} for RECO and GEN","M_{subleading L JJ}","subleadingLeptonThreeObjMass","","","subleadingLeptonThreeObjMass_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.15,0.53,0.3,0.68);
	makeAndSaveOverlayHisto(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts,"fourObjectMass>>fourObjectMassHist(50,0.,3100.)","fourObjectMassMatching>>fourObjectMassMatchingHist(50,0.,3100.)","fourObjectMassHist","fourObjectMassMatchingHist","M_{LLJJ} for RECO and GEN","M_{LLJJ}","fourObjectMass","","","fourObjectMass_matched_reco_and_gen_all_reco_cuts_applied.png",true,false,false,"RECO","GEN",0.15,0.53,0.3,0.68);


#endif


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///plot leading and subleading electron and jet pT and eta, dilepton mass, mWR, and all dR_ distributions (lead lepton to sublead lepton,
	///lead lepton to lead jet, lead lepton to sublead jet, sublead lepton to lead jet, sublead lepton to sublead jet, lead jet to sublead jet)
	///for different (MWR, MNu) pairs.  All (MWR, MNu) pairs should be shown on the same plot.

#ifdef StudyEffectOfMassPairs
	TChain * MWR2600MNu1300_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR2600MNu1300_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_2600_MNu_1300_using_GEN_SIM.root");
	TChain * MWR2600MNu520_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR2600MNu520_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_2600_MNu_520_using_GEN_SIM.root");
	TChain * MWR2600MNu2080_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR2600MNu2080_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_2600_MNu_2080_using_GEN_SIM.root");
	TChain * MWR1400MNu700_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR1400MNu700_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_1400_MNu_700_using_GEN_SIM.root");
	
	///low mass mNu = 300,200,100,and 50 GeV
	TChain * MWR2600MNu300_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR2600MNu300_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_2600_MNu_300_using_GEN_SIM.root");
	TChain * MWR2600MNu200_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR2600MNu200_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_2600_MNu_200_using_GEN_SIM.root");
	TChain * MWR2600MNu100_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR2600MNu100_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_2600_MNu_100_using_GEN_SIM.root");
	TChain * MWR2600MNu50_matchedGenNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	MWR2600MNu50_matchedGenNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_MWR_2600_MNu_50_using_GEN_SIM.root");
	
	//makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax)
	string branchNames[] = {"ptEle[0]","ptEle[1]","etaEle[0]","etaEle[1]","ptJet[0]","ptJet[1]","etaJet[0]","etaJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet"};
	string link=">>";
	string histoEndings[] = {"_leadLeptonPt(50,0.,1600.)","_subleadLeptonPt(50,0.,1000.)","_leadLeptonEta(50,-3.0,3.0)","_subleadLeptonEta(50,-3.0,3.0)","_leadJetPt(50,0.,1000.)","_subleadJetPt(50,0.,600.)","_leadJetEta(50,-3.0,3.0)","_subleadJetEta(50,-3.0,3.0)","_dileptonMass(50,0.,2600.)","_fourObjectMass(50,0.,3000.)","_dR_leadingLeptonLeadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingJet(50,0.,5.)","_dR_subleadingLeptonLeadingJet(50,0.,2.)","_dR_subleadingLeptonSubleadingJet(50,0.,2.5)","_dR_leadingLeptonSubleadingLepton(50,0.,5.)","_dR_leadingJetSubleadingJet(50,0.,2.5)"};
	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {"mWR2600mNu100","mWR2600mNu300","mWR2600mNu200","mWR2600mNu50"};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = MWR2600MNu100_matchedGenNoCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = MWR2600MNu300_matchedGenNoCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = MWR2600MNu200_matchedGenNoCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[3]+histoEndings[i]] = MWR2600MNu50_matchedGenNoCuts;
		string cName = "o"+to_string(i);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.75,0.6,0.98,0.95,false);
		placeHolderMap.clear();
	}///end loop over branchNames

#endif

#ifdef bkgndOverlaidOnMatchedSignal
	///declare and initialize TChains to analyzed bkgnd and matched signal files
	TChain * dyPlusJetsNoCuts = new TChain("bkgndRecoAnalyzerOne/bkgndRecoObjectsNoCuts");
	dyPlusJetsNoCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/dyPlusJets/*.root");
	TChain * ttBarNoCuts = new TChain("bkgndRecoAnalyzerOne/bkgndRecoObjectsNoCuts");
	ttBarNoCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/ttBar/*.root");
	
	TChain * dyPlusJetsAllCuts = new TChain("bkgndRecoAnalyzerFive/bkgndRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts");
	dyPlusJetsAllCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/dyPlusJets/analyzedPhys14Samples/*.root");
	TChain * ttBarAllCuts = new TChain("bkgndRecoAnalyzerFive/bkgndRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts");
	ttBarAllCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/ttBar/analyzedPhys14Samples/*.root");
	TChain * matchedRecoAllCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts");
	matchedRecoAllCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/analyzed_gen_matched_WRtoEEJJ_MWr2600_MNu1300.root");

	///setup inputs needed for makeAndSaveMultipleCurveOverlayHist() fxn, and call this fxn
	string branchNames[] = {"ptEle[0]","ptEle[1]","etaEle[0]","etaEle[1]","ptJet[0]","ptJet[1]","etaJet[0]","etaJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet","subleadingLeptonThreeObjMass","leadLeptonThreeObjMass"};
	string link=">>";
	string histoEndings[] = {"_leadLeptonPt(50,0.,1200.)","_subleadLeptonPt(50,0.,700.)","_leadLeptonEta(50,-3.0,3.0)","_subleadLeptonEta(50,-3.0,3.0)","_leadJetPt(70,0.,900.)","_subleadJetPt(50,0.,500.)","_leadJetEta(50,-3.0,3.0)","_subleadJetEta(50,-3.0,3.0)","_dileptonMass(50,0.,2500.)","_fourObjectMass(50,400.,3300.)","_dR_leadingLeptonLeadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingJet(50,0.,5.)","_dR_subleadingLeptonLeadingJet(50,0.,5.)","_dR_subleadingLeptonSubleadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingLepton(50,0.,5.)","_dR_leadingJetSubleadingJet(50,0.,5.)","_subleadingLeptonThreeObjMass(50,0.,1600.)","_leadLeptonThreeObjMass(50,0.,2800.)",};
	string histoTitles[]={"Lead Electron P_{T}  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","Sublead Electron P_{T}  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","Lead Electron #eta  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","Sublead Electron #eta  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","Lead Jet P_{T}  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","Sublead Jet P_{T}  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","Lead Jet #eta  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","Sublead Jet #eta  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","DiElectron Mass  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","M_{EEJJ}  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","#DeltaR Lead Electron Lead Jet  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","#DeltaR Lead Electron Sublead Jet  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","#DeltaR Sublead Electron Lead Jet  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","#DeltaR Sublead Electron Sublead Jet  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","#DeltaR Lead Electron Sublead Electron  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","#DeltaR Lead Jet Sublead Jet  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","M_{EJJ} using Sublead Electron  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV","M_{EJJ} using Lead Electron  #surds = 13 TeV M_{Nu} = 1.3 TeV M_{WR} = 2.6 TeV"};
	string histoXaxisLabels[]={"Lead Electron P_{T} [GeV]","Sublead Electron P_{T} [GeV]","Lead Electron #eta","Sublead Electron #eta","Lead Jet P_{T} [GeV]","Sublead Jet P_{T} [GeV]","Lead Jet #eta","Sublead Jet #eta","DiElectron Mass [GeV]","M_{EEJJ} [GeV]","#DeltaR Lead Electron Lead Jet","#DeltaR Lead Electron Sublead Jet","#DeltaR Sublead Electron Lead Jet","#DeltaR Sublead Electron Sublead Jet","#DeltaR Lead Electron Sublead Electron","#DeltaR Lead Jet Sublead Jet","M_{EJJ} using Sublead Electron [GeV]","M_{EJJ} using Lead Electron [GeV]"};
	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {"WR","TTBar","DY"};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = matchedRecoAllCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = ttBarAllCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = dyPlusJetsAllCuts;
		string cName = "o"+to_string(i);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.8,0.6,0.95,0.9,true,histoTitles[i],histoXaxisLabels[i]);
		placeHolderMap.clear();
	}///end loop over branchNames

	//makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea,string title,string xLabel){
	
#endif


#ifdef compareCentrallyProducedToPrivateWrSignal

	//compare WR kinematic distributions between centrally and privately produced WR datasets
	
	TString dirAndTreeName = "wrDecayChainAnalyzer/genAndMatchedRecoWrDecayNoCuts";
	
	//UPDATE these directory names
	TString privateDirName = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRGen/", centralDirName = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/WR_signal_MC_centralProduction/RunIISpring15/";
	
	TChain * genWRtoEEJJMWR800MNu400Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR800MNu400Private->Add(privateDirName + "PRIVATEGENERATIONFILE.root");
	TChain * genWRtoEEJJMWR800MNu400Central = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR800MNu400Central->Add(centralDirName + "analyzed_central_WREEJJ_MWR_XXX_MNuYYY.root");
	

	string branchNames[] = {"threeObjMassFromGenObjsFromScdHvyPtcl","ptGenLeptFromFstHvyPtcl","etaGenLeptFromFstHvyPtcl","phiGenLeptFromFstHvyPtcl","ptGenLeptFromScdHvyPtcl","etaGenLeptFromScdHvyPtcl","phiGenLeptFromScdHvyPtcl","ptGenQuarkOneFromScdHvyPtcl","etaGenQuarkOneFromScdHvyPtcl","phiGenQuarkOneFromScdHvyPtcl","ptGenQuarkTwoFromScdHvyPtcl","etaGenQuarkTwoFromScdHvyPtcl","phiGenQuarkTwoFromScdHvyPtcl","dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl","fourObjMassFromGenObjsFromFstAndScdHvyPtcl"};
	string link=">>";
	//string histoEndings[] = {"_massNu(50,395,405)","_genLeptPtFromWR(50,0.,450.)","_genLeptEtaFromWR(50,-4.,4.)","_genLeptPhiFromWR(50,-3.2,3.2)","_genLeptPtFromNu(50,0.,450.)","_genLeptEtaFromNu(50,-4.,4.)","_genLeptPhiFromNu(50,-3.2,3.2)","_genQrkOnePtFromNu(50,0.,450.)","_genQrkOneEtaFromNu(50,-4.,4.)","_genQrkOnePhiFromNu(50,-3.2,3.2)","_genQrkTwoPtFromNu(50,0.,450.)","_genQrkTwoEtaFromNu(50,-4.,4.)","_genQrkTwoPhiFromNu(50,-3.2,3.2)","_dileptonMassGenLeptsFromWRandNu(50,0.,850.)","_massWrGenLeptsAndQrksFromWRandNu(50,600.,1000.)"};	//for MWR 800 MNu 400
	//string histoEndings[] = {"_massNu(50,1250,1350)","_genLeptPtFromWR(50,0.,1450.)","_genLeptEtaFromWR(50,-4.,4.)","_genLeptPhiFromWR(50,-3.2,3.2)","_genLeptPtFromNu(50,0.,1450.)","_genLeptEtaFromNu(50,-4.,4.)","_genLeptPhiFromNu(50,-3.2,3.2)","_genQrkOnePtFromNu(50,0.,1450.)","_genQrkOneEtaFromNu(50,-4.,4.)","_genQrkOnePhiFromNu(50,-3.2,3.2)","_genQrkTwoPtFromNu(50,0.,1450.)","_genQrkTwoEtaFromNu(50,-4.,4.)","_genQrkTwoPhiFromNu(50,-3.2,3.2)","_dileptonMassGenLeptsFromWRandNu(50,0.,2600.)","_massWrGenLeptsAndQrksFromWRandNu(50,2200.,3000.)"};	//for MWR 2600 MNu 1300
	string histoEndings[] = {"_massNu(50,800,2500)","_genLeptPtFromWR(50,0.,2500.)","_genLeptEtaFromWR(50,-4.,4.)","_genLeptPhiFromWR(50,-3.2,3.2)","_genLeptPtFromNu(50,0.,2500.)","_genLeptEtaFromNu(50,-4.,4.)","_genLeptPhiFromNu(50,-3.2,3.2)","_genQrkOnePtFromNu(50,0.,2500.)","_genQrkOneEtaFromNu(50,-4.,4.)","_genQrkOnePhiFromNu(50,-3.2,3.2)","_genQrkTwoPtFromNu(50,0.,2500.)","_genQrkTwoEtaFromNu(50,-4.,4.)","_genQrkTwoPhiFromNu(50,-3.2,3.2)","_dileptonMassGenLeptsFromWRandNu(50,0.,5000.)","_massWrGenLeptsAndQrksFromWRandNu(100,1000.,6000.)"};	//for MWR 5000 MNu 2500, but modified to show masses as low as MWR 2000 MNu 1000 with poor binning
		
	string titles[] = {"CMS Private                     #surds = 13 TeV"};
	string xAxisLabels[] = {"M_{LJJ} of GEN Nu daughter lepton and quarks [GeV]","P_{T} of GEN WR daughter lepton [GeV]","#eta of GEN WR daughter lepton","#phi of GEN WR daughter lepton","P_{T} of GEN Nu daughter lepton [GeV]","#eta of GEN Nu daughter lepton","#phi of GEN Nu daughter lepton","P_{T} of GEN Nu daughter quark one [GeV]","#eta of GEN Nu daughter quark one","#phi of GEN Nu daughter quark one","P_{T} of GEN Nu daughter quark two [GeV]","#eta of GEN Nu daughter quark two","#phi of GEN Nu daughter quark two","M_{LL} of GEN Nu and WR daughter leptons [GeV]","M_{LLJJ} of GEN WR and Nu daughter leptons and quarks [GeV]"};
	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {"Central 50k events","Private 15k events"};	//compare private GENSIM to centrally produced miniAOD
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = genWRtoEEJJMWR800MNu400Central;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = genWRtoEEJJMWR800MNu400Private;
		
		string cName = "o"+to_string(i);
		
		//UPDATE the MWR_800_MNu_400 string to reflect the mass point which is used 
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.65,0.6,0.95,0.90,true,titles[0],xAxisLabels[i],"_MWR_800_MNu_400_centralVsPrivate",false,-1);
		placeHolderMap.clear();
	}///end loop over branchNames
	//use makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea,string title,string xLabel,string outputFileNameModifier,Bool_t specialGrouping,Float_t cutVal){



#endif
	//end compareCentrallyProducedToPrivateWrSignal

#ifdef studyGenWrKinematicsVsWrAndNuMasses

	//compare distributions of the gen WR mass using the WR particle itself between centrally produced WR->eejj datasets
	//and privately produced WR->eejj datasets
	//use makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea,string title,string xLabel,string outputFileNameModifier,Bool_t specialGrouping){

	TString dirAndTreeName = "wrDecayChainAnalyzer/genAndMatchedRecoWrDecayNoCuts";
	TString privateDirName = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRGen/";
	TChain * genWRtoEEJJMWR800MNu400Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR800MNu400Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_800_MNu_400.root");
	TChain * genWRtoEEJJMWR800MNu150Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR800MNu150Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_800_MNu_150.root");
	TChain * genWRtoEEJJMWR800MNu700Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR800MNu700Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_800_MNu_700.root");
	
	
	TChain * genWRtoEEJJMWR2600MNu1300Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR2600MNu1300Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_2600_MNu_1300.root");
	TChain * genWRtoEEJJMWR2600MNu300Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR2600MNu300Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_2600_MNu_300.root");
	TChain * genWRtoEEJJMWR2600MNu2500Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR2600MNu2500Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_2600_MNu_2500.root");


	TChain * genWRtoEEJJMWR5000MNu2500Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR5000MNu2500Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_5000_MNu_2500.root");
	TChain * genWRtoEEJJMWR5000MNu400Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR5000MNu400Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_5000_MNu_400.root");
	TChain * genWRtoEEJJMWR5000MNu4700Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR5000MNu4700Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_5000_MNu_4700.root");


	TChain * genWRtoEEJJMWR1600MNu800Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR1600MNu800Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_1600_MNu_800.root");
	TChain * genWRtoEEJJMWR2400MNu800Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR2400MNu800Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_2400_MNu_800.root");
	TChain * genWRtoEEJJMWR3400MNu800Private = new TChain(dirAndTreeName,"");
	genWRtoEEJJMWR3400MNu800Private->Add(privateDirName + "analyzed_private_WREEJJ_MWR_3400_MNu_800.root");

	string title = "CMS Private                          #surds = 13 TeV";
	string link=">>";
	
	//the arrays branchNames, cutVals, and xAxisLabels must have the same number of elements
	string branchNames[] = {"ptGenLeptFromFstHvyPtcl","ptGenLeptFromScdHvyPtcl","ptGenQuarkOneFromScdHvyPtcl","ptGenQuarkTwoFromScdHvyPtcl","dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl","dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl","dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl","subleadGenLeptonNotFromScdHvyPtcl"};
	vector<string> brNamesVect(branchNames,branchNames + sizeof(branchNames)/sizeof(string));
	Float_t cutVals[] = {-1, -1, -1, -1, -1, -1, -1, -1};	//cut values for different distributions, synched with branchNames array
	string xAxisLabels[] = {"P_{T} [GeV]","P_{T} [GeV]","P_{T} [GeV]","P_{T} [GeV]","M_{LL} [GeV]","#DeltaR lepton and quark","#DeltaR lepton and quark","1 = events in which subleading GEN lepton does not have Nu mother"};
	
	//string histoEndingsMaster[] = {};
	//TChain * chains[] = {genWRtoEEJJMWR800MNu150Private, genWRtoEEJJMWR800MNu400Private, genWRtoEEJJMWR800MNu700Private, genWRtoEEJJMWR2600MNu300Private, genWRtoEEJJMWR2600MNu1300Private, genWRtoEEJJMWR2600MNu2500Private, genWRtoEEJJMWR5000MNu400Private, genWRtoEEJJMWR5000MNu2500Private, genWRtoEEJJMWR5000MNu4700Private, genWRtoEEJJMWR1600MNu800Private, genWRtoEEJJMWR2400MNu800Private, genWRtoEEJJMWR3400MNu800Private};
	//string histoBeginnings[] = {"MWR 800 MNu 150","MWR 800 MNu 400","MWR 800 MNu 700","MWR 2600 MNu 300","MWR 2600 MNu 1300","MWR 2600 MNu 2500","MWR 5000 MNu 400","MWR 5000 MNu 2500","MWR 5000 MNu 4700","MWR 1600 MNu 800","MWR 2400 MNu 800","MWR 3400 MNu 800"};	//legend entries
	//string fileLabels[] = {"_MWR_800_several_MNu_private","_MWR_800_several_MNu_private","_MWR_800_several_MNu_private","_MWR_2600_several_MNu_private","_MWR_2600_several_MNu_private","_MWR_2600_several_MNu_private","_MWR_5000_several_MNu_private","_MWR_5000_several_MNu_private","_MWR_5000_several_MNu_private","_MNu_800_several_MWR_private","_MNu_800_several_MWR_private","_MNu_800_several_MWR_private"};	//output file labels

	//string histoEndings[] = {"_genLeptPtFromWR(70,0.,450.)","_genLeptPtFromNu(50,0.,450.)","_genQrkOnePtFromNu(50,0.,450.)","_genQrkTwoPtFromNu(50,0.,450.)","_dileptonMassGenLeptsFromWRandNu(100,0.,850.)","_deltaRgenLeptAndQrkOneFromNu(100,0.,4.5)","_deltaRgenLeptAndQrkTwoFromNu(100,0.,4.5)","_subleadGenLeptNotFromNu(2,0.,2.)"};	//low mass
	//string histoBeginnings[] = {"800 GeV WR  150 GeV Nu","800 GeV WR  400 GeV Nu","800 GeV WR  700 GeV Nu"};	//low mass
	
	string histoEndings[] = {"_genLeptPtFromWR(100,0.,1450.)","_genLeptPtFromNu(100,0.,1450.)","_genQrkOnePtFromNu(100,0.,1450.)","_genQrkTwoPtFromNu(100,0.,1450.)","_dileptonMassGenLeptsFromWRandNu(100,0.,2600.)","_deltaRgenLeptAndQrkOneFromNu(100,0.,4.5)","_deltaRgenLeptAndQrkTwoFromNu(100,0.,4.5)","_subleadGenLeptNotFromNu(2,0.,2.)"};	//medium mass
	string histoBeginnings[] = {"MWR 2600 MNu 300","MWR 2600 MNu 1300","MWR 2600 MNu 2500"};	//medium mass

	//string histoEndings[] = {"_genLeptPtFromWR(100,0.,2500.)","_genLeptPtFromNu(100,0.,2500.)","_genQrkOnePtFromNu(100,0.,2500.)","_genQrkTwoPtFromNu(100,0.,2500.)","_dileptonMassGenLeptsFromWRandNu(100,0.,5000.)","_deltaRgenLeptAndQrkOneFromNu(100,0.,4.5)","_deltaRgenLeptAndQrkTwoFromNu(100,0.,4.5)","_subleadGenLeptNotFromNu(2,0.,2.)"};	//high mass
	//string histoBeginnings[] = {"MWR 5000 MNu 400","MWR 5000 MNu 2500","MWR 5000 MNu 4700"};	//high mass

	//string histoEndings[] = {"_genLeptPtFromWR(80,0.,1800.)","_genLeptPtFromNu(50,0.,1800.)","_genQrkOnePtFromNu(50,0.,1800.)","_genQrkTwoPtFromNu(50,0.,1800.)","_dileptonMassGenLeptsFromWRandNu(100,0.,3400.)","_deltaRgenLeptAndQrkOneFromNu(100,0.,4.5)","_deltaRgenLeptAndQrkTwoFromNu(100,0.,4.5)","_subleadGenLeptNotFromNu(2,0.,2.)"};	//fixed Nu mass
	//string histoBeginnings[] = {"1600 GeV WR  800 GeV Nu","2400 GeV WR  800 GeV Nu","3400 GeV WR  800 GeV Nu"};	//fixed Nu mass
	
	//string histoEndings[] = {"_genLeptPtFromWR(100,0.,3000.)","_genLeptPtFromNu(70,0.,3000.)","_genQrkOnePtFromNu(70,0.,3000.)","_genQrkTwoPtFromNu(70,0.,3000.)","_dileptonMassGenLeptsFromWRandNu(70,0.,4800.)","_deltaRgenLeptAndQrkOneFromNu(100,0.,4.5)","_deltaRgenLeptAndQrkTwoFromNu(100,0.,4.5)","_subleadGenLeptNotFromNu(2,0.,2.)"};	//fixed Nu mass
	//string histoBeginnings[] = {"1600 GeV WR  800 GeV Nu","2600 GeV WR  1300 GeV Nu","5000 GeV WR  2500 GeV Nu"};	//several MWR up to 5000 GeV, MNu is half MWR
	
	vector<string> histoBeginningsVect(histoBeginnings,histoBeginnings + sizeof(histoBeginnings)/sizeof(string));
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = brNamesVect.size();
	for(unsigned int i=0; i<maxI; i++){
		string cName = "o"+to_string(i);
		//placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = genWRtoEEJJMWR800MNu150Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = genWRtoEEJJMWR800MNu400Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = genWRtoEEJJMWR800MNu700Private;

		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = genWRtoEEJJMWR2600MNu300Private;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = genWRtoEEJJMWR2600MNu1300Private;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = genWRtoEEJJMWR2600MNu2500Private;

		//placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = genWRtoEEJJMWR5000MNu400Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = genWRtoEEJJMWR5000MNu2500Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = genWRtoEEJJMWR5000MNu4700Private;

		//placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = genWRtoEEJJMWR1600MNu800Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = genWRtoEEJJMWR2400MNu800Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = genWRtoEEJJMWR3400MNu800Private;

		//placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = genWRtoEEJJMWR1600MNu800Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = genWRtoEEJJMWR2600MNu1300Private;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = genWRtoEEJJMWR5000MNu2500Private;

		makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.42,0.55,0.90,0.90,true,title,xAxisLabels[i],"_MWR_2600_several_MNu_private",true,cutVals[i]);
		placeHolderMap.clear();
	}///end loop over branchNames
	
	/*
	unsigned int maxIterations = histoBeginningsVect.size();
	unsigned int diffMassesPerPlot = 3;	///<user controlled
	for(unsigned int i=0; i<maxI; i++){
		string cName;
	
		for(unsigned int j=3; j<6; j+=diffMassesPerPlot){
			for(unsigned int n=j; n<(j+diffMassesPerPlot); n++){
				placeHolderMap[branchNames[i]+link+histoBeginnings[n]+histoEndingsMaster[n]] = chains[n];
				cName = "o"+to_string(n+i+j);
			}

			makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.65,0.6,0.95,0.90,true,title,xAxisLabels[i],fileLabels[j],true,cutVals[i]);
			placeHolderMap.clear();
		}//end loop over different mass points

	}///end loop over branchNames
	*/
	//use makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea,string title,string xLabel,string outputFileNameModifier,Bool_t specialGrouping,Float_t cutVal){



#endif
//end ifdef studyGenWrKinematicsVsWrAndNuMasses


#ifdef genPlotsUsingWRDecayProducts
	TChain * MWR2600MNu1300_matchedGenDecayProductsNoCuts = new TChain("genMatchedParticleAnalyzerOne/genLeptonsAndJetsNoCuts","");
	MWR2600MNu1300_matchedGenDecayProductsNoCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-2600_Nu_M-1300.root");
	TChain * MWR1400MNu700_matchedGenDecayProductsNoCuts = new TChain("genMatchedParticleAnalyzerOne/genLeptonsAndJetsNoCuts","");
	MWR1400MNu700_matchedGenDecayProductsNoCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-1400_Nu_M-700.root");
	TChain * MWR4400MNu2200_matchedGenDecayProductsNoCuts = new TChain("genMatchedParticleAnalyzerOne/genLeptonsAndJetsNoCuts","");
	MWR4400MNu2200_matchedGenDecayProductsNoCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-4400_Nu_M-2200.root");
	
	
	TChain * MWR2600MNu1300_matchedGenDecayProductsPtEtaDrCuts = new TChain("genMatchedParticleAnalyzerTwoPFive/genLeptonsAndJetsWithPtEtaDrCuts","");
	MWR2600MNu1300_matchedGenDecayProductsPtEtaDrCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-2600_Nu_M-1300.root");
	TChain * MWR1400MNu700_matchedGenDecayProductsPtEtaDrCuts = new TChain("genMatchedParticleAnalyzerTwoPFive/genLeptonsAndJetsWithPtEtaDrCuts","");
	MWR1400MNu700_matchedGenDecayProductsPtEtaDrCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-1400_Nu_M-700.root");
	TChain * MWR4400MNu2200_matchedGenDecayProductsPtEtaDrCuts = new TChain("genMatchedParticleAnalyzerTwoPFive/genLeptonsAndJetsWithPtEtaDrCuts","");
	MWR4400MNu2200_matchedGenDecayProductsPtEtaDrCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-4400_Nu_M-2200.root");


	TChain * MWR2600MNu1300_matchedGenDecayProductsAllCuts = new TChain("genMatchedParticleAnalyzerThree/genLeptonsAndJetsWithAllCuts","");
	MWR2600MNu1300_matchedGenDecayProductsAllCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-2600_Nu_M-1300.root");
	TChain * MWR1400MNu700_matchedGenDecayProductsAllCuts = new TChain("genMatchedParticleAnalyzerThree/genLeptonsAndJetsWithAllCuts","");
	MWR1400MNu700_matchedGenDecayProductsAllCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-1400_Nu_M-700.root");
	TChain * MWR4400MNu2200_matchedGenDecayProductsAllCuts = new TChain("genMatchedParticleAnalyzerThree/genLeptonsAndJetsWithAllCuts","");
	MWR4400MNu2200_matchedGenDecayProductsAllCuts->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-4400_Nu_M-2200.root");
	

	string branchNames[] = {"fourObjectMass","subleadingLeptonThreeObjMass","dileptonMass","dijetMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","ptEle[0]","etaEle[0]","phiEle[0]","ptEle[1]","etaEle[1]","phiEle[1]","ptJet[0]","etaJet[0]","phiJet[0]","ptJet[1]","etaJet[1]","phiJet[1]","subleadingLeptonThreeObjMass/fourObjectMass"};
	string link=">>";
	string histoEndings[] = {"_massWR(50,700,4900)","_subleadLeptonNuMass(50,300,2400)","_mLL(50,0,3700)","_mJJ(50,0,2500)","_leadLeptonLeadJetSeparation(30,0.,5)","_leadLeptonSubleadJetSeparation(30,0.,5)","_subleadLeptonLeadJetSeparation(30,0.,5)","_subleadLeptonSubleadJetSeparation(30,0.,5)","_leadLeptonPt(70,30,2100)","_leadLeptonEta(30,-2.5,2.5)","_leadLeptonPhi(30,-3.2,3.2)","_subleadLeptonPt(50,0,1600)","_subleadLeptonEta(30,-2.5,2.5)","_subleadLeptonPhi(30,-3.2,3.2)","_leadJetPt(60,0,1900)","_leadJetEta(30,-2.5,2.5)","_leadJetPhi(30,-3.2,3.2)","_subleadJetPt(50,0,1000)","_subleadJetEta(30,-2.5,2.5)","_subleadJetPhi(30,-3.2,3.2)","_threeObjMassDivByFourObjMass(50,0.,1.)"};
	string titles[] = {"GEN WR Mass using gen leptons and jets","GEN Nu Mass using gen leptons and jets","GEN M_{EE} using gen leptons","GEN M_{JJ} using gen jets","GEN #DeltaR(lead lepton,lead jet)","GEN #DeltaR(lead lepton,sublead jet)","GEN #DeltaR(sublead lepton,lead jet)","GEN #DeltaR(sublead lepton,sublead jet)","P_{T} of GEN lepton with WR mother","#eta of GEN lepton with WR mother","#phi of GEN lepton with WR mother","P_{T} of GEN lepton with Nu mother","#eta of GEN lepton with Nu mother","#phi of GEN lepton with Nu mother","P_{T} of lead GEN jet with Nu mother","#eta of lead GEN jet with Nu mother","#phi of lead GEN jet with Nu mother","P_{T} of sublead GEN jet with Nu mother","#eta of sublead GEN jet with Nu mother","#phi of sublead GEN jet with Nu mother","GEN Nu Mass/WR Mass using gen leptons and jets"};
	string xAxisLabels[] = {"WR mass [GeV]","Nu mass [GeV]","Dilepton mass [GeV]","Dijet mass [GeV]","#DeltaR(L,J)","#DeltaR(L,J)","#DeltaR(L,J)","#DeltaR(L,J)","pT [GeV]","#eta","#phi","pT [GeV]","#eta","#phi","pT [GeV]","#eta","#phi","pT [GeV]","#eta","#phi","Nu Mass/WR Mass"};
	string outputFileModifier[] = {"_MWR_4400_2600_1400_MNu_2200_1300_700_with_no_cuts","_MWR_4400_2600_1400_MNu_2200_1300_700_with_pt_eta_dR_cuts","_MWR_4400_2600_1400_MNu_2200_1300_700_with_all_cuts"};
	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {"MWR 1400 GeV","MWR 2600 GeV","MWR 4400 GeV"};
	map<string,TChain*> placeHolderMapOne, placeHolderMapTwo, placeHolderMapThree;
	unsigned int maxI = histoEndingVect.size();
	//unsigned int maxI = 2;
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMapOne[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = MWR1400MNu700_matchedGenDecayProductsNoCuts;
		placeHolderMapOne[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = MWR2600MNu1300_matchedGenDecayProductsNoCuts;
		placeHolderMapOne[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = MWR4400MNu2200_matchedGenDecayProductsNoCuts;
	
		placeHolderMapTwo[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = MWR1400MNu700_matchedGenDecayProductsPtEtaDrCuts;
		placeHolderMapTwo[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = MWR2600MNu1300_matchedGenDecayProductsPtEtaDrCuts;
		placeHolderMapTwo[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = MWR4400MNu2200_matchedGenDecayProductsPtEtaDrCuts;
		
		placeHolderMapThree[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = MWR1400MNu700_matchedGenDecayProductsAllCuts;
		placeHolderMapThree[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = MWR2600MNu1300_matchedGenDecayProductsAllCuts;
		placeHolderMapThree[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = MWR4400MNu2200_matchedGenDecayProductsAllCuts;
		
		string cNameOne = "o"+to_string(i), cNameTwo = "j"+to_string(i), cNameThree = "k"+to_string(i);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMapOne,cNameOne.c_str(),0.72,0.53,0.95,0.88,true,titles[i]+" after no cuts",xAxisLabels[i],outputFileModifier[0]);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMapTwo,cNameTwo.c_str(),0.72,0.53,0.95,0.88,true,titles[i]+" after pT, #eta, #DeltaR(L,J) cuts",xAxisLabels[i],outputFileModifier[1]);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMapThree,cNameThree.c_str(),0.72,0.53,0.95,0.88,true,titles[i]+" after all cuts",xAxisLabels[i],outputFileModifier[2]);
	
		placeHolderMapOne.clear(), placeHolderMapTwo.clear(), placeHolderMapThree.clear();
	
	}///end loop over branchNames


#endif


#ifdef twoDimPlotGenWrAcceptance
	///make a 2D plot with MNu as the vertical axis, MWR as the horizontal axis, and the fraction of GEN WR->eejj evts which
	///pass all offline cuts at each point in the (MWR, MNu) space
	///also print a table of these values

	///all input .root files should be in the same directory, and have file names which differ only in the WR and Nu mass values
	string dir= "/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/";
	string fileBegin = "analyzed_genWrToEEJJFullOfflineAnalysis_WR_";
	string fileEnd = "_1.root";
	string fileMiddle = "_NU_";
	gStyle->SetTitleOffset(1.6,"Y");
	gStyle->SetOptStat("");

	/*
	string cutEfficiencyVsMassFile = "offlineEfficienciesVsMasses.txt";
	ofstream writeToEfficiencyFile(cutEfficiencyVsMassFile.c_str(),ofstream::app);
	Float_t passingPercentage=-1;
	TH2F * twoDimAcceptanceHist = new TH2F("twoDimAccHist","Percentage of events passing GEN pT, eta, and #DeltaR cuts",58,700,3600,72,0,3600);

	int maxWrMass=3500, increment=50;
	//int maxWrMass=1850, increment=50;
	//starting value for wrMass equals the minimum wr mass
	for(int wrMass=800; wrMass<=maxWrMass ; wrMass+=increment){
		///loop over WR mass values
		
		for(int nuMass=50; nuMass<wrMass ; nuMass+=increment){
			///loop over Nu mass values

			///define input root file name
			string pfn = dir+fileBegin+to_string(wrMass)+fileMiddle+to_string(nuMass)+fileEnd;
			
			///define one TChain to count the number of events generated, and record the desired WR and Nu masses
			///define another TChain to count the number of evts passing all offline cuts at GEN lvl (no HLT or ID)
			TChain * genInfo = new TChain("genNuAnalyzerOne/genNu");
			genInfo->Add(pfn.c_str());
			//TChain * afterOfflineCuts = new TChain("genMatchedParticleAnalyzerThree/genLeptonsAndJetsWithAllCuts");
			TChain * afterOfflineCuts = new TChain("genMatchedParticleAnalyzerTwoPFive/genLeptonsAndJetsWithPtEtaDrCuts");
			afterOfflineCuts->Add(pfn.c_str());

			///calculate percentage of evts which pass GEN cuts, and store this percentage along with the nu and wr mass values in a txt file
			passingPercentage = (100)*((Float_t) afterOfflineCuts->GetEntries()/genInfo->GetEntries());
			//writeToEfficiencyFile << passingPercentage <<"\tpercent of events with WR mass=\t"<< wrMass <<"\tand Nu mass=\t"<< nuMass <<"\tpass all cuts"<< endl;
			
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
	r1->SaveAs("twoDimGenWrAcceptances_afterBasePtEtaDrCuts.png","recreate");
	r1->SaveAs("twoDimGenWrAcceptances_afterBasePtEtaDrCuts.pdf","recreate");
	*/

	/////////////////////////////////////////////////////////////////////////////////////////////
	///use a combination of privately produced GEN samples with centrally produced MINIAOD samples to make a TGraph of
	///the ratio of (RECO evts passing all WR signal region cuts|all GEN cuts passed)/(GEN evts passing all WR signal region cuts)
	///the RECO cuts include HEEP and jet ID requirements, whereas the GEN cuts do not
	///neither set of cuts requires that the trigger is fired (HLT_DoubleEle33)
	string recoDir = "/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/";
	string recoFileBegin = "all_analyzed_tree_allGenAndRecoOfflineCuts_eejjSignalRegion_WR_M-";
	string recoFileEnd = ".root";
	string recoFileMiddle = "_Nu_M-";
	string RecoOvrGenvsMassFile = "recoOvrGenScaleFactorsVsMass.txt";
	ofstream writeToRecoGenScaleFactorsFile(RecoOvrGenvsMassFile.c_str(),ofstream::app);
	Int_t nBins = 23;	///max is 23
	Float_t wrMassVals[nBins], recoGenScaleFactors[nBins];
	
	int wrMassArr[] = {800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3600,3800,4000,4200,4400,5000,5200,5600,5800,6000};
	//int wrMassArr[] = {800,1000};

	for(int i=0; i<nBins ; i++){
		string recoPfn = recoDir+recoFileBegin+to_string(wrMassArr[i])+recoFileMiddle+to_string(wrMassArr[i]/2)+recoFileEnd;
		wrMassVals[i] = (Float_t) wrMassArr[i];

		TChain * genAfterCuts = new TChain("genMatchedParticleAnalyzerAfterFullGenSelection/genLeptonsAndJetsWithAllCuts");
		genAfterCuts->Add(recoPfn.c_str());
		TChain * recoAfterCuts = new TChain("unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts");
		recoAfterCuts->Add(recoPfn.c_str());
		//TChain * genBeforeCuts = new TChain("genWRAnalyzer/genWR");
		//genBeforeCuts->Add(recoPfn.c_str());
	
		Float_t recoRatio = ((Float_t) recoAfterCuts->GetEntries()/genAfterCuts->GetEntries());
		//Float_t genRatio = ((Float_t) genAfterCuts->GetEntries()/genBeforeCuts->GetEntries());
		//recoGenScaleFactors[i] = recoRatio/genRatio;
		recoGenScaleFactors[i] = recoRatio;
		writeToRecoGenScaleFactorsFile << "for MWR=\t"<< wrMassArr[i] << "\tand MNu=\t"<< wrMassArr[i]/2 << "\tthe ratio of RECO evts passing all cuts (given evts which pass all GEN cuts) to GEN evts passing all cuts (RECO|GEN/GEN) is\t"<< recoGenScaleFactors[i] << endl;

		//genBeforeCuts->Delete();
		genAfterCuts->Delete();
		recoAfterCuts->Delete();

	}
	writeToRecoGenScaleFactorsFile.close();

	///make the TGraph, and save an image of it
	TGraph * recoGenFactorsGraph = new TGraph(nBins,wrMassVals,recoGenScaleFactors);
	recoGenFactorsGraph->GetXaxis()->SetTitle("WR Mass [GeV]");
	recoGenFactorsGraph->GetYaxis()->SetTitle("RECO/GEN scale factor");
	recoGenFactorsGraph->SetTitle("RECO to GEN scale factors vs MWR  MNu = MWR/2");
	recoGenFactorsGraph->SetMarkerStyle(20);
	recoGenFactorsGraph->SetMarkerColor(2);
	
	TCanvas * f1 = new TCanvas("f1","f1",800,700);
	f1->cd();
	recoGenFactorsGraph->Draw("AP");
	f1->SaveAs("recoOverGenScaleFactorVsWRMass.png","recreate");
	f1->SaveAs("recoOverGenScaleFactorVsWRMass.pdf","recreate");

#endif
	//end twoDimPlotGenWrAcceptance

#ifdef recoAndGenHLTEfficiency
	///similar to what is done in twoDimPlotGenWrAcceptance, but with HLT efficiency for evts passing:
	///		all offline cuts (pT, eta, dR, MLL, MLLJJ)
	///		gen pT (pT > 40 for gen leptons and jets matched to WR decay products), gen eta (|eta| < 2.5), and gen dR(L,J) cuts
	
	///all input .root files should be in the same directory, and have file names which differ only in the WR and Nu mass values
	string dir= "/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/";
	string fileBegin = "all_analyzed_tree_checkHLTEfficiencyRecoFullOfflineAndGenPtEtaDr_eejjSignalRegion_WR_M-";
	string fileEnd = ".root";
	string fileMiddle = "_Nu_M-";
	string hltEfficiencyVsMassFile = "hltEfficienciesVsMasses.txt";
	ofstream writeToHltEfficiencyFile(hltEfficiencyVsMassFile.c_str(),ofstream::app);
	gStyle->SetTitleOffset(1.4,"Y");
	Int_t nBins = 23;	///max is 23
	Float_t wrMassVals[nBins], genPtEtaDrHltEff[nBins], recoFullOfflineHltEff[nBins];
	gStyle->SetOptStat("");

	int wrMassArr[] = {800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3600,3800,4000,4200,4400,5000,5200,5600,5800,6000};
	//int wrMassArr[] = {800,1000};

	for(int i=0; i<nBins ; i++){
		///loop over WR mass values
		
		///define input root file name and add a Float_t to the x axis array used for later TGraph objects
		string pfn = dir+fileBegin+to_string(wrMassArr[i])+fileMiddle+to_string(wrMassArr[i]/2)+fileEnd;
		wrMassVals[i] = (Float_t) wrMassArr[i];

		///define one TChain to count the number of events which pass cuts before HLT
		///define another TChain to count the number of evts passing cuts and HLT
		TChain * genCutsBeforeHlt = new TChain("genMatchedParticleAnalyzerAfterGenPtEtaDrCuts/genLeptonsAndJetsWithPtEtaDrCuts");
		genCutsBeforeHlt->Add(pfn.c_str());
		TChain * genCutsAfterHlt = new TChain("genMatchedParticleAnalyzerAfterGenPtEtaDrCutsPostHLT/genLeptonsAndJetsWithPtEtaDrCutsPostHLT");
		genCutsAfterHlt->Add(pfn.c_str());


		TChain * recoOfflineCutsBeforeHlt = new TChain("unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts");
		recoOfflineCutsBeforeHlt->Add(pfn.c_str());
		TChain * recoOfflineCutsAfterHlt = new TChain("unmatchedSignalRecoAnalyzerFivePostHLT/signalRecoObjectsWithAllCutsAfterHLT");
		recoOfflineCutsAfterHlt->Add(pfn.c_str());


		///calculate percentage of evts which pass GEN cuts, and store this percentage along with the nu and wr mass values in a txt file
		//genHltPassingPercentage = (100)*((Float_t) genCutsAfterHlt->GetEntries()/genCutsBeforeHlt->GetEntries());
		//recoFullOfflineHltPassingPercentage = (100)*((Float_t) );
		genPtEtaDrHltEff[i] = (100)*((Float_t) genCutsAfterHlt->GetEntries()/genCutsBeforeHlt->GetEntries());
		recoFullOfflineHltEff[i] = (100)*((Float_t) recoOfflineCutsAfterHlt->GetEntries()/recoOfflineCutsBeforeHlt->GetEntries());

		writeToHltEfficiencyFile << genPtEtaDrHltEff[i] <<"\tpercent of events with WR mass=\t"<< wrMassArr[i] <<"\tand Nu mass=\t"<< wrMassArr[i]/2 <<"\tfire HLT_DoubleEle33 after passing GEN pT, eta, and dR(L,J) cuts"<< endl;
		writeToHltEfficiencyFile << recoFullOfflineHltEff[i] <<"\tpercent of events with WR mass=\t"<< wrMassArr[i] <<"\tand Nu mass=\t"<< wrMassArr[i]/2 <<"\tfire HLT_DoubleEle33 after passing all offline RECO cuts"<< endl;

		genCutsBeforeHlt->Delete();
		genCutsAfterHlt->Delete();
		recoOfflineCutsBeforeHlt->Delete();
		recoOfflineCutsAfterHlt->Delete();

	}///end loop over WR mass values	

	///close txt file
	writeToHltEfficiencyFile.close();


	///make TGraph objects after arrays of TGraph contents have been filled
	TGraph * genHltEffGraph = new TGraph(nBins,wrMassVals,genPtEtaDrHltEff);
	TGraph * recoFullOfflineHltEffGraph = new TGraph(nBins,wrMassVals,recoFullOfflineHltEff);
	
	///set axis labels and draw the histo, save an image of it (both pdf and png formats)
	genHltEffGraph->GetXaxis()->SetTitle("WR Mass [GeV]");
	genHltEffGraph->GetYaxis()->SetTitle("Percentage of events passing");
	genHltEffGraph->SetTitle("Percentage of WR signal events passing HLT after GEN pT, #eta, and #DeltaR cuts  MNu = MWR/2");
	genHltEffGraph->SetMarkerStyle(20);
	genHltEffGraph->SetMarkerColor(2);

	recoFullOfflineHltEffGraph->GetXaxis()->SetTitle("WR Mass [GeV]");
	recoFullOfflineHltEffGraph->GetYaxis()->SetTitle("Percentage of events passing");
	recoFullOfflineHltEffGraph->SetTitle("Percentage of WR signal events passing HLT after all RECO offline cuts  MNu = MWR/2");
	recoFullOfflineHltEffGraph->SetMarkerStyle(20);
	recoFullOfflineHltEffGraph->SetMarkerColor(2);

	
	TCanvas * r1 = new TCanvas("r1","r1",900,700);
	r1->cd();
	genHltEffGraph->Draw("AP");
	r1->SaveAs("hltEfficiencyAfterGenPtEtaDrCuts.png","recreate");
	r1->SaveAs("hltEfficiencyAfterGenPtEtaDrCuts.pdf","recreate");

	TCanvas * q1 = new TCanvas("q1","q1",900,700);
	q1->cd();
	recoFullOfflineHltEffGraph->Draw("AP");
	q1->SaveAs("hltEfficiencyAfterAllOfflineRecoCuts.png","recreate");
	q1->SaveAs("hltEfficiencyAfterAllOfflineRecoCuts.pdf","recreate");



#endif
	//end recoAndGenHLTEfficiency
	
	
#ifdef genAndRecoWrPlotsMinimalCuts
	///all input .root files should be in the same directory, and have file names which differ only in the WR and Nu mass values
	string dir= "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/WR_signal_MC_centralProduction/RunIISpring15/";
	string fileBegin = "analyzed_";
	string privOrCent = "central_";
	string fileWrMass = "WREEJJ_MWR_";
	string fileEnd = ".root";
	string fileMiddle = "_MNu_";
	string genCutEffVsMassFile = "genCutEfficienciesVsMasses.txt";
	ofstream writeToGenEfficiencyFile(genCutEffVsMassFile.c_str(),ofstream::trunc);
	string barrelAndEndcapFractionFile = "fractionOfEventsInBarrelAndEndcaps.txt";
	ofstream writeToBarrelEndcapFile(barrelAndEndcapFractionFile.c_str(), ofstream::trunc);
	gStyle->SetTitleOffset(1.4,"Y");
	Int_t nBins = 5;	///max is 24
	Float_t wrMassVals[nBins], genMatchedPtEtaEff[nBins], genLeadAndSubleadPtEtaEff[nBins];
	gStyle->SetOptStat("");

	int wrMassArr[] = {800,1600,2600,3400,5000};
	//int wrMassArr[] = {800,1000,1200,1400,1600,1800,2000,2400,2600,2800,3000,3200,3600,3800,4000,4200,4400,4600,4800,5000,5200,5600,5800,6000};
	//int wrMassArr[] = {2200};


	//element number i in plotArg is linked to element number i in plotCut   don't change the order
	//gen and reco branches
	//to make plots without cuts
	//string plotArg[] = {"numGenFstHvyPtcl","numGenScdHvyPtcl","numGenLeptons","numGenQuarks","leadGenLeptonNotFromFstHvyPtcl","subleadGenLeptonNotFromScdHvyPtcl","leadGenQuarkNotFromScdHvyPtcl","subleadGenQuarkNotFromScdHvyPtcl","etaGenFstHvyPtcl", "ptGenFstHvyPtcl", "massGenFstHvyPtcl", "etaGenScdHvyPtcl", "ptGenScdHvyPtcl", "massGenScdHvyPtcl", "etaGenLeptFromFstHvyPtcl", "ptGenLeptFromFstHvyPtcl", "phiGenLeptFromFstHvyPtcl", "etaGenLeptFromScdHvyPtcl", "ptGenLeptFromScdHvyPtcl", "phiGenLeptFromScdHvyPtcl", "etaGenQuarkOneFromScdHvyPtcl", "ptGenQuarkOneFromScdHvyPtcl", "phiGenQuarkOneFromScdHvyPtcl", "etaGenQuarkTwoFromScdHvyPtcl", "ptGenQuarkTwoFromScdHvyPtcl", "phiGenQuarkTwoFromScdHvyPtcl","ptLeadGenLepton", "etaLeadGenLepton", "phiLeadGenLepton", "ptSubleadGenLepton", "etaSubleadGenLepton", "phiSubleadGenLepton", "ptLeadGenQuark", "etaLeadGenQuark", "phiLeadGenQuark", "ptSubleadGenQuark", "etaSubleadGenQuark", "phiSubleadGenQuark","ptRecoLeptMatchedToWrDau", "etaRecoLeptMatchedToWrDau", "phiRecoLeptMatchedToWrDau", "ptRecoLeptMatchedToNuDau", "etaRecoLeptMatchedToNuDau", "phiRecoLeptMatchedToNuDau", "ptRecoJetOneMatchedToNuDau", "etaRecoJetOneMatchedToNuDau", "phiRecoJetOneMatchedToNuDau", "ptRecoJetTwoMatchedToNuDau", "etaRecoJetTwoMatchedToNuDau", "phiRecoJetTwoMatchedToNuDau", "ptGenJetFromMatchedRecoJetOne", "etaGenJetFromMatchedRecoJetOne", "phiGenJetFromMatchedRecoJetOne", "ptGenJetFromMatchedRecoJetTwo", "etaGenJetFromMatchedRecoJetTwo", "phiGenJetFromMatchedRecoJetTwo", "ptLeadRecoLept", "etaLeadRecoLept", "phiLeadRecoLept", "ptSubleadRecoLept", "etaSubleadRecoLept", "phiSubleadRecoLept", "ptLeadRecoJet", "etaLeadRecoJet", "phiLeadRecoJet", "ptSubleadRecoJet", "etaSubleadRecoJet", "phiSubleadRecoJet"};
	//string plotCut[] = {"numGenFstHvyPtcl>-1","numGenScdHvyPtcl>-1","numGenLeptons>-1","numGenQuarks>-1","leadGenLeptonNotFromFstHvyPtcl>-1","subleadGenLeptonNotFromScdHvyPtcl>-1","leadGenQuarkNotFromScdHvyPtcl>-1","subleadGenQuarkNotFromScdHvyPtcl>-1","etaGenFstHvyPtcl>-9", "ptGenFstHvyPtcl>-9", "massGenFstHvyPtcl>-9", "etaGenScdHvyPtcl>-9", "ptGenScdHvyPtcl>-9", "massGenScdHvyPtcl>-9", "etaGenLeptFromFstHvyPtcl>-9", "ptGenLeptFromFstHvyPtcl>-9", "phiGenLeptFromFstHvyPtcl>-9", "etaGenLeptFromScdHvyPtcl>-9", "ptGenLeptFromScdHvyPtcl>-9", "phiGenLeptFromScdHvyPtcl>-9", "etaGenQuarkOneFromScdHvyPtcl>-9", "ptGenQuarkOneFromScdHvyPtcl>-9", "phiGenQuarkOneFromScdHvyPtcl>-9", "etaGenQuarkTwoFromScdHvyPtcl>-9", "ptGenQuarkTwoFromScdHvyPtcl>-9", "phiGenQuarkTwoFromScdHvyPtcl>-9","ptLeadGenLepton>-9", "etaLeadGenLepton>-9", "phiLeadGenLepton>-9", "ptSubleadGenLepton>-9", "etaSubleadGenLepton>-9", "phiSubleadGenLepton>-9", "ptLeadGenQuark>-9", "etaLeadGenQuark>-9", "phiLeadGenQuark>-9", "ptSubleadGenQuark>-9", "etaSubleadGenQuark>-9", "phiSubleadGenQuark>-9","ptRecoLeptMatchedToWrDau>-9", "etaRecoLeptMatchedToWrDau>-9", "phiRecoLeptMatchedToWrDau>-9", "ptRecoLeptMatchedToNuDau>-9", "etaRecoLeptMatchedToNuDau>-9", "phiRecoLeptMatchedToNuDau>-9", "ptRecoJetOneMatchedToNuDau>-9", "etaRecoJetOneMatchedToNuDau>-9", "phiRecoJetOneMatchedToNuDau>-9", "ptRecoJetTwoMatchedToNuDau>-9", "etaRecoJetTwoMatchedToNuDau>-9", "phiRecoJetTwoMatchedToNuDau>-9", "ptGenJetFromMatchedRecoJetOne>-9", "etaGenJetFromMatchedRecoJetOne>-9", "phiGenJetFromMatchedRecoJetOne>-9", "ptGenJetFromMatchedRecoJetTwo>-9", "etaGenJetFromMatchedRecoJetTwo>-9", "phiGenJetFromMatchedRecoJetTwo>-9", "ptLeadRecoLept>-9", "etaLeadRecoLept>-9", "phiLeadRecoLept>-9", "ptSubleadRecoLept>-9", "etaSubleadRecoLept>-9", "phiSubleadRecoLept>-9", "ptLeadRecoJet>-9", "etaLeadRecoJet>-9", "phiLeadRecoJet>-9", "ptSubleadRecoJet>-9", "etaSubleadRecoJet>-9", "phiSubleadRecoJet>-9"};
	//string plotTitle[] = {"GEN W_{R} Multiplicity","GEN Nu Multiplicity","GEN Lepton Multiplicity","GEN Quark Multiplicity","Leading GEN Lepton is not W_{R} daughter","Subleading GEN Lepton is not Nu daughter","Leading GEN quark is not Nu daughter","Subleading GEN quark is not Nu daughter","GEN W_{R} #eta","GEN W_{R} P_{T}","GEN W_{R} Mass","GEN Nu #eta","GEN Nu P_{T}","GEN Nu Mass","#eta of GEN Lepton from W_{R}","P_{T} of GEN Lepton from W_{R}","#phi of GEN Lepton from W_{R}","#eta of GEN Lepton from Nu","P_{T} of GEN Lepton from Nu","#phi of GEN Lepton from Nu","#eta of first GEN Quark from Nu","P_{T} of first GEN Quark from Nu","#phi of first GEN Quark from Nu","#eta of second GEN Quark from Nu","P_{T} of second GEN Quark from Nu","#phi of second GEN Quark from Nu","Leading GEN Lepton P_{T}","Leading GEN Lepton #eta","Leading GEN Lepton #phi","Subleading GEN Lepton P_{T}","Subleading GEN Lepton #eta","Subleading GEN Lepton #phi","Leading GEN Quark P_{T}","Leading GEN Quark #eta","Leading GEN Quark #phi","Subleading GEN Quark P_{T}","Subleading GEN Quark #eta","Subleading GEN Quark #phi","P_{T} of RECO Lepton matched to W_{R} daughter","#eta of RECO Lepton matched to W_{R} daughter","#phi of RECO Lepton matched to W_{R} daughter","P_{T} of RECO Lepton matched to Nu daughter","#eta of RECO Lepton matched to Nu daughter","#phi of RECO Lepton matched to Nu daughter","P_{T} of first RECO Jet matched to Nu daughter","#eta of first RECO Jet matched to Nu daughter","#phi of first RECO Jet matched to Nu daughter","P_{T} of second RECO Jet matched to Nu daughter","#eta of second RECO Jet matched to Nu daughter","#phi of second RECO Jet matched to Nu daughter","P_{T} of first GEN Jet matched to Nu daughter","#eta of first GEN Jet matched to Nu daughter","#phi of first GEN Jet matched to Nu daughter","P_{T} of second GEN Jet matched to Nu daughter","#eta of second GEN Jet matched to Nu daughter","#phi of second GEN Jet matched to Nu daughter","Leading RECO Lepton P_{T}","Leading RECO Lepton #eta","Leading RECO Lepton #phi","Subleading RECO Lepton P_{T}","Subleading RECO Lepton #eta","Subleading RECO Lepton #phi","Leading RECO Jet P_{T}","Leading RECO Jet #eta","Leading RECO Jet #phi","Subleading RECO Jet P_{T}","Subleading RECO Jet #eta","Subleading RECO Jet #phi","GEN W_{R} Rapidity","GEN Nu Rapidity"};
	//string plotXaxisLabel[] = {"GEN W_{R} Multiplicity","GEN Nu Multiplicity","GEN Lepton Multiplicity","GEN Quark Multiplicity","1 = lead GEN lepton not W_{R} dau","1 = sublead GEN lepton not Nu dau","1 = lead GEN quark not Nu dau","1 = sublead GEN quark not Nu dau","GEN W_{R} #eta","GEN W_{R} P_{T} [GeV]","GEN W_{R} Mass [GeV]","GEN Nu #eta","GEN Nu P_{T} [GeV]","GEN Nu Mass [GeV]","#eta of GEN Lepton from W_{R}","P_{T} [GeV]","#phi","#eta","P_{T} [GeV]","#phi","#eta","P_{T} [GeV]","#phi","#eta","P_{T} [GeV]","#phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi","Lepton P_{T} [GeV]","Lepton #eta","Lepton #phi","Lepton P_{T} [GeV]","Lepton #eta","Lepton #phi","Lead Jet P_{T} [GeV]","Lead Jet #eta","Lead Jet #phi","Sublead Jet P_{T} [GeV]","Sublead Jet #eta","Sublead Jet #phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi","P_{T} [GeV]","#eta","#phi"};

	//default list of all gen branches
	string plotArg[] = {"numGenFstHvyPtcl","numGenScdHvyPtcl","numGenLeptons","numGenQuarks","leadGenLeptonNotFromFstHvyPtcl","subleadGenLeptonNotFromScdHvyPtcl","leadGenQuarkNotFromScdHvyPtcl","subleadGenQuarkNotFromScdHvyPtcl","etaGenFstHvyPtcl", "ptGenFstHvyPtcl", "massGenFstHvyPtcl", "etaGenScdHvyPtcl", "ptGenScdHvyPtcl", "massGenScdHvyPtcl", "etaGenLeptFromFstHvyPtcl", "ptGenLeptFromFstHvyPtcl", "phiGenLeptFromFstHvyPtcl", "etaGenLeptFromScdHvyPtcl", "ptGenLeptFromScdHvyPtcl", "phiGenLeptFromScdHvyPtcl", "etaGenQuarkOneFromScdHvyPtcl", "ptGenQuarkOneFromScdHvyPtcl", "phiGenQuarkOneFromScdHvyPtcl", "etaGenQuarkTwoFromScdHvyPtcl", "ptGenQuarkTwoFromScdHvyPtcl", "phiGenQuarkTwoFromScdHvyPtcl","ptLeadGenLepton", "etaLeadGenLepton", "phiLeadGenLepton", "ptSubleadGenLepton", "etaSubleadGenLepton", "phiSubleadGenLepton", "ptLeadGenQuark", "etaLeadGenQuark", "phiLeadGenQuark", "ptSubleadGenQuark", "etaSubleadGenQuark", "phiSubleadGenQuark", "fourObjMassFromGenObjsFromFstAndScdHvyPtcl", "dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl", "dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl", "dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl", "dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl", "dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl"};

	//standard plotting cuts to ensure existence of GEN particles
	string plotCut[] = {"numGenFstHvyPtcl>-1","numGenScdHvyPtcl>-1","numGenLeptons>-1","numGenQuarks>-1","leadGenLeptonNotFromFstHvyPtcl>-1","subleadGenLeptonNotFromScdHvyPtcl>-1","leadGenQuarkNotFromScdHvyPtcl>-1","subleadGenQuarkNotFromScdHvyPtcl>-1","etaGenFstHvyPtcl>-9", "ptGenFstHvyPtcl>-9", "massGenFstHvyPtcl>-9", "etaGenScdHvyPtcl>-9", "ptGenScdHvyPtcl>-9", "massGenScdHvyPtcl>-9", "etaGenLeptFromFstHvyPtcl>-9", "ptGenLeptFromFstHvyPtcl>-9", "phiGenLeptFromFstHvyPtcl>-9", "etaGenLeptFromScdHvyPtcl>-9", "ptGenLeptFromScdHvyPtcl>-9", "phiGenLeptFromScdHvyPtcl>-9", "etaGenQuarkOneFromScdHvyPtcl>-9", "ptGenQuarkOneFromScdHvyPtcl>-9", "phiGenQuarkOneFromScdHvyPtcl>-9", "etaGenQuarkTwoFromScdHvyPtcl>-9", "ptGenQuarkTwoFromScdHvyPtcl>-9", "phiGenQuarkTwoFromScdHvyPtcl>-9","ptLeadGenLepton>-9", "etaLeadGenLepton>-9", "phiLeadGenLepton>-9", "ptSubleadGenLepton>-9", "etaSubleadGenLepton>-9", "phiSubleadGenLepton>-9", "ptLeadGenQuark>-9", "etaLeadGenQuark>-9", "phiLeadGenQuark>-9", "ptSubleadGenQuark>-9", "etaSubleadGenQuark>-9", "phiSubleadGenQuark>-9", "fourObjMassFromGenObjsFromFstAndScdHvyPtcl>-9", "dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>-9", "dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>-9", "dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>-9", "dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>-9", "dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>-9"};

	//standard X axis labels when all GEN branches are being drawn
	string plotXaxisLabel[] = {"GEN W_{R} Multiplicity","GEN Nu Multiplicity","GEN Lepton Multiplicity","GEN Quark Multiplicity","1 = lead GEN lepton not W_{R} dau","1 = sublead GEN lepton not Nu dau","1 = lead GEN quark not Nu dau","1 = sublead GEN quark not Nu dau","GEN W_{R} #eta","GEN W_{R} P_{T} [GeV]","GEN W_{R} Mass [GeV]","GEN Nu #eta","GEN Nu P_{T} [GeV]","GEN Nu Mass [GeV]","#eta of GEN Lepton from W_{R}","P_{T} of GEN Lepton from W_{R} [GeV]","#phi of GEN Lepton from W_{R}","#eta of GEN Lepton from Nu","P_{T} of GEN Lepton from Nu [GeV]","#phi of GEN Lepton from Nu","#eta of GEN quark one from Nu","P_{T} of GEN quark one from Nu [GeV]","#phi of GEN quark one from Nu","#eta of GEN quark two from Nu","P_{T} of GEN quark two from Nu [GeV]","#phi of GEN quark two from Nu","P_{T} of leading GEN lepton [GeV]","#eta of leading GEN lepton","#phi of leading GEN lepton","P_{T} of subleading GEN lepton [GeV]","#eta of subleading GEN lepton","#phi of subleading GEN lepton","P_{T} of leading GEN quark [GeV]","#eta of leading GEN quark","#phi of leading GEN quark","P_{T} of subleading GEN quark [GeV]","#eta of subleading GEN quark","#phi of subleading GEN quark", "GEN M_{LLJJ} [GeV]", "GEN M_{LL} [GeV]", "#DeltaR[Gen Lept from WR, Gen Qrk One from Nu]", "#DeltaR[Gen Lept from WR, Gen Qrk Two from Nu]", "#DeltaR[Gen Lept from Nu, Gen Qrk One from Nu]", "#DeltaR[Gen Lept from Nu, Gen Qrk Two from Nu]"};

	//plotting cuts, arguments, and X axis labels for N-1 efficiencies of lepton and jet pT cuts, dilepton mass cut, dR lepton jet cuts, and four obj mass cut at GEN lvl
	//string cutEffOutputFilePath = "nMinusOneCutEfficienciesGENWR.txt";
	//string plotArg[] = {"ptGenLeptFromFstHvyPtcl", "ptGenLeptFromScdHvyPtcl", "ptGenQuarkOneFromScdHvyPtcl", "ptGenQuarkTwoFromScdHvyPtcl", "fourObjMassFromGenObjsFromFstAndScdHvyPtcl", "dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl", "dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl", "dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl", "dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl", "dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl"};

	//use this stdNMinusOneEffPlotCut in conjunction with other cuts when calculating GEN WR N-1 cut efficiencies
	//string stdNMinusOneEffPlotCut = " && fabs(etaGenLeptFromFstHvyPtcl)<2.4 && fabs(etaGenLeptFromScdHvyPtcl)<2.4 && fabs(etaGenQuarkOneFromScdHvyPtcl)<2.4 && fabs(etaGenQuarkTwoFromScdHvyPtcl)<2.4";
	
	//string plotCut[] = {"ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkTwoFromScdHvyPtcl>40 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4", "ptGenLeptFromFstHvyPtcl>60 && ptGenLeptFromScdHvyPtcl>50 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40 && fourObjMassFromGenObjsFromFstAndScdHvyPtcl>600 && dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl>200 && dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl>0.4 && dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl>0.4 && dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl>0.4"};

	//string plotXaxisLabel[] = {"P_{T} of GEN Lepton from W_{R} [GeV]","P_{T} of GEN Lepton from Nu [GeV]","P_{T} of GEN quark one from Nu [GeV]","P_{T} of GEN quark two from Nu [GeV]", "GEN M_{LLJJ} [GeV]", "GEN dilepton mass [GeV]", "#DeltaR[Gen Lept from WR, Gen Qrk One from Nu]", "#DeltaR[Gen Lept from WR, Gen Qrk Two from Nu]", "#DeltaR[Gen Lept from Nu, Gen Qrk One from Nu]", "#DeltaR[Gen Lept from Nu, Gen Qrk Two from Nu]"};


	string stdPlotTitle = "CMS Private                      #surds = 13 TeV";
	
	vector<string> vectForNEntries(plotCut,plotCut + sizeof(plotCut)/sizeof(string));
	int nBranches = vectForNEntries.size();

	for(int i=0; i<nBins ; i++){
		///loop over WR mass values
		
		///define input root file name and add a Float_t to the x axis array used for later TGraph objects
		string pfn = dir+fileBegin+privOrCent+fileWrMass+to_string(wrMassArr[i])+fileMiddle+to_string(wrMassArr[i]/2)+fileEnd;
		wrMassVals[i] = (Float_t) wrMassArr[i];
		string plotDir = "plotHolder/noCutsGenAndRecoWr_"+privOrCent+"_MWR-"+to_string(wrMassArr[i])+fileMiddle+to_string(wrMassArr[i]/2);
		//TF1 * fitcrv = new TF1("fitcrv","TMath::BreitWigner(x,[0],[1])",(0.7)*wrMassVals[i],(1.3)*wrMassVals[i]);
		//TF1 * fitcrv = new TF1("fitcrv","TMath::BreitWigner(x,[0],[1])");
		//fitcrv->FixParameter(0,wrMassVals[i]);

		TChain * wrChain = new TChain("wrDecayChainAnalyzer/genAndMatchedRecoWrDecayNoCuts", (to_string(wrMassArr[i]) ).c_str() );
		wrChain->Add(pfn.c_str());
	
		for(int j=0; j<nBranches ; j++){
			//make plots from every branch and save them to a directory with a unique name that identifies the WR and Nu masses
			
			//plot GEN WR mass with Breit Wigner fit overlaid on histo
			//if(plotArg[j] == "massGenFstHvyPtcl") makeAndSaveSingleHistoFromTreeWithFit(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(3000,"+to_string(wrMassArr[i]/2)+","+to_string((1.5)*wrMassArr[i]) +")",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_withBWFit", fitcrv);
		
			/*
			//determine efficiency of leading and subleading lepton and jet pt cuts
			if(plotArg[j] == "ptGenLeptFromFstHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],60.,true,cutEffOutputFilePath,"ptGenLeptFromFstHvyPtcl",false);
			if(plotArg[j] == "ptGenLeptFromScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],50.,true,cutEffOutputFilePath,"ptGenLeptFromScdHvyPtcl",false);
			if(plotArg[j] == "ptGenQuarkOneFromScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],40.,true,cutEffOutputFilePath,"ptGenQuarkOneFromScdHvyPtcl",false);
			if(plotArg[j] == "ptGenQuarkTwoFromScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],40.,true,cutEffOutputFilePath,"ptGenQuarkTwoFromScdHvyPtcl",false);

			//calculate efficiency of four object mass cut
			if(plotArg[j] == "fourObjMassFromGenObjsFromFstAndScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],600.,true,cutEffOutputFilePath,"fourObjMassFromGenObjsFromFstAndScdHvyPtcl",false);

			//efficiency of MLL cut
			if(plotArg[j] == "dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],200.,true,cutEffOutputFilePath,"dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl",false);

			//efficiency of dR lepton jet cuts
			if(plotArg[j] == "dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.4,true,cutEffOutputFilePath,"dRgenLeptonFromFstHvyPtclGenQuarkOneFromScdHvyPtcl",false);
			if(plotArg[j] == "dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.4,true,cutEffOutputFilePath,"dRgenLeptonFromFstHvyPtclGenQuarkTwoFromScdHvyPtcl",false);
			if(plotArg[j] == "dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.4,true,cutEffOutputFilePath,"dRgenLeptonFromScdHvyPtclGenQuarkOneFromScdHvyPtcl",false);
			if(plotArg[j] == "dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.4,true,cutEffOutputFilePath,"dRgenLeptonFromScdHvyPtclGenQuarkTwoFromScdHvyPtcl",false);
			*/
			
			//RETURN
			//make and save plots from each branch of the input tree
			if(plotArg[j] == "dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl" && wrMassArr[i] == 800) makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(50,0.,850.)",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, false);
			if(plotArg[j] == "dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl" && wrMassArr[i] == 2600) makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(50,0.,2600.)",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, false);
			if(plotArg[j] == "dileptonMassFromGenLeptonsFromFstAndScdHvyPtcl" && wrMassArr[i] == 5000) makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(50,0.,5000.)",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, false);
	
			//if(plotArg[j] == "ptRecoLeptMatchedToWrDau") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(50,0.,1100.)",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, false);
			//if(plotArg[j] == "ptRecoLeptMatchedToNuDau") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(50,0.,1100.)",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, false);
			//if(plotArg[j] == "ptRecoJetOneMatchedToNuDau") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(50,0.,1100.)",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, false);
			//if(plotArg[j] == "ptRecoJetTwoMatchedToNuDau") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist(40,0.,600.)",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, false);


			//else if(plotArg[j] != "ptRecoJetTwoMatchedToNuDau" && plotArg[j] != "ptRecoJetOneMatchedToNuDau" && plotArg[j] != "ptRecoLeptMatchedToNuDau" && plotArg[j] != "ptRecoLeptMatchedToWrDau") makeAndSaveSingleHistoFromTreeWithCuts(wrChain,"c"+to_string(j),plotCut[j],plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j],0.,false,"","",false, true);

		}//end loop over TChain branches


		/*
		//calculate the fraction of WR->LLJJ events in which the GEN leptons fall into the barrel or the endcaps or a mix
		//ignore the dead zone region between 1.444 and 1.566, and restrict the max eta to be less than 2.4
		Float_t eventsIgnoringDeadZoneAndHighEtaEndcap = (Float_t) wrChain->GetEntries("(fabs(etaGenLeptFromFstHvyPtcl) < 1.444 || (fabs(etaGenLeptFromFstHvyPtcl) > 1.566 && fabs(etaGenLeptFromFstHvyPtcl) < 2.4) ) && (fabs(etaGenLeptFromScdHvyPtcl) < 1.444 || (fabs(etaGenLeptFromScdHvyPtcl) > 1.566 && fabs(etaGenLeptFromScdHvyPtcl) < 2.4) )");
		Float_t eventsWithBothInBarrel = (Float_t) wrChain->GetEntries("fabs(etaGenLeptFromFstHvyPtcl) < 1.444 && fabs(etaGenLeptFromScdHvyPtcl) < 1.444");
		Float_t eventsWithOneBarrelOneEndcap = (Float_t) (wrChain->GetEntries("fabs(etaGenLeptFromFstHvyPtcl) < 1.444 && (fabs(etaGenLeptFromScdHvyPtcl) > 1.566 && fabs(etaGenLeptFromScdHvyPtcl) < 2.4)") + wrChain->GetEntries("fabs(etaGenLeptFromScdHvyPtcl) < 1.444 && (fabs(etaGenLeptFromFstHvyPtcl) > 1.566 && fabs(etaGenLeptFromFstHvyPtcl) < 2.4)"));
		Float_t eventsWithBothInEndcap = (Float_t) wrChain->GetEntries("(fabs(etaGenLeptFromFstHvyPtcl) > 1.566 && fabs(etaGenLeptFromFstHvyPtcl) < 2.4) && (fabs(etaGenLeptFromScdHvyPtcl) > 1.566 && fabs(etaGenLeptFromScdHvyPtcl) < 2.4)");
		Float_t bothBarrelEfficiency = eventsWithBothInBarrel/eventsIgnoringDeadZoneAndHighEtaEndcap;
		Float_t oneBarrelOneEndcapEfficiency = eventsWithOneBarrelOneEndcap/eventsIgnoringDeadZoneAndHighEtaEndcap;
		Float_t bothEndcapEfficiency = eventsWithBothInEndcap/eventsIgnoringDeadZoneAndHighEtaEndcap;
		
		if(i==0) writeToBarrelEndcapFile << "GEN acceptance requires both GEN leptons from the WR and Nu to have |eta| < 2.4 and outside the transition region 1.444 < |eta| < 1.566" << std::endl;;
		writeToBarrelEndcapFile << "" << std::endl;
		writeToBarrelEndcapFile << "WR mass=\t"<< wrMassArr[i] <<"\tNu mass=\t"<< wrMassArr[i]/2 << std::endl;
		writeToBarrelEndcapFile << "the percentage of events passing GEN acceptance with both leptons in the barrel=\t"<< 100*bothBarrelEfficiency << std::endl;
		writeToBarrelEndcapFile << "the percentage of events passing GEN acceptance with one lepton in the barrel and one lepton in the endcaps=\t"<< 100*oneBarrelOneEndcapEfficiency << std::endl;
		writeToBarrelEndcapFile << "the percentage of events passing GEN acceptance with both leptons in the endcaps=\t"<< 100*bothEndcapEfficiency << std::endl;
		*/

		/*
		TString baseGenMatchedCutString = "TMath::Abs(etaGenLeptFromFstHvyPtcl)<8 && TMath::Abs(etaGenLeptFromScdHvyPtcl)<8 && TMath::Abs(etaGenQuarkOneFromScdHvyPtcl)<8 && TMath::Abs(etaGenQuarkTwoFromScdHvyPtcl)<8";
		TString genMatchedEtaCutString = "TMath::Abs(etaGenLeptFromFstHvyPtcl)<2.4 && TMath::Abs(etaGenLeptFromScdHvyPtcl)<2.4 && TMath::Abs(etaGenQuarkOneFromScdHvyPtcl)<2.4 && TMath::Abs(etaGenQuarkTwoFromScdHvyPtcl)<2.4";
		TString genMatchedPtCutString = "ptGenLeptFromFstHvyPtcl>40 && ptGenLeptFromScdHvyPtcl>40 && ptGenQuarkOneFromScdHvyPtcl>40 && ptGenQuarkTwoFromScdHvyPtcl>40";
		
		TString baseGenLeadingCutString = "TMath::Abs(etaLeadGenLepton)<8 && TMath::Abs(etaSubleadGenLepton)<8 && TMath::Abs(etaLeadGenQuark)<8 && TMath::Abs(etaSubleadGenQuark)<8";
		TString genLeadingEtaCutString = "TMath::Abs(etaLeadGenLepton)<2.4 && TMath::Abs(etaSubleadGenLepton)<2.4 && TMath::Abs(etaLeadGenQuark)<2.4 && TMath::Abs(etaSubleadGenQuark)<2.4";
		TString genLeadingPtCutString = "ptLeadGenLepton>40 && ptSubleadGenLepton>40 && ptLeadGenQuark>40 && ptSubleadGenQuark>40";

		//genMatchedPtEtaEff[i] = (100)*(calculateEfficiencyWithOneChain(wrChain,"evWeight", baseGenMatchedCutString, genMatchedEtaCutString+" && "+genMatchedPtCutString));
		//genLeadAndSubleadPtEtaEff[i] = (100)*(calculateEfficiencyWithOneChain(wrChain,"evWeight", baseGenLeadingCutString, genLeadingEtaCutString+" && "+genLeadingPtCutString));

		writeToGenEfficiencyFile << genMatchedPtEtaEff[i] <<"\tpercent of events with WR mass=\t"<< wrMassArr[i] <<"\tand Nu mass=\t"<< wrMassArr[i]/2 <<"\tpass the requirement that the gen leptons and quarks from WR and Nu decays have |eta| < 2.5 and pt > 40"<< endl;
		writeToGenEfficiencyFile << genLeadAndSubleadPtEtaEff[i] <<"\tpercent of events with WR mass=\t"<< wrMassArr[i] <<"\tand Nu mass=\t"<< wrMassArr[i]/2 <<"\tpass the requirement that lead and sublead gen leptons and quarks have |eta| < 2.5 and pt > 40"<< endl;
		writeToGenEfficiencyFile <<" "<< endl;

		SaveTreePlots(wrChain, plotDir.c_str());
		*/

		wrChain->Delete();
		//fitcrv->Delete();

	}///end loop over WR mass values	

	///close txt file
	//writeToGenEfficiencyFile.close();
	//writeToBarrelEndcapFile.close();


#endif
	//end genAndRecoWrPlotsMinimalCuts


// /afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRGen/analyzedGen/withoutGenNuFilter/8TeV/analyzed_genWrToMuMuJJFullOfflineAnalysis_WR_2200_NU_1100_1.root	
#ifdef privateGenEff
	//use this ifdef to calculate the fraction of GEN evts in which the GEN WR daughter leptons and quarks pass a certain cut
	///all input .root files should be in the same directory, and have file names which differ only in the WR and Nu mass values
	//string dir= "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRGen/analyzedGen/withoutGenNuFilter/";
	string dir= "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/privateWRGen/analyzedGen/withoutGenNuFilter/8TeV/";
	string fileBegin = "analyzed_genWrToMuMuJJFullOfflineAnalysis_WR_";
	string fileEnd = "_1.root";
	string fileMiddle = "_NU_";
	string genCutEffVsMassFile = "genCutEfficienciesVsMasses.txt";
	ofstream writeToGenEfficiencyFile(genCutEffVsMassFile.c_str(),ofstream::trunc);
	gStyle->SetTitleOffset(1.4,"Y");
	Int_t nBins = 7;
	Float_t wrMassVals[nBins], genMatchedEtaEff[nBins];
	gStyle->SetOptStat("");

	int wrMassArr[] = {1000,2000,2200,2400,2600,2800,3000};

	string stdPlotTitle = "CMS Private                      #surds = 13 TeV";
	
	for(int i=0; i<nBins ; i++){
		///loop over WR mass values
		
		///define input root file name and add a Float_t to the x axis array used for later TGraph objects
		string pfn = dir+fileBegin+to_string(wrMassArr[i])+fileMiddle+to_string(wrMassArr[i]/2)+fileEnd;
		wrMassVals[i] = (Float_t) wrMassArr[i];
		string plotDir = "plotHolder/noCutsGenAndRecoWr_private_MWR-"+to_string(wrMassArr[i])+fileMiddle+to_string(wrMassArr[i]/2);

		//as long as the tree comes from genMatchedParticleAnalyzerOne or a later module defined in checkWRDecay_crabSafe_cfg.py, a cut which
		//requires the presence of two GEN leptons and quarks whose mothers are the WR and Nu is not needed
		TChain * wrChain = new TChain("genMatchedParticleAnalyzerOne/genLeptonsAndJetsNoCuts", (to_string(wrMassArr[i]) ).c_str() );
		wrChain->Add(pfn.c_str());
		
		TString baseGenMatchedCutString = "";
		TString genMatchedEtaCutString = "TMath::Abs(etaEle[0])<2.4 && TMath::Abs(etaEle[1])<2.4 && TMath::Abs(etaJet[0])<2.4 && TMath::Abs(etaJet[1])<2.4";
	
		Float_t eff = -1., effUnc = 0.;
		calculateEfficiencyWithOneChain(wrChain,"evWeight", baseGenMatchedCutString, genMatchedEtaCutString, effUnc, eff);
		genMatchedEtaEff[i] = 100*eff;

		writeToGenEfficiencyFile << genMatchedEtaEff[i] <<"\tpercent of events with WR mass=\t"<< wrMassArr[i] <<"\tand Nu mass=\t"<< wrMassArr[i]/2 <<"\tpass the requirement that the gen leptons and quarks from WR and Nu decays have |eta| < 2.4"<< endl;
		writeToGenEfficiencyFile <<" "<< endl;

		wrChain->Delete();

	}///end loop over WR mass values	

	///close txt file
	writeToGenEfficiencyFile.close();

#endif
	//end privateGenEff



#ifdef nMinusOneCutEffsGenAndReco
	///all input .root files should be in the same directory, and have file names which differ only in the WR and Nu mass values
	string dir= "../rootfilesfornminusoneefficiencies/";

	//inputFiles and processTags must have the same number of elements
	//string inputFiles[] = {"ttbarEEForNminusOne_withFlavorFilters.root", "dyjetsEEForNminusOne_withFlavorFilters.root", "wrEE_2000_1000_ForNminusOne_withFlavorFilters.root","ttbarMuMuForNminusOne_withFlavorFilters.root", "dyjetsMuMuForNminusOne_withFlavorFilters.root", "wrMuMu_2000_1000_ForNminusOne_withFlavorFilters.root","wrEE_1000_500_ForNminusOne_withFlavorFilters.root", "wrMuMu_1000_500_ForNminusOne_withFlavorFilters.root"};
	//string processTags[] = {"TTtoEE","DYtoEE","WRtoEE_2000_1000","TTtoMuMu","DYtoMuMu","WRtoMuMu_2000_1000","WRtoEE_1000_500","WRtoMuMu_1000_500"};
	string inputFiles[] = {"wrEE_800_400_ForNminusOne_withFlavorFilters.root", "wrMuMu_800_400_ForNminusOne_withFlavorFilters.root"};
	string processTags[] = {"WRtoEE_800_400","WRtoMuMu_800_400"};


	vector<string> inputFilesVect(inputFiles,inputFiles + sizeof(inputFiles)/sizeof(string));
	vector<string> processTagsVect(processTags,processTags + sizeof(processTags)/sizeof(string));
	int nBins = inputFilesVect.size();
	if(nBins != processTagsVect.size()){
		cout<<"inputFilesVect and processTagsVect do not have the same number of elements. exiting now"<<endl;
		exit(1);
	}
	gStyle->SetTitleOffset(1.4,"Y");
	gStyle->SetOptStat("");

	//chainNames contains the names of TChains stored in each input root file
	//chainNameVect and chainTagVect must have the same number of elements
	//string chainNames[] = {"analyzerTwo/genKinematicsUsingGenQuarksWithGenMotherRequirements","analyzerFour/genKinematicsUsingGenJetsWithGenMotherRequirements","analyzerFive/recoKinematics","analyzerSix/genKinematicsUsingGenQuarksWithoutGenMotherRequirementsWithGenFlavorReqs","analyzerSeven/genKinematicsUsingGenJetsWithoutGenMotherRequirementsWithGenFlavorReqs"};
	//string chainTags[] = {"genWithQuarksWithMatching","genWithJetsWithMatching","recoPassingJetLeptonId","genWithQuarksWithGenFlavorReqs","genWithJetsWithGenFlavorReqs"};
	string chainNames[] = {"analyzerFour/genKinematicsUsingGenJetsWithGenMotherRequirements","analyzerFive/recoKinematics",};
	string chainTags[] = {"genWithJetsWithMatching","recoPassingJetLeptonId"};
	
	vector<string> chainNameVect(chainNames,chainNames + sizeof(chainNames)/sizeof(string));
	int nChains = chainNameVect.size();
	vector<string> chainTagVect(chainTags,chainTags + sizeof(chainTags)/sizeof(string));
	if(nChains != chainTagVect.size()){
		cout<<"chainNameVect and chainTagVect do not have the same number of elements. exiting now"<<endl;
		exit(1);
	}
	

	//element number i in plotArg is tied to element number i in plotCut   don't change the order
	//plotting cuts, arguments, and X axis labels for N-1 efficiencies of lepton and jet pT cuts, dilepton mass cut, dR lepton jet cuts, and four obj mass cut at GEN lvl
	string cutEffOutputFilePath = "nMinusOneCutEfficienciesGENandRECO_DY_TTBar_WR.txt";
	string plotArg[] = {"etaLeptOne", "ptLeptOne", "etaLeptTwo", "ptLeptTwo", "etaQuarkOne", "ptQuarkOne", "etaQuarkTwo", "ptQuarkTwo", "dileptonMass", "dileptonPt", "dileptonEta", "dileptonPhi", "subleadLeptonBothHadronsMass", "subleadLeptonBothHadronsPt", "subleadLeptonBothHadronsEta", "subleadLeptonBothHadronsPhi", "leadLeptonBothHadronsMass", "leadLeptonBothHadronsPt", "leadLeptonBothHadronsEta", "leadLeptonBothHadronsPhi", "dileptonDihadronMass", "dileptonDihadronPt", "dileptonDihadronEta", "dileptonDihadronPhi", "dRleptonOneQuarkOne", "dRleptonOneQuarkTwo", "dRleptonTwoQuarkOne", "dRleptonTwoQuarkTwo"};

	//use this stdNMinusOneEffPlotCut in conjunction with other cuts when calculating N-1 cut efficiencies
	string stdNMinusOneEffPlotCut = " && fabs(etaLeptOne)<2.4 && fabs(etaLeptTwo)<2.4 && fabs(etaQuarkTwo)<2.4 && fabs(etaQuarkOne)<2.4";

	//use allCutsExceptStd to apply all offline cuts when plotting variables like dileptonPt or etaLeptTwo for which we are not trying to measure the N-1 efficiency
	string allCutsExceptStd = "ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4";
	string plotCut[] = {allCutsExceptStd,"ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4",allCutsExceptStd,"ptLeptOne>60 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4",allCutsExceptStd,"ptLeptOne>60 && ptLeptTwo>50 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4",allCutsExceptStd,"ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4","ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4",allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,"ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4",allCutsExceptStd,allCutsExceptStd,allCutsExceptStd,"ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4","ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4","ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkTwo>0.4","ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4"};

	string plotXaxisLabel[] = {"Lead lepton #eta", "Lead lepton P_{T} [GeV]","Sublead lepton #eta","Sublead lepton P_{T} [GeV]","Lead hadron #eta","Lead hadron P_{T} [GeV]","Sublead hadron #eta","Sublead hadron P_{T} [GeV]","Dilepton system mass [GeV]", "Dilepton system P_{T}","Dilepton system #eta","Dilepton system #phi","Sublead lepton plus dihadron system mass [GeV]","Sublead lepton plus dihadron system P_{T} [GeV]","Sublead lepton plus dihadron system #eta","Sublead lepton plus dihadron system #phi","Lead lepton plus dihadron system mass [GeV]","Lead lepton plus dihadron system P_{T} [GeV]","Lead lepton plus dihadron system #eta","Lead lepton plus dihadron system #phi","Dilepton Dihadron system mass M_{LLJJ} [GeV]","Dilepton Dihadron system P_{T} [GeV]","Dilepton Dihadron system #eta","Dilepton Dihadron system #phi", "#DeltaR[Lead lepton, Lead hadron]", "#DeltaR[Lead lepton, Sublead hadron]", "#DeltaR[Sublead lepton, Lead hadron]", "#DeltaR[Sublead lepton, Sublead hadron]"};


	string stdPlotTitle = "CMS Preliminary                 #surds = 13 TeV";
	string plotDir = "plotHolder/nMinusOneCutEfficiencies";
	
	vector<string> vectForNEntries(plotCut,plotCut + sizeof(plotCut)/sizeof(string));
	int nBranches = vectForNEntries.size();

	for(int i=0; i<nBins ; i++){
		///loop over input root files

		for(int r=0; r<nChains; r++){
			//loop over TChains within each input root file
	
			TChain * inputChain = new TChain(chainNameVect[r].c_str(),(processTagsVect[i]+"_"+chainTagVect[r]).c_str());
			inputChain->Add( (dir+inputFilesVect[i]).c_str() );

			for(int j=0; j<nBranches ; j++){
				//calculate N-1 efficiency of lepton and hadron pt cuts
				if(plotArg[j] == "ptLeptOne") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],60.,true,cutEffOutputFilePath,"ptLeptOne",false, true);
				if(plotArg[j] == "ptLeptTwo") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],50.,true,cutEffOutputFilePath,"ptLeptTwo",false, true);
				if(plotArg[j] == "ptQuarkOne") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],40.,true,cutEffOutputFilePath,"ptQuarkOne",false, true);
				if(plotArg[j] == "ptQuarkTwo") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],40.,true,cutEffOutputFilePath,"ptQuarkTwo",false, true);

				//calculate N-1 efficiency of four object mass cut
				if(plotArg[j] == "dileptonDihadronMass") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],600.,true,cutEffOutputFilePath,"dileptonDihadronMass",false, true);

				//N-1 efficiency of dilepton mass cut
				if(plotArg[j] == "dileptonMass") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],200.,true,cutEffOutputFilePath,"dileptonMass",false, true);

				//N-1 efficiency of dR lepton hadron cuts
				if(plotArg[j] == "dRleptonOneQuarkOne") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],0.4,true,cutEffOutputFilePath,"dRleptonOneQuarkOne",false, true);
				if(plotArg[j] == "dRleptonOneQuarkTwo") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],0.4,true,cutEffOutputFilePath,"dRleptonOneQuarkTwo",false, true);
				if(plotArg[j] == "dRleptonTwoQuarkOne") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],0.4,true,cutEffOutputFilePath,"dRleptonTwoQuarkOne",false, true);
				if(plotArg[j] == "dRleptonTwoQuarkTwo") makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],0.4,true,cutEffOutputFilePath,"dRleptonTwoQuarkTwo",false, true);

				//save plots from all other variables which are not used in offline cuts
				else makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r],0.,false,"","",false, true);

			}//end loop over TChain branches
			inputChain->Delete();
		}//end loop over TChains within each input root file

	}///end loop over input files	

#endif
	//end nMinusOneCutEffsGenAndReco


#ifdef multiStepCutEffsWRandBkgnds
	/**use this ifdef to calculate the efficiency of several cuts in a fixed order:
	 * the efficiency of events to pass the pT and eta cuts on leptons and jets, given events which have two jets and leptons passing ID
	 * the efficiency of evts to pass the dR(L,J) cuts, given events which pass the pT and eta cuts on leptons and jets
	 * the efficiency of evts to pass the combined MLL and MLLJJ cuts, given evts which pass the dR(L,J) and pT and eta cuts
	 **/

	///all input .root files should be in the same directory, and have file names which differ only in the WR and Nu mass values
	string dir= "../rootfilesfornminusoneefficiencies/";

	//inputFiles and processTags must have the same number of elements
	string inputFiles[] = {"ttbarEEForNminusOne_withFlavorFilters.root", "dyjetsEEForNminusOne_withFlavorFilters.root", "wrEE_2000_1000_ForNminusOne_withFlavorFilters.root","ttbarMuMuForNminusOne_withFlavorFilters.root", "dyjetsMuMuForNminusOne_withFlavorFilters.root", "wrMuMu_2000_1000_ForNminusOne_withFlavorFilters.root"};
	string processTags[] = {"TTtoEE","DYtoEE","WRtoEE_2000_1000","TTtoMuMu","DYtoMuMu","WRtoMuMu_2000_1000"};
	//string inputFiles[] = {"wrEE_2000_1000_ForNminusOne_withFlavorFilters.root", "wrMuMu_2000_1000_ForNminusOne_withFlavorFilters.root"};
	//string processTags[] = {"WRtoEE_2000_1000","WRtoMuMu_2000_1000"};
	//string inputFiles[] = {"wrMuMu_2000_1000_ForNminusOne_withFlavorFilters.root"};
	//string processTags[] = {"WRtoMuMu_2000_1000"};


	vector<string> inputFilesVect(inputFiles,inputFiles + sizeof(inputFiles)/sizeof(string));
	vector<string> processTagsVect(processTags,processTags + sizeof(processTags)/sizeof(string));
	int nBins = inputFilesVect.size();
	if(nBins != processTagsVect.size()){
		cout<<"inputFilesVect and processTagsVect do not have the same number of elements. exiting now"<<endl;
		exit(1);
	}
	gStyle->SetTitleOffset(1.4,"Y");
	gStyle->SetOptStat("");

	//chainNames contains the names of TChains stored in each input root file
	//chainNameVect and chainTagVect must have the same number of elements
	//string chainNames[] = {"analyzerTwo/genKinematicsUsingGenQuarksWithGenMotherRequirements","analyzerFour/genKinematicsUsingGenJetsWithGenMotherRequirements","analyzerFive/recoKinematics","analyzerSix/genKinematicsUsingGenQuarksWithoutGenMotherRequirementsWithGenFlavorReqs","analyzerSeven/genKinematicsUsingGenJetsWithoutGenMotherRequirementsWithGenFlavorReqs"};
	//string chainTags[] = {"genWithQuarksWithMatching","genWithJetsWithMatching","recoPassingJetLeptonId","genWithQuarksWithGenFlavorReqs","genWithJetsWithGenFlavorReqs"};
	string chainNames[] = {"analyzerFour/genKinematicsUsingGenJetsWithGenMotherRequirements","analyzerFive/recoKinematics",};
	string chainTags[] = {"genWithJetsWithMatching","recoPassingJetLeptonId"};
	
	vector<string> chainNameVect(chainNames,chainNames + sizeof(chainNames)/sizeof(string));
	int nChains = chainNameVect.size();
	vector<string> chainTagVect(chainTags,chainTags + sizeof(chainTags)/sizeof(string));
	if(nChains != chainTagVect.size()){
		cout<<"chainNameVect and chainTagVect do not have the same number of elements. exiting now"<<endl;
		exit(1);
	}
	

	//element number i in plotArg is tied to element number i in plotCut   don't change the order
	//plotting cuts, arguments, and X axis labels for N-1 efficiencies of lepton and jet pT cuts, dilepton mass cut, dR lepton jet cuts, and four obj mass cut at GEN lvl
	string cutEffOutputFilePath = "cutFlowEfficienciesGENandRECO_DY_TTBar_WR.txt";
	string plotArg[] = {"dileptonMass", "ptLeptOne", "ptLeptTwo"};

	//use this stdNMinusOneEffPlotCut in conjunction with other cuts when calculating N-1 cut efficiencies
	string stdNMinusOneEffPlotCut = " && fabs(etaLeptOne)<2.4 && fabs(etaLeptTwo)<2.4 && fabs(etaQuarkTwo)<2.4 && fabs(etaQuarkOne)<2.4";

	//use allCutsExceptStd to apply all offline cuts when plotting variables like dileptonPt or etaLeptTwo for which we are not trying to measure the N-1 efficiency
	string allCutsExceptStd = "ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dileptonMass>200 && dileptonDihadronMass>600 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4";
	
	string plotCut[] = {"","ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40","ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4"};

	//different variables can be plotted every time, the plots will not be used. only the efficiencies written to the txt file will be used
	string plotXaxisLabel[] = {"Dilepton Mass [GeV]", "Lead lepton P_{T} [GeV]", "Sublead lepton P_{T} [GeV]"};

	string stdPlotTitle = "CMS Private                 #surds = 13 TeV";
	string plotDir = "plotHolder/cutFlowEfficiencies";
	
	vector<string> vectForNEntries(plotCut,plotCut + sizeof(plotCut)/sizeof(string));
	int nBranches = vectForNEntries.size();

	for(int i=0; i<nBins ; i++){
		///loop over input root files

		for(int r=0; r<nChains; r++){
			//loop over TChains within each input root file
	
			TChain * inputChain = new TChain(chainNameVect[r].c_str(),(processTagsVect[i]+"_"+chainTagVect[r]).c_str());
			inputChain->Add( (dir+inputFilesVect[i]).c_str() );

			for(int j=0; j<nBranches ; j++){

				//need a cut which is always true in order for j equal 0 condition to make a plot
				if(j==0) makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),"2>1"+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r]+"_noCuts",0.,true,cutEffOutputFilePath,"ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40",true, true);
				if(j==1) makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r]+"_ptEtaCuts",0.,true,cutEffOutputFilePath,"ptLeptOne>60 && ptLeptTwo>50 && ptQuarkOne>40 && ptQuarkTwo>40 && dRleptonOneQuarkOne>0.4 && dRleptonOneQuarkTwo>0.4 && dRleptonTwoQuarkOne>0.4 && dRleptonTwoQuarkTwo>0.4",true, true);
				if(j==2) makeAndSaveSingleHistoFromTreeWithCuts(inputChain,"c"+to_string(j),plotCut[j]+stdNMinusOneEffPlotCut,plotArg[j]+">>"+plotArg[j]+"Hist()",plotArg[j]+"Hist",stdPlotTitle,plotXaxisLabel[j],plotDir+"_"+plotArg[j]+"_"+processTagsVect[i]+"_"+chainTagVect[r]+"_ptEtaDrCuts",0.,true,cutEffOutputFilePath,allCutsExceptStd,true, true);
		
		
			}//end loop over TChain branches
			inputChain->Delete();
		}//end loop over TChains within each input root file

	}///end loop over input files	

#endif
	//end multiStepCutEffsWRandBkgnds

#ifdef sigRgnWrAndBkgnds

	///use this ifdef to make kinematic plots with three or more curves
	///one curve for ttbar, one for dyjets, and one or more curves for WR signal


#endif
	//end sigRgnWrAndBkgnds

#ifdef makeShiftedMCPileupFiles
	//read the MCPileup.root file in. this file contains one object - a TH1F histo with name pileup (different title)
	//make two copies of the file with different names. In one file, shift the x bin values (nPU) up by 5 percent.
	//In the other file, shift the x bin values (nPU) down by 5 percent.  The names of the histo in these two new
	//files will be pileup, but the file names will be different.
	/*
	TFile * existingPileupFile, upPileupFile, downPileupFile;
	existingPileupFile = new TFile("../data/MCPileup.root","READ");

	//read the pileup histo from existingPileupFile
	TH1F * puShiftUp = ((TH1F*) existingPileupFile->FindObject("pileup"))->Clone("pileup");
	puShiftUp->SetTitle("nPUMCShiftUp");
	Int_t nbins = puShiftUp->GetNbinsX();
	//shift all the nPU x axis values up by 5 percent
	for(Int_t i=0){

	}//end loop over all bins of puShiftUp

	upPileupFile = new TFile("../data/MCPileupShiftUp.root","RECREATE");


	delete upPileupFile;
	
	downPileupFile = new TFile("../data/MCPileupShiftDown.root","RECREATE");


	delete existingPileupFile;
	delete downPileupFile;
*/


#endif
	//end makeShiftedMCPileupFiles
	
	
#ifdef sOverBsensitivity
	//calculate S over B using a WR signal sample and DY plus TT after all signal region selections are applied
	//default SoverB is calculated using default signal region cuts
	TString treeName = "central_value_tree";
	int massWindowIndex = 7;	//linked to a specific mass window. change this index if different MWR hypothesis is used
	Int_t numBkgnds = 2;


	//this vector contains all TChains declared below in a specific order. one signal TChain followed by
	//all relevant bkgnd chains for that channel and set of cuts
	std::vector<TChain*> sigAndBkgnds;
	std::vector<std::string> comments;
	
	//electron channel
	TString defDir = "../analysisCppOutputRootFiles/";
	TChain * defWrSigEE = new TChain(treeName);
	defWrSigEE->Add(defDir+"selected_tree_WRtoEEJJ_2200_1100_signal_eeEE.root");
	TChain * defTTEE = new TChain(treeName);
	defTTEE->Add(defDir+"selected_tree_TT_signal_eeEE.root");
	TChain * defDYEE = new TChain(treeName);
	defDYEE->Add(defDir+"selected_tree_DYAMC_signal_eeEE.root");
	sigAndBkgnds.push_back(defWrSigEE), sigAndBkgnds.push_back(defTTEE), sigAndBkgnds.push_back(defDYEE);
	comments.push_back("default cuts e chnl");

	TString lowJetDir = "../analysisCppOutputRootFiles/withJetPtThirtyCut/";
	TChain * lowJetWrSigEE = new TChain(treeName);
	lowJetWrSigEE->Add(lowJetDir+"selected_tree_WRtoEEJJ_2200_1100_signal_eeEE.root");
	TChain * lowJetTTEE = new TChain(treeName);
	lowJetTTEE->Add(lowJetDir+"selected_tree_TT_signal_eeEE.root");
	TChain * lowJetDYEE = new TChain(treeName);
	lowJetDYEE->Add(lowJetDir+"selected_tree_DYAMC_signal_eeEE.root");
	sigAndBkgnds.push_back(lowJetWrSigEE), sigAndBkgnds.push_back(lowJetTTEE), sigAndBkgnds.push_back(lowJetDYEE);
	comments.push_back("lower jet pT cut e chnl");

	TString lowSubLeptDir = "../analysisCppOutputRootFiles/withSubleadLeptonPtFortyCut/";
	TChain * lowSubLeptWrSigEE = new TChain(treeName);
	lowSubLeptWrSigEE->Add(lowSubLeptDir+"selected_tree_WRtoEEJJ_2200_1100_signal_eeEE.root");
	TChain * lowSubLeptTTEE = new TChain(treeName);
	lowSubLeptTTEE->Add(lowSubLeptDir+"selected_tree_TT_signal_eeEE.root");
	TChain * lowSubLeptDYEE = new TChain(treeName);
	lowSubLeptDYEE->Add(lowSubLeptDir+"selected_tree_DYAMC_signal_eeEE.root");
	sigAndBkgnds.push_back(lowSubLeptWrSigEE), sigAndBkgnds.push_back(lowSubLeptTTEE), sigAndBkgnds.push_back(lowSubLeptDYEE);
	comments.push_back("lower sublead lept pT cut e chnl");

	TString lowLeadLeptDir = "../analysisCppOutputRootFiles/withLeadLeptonPtFiftyCut/";
	TChain * lowLeadLeptWrSigEE = new TChain(treeName);
	lowLeadLeptWrSigEE->Add(lowLeadLeptDir+"selected_tree_WRtoEEJJ_2200_1100_signal_eeEE.root");
	TChain * lowLeadLeptTTEE = new TChain(treeName);
	lowLeadLeptTTEE->Add(lowLeadLeptDir+"selected_tree_TT_signal_eeEE.root");
	TChain * lowLeadLeptDYEE = new TChain(treeName);
	lowLeadLeptDYEE->Add(lowLeadLeptDir+"selected_tree_DYAMC_signal_eeEE.root");
	sigAndBkgnds.push_back(lowLeadLeptWrSigEE), sigAndBkgnds.push_back(lowLeadLeptTTEE), sigAndBkgnds.push_back(lowLeadLeptDYEE);
	comments.push_back("lower lead lept pT cut e chnl");


	//muon channel
	TChain * defWrSigMuMu = new TChain(treeName);
	defWrSigMuMu->Add(defDir+"selected_tree_WRtoMuMuJJ_2200_1100_signal_mumuMuMu.root");
	TChain * defTTMuMu = new TChain(treeName);
	defTTMuMu->Add(defDir+"selected_tree_TT_signal_mumuMuMu.root");
	TChain * defDYMuMu = new TChain(treeName);
	defDYMuMu->Add(defDir+"selected_tree_DYAMC_signal_mumuMuMu.root");
	sigAndBkgnds.push_back(defWrSigMuMu), sigAndBkgnds.push_back(defTTMuMu), sigAndBkgnds.push_back(defDYMuMu);
	comments.push_back("default cuts mu chnl");

	TChain * lowJetWrSigMuMu = new TChain(treeName);
	lowJetWrSigMuMu->Add(lowJetDir+"selected_tree_WRtoMuMuJJ_2200_1100_signal_mumuMuMu.root");
	TChain * lowJetTTMuMu = new TChain(treeName);
	lowJetTTMuMu->Add(lowJetDir+"selected_tree_TT_signal_mumuMuMu.root");
	TChain * lowJetDYMuMu = new TChain(treeName);
	lowJetDYMuMu->Add(lowJetDir+"selected_tree_DYAMC_signal_mumuMuMu.root");
	sigAndBkgnds.push_back(lowJetWrSigMuMu), sigAndBkgnds.push_back(lowJetTTMuMu), sigAndBkgnds.push_back(lowJetDYMuMu);
	comments.push_back("lower jet pT cut mu chnl");

	TChain * lowSubLeptWrSigMuMu = new TChain(treeName);
	lowSubLeptWrSigMuMu->Add(lowSubLeptDir+"selected_tree_WRtoMuMuJJ_2200_1100_signal_mumuMuMu.root");
	TChain * lowSubLeptTTMuMu = new TChain(treeName);
	lowSubLeptTTMuMu->Add(lowSubLeptDir+"selected_tree_TT_signal_mumuMuMu.root");
	TChain * lowSubLeptDYMuMu = new TChain(treeName);
	lowSubLeptDYMuMu->Add(lowSubLeptDir+"selected_tree_DYAMC_signal_mumuMuMu.root");
	sigAndBkgnds.push_back(lowSubLeptWrSigMuMu), sigAndBkgnds.push_back(lowSubLeptTTMuMu), sigAndBkgnds.push_back(lowSubLeptDYMuMu);
	comments.push_back("lower sublead lept pT cut mu chnl");

	TChain * lowLeadLeptWrSigMuMu = new TChain(treeName);
	lowLeadLeptWrSigMuMu->Add(lowLeadLeptDir+"selected_tree_WRtoMuMuJJ_2200_1100_signal_mumuMuMu.root");
	TChain * lowLeadLeptTTMuMu = new TChain(treeName);
	lowLeadLeptTTMuMu->Add(lowLeadLeptDir+"selected_tree_TT_signal_mumuMuMu.root");
	TChain * lowLeadLeptDYMuMu = new TChain(treeName);
	lowLeadLeptDYMuMu->Add(lowLeadLeptDir+"selected_tree_DYAMC_signal_mumuMuMu.root");
	sigAndBkgnds.push_back(lowLeadLeptWrSigMuMu), sigAndBkgnds.push_back(lowLeadLeptTTMuMu), sigAndBkgnds.push_back(lowLeadLeptDYMuMu);
	comments.push_back("lower lead lept pT cut mu chnl");
	/**/

	//loop over elements in vector of TChains, calculate S over B, and write S over B plus a comment about the cuts to a file
	Int_t vSize = sigAndBkgnds.size();
	std::string brName = "NEventsInRange";
	std::string errBrName = "ErrorEventsInRange";
	//only one entry in every branch of this tree
	for(Int_t i=0; i<vSize; i+= (numBkgnds+1) ){
		Float_t sigEvts = 0, bkgndEvts = 0, sOverB = 0;
		Float_t sOverSqrtB = 0;
	
		//load tree contents
		for(Int_t j=i; j<(i+numBkgnds+1) ; j++){
			Float_t tempHolder[64];
			Float_t tempErrHolder[64];
			sigAndBkgnds[j]->SetBranchAddress(brName.c_str(), &tempHolder);
			sigAndBkgnds[j]->SetBranchAddress(errBrName.c_str(), &tempErrHolder);
			sigAndBkgnds[j]->GetEntry(0);
			if(tempHolder[massWindowIndex] >= 0){
				if(j==i) sigEvts += tempHolder[massWindowIndex];
				else bkgndEvts += tempHolder[massWindowIndex];
			}
			if(tempHolder[massWindowIndex] < 0){//the error on the events in range is always positive
				if(j==i) sigEvts += tempErrHolder[massWindowIndex];
				else bkgndEvts += tempErrHolder[massWindowIndex];
			}
		}//load tree contents, only one entry in each tree

		//calculate S over B
		if(bkgndEvts > 0) sOverB = sigEvts/bkgndEvts;
		if(bkgndEvts > 0) sOverSqrtB = sigEvts/TMath::Sqrt(bkgndEvts);
		std::cout<< comments[(Int_t) (i/3)] << std::endl;
		std::cout<<"signal =\t"<< sigEvts << std::endl;
		std::cout<<"bkgnd =\t"<< bkgndEvts << std::endl;
		std::cout<<"S/B =\t"<< sOverB << std::endl;
		std::cout<<"S/sqrt(B) =\t"<< sOverSqrtB << std::endl;
		std::cout<<" "<< std::endl;

	}//end loop over vector of TChains


#endif
	//end sOverBsensitivity


#ifdef showMassWindows
	//make two sets of plots
	//1 - show MLLJJ stacked backgrounds with different colors assigned to different bkgnds (no error bars or points)
	//2 - show each MLLJJ background as data points with a second set of points showing only stat unc, and a third set of pts showing only syst unc
	//    for each bkgnd there should be 5 data points in each bin: 2 showing syst uncert, 2 showing stat unc, and 1 showing nominal value
	//    for N bkgnds, there should be N bins at each mass point, don't stack the bkgnds

	//identify mass windows to show, and declare labels for them
	//changes to the mass window indices should be propagated to the mass window labels
	int massIndices[] = {2,8,15};	//mass windows which will be shown
	TString massWindowLabel[] = {"1.0 TeV","2.2 TeV","3.6 TeV"};	//labels for mass windows
	std::vector<TString> massWindowLabelVect(massWindowLabel,massWindowLabel + sizeof(massWindowLabel)/sizeof(TString) );
	Int_t numMassWindows = massWindowLabelVect.size();

	//declare pointers to input TTrees and emu weights which should be applied to EMu data to estimate TTBar
	TString treeName = "syst_tree";
  	TString localDir = "/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysAllSyst/";
	Double_t EEWgt = 0.414, MuMuWgt = 0.657;	//weights to apply to EMu data to estimate TTBar

   	//declare names of bkgnds which will be shown, and histos used in stacking or overlaying
	TString bkgndNames[] = {"TT Data Driven","DYAMCNLO"};
	Int_t numBkgnds = 2;
	TH1D * EEBkgndOneForStack = new TH1D("EEBkgndOneForStack","",numMassWindows,0,numMassWindows-1);
	TH1D * EEBkgndTwoForStack = new TH1D("EEBkgndTwoForStack","",numMassWindows,0,numMassWindows-1);
	TH1D * MuMuBkgndOneForStack = new TH1D("MuMuBkgndOneForStack","",numMassWindows,0,numMassWindows-1);
	TH1D * MuMuBkgndTwoForStack = new TH1D("MuMuBkgndTwoForStack","",numMassWindows,0,numMassWindows-1);

	//syst uncert of each bkgnd
	TH1D * EEBkgndOneSystUnc = new TH1D("EEBkgndOneSystUnc","",numMassWindows,0,numMassWindows-1);
	TH1D * EEBkgndTwoSystUnc = new TH1D("EEBkgndTwoSystUnc","",numMassWindows,0,numMassWindows-1);
	TH1D * MuMuBkgndOneSystUnc = new TH1D("MuMuBkgndOneSystUnc","",numMassWindows,0,numMassWindows-1);
	TH1D * MuMuBkgndTwoSystUnc = new TH1D("MuMuBkgndTwoSystUnc","",numMassWindows,0,numMassWindows-1);

	//stat uncert of each bkgnd
	TH1D * EEBkgndOneStatUnc = new TH1D("EEBkgndOneStatUnc","",numMassWindows,0,numMassWindows-1);
	TH1D * EEBkgndTwoStatUnc = new TH1D("EEBkgndTwoStatUnc","",numMassWindows,0,numMassWindows-1);
	TH1D * MuMuBkgndOneStatUnc = new TH1D("MuMuBkgndOneStatUnc","",numMassWindows,0,numMassWindows-1);
	TH1D * MuMuBkgndTwoStatUnc = new TH1D("MuMuBkgndTwoStatUnc","",numMassWindows,0,numMassWindows-1);

	//stat and syst and nominal for all bkgnds at every mass window
	TH1D * EEAllBkgndsNoUncs = new TH1D("EEAllBkgndsNoUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	TH1D * MuMuAllBkgndsNoUncs = new TH1D("MuMuAllBkgndsNoUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	
	TH1D * EEAllBkgndsPlusSystUncs = new TH1D("EEAllBkgndsPlusSystUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	TH1D * MuMuAllBkgndsPlusSystUncs = new TH1D("MuMuAllBkgndsPlusSystUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	TH1D * EEAllBkgndsMinusSystUncs = new TH1D("EEAllBkgndsMinusSystUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	TH1D * MuMuAllBkgndsMinusSystUncs = new TH1D("MuMuAllBkgndsMinusSystUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	
	TH1D * EEAllBkgndsPlusStatUncs = new TH1D("EEAllBkgndsPlusStatUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	TH1D * MuMuAllBkgndsPlusStatUncs = new TH1D("MuMuAllBkgndsPlusStatUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	TH1D * EEAllBkgndsMinusStatUncs = new TH1D("EEAllBkgndsMinusStatUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	TH1D * MuMuAllBkgndsMinusStatUncs = new TH1D("MuMuAllBkgndsMinusStatUncs","",numBkgnds*numMassWindows,0,(numBkgnds*numMassWindows)-1);
	


	//electron channel
	TChain * TTEE = new TChain(treeName);
	TTEE->Add(localDir+"selected_tree_data_flavoursidebandEMu.root");
	//TTEE->SetWeight(EEWgt,"global");	//trick to get EMu data evts to be counted as EE TTBar evts
	TChain * DYEE = new TChain(treeName);
	DYEE->Add(localDir+"selected_tree_DYAMC_signal_eeEE_withMllWeight.root");
	//TChain * WEE = new TChain(treeName);
	//WEE->Add(localDir+"selected_tree_W_signal_eeEE.root");
	//TChain * WZEE = new TChain(treeName);
	//WZEE->Add(localDir+"selected_tree_WZ_signal_eeEE.root");
	//TChain * ZZEE = new TChain(treeName);
	//ZZEE->Add(localDir+"selected_tree_ZZ_signal_eeEE.root");
	
	//muon channel
	TChain * TTMuMu = new TChain(treeName);
	TTMuMu->Add(localDir+"selected_tree_data_flavoursidebandEMuCopy.root");
	//TTMuMu->SetWeight(MuMuWgt,"global");	//trick to get EMu data evts to be counted as MuMu TTBar evts
	TChain * DYMuMu = new TChain(treeName);
	DYMuMu->Add(localDir+"selected_tree_DYAMC_signal_mumuMuMu_withMllWeight.root");
	//TChain * WMuMu = new TChain(treeName);
	//WMuMu->Add(localDir+"selected_tree_W_signal_mumuMuMu.root");
	//TChain * WZMuMu = new TChain(treeName);
	//WZMuMu->Add(localDir+"selected_tree_WZ_signal_mumuMuMu.root");
	//TChain * ZZMuMu = new TChain(treeName);
	//ZZMuMu->Add(localDir+"selected_tree_ZZ_signal_mumuMuMu.root");

	
	//fill histos corresponding to individual bkgnds
	//loop over different mass windows
	for(Int_t i=0; i<numMassWindows; i++){
		//TTBar ele chnl
		EEBkgndOneForStack->SetBinContent(i+1, EEWgt*(calculateBranchMean(TTEE, "NEventsInRange["+to_string(massIndices[i])+"]", 1)) );
		//std::cout<<"in the mass window "<< massWindowLabel[i] <<" there are this many TT EE evts\t" << EEBkgndOneForStack->GetBinContent(i+1) <<std::endl;
		EEBkgndOneForStack->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		EEBkgndOneForStack->SetFillColor(kGreen);
		EEBkgndOneSystUnc->SetBinContent(i+1, EEWgt*(calculateBranchStdDev(TTEE, "NEventsInRange["+to_string(massIndices[i])+"]", 1)) );
		EEBkgndOneSystUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		EEBkgndOneStatUnc->SetBinContent(i+1, EEWgt*(calculateBranchMean(TTEE, "ErrorEventsInRange["+to_string(massIndices[i])+"]", 1)) );
		EEBkgndOneStatUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
	
		//DY ele chnl
		EEBkgndTwoForStack->SetBinContent(i+1,calculateBranchMean(DYEE, "NEventsInRange["+to_string(massIndices[i])+"]", 1));
		EEBkgndTwoForStack->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		EEBkgndTwoForStack->SetFillColor(kYellow);
		EEBkgndTwoSystUnc->SetBinContent(i+1,calculateBranchStdDev(DYEE, "NEventsInRange["+to_string(massIndices[i])+"]", 1));
		EEBkgndTwoSystUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		EEBkgndTwoStatUnc->SetBinContent(i+1,calculateBranchMean(DYEE, "ErrorEventsInRange["+to_string(massIndices[i])+"]", 1));
		EEBkgndTwoStatUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
	
		//TTBar mu chnl
		MuMuBkgndOneForStack->SetBinContent(i+1, MuMuWgt*(calculateBranchMean(TTMuMu, "NEventsInRange["+to_string(massIndices[i])+"]", 1)) );
		//std::cout<<"in the mass window "<< massWindowLabel[i] <<" there are this many TT MuMu evts\t" << MuMuBkgndOneForStack->GetBinContent(i+1) <<std::endl;
		MuMuBkgndOneForStack->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		MuMuBkgndOneForStack->SetFillColor(kGreen);
		MuMuBkgndOneSystUnc->SetBinContent(i+1, MuMuWgt*(calculateBranchStdDev(TTMuMu, "NEventsInRange["+to_string(massIndices[i])+"]", 1)) );
		MuMuBkgndOneSystUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		MuMuBkgndOneStatUnc->SetBinContent(i+1, MuMuWgt*(calculateBranchMean(TTMuMu, "ErrorEventsInRange["+to_string(massIndices[i])+"]", 1)) );
		MuMuBkgndOneStatUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
	
		//DY mu chnl
		MuMuBkgndTwoForStack->SetBinContent(i+1,calculateBranchMean(DYMuMu, "NEventsInRange["+to_string(massIndices[i])+"]", 1));
		MuMuBkgndTwoForStack->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		MuMuBkgndTwoForStack->SetFillColor(kYellow);
		MuMuBkgndTwoSystUnc->SetBinContent(i+1,calculateBranchStdDev(DYMuMu, "NEventsInRange["+to_string(massIndices[i])+"]", 1));
		MuMuBkgndTwoSystUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);
		MuMuBkgndTwoStatUnc->SetBinContent(i+1,calculateBranchMean(DYMuMu, "ErrorEventsInRange["+to_string(massIndices[i])+"]", 1));
		MuMuBkgndTwoStatUnc->GetXaxis()->SetBinLabel(i+1, massWindowLabel[i]);

		//combined DY and TT with stat and syst uncs
		//from left to right: TT first, then DY for each mass window
		//overlay all of these histos, draw them as points without errors
		Int_t binNum = i*numBkgnds;
		EEAllBkgndsNoUncs->SetBinContent(binNum+1, EEBkgndOneForStack->GetBinContent(i+1) );//TT
		EEAllBkgndsNoUncs->SetBinContent(binNum+2, EEBkgndTwoForStack->GetBinContent(i+1) );//DY
		EEAllBkgndsNoUncs->GetXaxis()->SetBinLabel(binNum+1,"TT " + massWindowLabel[i]);
		EEAllBkgndsNoUncs->GetXaxis()->SetBinLabel(binNum+2,"DY " + massWindowLabel[i]);
	
		EEAllBkgndsPlusSystUncs->SetBinContent(i*numBkgnds+1, EEBkgndOneForStack->GetBinContent(i+1) + EEBkgndOneSystUnc->GetBinContent(i+1) );//TT+syst
		EEAllBkgndsPlusSystUncs->SetBinContent(i*numBkgnds+2, EEBkgndTwoForStack->GetBinContent(i+1) + EEBkgndTwoSystUnc->GetBinContent(i+1) );//DY+syst
		EEAllBkgndsPlusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		EEAllBkgndsPlusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);
	
		EEAllBkgndsMinusSystUncs->SetBinContent(i*numBkgnds+1, EEBkgndOneForStack->GetBinContent(i+1) - EEBkgndOneSystUnc->GetBinContent(i+1) );//TT-syst
		EEAllBkgndsMinusSystUncs->SetBinContent(i*numBkgnds+2, EEBkgndTwoForStack->GetBinContent(i+1) - EEBkgndTwoSystUnc->GetBinContent(i+1) );//DY-syst
		EEAllBkgndsMinusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		EEAllBkgndsMinusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);
	
	
		EEAllBkgndsPlusStatUncs->SetBinContent(i*numBkgnds+1, EEBkgndOneForStack->GetBinContent(i+1) + EEBkgndOneStatUnc->GetBinContent(i+1) );//TT+stat
		EEAllBkgndsPlusStatUncs->SetBinContent(i*numBkgnds+2, EEBkgndTwoForStack->GetBinContent(i+1) + EEBkgndTwoStatUnc->GetBinContent(i+1) );//DY+stat
		EEAllBkgndsPlusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		EEAllBkgndsPlusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);
		
		
		EEAllBkgndsMinusStatUncs->SetBinContent(i*numBkgnds+1, EEBkgndOneForStack->GetBinContent(i+1) - EEBkgndOneStatUnc->GetBinContent(i+1) );//TT-stat
		EEAllBkgndsMinusStatUncs->SetBinContent(i*numBkgnds+2, EEBkgndTwoForStack->GetBinContent(i+1) - EEBkgndTwoStatUnc->GetBinContent(i+1) );//DY-stat
		EEAllBkgndsMinusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		EEAllBkgndsMinusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);


		MuMuAllBkgndsNoUncs->SetBinContent(binNum+1, MuMuBkgndOneForStack->GetBinContent(i+1) );//TT
		MuMuAllBkgndsNoUncs->SetBinContent(binNum+2, MuMuBkgndTwoForStack->GetBinContent(i+1) );//DY
		MuMuAllBkgndsNoUncs->GetXaxis()->SetBinLabel(binNum+1,"TT " + massWindowLabel[i]);
		MuMuAllBkgndsNoUncs->GetXaxis()->SetBinLabel(binNum+2,"DY " + massWindowLabel[i]);
	
		MuMuAllBkgndsPlusSystUncs->SetBinContent(i*numBkgnds+1, MuMuBkgndOneForStack->GetBinContent(i+1) + MuMuBkgndOneSystUnc->GetBinContent(i+1) );//TT+syst
		MuMuAllBkgndsPlusSystUncs->SetBinContent(i*numBkgnds+2, MuMuBkgndTwoForStack->GetBinContent(i+1) + MuMuBkgndTwoSystUnc->GetBinContent(i+1) );//DY+syst
		MuMuAllBkgndsPlusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		MuMuAllBkgndsPlusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);
	
		MuMuAllBkgndsMinusSystUncs->SetBinContent(i*numBkgnds+1, MuMuBkgndOneForStack->GetBinContent(i+1) - MuMuBkgndOneSystUnc->GetBinContent(i+1) );//TT-syst
		MuMuAllBkgndsMinusSystUncs->SetBinContent(i*numBkgnds+2, MuMuBkgndTwoForStack->GetBinContent(i+1) - MuMuBkgndTwoSystUnc->GetBinContent(i+1) );//DY-syst
		MuMuAllBkgndsMinusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		MuMuAllBkgndsMinusSystUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);
	
	
		MuMuAllBkgndsPlusStatUncs->SetBinContent(i*numBkgnds+1, MuMuBkgndOneForStack->GetBinContent(i+1) + MuMuBkgndOneStatUnc->GetBinContent(i+1) );//TT+stat
		MuMuAllBkgndsPlusStatUncs->SetBinContent(i*numBkgnds+2, MuMuBkgndTwoForStack->GetBinContent(i+1) + MuMuBkgndTwoStatUnc->GetBinContent(i+1) );//DY+stat
		MuMuAllBkgndsPlusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		MuMuAllBkgndsPlusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);
		
		
		MuMuAllBkgndsMinusStatUncs->SetBinContent(i*numBkgnds+1, MuMuBkgndOneForStack->GetBinContent(i+1) - MuMuBkgndOneStatUnc->GetBinContent(i+1) );//TT-stat
		MuMuAllBkgndsMinusStatUncs->SetBinContent(i*numBkgnds+2, MuMuBkgndTwoForStack->GetBinContent(i+1) - MuMuBkgndTwoStatUnc->GetBinContent(i+1) );//DY-stat
		MuMuAllBkgndsMinusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+1, "TT " + massWindowLabel[i]);
		MuMuAllBkgndsMinusStatUncs->GetXaxis()->SetBinLabel(i*numBkgnds+2,"DY " + massWindowLabel[i]);
	
	}//end loop over different mass windows
	

	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	//set 1  stacked backgrounds
	TLegend *legEE = new TLegend(0.6, 0.6, 0.9, 0.9);
	legEE->AddEntry(EEBkgndOneForStack, "TTBar Data Driven");
	legEE->AddEntry(EEBkgndTwoForStack, "DY AMCNLO");
	TCanvas* canvasEE = new TCanvas("canvasEE","canvasEE",0,0,600,600);
	canvasEE->cd();
	THStack *EEStack = new THStack();	//add smallest bkgnds first, and largest bkgnds last
	EEStack->SetMinimum(-10);
	EEStack->Add(EEBkgndTwoForStack);
	EEStack->Add(EEBkgndOneForStack);
	EEStack->Draw("histo");
	EEStack->GetYaxis()->SetTitle("Ele Events per W_{R} mass window");
	EEStack->GetYaxis()->SetTitleOffset(1.35);
	EEStack->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
	legEE->Draw();
	canvasEE->cd();
	TString outFileNameStackedBkgndsEE = "EEChnlStackedBkgnds";
	canvasEE->Print(outFileNameStackedBkgndsEE + ".pdf");
	canvasEE->Print(outFileNameStackedBkgndsEE + ".png");
	EEStack->SetMinimum(0.1);
	canvasEE->SetLogy();
	canvasEE->Print(outFileNameStackedBkgndsEE + "_log.pdf");
	canvasEE->Print(outFileNameStackedBkgndsEE + "_log.png");
	canvasEE->Close();

	TLegend *legMuMu = new TLegend(0.6, 0.6, 0.9, 0.9);
	legMuMu->AddEntry(MuMuBkgndOneForStack, "TTBar Data Driven");
	legMuMu->AddEntry(MuMuBkgndTwoForStack, "DY AMCNLO");
	TCanvas* canvasMuMu = new TCanvas("canvasMuMu","canvasMuMu",0,0,600,600);
	canvasMuMu->cd();
	THStack *MuMuStack = new THStack();	//add smallest bkgnds first, and largest bkgnds last
	MuMuStack->SetMinimum(-10);
	MuMuStack->Add(MuMuBkgndTwoForStack);
	MuMuStack->Add(MuMuBkgndOneForStack);
	MuMuStack->Draw("histo");
	MuMuStack->GetYaxis()->SetTitle("Muon Events per W_{R} mass window");
	MuMuStack->GetYaxis()->SetTitleOffset(1.35);
	MuMuStack->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
	legMuMu->Draw();
	canvasMuMu->cd();
	TString outFileNameStackedBkgndsMuMu = "MuMuChnlStackedBkgnds";
	canvasMuMu->Print(outFileNameStackedBkgndsMuMu + ".pdf");
	canvasMuMu->Print(outFileNameStackedBkgndsMuMu + ".png");
	MuMuStack->SetMinimum(0.1);
	canvasMuMu->SetLogy();
	canvasMuMu->Print(outFileNameStackedBkgndsMuMu + "_log.pdf");
	canvasMuMu->Print(outFileNameStackedBkgndsMuMu + "_log.png");
	canvasMuMu->Close();


	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	//set 2 individual bkgnds shown as points with a second set of pts representing stat unc, and a third set of pts representing syst unc
	//each mass bin is shown N times, where N = number of different bkgnds. 5 data points per bin
	gStyle->SetOptStat("");
	TCanvas* canvasTotEE = new TCanvas("canvasTotEE","canvasTotEE",0,0,600,600);
	canvasTotEE->cd();
	EEAllBkgndsNoUncs->SetLineColor(kBlack);
	EEAllBkgndsNoUncs->SetLineWidth(2);


	EEAllBkgndsPlusSystUncs->SetMarkerColor(kRed);
	EEAllBkgndsPlusSystUncs->SetMarkerStyle(20);
	EEAllBkgndsPlusSystUncs->SetMarkerSize(1);
	
	EEAllBkgndsMinusSystUncs->SetMarkerColor(kRed);
	EEAllBkgndsMinusSystUncs->SetMarkerStyle(20);
	EEAllBkgndsMinusSystUncs->SetMarkerSize(1);
	
	EEAllBkgndsPlusStatUncs->SetMarkerColor(kBlue);
	EEAllBkgndsPlusStatUncs->SetMarkerStyle(20);
	EEAllBkgndsPlusStatUncs->SetMarkerSize(1);
	
	EEAllBkgndsMinusStatUncs->SetMarkerColor(kBlue);
	EEAllBkgndsMinusStatUncs->SetMarkerStyle(20);
	EEAllBkgndsMinusStatUncs->SetMarkerSize(1);
	
	TLegend *legTotEE = new TLegend(0.6, 0.6, 0.9, 0.9);
	legTotEE->AddEntry(EEAllBkgndsNoUncs, "Nominal");
	legTotEE->AddEntry(EEAllBkgndsPlusSystUncs, "Nominal+Syst Unc","p");
	legTotEE->AddEntry(EEAllBkgndsMinusSystUncs, "Nominal-Syst Unc","p");
	legTotEE->AddEntry(EEAllBkgndsPlusStatUncs, "Nominal+Stat Unc","p");
	legTotEE->AddEntry(EEAllBkgndsMinusStatUncs, "Nominal-Stat Unc","p");


	EEAllBkgndsNoUncs->Draw("histo");
	EEAllBkgndsPlusSystUncs->Draw("Psame");
	EEAllBkgndsMinusSystUncs->Draw("Psame");
	EEAllBkgndsPlusStatUncs->Draw("Psame");
	EEAllBkgndsMinusStatUncs->Draw("Psame");
	EEAllBkgndsNoUncs->GetYaxis()->SetTitle("Ele Events per W_{R} mass window");
	EEAllBkgndsNoUncs->GetYaxis()->SetTitleOffset(1.35);
	EEAllBkgndsNoUncs->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
	canvasTotEE->Update();
	legTotEE->Draw();
	canvasTotEE->cd();
	TString outFileNameIndivBkgndsEE = "EEChnlIndivBkgndsWithUncs";
	canvasTotEE->Print(outFileNameIndivBkgndsEE + ".pdf");
	canvasTotEE->Print(outFileNameIndivBkgndsEE + ".png");
	canvasTotEE->Print(outFileNameIndivBkgndsEE + ".C");
	EEAllBkgndsNoUncs->SetMinimum(0.1);
	canvasTotEE->SetLogy();
	canvasTotEE->Print(outFileNameIndivBkgndsEE + "_log.pdf");
	canvasTotEE->Print(outFileNameIndivBkgndsEE + "_log.png");
	canvasTotEE->Close();


	TCanvas* canvasTotMuMu = new TCanvas("canvasTotMuMu","canvasTotMuMu",0,0,600,600);
	canvasTotMuMu->cd();
	MuMuAllBkgndsNoUncs->SetLineColor(kBlack);
	MuMuAllBkgndsNoUncs->SetLineWidth(2);


	MuMuAllBkgndsPlusSystUncs->SetMarkerColor(kRed);
	MuMuAllBkgndsPlusSystUncs->SetMarkerStyle(20);
	MuMuAllBkgndsPlusSystUncs->SetMarkerSize(1);
	
	MuMuAllBkgndsMinusSystUncs->SetMarkerColor(kRed);
	MuMuAllBkgndsMinusSystUncs->SetMarkerStyle(20);
	MuMuAllBkgndsMinusSystUncs->SetMarkerSize(1);
	
	MuMuAllBkgndsPlusStatUncs->SetMarkerColor(kBlue);
	MuMuAllBkgndsPlusStatUncs->SetMarkerStyle(20);
	MuMuAllBkgndsPlusStatUncs->SetMarkerSize(1);
	
	MuMuAllBkgndsMinusStatUncs->SetMarkerColor(kBlue);
	MuMuAllBkgndsMinusStatUncs->SetMarkerStyle(20);
	MuMuAllBkgndsMinusStatUncs->SetMarkerSize(1);
	
	TLegend *legTotMuMu = new TLegend(0.6, 0.6, 0.9, 0.9);
	legTotMuMu->AddEntry(MuMuAllBkgndsNoUncs, "Nominal");
	legTotMuMu->AddEntry(MuMuAllBkgndsPlusSystUncs, "Nominal+Syst Unc","p");
	legTotMuMu->AddEntry(MuMuAllBkgndsMinusSystUncs, "Nominal-Syst Unc","p");
	legTotMuMu->AddEntry(MuMuAllBkgndsPlusStatUncs, "Nominal+Stat Unc","p");
	legTotMuMu->AddEntry(MuMuAllBkgndsMinusStatUncs, "Nominal-Stat Unc","p");


	MuMuAllBkgndsNoUncs->Draw("histo");
	MuMuAllBkgndsPlusSystUncs->Draw("Psame");
	MuMuAllBkgndsMinusSystUncs->Draw("Psame");
	MuMuAllBkgndsPlusStatUncs->Draw("Psame");
	MuMuAllBkgndsMinusStatUncs->Draw("Psame");
	MuMuAllBkgndsNoUncs->GetYaxis()->SetTitle("Mu Events per W_{R} mass window");
	MuMuAllBkgndsNoUncs->GetYaxis()->SetTitleOffset(1.35);
	MuMuAllBkgndsNoUncs->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
	canvasTotMuMu->Update();
	legTotMuMu->Draw();
	canvasTotMuMu->cd();
	TString outFileNameIndivBkgndsMuMu = "MuMuChnlIndivBkgndsWithUncs";
	canvasTotMuMu->Print(outFileNameIndivBkgndsMuMu + ".pdf");
	canvasTotMuMu->Print(outFileNameIndivBkgndsMuMu + ".png");
	canvasTotMuMu->Print(outFileNameIndivBkgndsMuMu + ".C");
	MuMuAllBkgndsNoUncs->SetMinimum(0.1);
	canvasTotMuMu->SetLogy();
	canvasTotMuMu->Print(outFileNameIndivBkgndsMuMu + "_log.pdf");
	canvasTotMuMu->Print(outFileNameIndivBkgndsMuMu + "_log.png");
	canvasTotMuMu->Close();


#endif
	//end showMassWindows

#ifdef printNewDySyst
	///open the DYMADHT and DYAMC signal region files and print evt count difference btwn them in one
	///mass window divided by 2 times the DYAMC evt count
	string pathToMassWindowFile = "../configs/mass_cuts.txt";

	int wrMassPoints[] = {1000,1600,2200,2800,3600};
	vector<int> wrMassVect(wrMassPoints, wrMassPoints + sizeof(wrMassPoints)/sizeof(int) );
	
	string pathToDySystFile = "dyAmcMinusMadHtSyst.txt";
	ofstream writeToDySystFile(pathToDySystFile.c_str(),ofstream::trunc);

	string treeName = "central_value_tree";
	TChain * dyMuMuAMC = new TChain(treeName.c_str(),"DYMuMuAMC");
	dyMuMuAMC->Add("../analysisCppOutputRootFiles/selected_tree_DYAMC_signal_mumuMuMu_withMllWeight.root" );
	TChain * dyMuMuMADHT = new TChain(treeName.c_str(),"DYMuMuMADHT");
	dyMuMuMADHT->Add("../analysisCppOutputRootFiles/selected_tree_DYMADHT_signal_mumuMuMu_withMllWeight.root" );

	TChain * dyEEAMC = new TChain(treeName.c_str(),"DYEEAMC");
	dyEEAMC->Add("../analysisCppOutputRootFiles/selected_tree_DYAMC_signal_eeEE_withMllWeight.root" );
	TChain * dyEEMADHT = new TChain(treeName.c_str(),"DYEEMADHT");
	dyEEMADHT->Add("../analysisCppOutputRootFiles/selected_tree_DYMADHT_signal_eeEE_withMllWeight.root" );

	map<string,TChain*> dyChainMap;
	dyChainMap["dyMuMuAMC"] = dyMuMuAMC;
	dyChainMap["dyMuMuMADHT"] = dyMuMuMADHT;
	dyChainMap["dyEEAMC"] = dyEEAMC;
	dyChainMap["dyEEMADHT"] = dyEEMADHT;

	int numWrMasses = wrMassVect.size();

	for(int m=0; m<numWrMasses; m++){
		//loop over each WR mass window
		map<string,Double_t> means;

		for(map<string,TChain*>::const_iterator chMapIt=dyChainMap.begin(); chMapIt!=dyChainMap.end(); chMapIt++){
			//loop over each TChain in dyChainMap

			string nEventsIndex = to_string( getIndexForSystematicUncertainty( to_string(wrMassVect[m]), chMapIt->first, pathToMassWindowFile) );
			string drawArg = "NEventsInRange["+nEventsIndex+"]>>";
			string histName = (chMapIt->first) + "tempHist";
			(chMapIt->second)->Draw( (drawArg + histName + "()").c_str() );
			TH1F* tempHist = (TH1F*) gROOT->FindObject(histName.c_str());
			if(tempHist->GetMean() > 0.) means[chMapIt->first] = tempHist->GetMean();
			else{
				//num weighted entries in dyAMC is negative or zero
				string statUncDrawArg = "ErrorEventsInRange["+nEventsIndex+"]>>";
				string statUncHistName = (chMapIt->first) + "StatTempHist";
				(chMapIt->second)->Draw( (statUncDrawArg + statUncHistName + "()").c_str() );
				TH1F* statUncTempHist = (TH1F*) gROOT->FindObject(statUncHistName.c_str());
				means[chMapIt->first] = statUncTempHist->GetMean();
			}
		}//end loop over different TChains
		writeToDySystFile << "EE\t" << wrMassVect[m] << "\t" << (fabs(means["dyEEAMC"] - means["dyEEMADHT"])/(2*means["dyEEAMC"])) << endl;
		writeToDySystFile << "MuMu\t" << wrMassVect[m] << "\t" << (fabs(means["dyMuMuAMC"] - means["dyMuMuMADHT"])/(2*means["dyMuMuAMC"])) << endl;

	}//end loop over different WR mass points

	writeToDySystFile.close();

#endif
	//end printNewDySyst

	// sumPt
#ifdef DYHTPlot
	//make a plot of the sum pT of all GEN quarks and gluons leaving the hard pp->Z->ll interaction
	//the plot will show one curve for the DYMad inclusive dataset, and another curve for the DYMadHT100-200 dataset
	
	/**/
	//chains for GEN hadron sumPt
	TString treeName = "genMatchedAnalyzerOne/genGluonsAndQuarksNoCuts";
	TChain * dyMadIncl = new TChain(treeName,"dyMadIncl");
	dyMadIncl->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadIncl.root");
	TChain * dyMadHT100to200 = new TChain(treeName,"dyMadHT100to200");
	dyMadHT100to200->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT100to200.root");
	TChain * dyMadHT200to400 = new TChain(treeName,"dyMadHT200to400");
	dyMadHT200to400->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT200to400.root");
	TChain * dyMadHT400to600 = new TChain(treeName,"dyMadHT400to600");
	dyMadHT400to600->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT400to600.root");
	TChain * dyMadHT600toInf = new TChain(treeName,"dyMadHT600toInf");
	dyMadHT600toInf->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT600toInf.root");


	//chains for GEN lepton Z pT
	TString treeNameLepton = "genMatchedAnalyzerTwo/genElectronsNoCuts";
	TChain * dyMadLeptIncl = new TChain(treeNameLepton,"dyMadLeptIncl");
	dyMadLeptIncl->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadIncl.root");
	TChain * dyMadLeptHT100to200 = new TChain(treeNameLepton,"dyMadLeptHT100to200");
	dyMadLeptHT100to200->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT100to200.root");
	TChain * dyMadLeptHT200to400 = new TChain(treeNameLepton,"dyMadLeptHT200to400");
	dyMadLeptHT200to400->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT200to400.root");
	TChain * dyMadLeptHT400to600 = new TChain(treeNameLepton,"dyMadLeptHT400to600");
	dyMadLeptHT400to600->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT400to600.root");
	TChain * dyMadLeptHT600toInf = new TChain(treeNameLepton,"dyMadLeptHT600toInf");
	dyMadLeptHT600toInf->Add("../filesForDYMadStudies/hadronAndEleKinematicsGenDYJetsMadHT600toInf.root");


	//add hadron sumPt chains as friends to lepton Z pT chains
	dyMadLeptIncl->AddFriend(dyMadIncl, "MadInclHadronChain");
	dyMadLeptHT100to200->AddFriend(dyMadHT100to200, "MadHTBinnedHadronChain");
	dyMadLeptHT200to400->AddFriend(dyMadHT200to400, "MadHTBinnedHadronChain");
	dyMadLeptHT400to600->AddFriend(dyMadHT400to600, "MadHTBinnedHadronChain");
	dyMadLeptHT600toInf->AddFriend(dyMadHT600toInf, "MadHTBinnedHadronChain");

	gStyle->SetOptStat("");
	
	//draw 1D plot of GEN HT for inclusive DYMadIncl and DYMadHT100to200
	dyMadIncl->Draw("sumPt>>genSumPtMadInclHist(42,0.,210.)");
	TH1F * sumPtMadInclHist = (TH1F*) gROOT->FindObject("genSumPtMadInclHist");
	dyMadHT100to200->Draw("sumPt>>genSumPtMadHT100to200Hist(42,0.,210.)");
	TH1F * sumPtMadHT100to200Hist = (TH1F*) gROOT->FindObject("genSumPtMadHT100to200Hist");
	TLegend * sumPtLeg = new TLegend( 0.6, 0.60, 0.90, 0.90);
	sumPtLeg->AddEntry(sumPtMadInclHist, "DYMadIncl HT<100 filter");
	sumPtLeg->AddEntry(sumPtMadHT100to200Hist, "DYMadHT100to200");
	TCanvas * canvMadInclAndHT100to200 = new TCanvas("MadInclAndHT100to200","MadInclAndHT100to200",800,800);
	canvMadInclAndHT100to200->cd();
	sumPtMadInclHist->SetLineWidth(2);
	sumPtMadInclHist->SetLineColor(kBlack);
	sumPtMadHT100to200Hist->SetLineWidth(2);
	sumPtMadHT100to200Hist->SetLineColor(kRed);
	sumPtMadInclHist->GetXaxis()->SetTitle("GEN H_{T} (GeV)");
	sumPtMadInclHist->GetYaxis()->SetTitle("Events");
	sumPtMadInclHist->SetTitle("GEN H_{T}");
	sumPtMadInclHist->Draw("histo");
	sumPtMadHT100to200Hist->Draw("histo same");
	sumPtLeg->Draw();
	canvMadInclAndHT100to200->Update();
	TString gensumPtFileName = "genHT_DYMadInclAndMadHT100to200";
	canvMadInclAndHT100to200->Print(gensumPtFileName+".png");
	canvMadInclAndHT100to200->Print(gensumPtFileName+".pdf");
	canvMadInclAndHT100to200->SetLogy();
	canvMadInclAndHT100to200->Print(gensumPtFileName+"_log.png");
	canvMadInclAndHT100to200->Print(gensumPtFileName+"_log.pdf");
	canvMadInclAndHT100to200->Close();

	//draw 2D plot of GEN HT and GEN Z pT   yAxis:xAxis
	dyMadLeptIncl->Draw("sumPt:MadInclHadronChain.zPt>>twoDimMadInclHist()");
	TH2F * twoDimMadInclHist = (TH2F*) gROOT->FindObject("twoDimMadInclHist");
	gStyle->SetOptStat("");
	TCanvas * canvTwoDimMadIncl = new TCanvas("twoDimMadIncl","twoDimMadIncl",800,800);
	canvTwoDimMadIncl->cd();
	twoDimMadInclHist->GetXaxis()->SetTitle("GEN H_{T} (GeV)");
	twoDimMadInclHist->GetYaxis()->SetTitle("GEN Z P_{T} (GeV)");
	twoDimMadInclHist->GetYaxis()->SetTitleOffset(1.4);
	twoDimMadInclHist->SetTitle("Mad Inclusive GEN Z P_{T} vs H_{T}");
	twoDimMadInclHist->Draw("colz");
	TString twoDimMadInclOutFileName = "twoDimMadInclGenHTvsGenZpT";
	canvTwoDimMadIncl->Print(twoDimMadInclOutFileName + ".pdf");
	canvTwoDimMadIncl->Print(twoDimMadInclOutFileName + ".png");
	canvTwoDimMadIncl->Close();

	dyMadLeptHT600toInf->Draw("sumPt:MadHTBinnedHadronChain.zPt>>twoDimMadHTBinnedHist()");
	TH2F * twoDimMadHTBinnedHist = (TH2F*) gROOT->FindObject("twoDimMadHTBinnedHist");
	gStyle->SetOptStat("");
	TCanvas * canvTwoDimMadHTBinned = new TCanvas("twoDimMadHTBinned","twoDimMadHTBinned",800,800);
	canvTwoDimMadHTBinned->cd();
	twoDimMadHTBinnedHist->GetXaxis()->SetTitle("GEN H_{T} (GeV)");
	twoDimMadHTBinnedHist->GetYaxis()->SetTitle("GEN Z P_{T} (GeV)");
	twoDimMadHTBinnedHist->GetYaxis()->SetTitleOffset(1.4);
	twoDimMadHTBinnedHist->SetTitle("MadHT600toInf GEN Z P_{T} vs H_{T}");
	twoDimMadHTBinnedHist->Draw("colz");
	TString twoDimMadHTBinnedOutFileName = "twoDimMadHT600toInfBinnedGenHTvsGenZpT";
	canvTwoDimMadHTBinned->Print(twoDimMadHTBinnedOutFileName + ".pdf");
	canvTwoDimMadHTBinned->Print(twoDimMadHTBinnedOutFileName + ".png");
	canvTwoDimMadHTBinned->Close();

	/**/

	/*
	//chains for RECO hadron sumPt
	TString treeName = "recoAnalyzerTwo/recoJetsAndLeptonsNoCuts";
	TChain * dyMadRecoIncl = new TChain(treeName,"dyMadRecoIncl");
	dyMadRecoIncl->Add("../filesForDYMadStudies/recoHadronAndLeptonKinematicsGenDYJetsMadInclAllSkims.root");
	TChain * dyMadRecoHT100to200 = new TChain(treeName,"dyMadRecoHT100to200");
	dyMadRecoHT100to200->Add("../filesForDYMadStudies/recoHadronAndLeptonKinematicsGenDYJetsMadHT100to200AllSkims.root");
	TChain * dyMadRecoHT200to400 = new TChain(treeName,"dyMadRecoHT200to400");
	dyMadRecoHT200to400->Add("../filesForDYMadStudies/recoHadronAndLeptonKinematicsGenDYJetsMadHT200to400AllSkims.root");
	TChain * dyMadRecoHT400to600 = new TChain(treeName,"dyMadRecoHT400to600");
	dyMadRecoHT400to600->Add("../filesForDYMadStudies/recoHadronAndLeptonKinematicsGenDYJetsMadHT400to600AllSkims.root");
	TChain * dyMadRecoHT600toInf = new TChain(treeName,"dyMadRecoHT600toInf");
	dyMadRecoHT600toInf->Add("../filesForDYMadStudies/recoHadronAndLeptonKinematicsGenDYJetsMadHT600toInfAllSkims.root");

	//draw the sumPt of the two leading RECO jets passing tight ID
	//crossSxn norm for different DYMad samples:
	Double_t lumi = 2640.523267;
	//lumi: 2640.523267
	//Incl: (5991*lumi)/9042031
	//ht100to200: (181.302*lumi)/2725655
	//ht200to400: (50.4177*lumi)/973937
	//ht400to600: (6.98394*lumi)/1067758
	//ht600toInf: (2.70354*lumi)/998912

	//single plots
	//dyMadRecoIncl->Draw("sumPt>>recoJetMadInclHist()");
	//TH1F * recoJetMadInclHist = (TH1F*) gROOT->FindObject("recoJetMadInclHist");
	//gStyle->SetOptStat("");
	//TCanvas * canvRecoJetMadIncl = new TCanvas("recoJetMadIncl","recoJetMadIncl",800,800);
	//canvRecoJetMadIncl->cd();
	//recoJetMadInclHist->GetXaxis()->SetTitle("RECO #Sigma P_{T} Two Lead Jets (GeV)");
	//recoJetMadInclHist->GetYaxis()->SetTitle("Events");
	//recoJetMadInclHist->GetYaxis()->SetTitleOffset(1.3);
	//recoJetMadInclHist->SetTitle("MadIncl Reco H_{T}");
	//recoJetMadInclHist->SetLineWidth(2);
	//recoJetMadInclHist->Scale((5991*lumi)/9042031);
	//recoJetMadInclHist->Draw("histo");
	//TString recoJetMadInclOutFileName = "MadInclRecoHTAllSkims";
	//canvRecoJetMadIncl->Print(recoJetMadInclOutFileName + ".pdf");
	//canvRecoJetMadIncl->Print(recoJetMadInclOutFileName + ".png");
	//canvRecoJetMadIncl->Close();

#ifdef DEBUG
	std::cout<<"about to call TTree Draw to make histos from all DYMadIncl and HTBinned evts"<<std::endl;
#endif

	//one stacked plot showing DYMadIncl (with GEN HT < 100 filter) + all DYMadHT binned samples
	gStyle->SetOptStat("");
	//TString xAxisLabel = "RECO #Sigma P_{T} Two Lead Jets (GeV)";
	TString xAxisLabel = "RECO M_{EEJJ} (GeV)";
	TString yAxisLabel = "Events";
	dyMadRecoIncl->Draw("dileptonMass>>recoJetMadInclHist(75,0.,1200.)");
	TH1F * recoJetMadInclHist = (TH1F*) gROOT->FindObject("recoJetMadInclHist");
	recoJetMadInclHist->Scale((5991*lumi)/9042031);
	recoJetMadInclHist->SetFillColor(kYellow);
	recoJetMadInclHist->GetXaxis()->SetTitle(xAxisLabel);
	recoJetMadInclHist->GetYaxis()->SetTitle(yAxisLabel);
	dyMadRecoHT100to200->Draw("dileptonMass>>recoJetMadHT100to200Hist(75,0.,1200.)");
	TH1F * recoJetMadHT100to200Hist = (TH1F*) gROOT->FindObject("recoJetMadHT100to200Hist");
	recoJetMadHT100to200Hist->Scale((181.302*lumi)/2725655);
	recoJetMadHT100to200Hist->SetFillColor(kGreen);
	dyMadRecoHT200to400->Draw("dileptonMass>>recoJetMadHT200to400Hist(75,0.,1200.)");
	TH1F * recoJetMadHT200to400Hist = (TH1F*) gROOT->FindObject("recoJetMadHT200to400Hist");
	recoJetMadHT200to400Hist->Scale((50.4177*lumi)/973937);
	recoJetMadHT200to400Hist->SetFillColor(kBlue);
	dyMadRecoHT400to600->Draw("dileptonMass>>recoJetMadHT400to600Hist(75,0.,1200.)");
	TH1F * recoJetMadHT400to600Hist = (TH1F*) gROOT->FindObject("recoJetMadHT400to600Hist");
	recoJetMadHT400to600Hist->Scale((6.98394*lumi)/1067758);
	recoJetMadHT400to600Hist->SetFillColor(kCyan);
	dyMadRecoHT600toInf->Draw("dileptonMass>>recoJetMadHT600toInfHist(75,0.,1200.)");
	TH1F * recoJetMadHT600toInfHist = (TH1F*) gROOT->FindObject("recoJetMadHT600toInfHist");
	recoJetMadHT600toInfHist->Scale((2.70354*lumi)/998912);
	recoJetMadHT600toInfHist->SetFillColor(kMagenta);
	recoJetMadHT600toInfHist->GetXaxis()->SetTitle(xAxisLabel);
	recoJetMadHT600toInfHist->GetYaxis()->SetTitle(yAxisLabel);


#ifdef DEBUG
	std::cout<<"made individual histos, applied cross sxn scaling, and set fill color"<<std::endl;
#endif

	//now make a stacked hist, showing all HTBins and the inclusive Mad distribution for RECO HT
	TLegend *legStackedMadInclAndHTBinned = new TLegend( 0.6, 0.60, 0.90, 0.90 ) ; 
	legStackedMadInclAndHTBinned->AddEntry(recoJetMadHT600toInfHist , "MadHT600toInf" ) ; 
	legStackedMadInclAndHTBinned->AddEntry(recoJetMadHT400to600Hist , "MadHT400to600" ) ; 
	legStackedMadInclAndHTBinned->AddEntry(recoJetMadHT200to400Hist , "MadHT200to400" ) ; 
	legStackedMadInclAndHTBinned->AddEntry(recoJetMadHT100to200Hist , "MadHT100to200" ) ; 
	legStackedMadInclAndHTBinned->AddEntry(recoJetMadInclHist , "MadIncl" ) ; 
	legStackedMadInclAndHTBinned->SetFillColor( kWhite ) ; 

#ifdef DEBUG
	std::cout<<"made legend"<<std::endl;
#endif

	TCanvas * canvStackMadInclAndHT = new TCanvas("canvStackMadInclAndHT","canvStackMadInclAndHT",800,800);
	canvStackMadInclAndHT->cd();
#ifdef DEBUG
	std::cout<<"made canvas for stacked bkgnds"<<std::endl;
#endif
	
	THStack* stackedDYMadRecoHT = new THStack();
	stackedDYMadRecoHT->Add(recoJetMadHT600toInfHist);
	stackedDYMadRecoHT->Add(recoJetMadHT400to600Hist);
	stackedDYMadRecoHT->Add(recoJetMadHT200to400Hist);
	stackedDYMadRecoHT->Add(recoJetMadHT100to200Hist);
	stackedDYMadRecoHT->Add(recoJetMadInclHist);
	stackedDYMadRecoHT->SetMinimum(5);
	stackedDYMadRecoHT->SetTitle("CMS Private   #surds = 13 TeV #int lumi = 2.6 fb^{-1}");

#ifdef DEBUG
	std::cout<<"about to draw stacked histo"<<std::endl;
#endif

	stackedDYMadRecoHT->Draw("histo");
	legStackedMadInclAndHTBinned->Draw();

	TString stackedMadInclAndHTOutFileName = "stackedMadInclAndHTBinnedRecoDiElectronMassAllSkims";
	canvStackMadInclAndHT->Print(stackedMadInclAndHTOutFileName+".pdf");
	canvStackMadInclAndHT->Print(stackedMadInclAndHTOutFileName+".png");
	canvStackMadInclAndHT->SetLogy();
	canvStackMadInclAndHT->Print(stackedMadInclAndHTOutFileName+"_log.pdf");
	canvStackMadInclAndHT->Print(stackedMadInclAndHTOutFileName+"_log.png");
	*/


#endif
	//end DYHTPlot
}///end macroSandBox()

