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

using namespace std;

#define twoDimPlotGenWrAcceptance
//#define recoAndGenHLTEfficiency
//#define genPlotsUsingWRDecayProducts
//#define compareCentrallyProducedToPrivateWrSignal
//#define genWrAndNuMass
//#define signalRegionEEJJBkgnds
//#define lowMassSkimmedBkgndOnRealData
//#define lowMassFlavorSidebandBkgndOnData
//#define checkWellSeparatedGenPtBins
//#define PtRatioProfiles
//#define DEBUG
//#define RecoGenOverlays
//#define StudyEffectOfMassPairs
//#define bkgndOverlaidOnMatchedSignal
//#define DEBUGEVTWEIGHTMTHD
//#define DEBUGVECTOR



///use this fxn to plot and save one histogram using one branch from one TTree (or TChain)
void makeAndSaveSingleHistoFromTree(TChain * chain,string canvName,string treeDrawArgs,string histName,string histTitle,string xAxisTitle,string outputFileName){
	gStyle->SetOptStat("");
	TCanvas * cc = new TCanvas(canvName.c_str(),canvName.c_str(),750,700);
	cc->cd();
	chain->Draw(treeDrawArgs.c_str());
	TH1 * tempHist = (TH1*) gROOT->FindObject(histName.c_str());
	tempHist->SetTitle(histTitle.c_str());
	tempHist->GetXaxis()->SetTitle(xAxisTitle.c_str());
	tempHist->SetLineWidth(3);
	tempHist->SetLineColor(1);

	//adjust the horizontal axis range
	Double_t mean = tempHist->GetMean(1);
	tempHist->GetXaxis()->SetRangeUser((mean)/2,1.3*(mean));
	tempHist->Draw("");
	cc->SaveAs(outputFileName.c_str(),"recreate");
	tempHist->Delete();

}///end makeAndSaveSingleHistoFromTree()


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
void makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea,string title,string xLabel,string outputFileNameModifier,Bool_t specialGrouping){
	gStyle->SetOptStat("");
	TCanvas * canv = new TCanvas(canvName,canvName,750,700);
	canv->cd();
	TLegend * leg = new TLegend(legXmin,legYmin,legXmax,legYmax);
	map<string,TH1F*> overlayHistoMap;	///< links string keys to TH1F histos which will ultimately be overlaid
	for(map<string,TChain*>::const_iterator chMapIt=inputChainMap.begin(); chMapIt!=inputChainMap.end(); chMapIt++){
		size_t openParenth = (chMapIt->first).find_first_of('('), lastChevron = (chMapIt->first).find_last_of('>');
		string uncutHistoName(chMapIt->first);
		///now initialize a new string, get rid of the content in uncutHistoName before '>>' and after '(',
		///and store the substring in the new string object
		string oneHistoName( uncutHistoName.substr(lastChevron+1,openParenth-lastChevron-1) );
		(chMapIt->second)->Draw((chMapIt->first).c_str());
		
		///save pointers to all histograms into overlayHistoMap
		overlayHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
	}///end loop over elements in inputChainMap

	///now overlay all TH1F objects in overlayHistoMap onto one TCanvas
	int colors[] = {1,2,4,8,12,25,30,40,45};
	vector<int> colorVect(colors,colors + sizeof(colors)/sizeof(int) );
	//short colorsSpecial[] = {kBlack,kBlack,kBlack,kRed,kRed,kRed,kBlue,kBlue,kBlue};
	short colorsSpecial[] = {kRed,kBlack,kBlue};
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
	
		size_t lastChevron = (histIt->first).find_last_of('>');
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
	leg->Draw();
	canv->SaveAs(outputFile.c_str(),"recreate");
	canv->SaveAs(outputFilePdf.c_str(),"recreate");

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

#ifdef signalRegionEEJJBkgnds

	TString treeName = "unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts";
	TString dirName = "/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/";
	TChain * ttBarTree = new TChain(treeName,"");
	ttBarTree->Add(dirName+"analyzed_TTOnly_PowhegPythia_25ns_eejj_signal_region.root");
	TChain * dyJetsTree = new TChain(treeName,"");
	dyJetsTree->Add(dirName+"analyzed_DYJets_Madgraph_25ns_eejj_signal_region.root");
	TChain * wJetsTree = new TChain(treeName,"");
	wJetsTree->Add(dirName+"analyzed_WJets_Madgraph_25ns_eejj_signal_region.root");
	//wJetsTree->Add(dirName+"analyzed_WJets_25ns_eejj_signal_region.root");	///< wJets aMCNLO
	TChain * wzTree = new TChain(treeName,"");
	wzTree->Add(dirName+"analyzed_WZ_25ns_eejj_signal_region.root");
	//wzTree->Add(dirName+"analyzed_WZPlusJets_25ns_eejj_signal_region.root");	///< wzPlusJets aMCNLO
	TChain * zzTree = new TChain(treeName,"");
	zzTree->Add(dirName+"analyzed_ZZ_25ns_eejj_signal_region.root");

	Float_t integratedLumi = 1000.;	///inverse picobarns
	string branchNames[] = {"fourObjectMass"};
	string histoEndings[] = {"_fourObjectMass(20,0.,3600.)"};
	string link = ">>";

	TString evWeightCut = "(evWeightSign < 0 ? -1. : 1.)";

	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	//string histoBeginnings[] = {zzKey,wzKey,ttBarKey,dyPlusJetsKey};
	string histoBeginnings[] = {zzKey,wzKey,dyPlusJetsKey};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = dyJetsTree;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = ttBarTree;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = wzTree;
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = zzTree;
		string cName = "o"+to_string(i);
		makeStackedHisto(placeHolderMap,cName.c_str(),0.5,0.75,0.9,0.89,xSxnsFiftyNs,integratedLumi,numEvtsTwentyFiveNs,evWeightCut,false);
		placeHolderMap.clear();
	}///end loop over branchNames
	//makeStackedHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,map<string,Float_t> crossSxnsMap,Float_t intLumi,map<string,Float_t> nEvtsMap,TString treeCuts,Bool_t doLogYaxis){
	

#endif
///end signalRegionEEJJBkgnds


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


#ifdef genWrAndNuMass
	//plot the gen WR and Nu mass distributions using the GEN lvl WR and Nu, not their decay products
	string dir= "/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/";
	string fileBegin = "all_genWrNuAndDecayKinematicsNoMatchingInfo_WR_M-";
	string fileEnd = ".root";
	string fileMiddle = "_Nu_M-";
	string wrMassOutputBegin = "genWRMass_WR_M-";
	string nuMassOutputBegin = "genNuMass_WR_M-";
	string nuWrMassRatioOutputBegin = "genNuMassOverWrMass_WR_M-";
	string outputMiddle = "_Nu_M-";
	string outputEnd = ".png";
	int wrMass[] = {800,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3200,3600,3800,4400,5000,5200,5600,5800,6000};
	//int wrMass[] = {800,1000};
	vector<int> wrMassVect(wrMass,wrMass + sizeof(wrMass)/sizeof(int));

	string massDrawArgs="massWr>>+particleMassHisto";
	string massHistoName="particleMassHisto";

	//loop over all root files, plot the wr mass, and save png image of the plot
	unsigned int maxI = wrMassVect.size();
	for(unsigned int i=0;i<maxI ;i++){
		string pfn = dir + fileBegin + to_string(wrMassVect[i]) + fileMiddle + to_string(wrMassVect[i]/2) + fileEnd;
		TChain * wrChain = new TChain("genWRAnalyzerOne/genWR");
		wrChain->Add(pfn.c_str());
		TChain * nuChain = new TChain("genNuAnalyzerOne/genNu");
		nuChain->Add(pfn.c_str());
		makeAndSaveSingleHistoFromTree(wrChain,"c"+to_string(i),massDrawArgs,massHistoName,"GEN WR mass when Nu mass = "+to_string(wrMassVect[i]/2)+" GeV","WR Mass [GeV]",wrMassOutputBegin+to_string(wrMassVect[i])+outputMiddle+to_string(wrMassVect[i]/2)+outputEnd);
		makeAndSaveSingleHistoFromTree(nuChain,"d"+to_string(i),massDrawArgs,massHistoName,"GEN Nu mass when WR mass = "+to_string(wrMassVect[i])+" GeV","Nu Mass [GeV]",nuMassOutputBegin+to_string(wrMassVect[i])+outputMiddle+to_string(wrMassVect[i]/2)+outputEnd);

		//use the friend tree fxn to plot Nu mass / WR mass
		string chAlias="wrTree", massRatioDrawArgs="(massWr/"+chAlias+".massWr)>>+nuMassOverWrMassHisto", massRatioHistoName="nuMassOverWrMassHisto";
		makeAndSaveSingleHistoFromTreeWithFriend(nuChain,wrChain,chAlias,"e"+to_string(i),massRatioDrawArgs,massRatioHistoName,"GEN Nu Mass / WR Mass when WR mass = "+to_string(wrMassVect[i])+" GeV","Nu Mass / WR Mass",nuWrMassRatioOutputBegin+to_string(wrMassVect[i])+outputMiddle+to_string(wrMassVect[i]/2)+outputEnd);

		nuChain->Delete(), wrChain->Delete();

	}//end loop over root files
	

#endif


#ifdef compareCentrallyProducedToPrivateWrSignal

	//compare distributions of the gen WR mass using the WR particle itself between centrally produced WR->eejj datasets
	//and privately produced WR->eejj datasets
	//use makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea,string title,string xLabel)

	TChain * Gen1000MWR500MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen1000MWR500MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_1000_NU_500_1.root");
	TChain * Gen1000MWR150MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen1000MWR150MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_1000_NU_150_1.root");
	TChain * Gen1000MWR900MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen1000MWR900MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_1000_NU_900_1.root");
	
	TChain * Gen2000MWR1000MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen2000MWR1000MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_2000_NU_1000_1.root");
	TChain * Gen2000MWR200MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen2000MWR200MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_2000_NU_200_1.root");
	TChain * Gen2000MWR1900MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen2000MWR1900MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_2000_NU_1900_1.root");
	
	TChain * Gen3200MWR1600MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen3200MWR1600MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_3200_NU_1600_1.root");
	TChain * Gen3200MWR300MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen3200MWR300MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_3200_NU_300_1.root");
	TChain * Gen3200MWR3100MNu = new TChain("genWRAnalyzerOne/genWR","");
	Gen3200MWR3100MNu->Add("/eos/uscms/store/user/skalafut/WR/13TeV/analyzed_GEN_WRSignal_grid/analyzed_genWrToEEJJFullOfflineAnalysis_WR_3200_NU_3100_1.root");
	
	//TChain * GenMWR1300_matchedGenCentralProduction = new TChain("genWRAnalyzerOne/genWR","");
	//MWR2600MNu1300_matchedGenCentralProduction->Add("/eos/uscms/store/user/skalafut/analyzed_25ns_WR_MC_check_WR_mass/all_genWrKinematics_WR_M-2600_Nu_M-1300.root");

	string branchNames[] = {"massWr"};
	string link=">>";
	//string histoEndings[] = {"_massWR(50,700,1300)"};	//compare privately produced GENSIM with different Nu masses and the same WR mass
	//string histoEndings[] = {"_massWR(50,1600,2400)"};	//compare privately produced GENSIM with different Nu masses and the same WR mass
	string histoEndings[] = {"_massWR(50,2800,3800)"};	//compare privately produced GENSIM with different Nu masses and the same WR mass
		
	//string histoEndings[] = {"_numWrs(4,0.,3.)","_etaWR(60,-6.,6.)","_ptWR(20,0.,200.)","_phiWR(30,-3.2,3.2)","_massWR(50,1800,3000)"};	//compare private GENSIM to centrally produced miniAOD
	string titles[] = {"GEN WR Mass"};
	string xAxisLabels[] = {"WR mass [GeV]"};
	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	//string histoBeginnings[] = {"MWR 1000 MNU 150","MWR 1000 MNU 500","MWR 1000 MNU 900"};	//compare privately produced GENSIM with different Nu masses and the same WR mass
	//string histoBeginnings[] = {"MWR 2000 MNU 200","MWR 2000 MNU 1000","MWR 2000 MNU 1900"};	//compare privately produced GENSIM with different Nu masses and the same WR mass
	string histoBeginnings[] = {"MWR 3200 MNU 300","MWR 3200 MNU 1600","MWR 3200 MNU 3100"};	//compare privately produced GENSIM with different Nu masses and the same WR mass
	//string histoBeginnings[] = {"MWR 1000 MNU 150","MWR 1000 MNU 500","MWR 1000 MNU 900","MWR 2000 MNU 200","MWR 2000 MNU 1000","MWR 2000 MNU 1900","MWR 3200 MNU 300","MWR 3200 MNU 1600","MWR 3200 MNU 3100"};	//compare privately produced GENSIM with different Nu masses and the same WR mass
	
	//string histoBeginnings[] = {"Central","Private"};	//compare private GENSIM to centrally produced miniAOD
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		//placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = Gen1000MWR150MNu;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = Gen1000MWR500MNu;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = Gen1000MWR900MNu;
	
		//placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = Gen2000MWR200MNu;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = Gen2000MWR1000MNu;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = Gen2000MWR1900MNu;
	
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = Gen3200MWR300MNu;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = Gen3200MWR1600MNu;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = Gen3200MWR3100MNu;
		
		//placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = MWR2600MNu1300_matchedGenCentralProduction;
		//placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = MWR2600MNu1300;
		
		string cName = "o"+to_string(i);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.65,0.45,0.95,0.90,true,titles[i]+" at different GEN Nu Masses",xAxisLabels[i],"_GEN_WR_mass_3200_Nu_masses_300_1600_3100",true);
		//makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.75,0.6,0.98,0.95,true,titles[i]+" Central vs Private Production",xAxisLabels[i],"_GEN_WR_mass_multipleNuAndWRMasses");
		placeHolderMap.clear();
	}///end loop over branchNames


#endif

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
	gStyle->SetTitleOffset(1.5,"Y");
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
	///the ratio of (RECO evts passing all WR signal region cuts)/(GEN evts passing all WR signal region cuts)
	///the RECO cuts will include HEEP and jet ID requirements, whereas the GEN cuts will not
	///neither set of cuts requires that the trigger is fired (HLT_DoubleEle33)
	string recoDir = "/eos/uscms/store/user/skalafut/analyzed_25ns_eejj_signal_region/";
	string recoFileBegin = "all_analyzed_tree_checkHLTEfficiencyRecoFullOfflineAndGenPtEtaDr_eejjSignalRegion_WR_M-";
	string recoFileEnd = ".root";
	string recoFileMiddle = "_Nu_M-";
	string RecoOvrGenvsMassFile = "recoOvrGenScaleFactorsVsMass.txt";
	ofstream writeToRecoGenScaleFactorsFile(RecoOvrGenvsMassFile.c_str(),ofstream::app);
	Int_t nBins = 12;	///max is 12
	Float_t wrMassVals[nBins], recoGenScaleFactors[nBins];

	int wrMassArr[] = {800,1000,1200,1400,1600,2000,2200,2400,2600,2800,3000,3200};
	//int wrMassArr[] = {800,1000};

	for(int i=0; i<nBins ; i++){
		string recoPfn = recoDir+recoFileBegin+to_string(wrMassArr[i])+recoFileMiddle+to_string(wrMassArr[i]/2)+recoFileEnd;
		string genPfn = dir+fileBegin+to_string(wrMassArr[i])+fileMiddle+to_string(wrMassArr[i]/2)+fileEnd;
		wrMassVals[i] = (Float_t) wrMassArr[i];
	
		///for simplicity, assume 50k MINIAOD evts
		TChain * genAfterCuts = new TChain("genMatchedParticleAnalyzerThree/genLeptonsAndJetsWithAllCuts");
		genAfterCuts->Add(genPfn.c_str());
		TChain * recoAfterCuts = new TChain("unmatchedSignalRecoAnalyzerFive/signalRecoObjectsWithAllCuts");
		recoAfterCuts->Add(recoPfn.c_str());
		TChain * genBeforeCuts = new TChain("genNuAnalyzerOne/genNu");
		genBeforeCuts->Add(genPfn.c_str());
	
		Float_t recoRatio = ((Float_t) recoAfterCuts->GetEntries()/50000);
		Float_t genRatio = ((Float_t) genAfterCuts->GetEntries()/genBeforeCuts->GetEntries());
		recoGenScaleFactors[i] = recoRatio/genRatio;
		writeToRecoGenScaleFactorsFile << "for MWR=\t"<< wrMassArr[i] << "\tand MNu=\t"<< wrMassArr[i]/2 << "\tthe ratio of RECO evts passing all cuts to GEN evts passing all cuts (RECO/GEN) is\t"<< recoGenScaleFactors[i] << endl;

		genAfterCuts->Delete();
		recoAfterCuts->Delete();

	}
	writeToRecoGenScaleFactorsFile.close();

	///make the TGraph, and save an image of it
	TGraph * recoGenFactorsGraph = new TGraph(nBins,wrMassVals,recoGenScaleFactors);
	recoGenFactorsGraph->GetXaxis()->SetTitle("WR Mass [GeV]");
	recoGenFactorsGraph->GetYaxis()->SetTitle("RECO/GEN");
	recoGenFactorsGraph->SetTitle("RECO evts passing / GEN evts passing  MNu = MWR/2");
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
	
}///end macroSandBox()

