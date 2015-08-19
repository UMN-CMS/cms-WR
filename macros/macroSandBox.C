#include <TFile.h>
#include <TList.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
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
#include <vector>
#include <map>
#include <TEventList.h>
#include <TEntryList.h>
#include <TEntryListArray.h>
#include "dumpTreePlots.C"

using namespace std;

#define lowMassSkimmedBkgndOnRealData
//#define checkWellSeparatedGenPtBins
//#define PtRatioProfiles
//#define DEBUG
//#define RecoGenOverlays
//#define StudyEffectOfMassPairs
//#define bkgndOverlaidOnMatchedSignal
//#define DEBUGGETEVTWGT

/**
 * use this fxn to compute the evt weight for one evt in a TChain
 * return the weight as a Float_t value
 */
Float_t getEvtWeight(Float_t mcEvWgtSign, Int_t numVertices, map<string,vector<Float_t> > mcPuWeights, string mcName){
	Float_t weight = 0;
	for(map<string,vector<Float_t> >::const_iterator puWgtIt=mcPuWeights.begin(); puWgtIt!=mcPuWeights.end(); puWgtIt++){
#ifdef DEBUGGETEVTWGT
		cout<<"pu weights map has string \t"<< puWgtIt->first <<endl;
		cout<<"mc process name is \t"<< mcName <<endl;
#endif
		///loop over the unique keys in mcPuWeights and find the one which matches the MC process name (ttBar, WZ, etc)
		if( (puWgtIt->second).size() < numVertices ) return 0;
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
	string val="45";

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
 * use this fxn to overlay kinematic distributions from real data (points) onto MC (filled histos) using a THStack object and an overlaid TH1F object
 * the keys of crossSxnsMap and nEvtsMap are names of physics processes, like ttBar and dyPlusJets
 * the keys of inputChainMap contains the histogram plotting argument to use with TChain->Draw("plottingArgs")
 * the areas of the MC histograms are normalized to the integrated luminosity of real data
 * this is done by calling Scale(crossSxn * integratedLumi / numEvts) on each of the histos from MC data
 * if there are N unique keys in inputChainMap, then there are N-1 unique keys in the maps called crossSxnsMap and nEvtsMap
 *
 * the TChain to real data should contain the phrase "ealData" in the inputChainMap key
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
		//if( (chMapIt->first).find("ealData") != string::npos) continue;
		size_t openParenth = (chMapIt->first).find_first_of('('), lastChevron = (chMapIt->first).find_last_of('>');
		string uncutHistoName(chMapIt->first);
		///now initialize a new string, get rid of the content in uncutHistoName before '>>' and after '(',
		///and store the substring in the new string object
		string oneHistoName( uncutHistoName.substr(lastChevron+1,openParenth-lastChevron-1) );
		string afterLastChevron( uncutHistoName.substr(lastChevron+1) );
		if((chMapIt->first).find("ealData") == string::npos ){
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
				(chMapIt->second)->SetBranchAddress("nVertices",&vertices);
				(chMapIt->second)->SetBranchAddress("evWeightSign",&evWgtSign);
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
						wgt = getEvtWeight(evWgtSign, vertices, puWeights, bkgndString);
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
						wgt = getEvtWeight(evWgtSign, vertices, puWeights, bkgndString);
						hTemp->Fill(desiredArray[stoi(arrElementString)], wgt);
					}

				}///end Float_t[] array branch filter
			
			}///end if(doPuReweighting == true)
			else if(doPuReweighting == false) (chMapIt->second)->Draw((chMapIt->first).c_str(), treeCuts);


			stackedHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
		}
		else if((chMapIt->first).find("ealData") != string::npos ){
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
		if((chMapIt->first).find("ealData") == string::npos ){
			stackedHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
		}///end search for histos made from MC datasets
		
		else if((chMapIt->first).find("ealData") != string::npos){
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

	///now add all TH1F objects in stackedHistoMap into one THStack object, and overlay the
	///histo in overlaidHistoMap onto the THStack object
	///draw the THStack object first!
	///1 = black for colors
	int colors[] = {2,4,5,8,12,30,40,45};
	vector<int> colorVect(colors,colors + sizeof(colors)/sizeof(int) );
	Int_t i=0;
	THStack * histoStack = new THStack("","");
	if(stackedHistoMap.size() > colorVect.size() ) cout<<"not enough unique colors in MultipleCurveOverlayHisto fxn!"<<endl;
#ifdef DEBUG
	std::cout<<"made a THStack object with null name and title"<<std::endl;
#endif

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

#ifdef DEBUG
		std::cout<<"rescaled the histo \t"<< histIt->first <<std::endl;
		std::cout<<"after rescaling there are "<< (histIt->second)->Integral() <<" entries in the histo mentioned immediately above"<<std::endl;
#endif
	
		histoStack->Add((histIt->second));	///< add each histo in stackedHistoMap to the THStack object

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

	Double_t originalMax = histoStack->GetMaximum();
	histoStack->SetMaximum(20*originalMax);
	histoStack->SetMinimum(0.1);
	histoStack->Draw("hist");

#ifdef DEBUG
		std::cout<<"stacked histo has been drawn"<<std::endl;
#endif

	string outputFile;
	for(map<string,TH1F*>::const_iterator hIt=overlaidHistoMap.begin(); hIt!=overlaidHistoMap.end(); hIt++){
		size_t lastChevron = (hIt->first).find_last_of('>');
		size_t underscorePos = (hIt->first).find_first_of("_",lastChevron);
		string legEntryName = (hIt->first).substr(lastChevron+1,underscorePos-lastChevron-1);
		leg->AddEntry(hIt->second,legEntryName.c_str(),"ep");
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

}///end overlayPointsOnStackedHistos


/** the TString key in inputChainMap contains the histogram plotting argument to use with TChain->Draw("plottingArgs")
 * the histogram name does not need to be passed into the fxn as an argument, it will be pulled from the map key
 * if doNormalizationByArea is true, then normalize the plotted histos so that the area under each curve = 1
 *
 *
 */
void makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,Bool_t doNormalizationByArea){
	TCanvas * canv = new TCanvas(canvName,canvName,750,700);
	canv->cd();
	TLegend * leg = new TLegend(legXmin,legYmin,legXmax,legYmax);
	map<string,TH1F*> overlayHistoMap;	///< links TString keys to TH1F histos which will ultimately be overlaid
	for(map<string,TChain*>::const_iterator chMapIt=inputChainMap.begin(); chMapIt!=inputChainMap.end(); chMapIt++){
		size_t openParenth = (chMapIt->first).find_first_of('('), lastChevron = (chMapIt->first).find_last_of('>');
		string uncutHistoName(chMapIt->first);
		///now initialize a new string, get rid of the content in uncutHistoName before '>>' and after '(',
		///and store the substring in the new string object
		string oneHistoName( uncutHistoName.substr(lastChevron+1,openParenth-lastChevron-1) );
		(chMapIt->second)->Draw((chMapIt->first).c_str());
		overlayHistoMap[chMapIt->first]= (TH1F*) gROOT->FindObject(oneHistoName.c_str());
	}///end loop over elements in inputChainMap
	///now overlay all TH1F objects in overlayHistoMap onto one TCanvas
	int colors[] = {1,2,4,8,12,30,40,45};
	vector<int> colorVect(colors,colors + sizeof(colors)/sizeof(int) );
	Int_t i=0;
	if(overlayHistoMap.size() > colorVect.size() ) cout<<"not enough unique colors in MultipleCurveOverlayHisto fxn!"<<endl;
	for(map<string,TH1F*>::const_iterator histIt=overlayHistoMap.begin(); histIt!=overlayHistoMap.end(); histIt++){
		(histIt->second)->SetLineColor(colorVect[i]);
		(histIt->second)->SetLineWidth(3);
		if(doNormalizationByArea){
			Double_t oldIntegral = (histIt->second)->Integral();
			(histIt->second)->Scale(1.0/oldIntegral);
		}///end filter to normalize histo area to 1
		if(histIt==overlayHistoMap.begin()){
			string xLabel="";
			if((histIt->first).find("Mass") !=string::npos || (histIt->first).find("ptEle") !=string::npos || (histIt->first).find("ptJet") !=string::npos){
				xLabel += "GeV";
			}///end if(histo is plotting a variable with dimension of energy)
			(histIt->second)->GetXaxis()->SetTitle(xLabel.c_str());
			Double_t oldMax = (histIt->second)->GetBinContent((histIt->second)->GetMaximumBin());
			(histIt->second)->SetMaximum(3*oldMax);
		}///end filter to set histogram X axis label
		size_t lastChevron = (histIt->first).find_last_of('>');
		size_t underscorePos = (histIt->first).find_first_of("_",lastChevron);
		string legEntryName = (histIt->first).substr(lastChevron+1,underscorePos-lastChevron);
#ifdef DEBUG
		cout<<"histo name = \t"<< histIt->first <<endl;
		cout<<"legEntry has name = \t"<< legEntryName << endl;
#endif
		leg->AddEntry(histIt->second,legEntryName.c_str());
		i++;
	}///end loop to set line colors of histos, and add entries to TLegend

	string outputFile;
	for(map<string,TH1F*>::const_iterator hIt=overlayHistoMap.begin(); hIt!=overlayHistoMap.end(); hIt++){
		if(hIt==overlayHistoMap.begin()){
			(hIt->second)->Draw();
			size_t firstChevron = (hIt->first).find_first_of('>');
			outputFile = (hIt->first).substr(0,firstChevron);
			outputFile += ".png";
#ifdef DEBUG
			cout<<"outputFile = \t"<< outputFile <<endl;
#endif
		}
		else (hIt->second)->Draw("same");
	}///end loop which draws histograms
	leg->Draw();
	canv->SaveAs(outputFile.c_str(),"recreate");

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
	xSxnsFiftyNs[ttBarKey]=815.9;
	xSxnsFiftyNs[dyPlusJetsKey]=6025.2;
	xSxnsFiftyNs[wJetsKey]=61500;
	xSxnsFiftyNs[wzKey]=66.1;
	xSxnsFiftyNs[zzKey]=15.4;
	map<string,Float_t> numEvtsFiftyNs;
	numEvtsFiftyNs[ttBarKey]=4994250;
	numEvtsFiftyNs[dyPlusJetsKey]=19925500;
	numEvtsFiftyNs[wJetsKey]=24089991;
	numEvtsFiftyNs[wzKey]=996920;
	numEvtsFiftyNs[zzKey]=998848;
	


#ifdef lowMassSkimmedBkgndOnRealData
	//overlayPointsOnStackedHistos(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax,map<string,Float_t> crossSxnsMap,Float_t intLumi,map<string,Float_t> nEvtsMap, TString treeCuts, Bool_t doPuReweighting, map<string,vector<Float_t> > puWeights,Bool_t doLogYaxis)
	string moduleAndTreeName = "recoAnalyzerTwo/recoObjectsWithPtEtaCuts";
	TChain * dyPlusJetsEEJJLowMassSkim = new TChain(moduleAndTreeName.c_str());
	dyPlusJetsEEJJLowMassSkim->Add("/eos/uscms/store/user/skalafut/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/analyzed_DYJets_50ns_skim_low_mass_region_eejj.root");
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


	Float_t integratedLumi = 40.001/0.962;	///< in picobarns

	//string branchNames[] = {"ptEle[0]","ptEle[1]","etaEle[0]","etaEle[1]","ptJet[0]","ptJet[1]","etaJet[0]","etaJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet","leadLeptonThreeObjMass","subleadingLeptonThreeObjMass"};
	string link=">>";
	//string histoEndings[] = {"_leadLeptonPt(40,0.,200.)","_subleadLeptonPt(20,0.,100.)","_leadLeptonEta(50,-3.0,3.0)","_subleadLeptonEta(50,-3.0,3.0)","_leadJetPt(50,0.,200.)","_subleadJetPt(50,0.,100.)","_leadJetEta(50,-3.0,3.0)","_subleadJetEta(50,-3.0,3.0)","_dileptonMass(50,0.,400.)","_fourObjectMass(50,0.,600.)","_dR_leadingLeptonLeadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingJet(50,0.,5.)","_dR_subleadingLeptonLeadingJet(50,0.,5.)","_dR_subleadingLeptonSubleadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingLepton(50,0.,5.)","_dR_leadingJetSubleadingJet(50,0.,5.)","_leadLeptonThreeObjMass(50,0.,500.)","_subleadingLeptonThreeObjMass(50,0.,500.)"};
	
	string branchNames[] = {"ptEle[0]","ptEle[1]","ptJet[0]","ptJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet"};
	string histoEndings[] = {"_leadLeptonPt(20,30.,200.)","_subleadLeptonPt(20,30.,110.)","_leadJetPt(25,30.,200.)","_subleadJetPt(20,30.,110.)","_dileptonMass(25,0.,200.)","_fourObjectMass(20,0.,600.)","_dR_leadingLeptonLeadingJet(30,0.,5.)","_dR_leadingLeptonSubleadingJet(30,0.,5.)","_dR_subleadingLeptonLeadingJet(30,0.,5.)","_dR_subleadingLeptonSubleadingJet(30,0.,5.)","_dR_leadingLeptonSubleadingLepton(30,0.,5.)","_dR_leadingJetSubleadingJet(30,0.,5.)"};
	
	//string branchNames[] = {"leadLeptonThreeObjMass","subleadingLeptonThreeObjMass"};
	//string histoEndings[] = {"_leadLeptonThreeObjMass(25,50.,600.)","_subleadingLeptonThreeObjMass(25,50.,600.)"};
	
	//string branchNames[] = {"etaEle[0]","etaEle[1]","etaJet[0]","etaJet[1]"};
	//string histoEndings[] = {"_leadLeptonEta(15,-3.0,3.0)","_subleadLeptonEta(15,-3.0,3.0)","_leadJetEta(15,-3.0,3.0)","_subleadJetEta(15,-3.0,3.0)"};
	
	//string branchNames[] = {"nJets","nLeptons","nVertices"};
	//string histoEndings[] = {"_nJets(12,0.,12.)","_nLeptons(6,0.,6.)","_nVertices(25,0.,50.)"};
	
	TString evWeightCut = "(evWeightSign < 0 ? -1. : 1.)";

	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {"doubleEGRealData",ttBarKey,dyPlusJetsKey,wJetsKey,wzKey,zzKey};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=0; i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = doubleEGEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = ttBarEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = dyPlusJetsEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[3]+histoEndings[i]] = wJetsEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[4]+histoEndings[i]] = wzEEJJLowMassSkim;
		placeHolderMap[branchNames[i]+link+histoBeginnings[5]+histoEndings[i]] = zzEEJJLowMassSkim;
		string cName = "o"+to_string(i);
		overlayPointsOnStackedHistos(placeHolderMap,cName.c_str(),0.5,0.67,0.98,0.9,xSxnsFiftyNs,integratedLumi,numEvtsFiftyNs,evWeightCut,false,pileupWeights,true);
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
	dyPlusJetsAllCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/dyPlusJets/*.root");
	TChain * ttBarAllCuts = new TChain("bkgndRecoAnalyzerFive/bkgndRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts");
	ttBarAllCuts->Add("/eos/uscms/store/user/skalafut/WR/13TeV/bkgnds/ttBar/*.root");
	TChain * matchedRecoNoCuts = new TChain("matchedRecoAnalyzerOne/matchedRecoObjectsNoCuts");
	matchedRecoNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	
	///get # of evts before all cuts, and after all cuts are applied
	//cout<<"ttBar evts before cuts = \t"<< ttBarNoCuts->GetEntriesFast() <<"\t after cuts = \t"<< ttBarAllCuts->GetEntriesFast() << endl;
	//cout<<"dyPlusJets evts before cuts = \t"<< dyPlusJetsNoCuts->GetEntriesFast() <<"\t after cuts = \t"<< dyPlusJetsAllCuts->GetEntriesFast() << endl;
	//cout<<"signal evts before cuts = \t"<< matchedRecoNoCuts->GetEntriesFast() <<"\t after cuts = \t"<< matchedRecoPtEtaDileptonMassDrFourObjMassCuts->GetEntriesFast() << endl;
	//ttBarNoCuts->Draw("evtNumber","");
	//ttBarAllCuts->Draw("evtNumber","");
	//dyPlusJetsNoCuts->Draw("evtNumber","");
	//dyPlusJetsAllCuts->Draw("evtNumber","");
	matchedRecoNoCuts->Draw("evtNumber","");


	///setup inputs needed for makeAndSaveMultipleCurveOverlayHist() fxn, and call this fxn
	string branchNames[] = {"ptEle[0]","ptEle[1]","etaEle[0]","etaEle[1]","ptJet[0]","ptJet[1]","etaJet[0]","etaJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet","subleadingLeptonThreeObjMass","leadLeptonThreeObjMass"};
	string link=">>";
	string histoEndings[] = {"_leadLeptonPt(50,0.,1200.)","_subleadLeptonPt(50,0.,700.)","_leadLeptonEta(50,-3.0,3.0)","_subleadLeptonEta(50,-3.0,3.0)","_leadJetPt(70,0.,900.)","_subleadJetPt(50,0.,500.)","_leadJetEta(50,-3.0,3.0)","_subleadJetEta(50,-3.0,3.0)","_dileptonMass(50,0.,2500.)","_fourObjectMass(50,400.,3300.)","_dR_leadingLeptonLeadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingJet(50,0.,5.)","_dR_subleadingLeptonLeadingJet(50,0.,5.)","_dR_subleadingLeptonSubleadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingLepton(50,0.,5.)","_dR_leadingJetSubleadingJet(50,0.,5.)","_subleadingLeptonThreeObjMass(50,0.,1600.)","_leadLeptonThreeObjMass(50,0.,2800.)",};
	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {"mWR2600mNu1300","ttBar","dyPlusJets"};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	/*
	for(unsigned int i=0; i<1; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = matchedRecoPtEtaDileptonMassDrFourObjMassCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = ttBarAllCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = dyPlusJetsAllCuts;
		string cName = "o"+to_string(i);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.75,0.6,0.98,0.95,true);
		placeHolderMap.clear();
	}///end loop over branchNames
	*/

#endif

}///end macroSandBox()

