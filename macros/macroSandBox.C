#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
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

#define checkWellSeparatedGenPtBins
//#define PtRatioProfiles
//#define DEBUG
//#define RecoGenOverlays
//#define StudyEffectOfMassPairs
//#define bkgndOverlaidOnMatchedSignal


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
	string histoBeginnings[] = {"mWR2600mNu1300","TTBar","DYPlusJets"};
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

