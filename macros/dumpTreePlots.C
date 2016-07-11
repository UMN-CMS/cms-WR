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
#include <TObjArray.h>
#include "TCollection.h"
#include <TCut.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <TEventList.h>
#include <TEntryList.h>
#include <TEntryListArray.h>

//using namespace std;

//#define DEBUG
//#define GENLVL
#define GENMATCHEDRECO

///calculate the x axis limits of a plot such that 99.9 percent of the histogram integral is shown
///return a pair of values representing the low and high x axis values
std::pair<Double_t,Double_t> determinePlotDomain(TH1F * histo){
#ifdef DEBUG
	std::cout<<"in determinePlotDomain fxn"<<std::endl;
#endif
	Double_t trueIntegral = histo->Integral();
	Double_t integralFrxn, tempIntegral;
	Int_t nBins = histo->GetNbinsX(), lowBin, highBin;
	Double_t lowVal = -10, highVal = -10;
	std::string hName(histo->GetName());
	///IF the plot is not a GEN mass plot
	///start from a bin half way between the furthest left and furthest right bin, and expand the integral range
	///for each loop iteration
	if(hName.find("mass") == string::npos){
#ifdef DEBUG
		std::cout<<"the histo named\t"<< hName <<"\tdoes not contain the word mass"<<std::endl;
#endif
		for(Int_t i=0; i<((Int_t) nBins/2)-1; i++){
			lowBin = ((Int_t) nBins/2)-i, highBin = ((Int_t) nBins/2)+i;
			tempIntegral = histo->Integral(lowBin,highBin);
			if(tempIntegral/trueIntegral >= 0.999){
				lowVal = histo->GetXaxis()->GetBinLowEdge(lowBin);
				highVal = (histo->GetXaxis()->GetBinCenter(highBin)) + ( (histo->GetXaxis()->GetBinWidth(2))/2 );
				break;
			}//end if(tempIntegral is within 99.9 percent of total histo integral)
		}//end loop over histo bins
	}//end if(histo name does not contain the word mass)

	///ELSE (the plot is a GEN mass plot) start at the bin which has the greatest number of entries, and expand
	///the integral range above and below this bin
	else{
		Int_t maxContentBin = histo->GetMaximumBin();	///<expand Integral() range symmetrically about this bin
#ifdef DEBUG
		std::cout<<"the histo named\t"<< hName <<"\tdoes contain the word mass"<<std::endl;
		std::cout<<"the bin number with the greatest number of entries is\t"<< maxContentBin <<std::endl;
#endif
		Double_t maxContent = histo->GetBinContent(maxContentBin);
		Double_t thresholdBinContents = maxContent/5;	///<the number of entries in lowBin and highBin must be less than this threshold
		for(Int_t i=0; i<((Int_t) nBins/2)-1; i++){
			lowBin = maxContentBin-i, highBin = maxContentBin+i;
			if(lowBin < 1) lowBin = 1;
			if(highBin > nBins) highBin = nBins;
			tempIntegral = histo->Integral(lowBin,highBin);
			if(tempIntegral/trueIntegral >= 0.999 && histo->GetBinContent(lowBin) <= thresholdBinContents && histo->GetBinContent(highBin) <= thresholdBinContents){
				lowVal = histo->GetXaxis()->GetBinLowEdge(lowBin);
				highVal = (histo->GetXaxis()->GetBinCenter(highBin)) + ( (histo->GetXaxis()->GetBinWidth(2))/2 );
				break;
			}//end if(tempIntegral is within 99.9 percent of total histo integral)
		}//end loop over histo bins
	}//end else(plot shows a GEN mass distribution

#ifdef DEBUG
	std::cout<<"min domain value=\t"<<lowVal<<std::endl;
	std::cout<<"max domain value=\t"<<highVal<<std::endl;
#endif
	std::pair<Double_t,Double_t> domainExtrema = std::make_pair(lowVal,highVal);
	return domainExtrema;

}///end determinePlotDomain()


///adjust the lower and upper limit values for the horizontal axis of the histo
void resetXaxisLimits(TH1F * histo){
	string hName(histo->GetName());
	Double_t mean = histo->GetMean(1);
	Double_t rms = histo->GetRMS(1);
	Double_t peakVal = histo->GetXaxis()->GetBinCenter(histo->GetMaximumBin());
	Double_t min = 0, max = 0;

	///easy to compute limits based on histo name
	if(hName.find("phi") != string::npos){
		histo->GetXaxis()->SetRangeUser(-3.3,3.3);
	}///end phi
	else{
		std::pair<Double_t,Double_t> plotMinAndMax = determinePlotDomain(histo);
		histo->GetXaxis()->SetRangeUser(plotMinAndMax.first, plotMinAndMax.second);
	}

}///end resetXaxisLimits()

void resetNumXaxisBins(TH1F * histo){
	Int_t oldNbins = histo->GetNbinsX();
	Double_t min = histo->GetXaxis()->GetBinCenter(1);
	Double_t max = histo->GetXaxis()->GetBinCenter(oldNbins);
	histo->SetBins(30,min,max);

}///end resetNumXaxisBins()

//use this fxn to compare distributions from two TChains 
//this function is essentially two copies of makeAndSaveHistoUsingEntryList()
//two TChains, two listFillArgs, two list names, two histPlotArgs, two histNames, one histTitle,
//one xAxisTitle, one canvName, on TCut object, one outputFile, and two Bool_t args are given
//to this fxn as inputs
void makeAndSaveOverlayHistoUsingEntryLists(TChain * sigChain,TChain * bkgndChain,TString sigListFillArgs,TString sigListName,TString bkgndListFillArgs,TString bkgndListName,TString sigHistPlotArg,TString bkgndHistPlotArg,TString sigHistName,TString bkgndHistName,TString histTitle,TString xAxisTitle,TString canvName,TCut sigFilt,TCut bkgndFilt,TString outputFile,Bool_t isPlottingEnergy,Bool_t isPlottingInverseEnergy,Bool_t normalizeHistos){
	sigChain->Draw(sigListFillArgs,sigFilt,"entrylistarray");
	sigChain->SetEntryList((TEntryListArray*) gROOT->FindObject(sigListName) );
	bkgndChain->Draw(bkgndListFillArgs,bkgndFilt,"entrylistarray");
	bkgndChain->SetEntryList((TEntryListArray*) gROOT->FindObject(bkgndListName) );
	gStyle->SetOptStat("emriou");
	
	TCanvas * canv = new TCanvas(canvName,canvName,700,700);
	canv->cd();
	sigChain->Draw(sigHistPlotArg);
	TH1F * sigHist = (TH1F*) gROOT->FindObject(sigHistName);
	bkgndChain->Draw(bkgndHistPlotArg);
	TH1F * bkgndHist = (TH1F*) gROOT->FindObject(bkgndHistName);

	if(normalizeHistos){
		Double_t sigIntegral = sigHist->Integral();
		sigHist->Scale(1/sigIntegral);
		Double_t bkgndIntegral = bkgndHist->Integral();
		bkgndHist->Scale(1/bkgndIntegral);
	}
	
	sigHist->SetLineColor(1);	//black
	bkgndHist->SetLineColor(2);	//red

	//sigHist will be drawn first, bkgndHist overlaid on top.  If the largest bin in bkgndHist > the largest bin in sigHist,
	//then increase the y axis max on sigHist to accommodate the peak in bkgndHist
	if(sigHist->GetBinContent(sigHist->GetMaximumBin()) < bkgndHist->GetBinContent(bkgndHist->GetMaximumBin()) ){
		sigHist->SetMaximum((1.1)*( bkgndHist->GetBinContent(bkgndHist->GetMaximumBin()) ) );
	}
	sigHist->SetTitle(histTitle);
	//if isPlottingEnergy or isPlottingInverseEnergy is true, then append units to the x axis label
	TString enrg = " (GeV)";
	TString invEnrg = " (1/GeV)";
	TString completeXaxisTitle;
	if(isPlottingEnergy) completeXaxisTitle = xAxisTitle+enrg;
	if(isPlottingInverseEnergy) completeXaxisTitle = xAxisTitle+invEnrg;
	if(!isPlottingEnergy && !isPlottingInverseEnergy) completeXaxisTitle = xAxisTitle;
	sigHist->GetXaxis()->SetTitle(completeXaxisTitle);
	if(histTitle.Contains("Iso") ){
		canv->SetLogy(1);
		sigHist->SetMinimum(1);
	}
	char temp[130];
	if(isPlottingInverseEnergy && sigHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f / GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if(isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f / GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.001){
		sprintf(temp,"Events / %.4f / GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f ", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f ", sigHist->GetXaxis()->GetBinWidth(1));
	}
	sigHist->GetYaxis()->SetTitle(temp);
	sigHist->Draw();
	bkgndHist->Draw("same");
	canv->SaveAs(outputFile,"recreate");
	
}//end makeAndSaveOverlayHistoUsingEntryLists()


///make an overlay plot without using entrylistarray cuts
void makeAndSaveOverlayHisto(TChain * sigChain,TChain * bkgndChain,TString sigHistPlotArg,TString bkgndHistPlotArg,TString sigHistName,TString bkgndHistName,TString histTitle,TString xAxisTitle,TString canvName,TCut sigFilt,TCut bkgndFilt,TString outputFile,Bool_t isPlottingEnergy,Bool_t isPlottingInverseEnergy,Bool_t normalizeHistos,TString sigLegendLabel,TString bkgndLegendLabel,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax){
	TCanvas * canv = new TCanvas(canvName,canvName,750,700);
	canv->cd();
	sigChain->Draw(sigHistPlotArg);
	TH1F * sigHist = (TH1F*) gROOT->FindObject(sigHistName);
	bkgndChain->Draw(bkgndHistPlotArg);
	TH1F * bkgndHist = (TH1F*) gROOT->FindObject(bkgndHistName);
	gStyle->SetOptStat("emriou");
	
	if(normalizeHistos){
		Double_t sigIntegral = sigHist->Integral();
		sigHist->Scale(1/sigIntegral);
		Double_t bkgndIntegral = bkgndHist->Integral();
		bkgndHist->Scale(1/bkgndIntegral);
	}
	
	sigHist->SetLineColor(1);	//black
	bkgndHist->SetLineColor(2);	//red
	TLegend * leg = new TLegend(legXmin,legYmin,legXmax,legYmax);
	leg->AddEntry(sigHist,sigLegendLabel);
	leg->AddEntry(bkgndHist,bkgndLegendLabel);


	//sigHist will be drawn first, bkgndHist overlaid on top.  If the largest bin in bkgndHist > the largest bin in sigHist,
	//then increase the y axis max on sigHist to accommodate the peak in bkgndHist
	if(sigHist->GetBinContent(sigHist->GetMaximumBin()) < bkgndHist->GetBinContent(bkgndHist->GetMaximumBin()) ){
		sigHist->SetMaximum((1.1)*( bkgndHist->GetBinContent(bkgndHist->GetMaximumBin()) ) );
	}
	sigHist->SetTitle(histTitle);
	//if isPlottingEnergy or isPlottingInverseEnergy is true, then append units to the x axis label
	TString enrg = " (GeV)";
	TString invEnrg = " (1/GeV)";
	TString completeXaxisTitle;
	if(isPlottingEnergy) completeXaxisTitle = xAxisTitle+enrg;
	if(isPlottingInverseEnergy) completeXaxisTitle = xAxisTitle+invEnrg;
	if(!isPlottingEnergy && !isPlottingInverseEnergy) completeXaxisTitle = xAxisTitle;
	sigHist->GetXaxis()->SetTitle(completeXaxisTitle);
	char temp[130];
	if(isPlottingInverseEnergy && sigHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f / GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if(isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f / GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.001){
		sprintf(temp,"Events / %.4f / GeV", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f ", sigHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && sigHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f ", sigHist->GetXaxis()->GetBinWidth(1));
	}
	sigHist->GetYaxis()->SetTitle(temp);
	sigHist->GetYaxis()->SetTitleOffset(1.5);
	sigHist->Draw();
	bkgndHist->Draw("same");
	leg->Draw();
	canv->SaveAs(outputFile,"recreate");
	
}//end makeAndSaveOverlayHisto()



void makeAndSaveHistoUsingEntryList(TChain * chain,TString listFillArgs,TString listName,TString histPlotArgs,TString histName,TString histTitle,TString xAxisTitle,TString canvName,TCut filters,TString outputFile, Bool_t isPlottingEnergy, Bool_t isPlottingInverseEnergy){
	//fill the TEntryList named listName, and apply the filters when the list is made
	chain->Draw(listFillArgs,filters,"entrylistarray");
	//tell the chain to only use entries in the object named listName when calling TTree::Draw() in the future
	chain->SetEntryList((TEntryListArray*) gROOT->FindObject(listName) );
	gStyle->SetOptStat("emriou");
	
	//run the code in makeAndSaveSingleTreeHisto() to make comprehendible plots with useful labels and title 
	TCanvas * canv = new TCanvas(canvName,canvName,500,500);
	canv->cd();
	chain->Draw(histPlotArgs);
	TH1F * pHist = (TH1F*) gROOT->FindObject(histName);
	pHist->SetTitle(histTitle);
	pHist->GetXaxis()->SetTitle(xAxisTitle);
	if(histTitle.Contains("Iso") ){
		canv->SetLogy(1);
		pHist->SetMinimum(1);
	}
	//every histo should have at least three bins.
	//if a histo is created with numBins = 1, then bin #1 is the bin which is plotted 
	char temp[130];
	if(isPlottingInverseEnergy && pHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f / GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if(isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f / GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.001){
		sprintf(temp,"Events / %.4f / GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f ", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f ", pHist->GetXaxis()->GetBinWidth(1));
	}
	pHist->GetYaxis()->SetTitle(temp);
	pHist->Draw();
	canv->SaveAs(outputFile,"recreate");

}//end makeAndSaveHistoUsingEntryList()


//use this to make and save a histogram from a single TTree branch
//plotArgs is used in TTree::Draw, and could be something like "etaGenEle[0]>>leadingEta(100,-3.0,3.0)"
//histName is the name of the histogram object, and is contained in plotArgs
//histTitle and xAxisTitle will be used with pHist
void makeAndSaveSingleTreeHisto(TTree * tree,TString plotArgs,TString histName,TString histTitle,TString xAxisTitle,TString canvName,TCut filters,TString outputFile, Bool_t isPlottingEnergy, Bool_t isPlottingInverseEnergy){
	TCanvas * canv = new TCanvas(canvName,canvName,500,500);
	canv->cd();
	tree->Draw(plotArgs,filters);
	TH1F * pHist = (TH1F*) gROOT->FindObject(histName);
	pHist->SetTitle(histTitle);
	pHist->GetXaxis()->SetTitle(xAxisTitle);
	if(histTitle.Contains("Iso") ){
		canv->SetLogy(1);
		pHist->SetMinimum(1);
	}
	gStyle->SetOptStat("emriou");
	
	//every histo should have at least three bins.
	//if a histo is created with numBins = 1, then bin #1 is the bin which is plotted 
	char temp[130];
	if(isPlottingInverseEnergy && pHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f / GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if(isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f / GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( isPlottingInverseEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.001){
		sprintf(temp,"Events / %.4f / GeV", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) <= 0.01){
		sprintf(temp,"Events / %.3f ", pHist->GetXaxis()->GetBinWidth(1));
	}
	if( !isPlottingInverseEnergy && !isPlottingEnergy && pHist->GetXaxis()->GetBinWidth(1) > 0.01){
		sprintf(temp,"Events / %.2f ", pHist->GetXaxis()->GetBinWidth(1));
	}
	pHist->GetYaxis()->SetTitle(temp);
	pHist->Draw();
	canv->SaveAs(outputFile,"recreate");

}//end makeAndSaveSingleTreeHisto()


///use this fxn to grab a histogram made in a TTree->Draw() call, draw this histo
///with a thicker line, and save it to a .png file 
void saveSingleHisto(TString canvName,TString histName,TString histTitle,TString outFile){
	Bool_t isPlottingGeV = false;
	
	if(histTitle.Contains("pt") || histTitle.Contains("Mass")) isPlottingGeV = true;

	gStyle->SetOptStat("emrou");
	TCanvas * canv = new TCanvas(canvName,canvName,600,600);
	canv->cd();
	TH1F * hist = (TH1F*) gROOT->FindObject(histName);
	hist->SetTitle(histTitle);
	hist->SetLineColor(1);
	hist->SetLineWidth(3);
	resetXaxisLimits(hist);
	resetNumXaxisBins(hist);
	
	//char temp[130];
	//if(isPlottingGeV) sprintf(temp,"Events / %.2f GeV",hist->GetXaxis()->GetBinWidth(1));
	//else sprintf(temp,"Events / %.2f ",hist->GetXaxis()->GetBinWidth(1));
	//hist->GetYaxis()->SetTitle(temp);

	hist->Draw();
	canv->Update();
	canv->SaveAs(outFile,"recreate");
	canv->Close();
	canv->Delete();
}///end saveSingleHisto()

///use this fxn to take one tuple and dump all branches into plots with unique names
void SaveTreePlots(TChain * chain, TString outputFileName){
	TObjArray * branches = chain->GetListOfBranches();
	TIter brItr(branches);
	for(TBranch * aBranch=(TBranch*) brItr.Next(); aBranch!=NULL; aBranch=(TBranch*) brItr.Next()){
		TString brName = aBranch->GetName();
#ifdef DEBUG
		std::cout<<"brName = \t"<< brName <<std::endl;
#endif
		string tempBrName(brName);
		if(tempBrName.find("Weight") != string::npos || brName.CompareTo("evtNumber")==0 || brName.CompareTo("runNumber")==0) continue;	///<ignore branches with miscellaneous information
		if(brName.CompareTo("etaEle")==0 || brName.CompareTo("phiEle")==0 || brName.CompareTo("ptEle")==0 || 
				brName.CompareTo("ptJet")==0 || brName.CompareTo("etaJet")==0 || brName.CompareTo("phiJet")==0 ){
			///NOTE the max value of 2 is set knowing that the max number of array elements per entry is 2
			///if the tree structure changes to allow arrays with more than 2 entries, this loop upper limit value
			///will have to change!
			for(unsigned int i=0; i<2; i++){
				TString num = std::to_string(i);
				TString updatedBrName = brName+"["+num+"]"+">>"+brName+num+"Hist";
				chain->Draw(updatedBrName);

				///now updatedBrName has the array name and element of interest, either [0] or [1]
				saveSingleHisto(brName+num,brName+num+"Hist",brName+"["+num+"]",outputFileName+"_"+brName+"["+num+"]"+".png");
				saveSingleHisto(brName+num,brName+num+"Hist",brName+"["+num+"]",outputFileName+"_"+brName+"["+num+"]"+".pdf");
			}///end for loop

		}///end filter to do something different if the branch name corresponds to an array branch
		
		else{
			TString updatedBrName = brName+">>"+brName+"Hist";
			TString histName = brName+"Hist";
#ifdef DEBUG
			std::cout<<"updatedBrName = \t"<< updatedBrName <<std::endl;
#endif

			chain->Draw(updatedBrName);
			
			TString outFilePath = outputFileName+"_"+brName+".png";
			TString outFilePathPdf = outputFileName+"_"+brName+".pdf";
#ifdef DEBUG
			std::cout<<"outFilePath = \t"<< outFilePath <<std::endl;
#endif
			saveSingleHisto(brName,histName,brName,outFilePath);
			saveSingleHisto(brName,histName,brName,outFilePathPdf);
		
		}

	}///end loop over branches in TChain pointer named chain

}///end SaveTreePlots()


///use this macro to take a file with multiple directories, one tree per directory, and dump
///all of the plots in all trees into pdf files with unique names
void dumpTreePlots(){
	///chains made without pdgId matching
	/*
	TChain * noCuts = new TChain("genAnalyzerOne/genObjectsNoCuts","");
	noCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	TChain * ptEtaCuts = new TChain("genAnalyzerTwo/genObjectsWithPtEtaCuts","");
	ptEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	TChain * ptEtaDileptonMassCuts = new TChain("genAnalyzerThree/genObjectsWithPtEtaAndDileptonMassCuts","");
	ptEtaDileptonMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel.root");
	*/

#ifdef GENLVL
	///chains made with pdgId matching
	TChain * matchedNoCuts = new TChain("matchedGenAnalyzerOne/matchedGenObjectsNoCuts","");
	matchedNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaCuts = new TChain("matchedGenAnalyzerTwo/matchedGenObjectsWithPtEtaCuts","");
	matchedPtEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaDileptonMassCuts = new TChain("matchedGenAnalyzerThree/matchedGenObjectsWithPtEtaAndDileptonMassCuts","");
	matchedPtEtaDileptonMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaDileptonMassDrCuts = new TChain("matchedGenAnalyzerFour/matchedGenObjectsWithPtEtaDileptonMassAndDrCuts","");
	matchedPtEtaDileptonMassDrCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedGenAnalyzerFive/matchedGenObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	matchedPtEtaDileptonMassDrFourObjMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");


	TString plotDir_noCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/noCuts";
	TString plotDir_withPtEtaCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/withPtEtaCuts";
	TString plotDir_withPtEtaDileptonMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/withPtEtaDileptonMassCuts";

	TString plotDir_matched_noCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_noCuts/noCuts";
	TString plotDir_matched_withPtEtaCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaCuts/withPtEtaCuts";
	TString plotDir_matched_withPtEtaDileptonMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaDileptonMassCuts/withPtEtaDileptonMassCuts";
	TString plotDir_matched_withPtEtaDileptonMassDrCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaDileptonMassDrCuts/withPtEtaDileptonMassDrCuts";
	TString plotDir_matched_withPtEtaDileptonMassDrFourObjMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/GEN/matched_ptEtaDileptonMassDrFourObjMassCuts/withPtEtaDileptonMassDrFourObjMassCuts";




	//SaveTreePlots(noCuts, plotDir_noCuts);
	//SaveTreePlots(ptEtaCuts, plotDir_withPtEtaCuts);
	//SaveTreePlots(ptEtaDileptonMassCuts, plotDir_withPtEtaDileptonMassCuts);
	
	SaveTreePlots(matchedNoCuts, plotDir_matched_noCuts);
	SaveTreePlots(matchedPtEtaCuts, plotDir_matched_withPtEtaCuts);
	SaveTreePlots(matchedPtEtaDileptonMassCuts, plotDir_matched_withPtEtaDileptonMassCuts);
	SaveTreePlots(matchedPtEtaDileptonMassDrCuts, plotDir_matched_withPtEtaDileptonMassDrCuts);
	SaveTreePlots(matchedPtEtaDileptonMassDrFourObjMassCuts, plotDir_matched_withPtEtaDileptonMassDrFourObjMassCuts);

#endif	



#ifdef GENMATCHEDRECO

	///chains made with deltaR matching btwn reco and gen
	/*
	TChain * matchedGenJetsToGenQuarksNoCuts = new TChain("matchGenJetsToGenQuarksNoCutsNewPath/matchedGenJetsNoCutsTree","");
	matchedGenJetsToGenQuarksNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	
	TChain * matchedRecoJetsToGenJetsNoCuts = new TChain("matchRecoJetsToGenJetsNoCutsNewPath/matchedRecoJetsNoCutsTree","");
	matchedRecoJetsToGenJetsNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	
	TChain * matchedRecoEleToLeadingGenEleNoCuts = new TChain("matchRecoEleToLeadingGenEleNoCutsNewPath/matchedLeadingRecoEleNoCutsTree","");
	matchedRecoEleToLeadingGenEleNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	TChain * matchedRecoEleToSubleadingGenEleNoCuts = new TChain("matchRecoEleToSubleadingGenEleNoCutsNewPath/matchedSubleadingRecoEleNoCutsTree","");
	matchedRecoEleToSubleadingGenEleNoCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	
	TString plotDir_reco_matched_noCuts_dR_genJetsToGenQuarks = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_genJetsToGenQuarks/noCuts_genJetsToGenQuarks";
	TString plotDir_reco_matched_noCuts_dR_recoElesToGenLeadingEles = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_recoElesToGenLeadingEles/noCuts_recoElesToGenLeadingEles";
	TString plotDir_reco_matched_noCuts_dR_recoElesToGenSubleadingEles = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_recoElesToGenSubleadingEles/noCuts_recoElesToGenSubleadingEles";
	TString plotDir_reco_matched_noCuts_dR_recoJetsToGenJets = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts_dR_recoJetsToGenJets/noCuts_recoJetsToGenJets";
	
	///make plots of dR btwn gen jets and gen quarks, reco jets and gen jets, and reco eles and gen eles
	SaveTreePlots(matchedGenJetsToGenQuarksNoCuts,plotDir_reco_matched_noCuts_dR_genJetsToGenQuarks);
	SaveTreePlots(matchedRecoJetsToGenJetsNoCuts,plotDir_reco_matched_noCuts_dR_recoJetsToGenJets);
	SaveTreePlots(matchedRecoEleToLeadingGenEleNoCuts,plotDir_reco_matched_noCuts_dR_recoElesToGenLeadingEles);
	SaveTreePlots(matchedRecoEleToSubleadingGenEleNoCuts,plotDir_reco_matched_noCuts_dR_recoElesToGenSubleadingEles);
	
	*/

	TString inputFile = "/eos/uscms/store/user/skalafut/WR/13TeV/RunIISpring15_MiniAODSignalSamples/analyzed_gen_matched_WRtoEEJJ_MWr800_MNu400.root";
	TChain * matchedRecoNoCuts = new TChain("matchedRecoAnalyzerOne/matchedRecoObjectsNoCuts","");
	matchedRecoNoCuts->Add(inputFile);
	
	TChain * matchedRecoPtEtaCuts = new TChain("matchedRecoAnalyzerTwo/matchedRecoObjectsWithPtEtaCuts","");
	matchedRecoPtEtaCuts->Add(inputFile);
	TChain * matchedRecoPtEtaDileptonMassCuts = new TChain("matchedRecoAnalyzerThree/matchedRecoObjectsWithPtEtaAndDileptonMassCuts","");
	matchedRecoPtEtaDileptonMassCuts->Add(inputFile);
	TChain * matchedRecoPtEtaDileptonMassDrCuts = new TChain("matchedRecoAnalyzerFour/matchedRecoObjectsWithPtEtaDileptonMassAndDrCuts","");
	matchedRecoPtEtaDileptonMassDrCuts->Add(inputFile);
	TChain * matchedRecoPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Add(inputFile);


	/*
	TString plotDir_reco_matched_noCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_noCuts/noCuts";
	TString plotDir_reco_matched_withPtEtaCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaCuts/withPtEtaCuts";
	TString plotDir_reco_matched_withPtEtaDileptonMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaDileptonMassCuts/withPtEtaDileptonMassCuts";
	TString plotDir_reco_matched_withPtEtaDileptonMassDrCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaDileptonMassDrCuts/withPtEtaDileptonMassDrCuts";
	TString plotDir_reco_matched_withPtEtaDileptonMassDrFourObjMassCuts = "/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/plots/RECO/matched_ptEtaDileptonMassDrFourObjMassCuts/withPtEtaDileptonMassDrFourObjMassCuts";
	*/

	///make plots of reco jets and eles matched to GEN objects after different levels of cuts
	///the matching is done before any cuts are applied at GEN or reco lvl
	SaveTreePlots(matchedRecoNoCuts, "noCuts");
	//SaveTreePlots(matchedRecoPtEtaCuts, plotDir_reco_matched_withPtEtaCuts);
	//SaveTreePlots(matchedRecoPtEtaDileptonMassCuts, plotDir_reco_matched_withPtEtaDileptonMassCuts);
	//SaveTreePlots(matchedRecoPtEtaDileptonMassDrCuts, plotDir_reco_matched_withPtEtaDileptonMassDrCuts);
	//SaveTreePlots(matchedRecoPtEtaDileptonMassDrFourObjMassCuts, plotDir_reco_matched_withPtEtaDileptonMassDrFourObjMassCuts);

#endif


}///end dumpTreePlots()

