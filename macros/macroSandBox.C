#include <TFile.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
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

//#define PtRatioProfiles
#define DEBUG
//#define RecoGenOverlays
#define StudyEffectOfMassPairs


/** the TString key in inputChainMap contains the histogram plotting argument to use with TChain->Draw("plottingArgs")
 * the histogram name does not need to be passed into the fxn as an argument, it will be pulled from the map key
 * aa
 *
 *
 */
void makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax){
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
	int colors[] = {1,2,4,9,12,30,40,45};
	vector<int> colorVect(colors,colors + sizeof(colors)/sizeof(int) );
	Int_t i=0;
	if(overlayHistoMap.size() > colorVect.size() ) cout<<"not enough unique colors in MultipleCurveOverlayHisto fxn!"<<endl;
	for(map<string,TH1F*>::const_iterator histIt=overlayHistoMap.begin(); histIt!=overlayHistoMap.end(); histIt++){
		(histIt->second)->SetLineColor(colorVect[i]);
		(histIt->second)->SetLineWidth(3);
		if(histIt==overlayHistoMap.begin()){
			string xLabel="";
			if((histIt->first).find("Mass") !=string::npos || (histIt->first).find("ptEle") !=string::npos || (histIt->first).find("ptJet") !=string::npos){
				xLabel += "GeV";
			}///end if(histo is plotting a variable with dimension of energy)
			(histIt->second)->GetXaxis()->SetTitle(xLabel.c_str());
			Double_t oldMax = (histIt->second)->GetBinContent((histIt->second)->GetMaximumBin());
			(histIt->second)->SetMaximum(2*oldMax);
		}///end filter to set histogram X axis label
		size_t mNuPosition = (histIt->first).find_first_of("mNu"), lastChevron = (histIt->first).find_last_of('>');
		string legEntryName;
		if(mNuPosition!=string::npos) legEntryName = (histIt->first).substr(lastChevron+1,mNuPosition+14-lastChevron);
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
 */
void calcAndPlotPtRatioProfile(TChain * chain, Int_t arrSize, Int_t arrIndex, TString profileName, TString profileTitle, Int_t nProfileBins, Float_t minX, Float_t maxX, Float_t minY, Float_t maxY, TString matchedGenBranchName, TString recoBranchName, TString otherMatchedGenBranchName,TString canvName, TString drawOptions, TString xAxisLabel, TString outputFile){
	TCanvas * g0 = new TCanvas(canvName,canvName,600,600);
	g0->cd();
	gStyle->SetOptStat("emriou");
	///declare two arrays which will store info from the TChain
	Float_t recoArr[arrSize], matchedGenArr[arrSize], otherMatchedGenArr[arrSize];
	///initialize the TProfile object
	TProfile * profilePtr = new TProfile(profileName, profileTitle, nProfileBins, minX, maxX, minY, maxY);
	///link the two arrays to the TChain branches which should be read out
	chain->SetBranchAddress(recoBranchName, &recoArr);
	chain->SetBranchAddress(matchedGenBranchName, &matchedGenArr);
	///do not use otherMatchedGenArr if otherMatchedGenBranchName is identical to matchedGenBranchName 
	if(otherMatchedGenBranchName.CompareTo(matchedGenBranchName)!=0) chain->SetBranchAddress(otherMatchedGenBranchName, &otherMatchedGenArr);
	Long64_t nEntries = chain->GetEntriesFast();
	for(Long64_t i=0;i<nEntries; i++){
		chain->GetEntry(i);
		if(otherMatchedGenBranchName.CompareTo(matchedGenBranchName)!=0){
			profilePtr->Fill(otherMatchedGenArr[arrIndex],recoArr[arrIndex]/matchedGenArr[arrIndex]);
		}
		else profilePtr->Fill(matchedGenArr[arrIndex],recoArr[arrIndex]/matchedGenArr[arrIndex]);
	}///end loop over entries in TChain
	
	///now draw the TProfile and save an image of the plot
	profilePtr->SetMarkerColor(kRed);
	profilePtr->GetXaxis()->SetTitle(xAxisLabel);
	profilePtr->Draw(drawOptions);
	g0->SaveAs(outputFile,"recreate");

}///end calcAndPlotPtRatioProfile()


///use this macro to develop and run new macros which don't have a central theme, other than being useful to the WR analysis
///the first use of this macro was to plot reco pT/gen pT for reco jets and leptons matched to GEN counterparts

void macroSandBox(){
	TChain * matchedRecoPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	TChain * copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts = new TChain("matchedRecoAnalyzerFive/matchedRecoObjectsWithPtEtaDileptonMassDrAndFourObjMassCuts","");
	copy_matchedRecoPtEtaDileptonMassDrFourObjMassCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");
	

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///plot reco pT/matched gen pT for both leptons and jets as a function of gen pT, eta, and phi
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef PtRatioProfiles	
	///leading lepton
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingLeptonPtRatioVsGenPt","leading electron P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1100.,0.,1.1,"ptMatchingEle","ptEle","ptMatchingEle","leadingLeptonGenPt","*H","P_{T,GEN} (GeV)","leadingLeptonPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingLeptonPtRatioVsGenEta","leading electron P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.0,3.0,0.,1.1,"ptMatchingEle","ptEle","etaMatchingEle","leadingLeptonGenEta","*H","#eta_{GEN}","leadingLeptonPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingLeptonPtRatioVsGenPhi","leading electron P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingEle","ptEle","phiMatchingEle","leadingLeptonGenPhi","*H","#phi_{GEN}","leadingLeptonPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png");

	///subsubleading lepton
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingLeptonPtRatioVsGenPt","subleading electron P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1000.,0.,1.1,"ptMatchingEle","ptEle","ptMatchingEle","subleadingLeptonGenPt","*H","P_{T,GEN} (GeV)","subleadingLeptonPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingLeptonPtRatioVsGenEta","subleading electron P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.2,3.2,0.,1.1,"ptMatchingEle","ptEle","etaMatchingEle","subleadingLeptonGenEta","*H","#eta_{GEN}","subleadingLeptonPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingLeptonPtRatioVsGenPhi","subleading electron P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingEle","ptEle","phiMatchingEle","subleadingLeptonGenPhi","*H","#phi_{GEN}","subleadingLeptonPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png");


	///leading jet
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingJetPtRatioVsGenPt","leading jet P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,1000.,0.,1.1,"ptMatchingJet","ptJet","ptMatchingJet","leadingJetGenPt","*H","P_{T,GEN} (GeV)","leadingJetPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingJetPtRatioVsGenEta","leading jet P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.0,3.0,0.,1.1,"ptMatchingJet","ptJet","etaMatchingJet","leadingJetGenEta","*H","#eta_{GEN}","leadingJetPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,0,"leadingJetPtRatioVsGenPhi","leading jet P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingJet","ptJet","phiMatchingJet","leadingJetGenPhi","*H","#phi_{GEN}","leadingJetPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png");


	///subleading jet
	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingJetPtRatioVsGenPt","subleading jet P_{T,RECO}/P_{T,GEN} vs P_{T,GEN}",20,0.,600.,0.,1.1,"ptMatchingJet","ptJet","ptMatchingJet","subleadingJetGenPt","*H","P_{T,GEN} (GeV)","subleadingJetPtRatio_vs_GenPt_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingJetPtRatioVsGenEta","subleading jet P_{T,RECO}/P_{T,GEN} vs #eta_{GEN}",20,-3.2,3.2,0.,1.1,"ptMatchingJet","ptJet","etaMatchingJet","subleadingJetGenEta","*H","#eta_{GEN}","subleadingJetPtRatio_vs_GenEta_matched_signal_eejj_after_all_reco_cuts.png");

	calcAndPlotPtRatioProfile(matchedRecoPtEtaDileptonMassDrFourObjMassCuts,2,1,"subleadingJetPtRatioVsGenPhi","subleading jet P_{T,RECO}/P_{T,GEN} vs #phi_{GEN}",20,-3.5,3.5,0.,1.1,"ptMatchingJet","ptJet","phiMatchingJet","subleadingJetGenPhi","*H","#phi_{GEN}","subleadingJetPtRatio_vs_GenPhi_matched_signal_eejj_after_all_reco_cuts.png");

#endif


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///plot pT, eta, and phi for matched reco objects (leptons and jets) on top of the gen objects which they are matched to
	///these plots show the distributions BEFORE any cuts are applied at reco level
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
	
	//makeAndSaveMultipleCurveOverlayHisto(map<string,TChain *> inputChainMap,TString canvName,Float_t legXmin,Float_t legYmin,Float_t legXmax,Float_t legYmax)
	string branchNames[] = {"ptEle[0]","ptEle[1]","etaEle[0]","etaEle[1]","ptJet[0]","ptJet[1]","etaJet[0]","etaJet[1]","dileptonMass","fourObjectMass","dR_leadingLeptonLeadingJet","dR_leadingLeptonSubleadingJet","dR_subleadingLeptonLeadingJet","dR_subleadingLeptonSubleadingJet","dR_leadingLeptonSubleadingLepton","dR_leadingJetSubleadingJet"};
	string link=">>";
	string histoEndings[] = {"_leadLeptonPt(50,0.,1600.)","_subleadLeptonPt(50,0.,1200.)","_leadLeptonEta(50,-3.0,3.0)","_subleadLeptonEta(50,-3.0,3.0)","_leadJetPt(50,0.,1000.)","_subleadJetPt(50,0.,1000.)","_leadJetEta(50,-3.0,3.0)","_subleadJetEta(50,-3.0,3.0)","_dileptonMass(50,0.,2600.)","_fourObjectMass(50,0.,3000.)","_dR_leadingLeptonLeadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingJet(50,0.,5.)","_dR_subleadingLeptonLeadingJet(50,0.,5.)","_dR_subleadingLeptonSubleadingJet(50,0.,5.)","_dR_leadingLeptonSubleadingLepton(50,0.,5.)","_dR_leadingJetSubleadingJet(50,0.,5.)"};
	vector<string> histoEndingVect(histoEndings,histoEndings + sizeof(histoEndings)/sizeof(string));
	string histoBeginnings[] = {"mWR2600_mNu2080","mWR2600_mNu1300","mWR2600_mNu520","mWR1400_mNu700"};
	map<string,TChain*> placeHolderMap;
	unsigned int maxI = histoEndingVect.size();
	for(unsigned int i=(maxI-1); i<maxI; i++){
		placeHolderMap[branchNames[i]+link+histoBeginnings[0]+histoEndings[i]] = MWR2600MNu2080_matchedGenNoCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[1]+histoEndings[i]] = MWR2600MNu1300_matchedGenNoCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[2]+histoEndings[i]] = MWR2600MNu520_matchedGenNoCuts;
		placeHolderMap[branchNames[i]+link+histoBeginnings[3]+histoEndings[i]] = MWR1400MNu700_matchedGenNoCuts;
		string cName = "o"+to_string(i);
		makeAndSaveMultipleCurveOverlayHisto(placeHolderMap,cName.c_str(),0.75,0.6,0.98,0.95);
		placeHolderMap.clear();
	}///end loop over branchNames

#endif

}///end macroSandBox()

