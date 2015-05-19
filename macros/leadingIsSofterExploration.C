#include <TEventList.h>
#include <TEntryList.h>
#include <TEntryListArray.h>
#include <TFile.h>
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
#include "dumpTreePlots.C" 

using namespace std;

//#define DEBUG

void leadingIsSofterExploration(){
	TChain * matchedGenPtEtaCuts = new TChain("matchedGenAnalyzerTwo/matchedGenObjectsWithPtEtaCuts","");
	matchedGenPtEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_genElectronChannel_using_miniAOD.root");
	TChain * matchedRecoPtEtaCuts = new TChain("matchedRecoAnalyzerTwo/matchedRecoObjectsWithPtEtaCuts","");
	matchedRecoPtEtaCuts->Add("/uscms/home/skalafut/WR/CMSSW_7_4_0_pre9/src/ExoAnalysis/cmsWR/analysis_recoElectronChannel_two_stage_matching_for_jets.root");

	///now use makeAndSaveSingleTreeHisto() from dumpTreePlots.C to make kinematics plots in gen and reco evts
	///where the higher pT electron comes from the decay of the heavy Nu
	TCut softLeadLept = "leadingIsHardest<1";
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"etaGenEle[1]>>subleadingGenEleeta","subleadingGenEleeta","#eta subleading Gen lepton","#eta","g1",softLeadLept,"subleadingGenEleeta_leadingLeptIsSoftest.png",false,false); 
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"etaGenEle[0]>>leadingGenEleeta","leadingGenEleeta","#eta leading Gen lepton","#eta","g2",softLeadLept,"leadingGenEleeta_leadingLeptIsSoftest.png",false,false); 
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"cosh(etaGenEle[1])*ptGenEle[1]>>subleadingGenEleenergy","subleadingGenEleenergy","Energy subleading Gen lepton","E (GeV)","g3",softLeadLept,"subleadingGenEleenergy_leadingLeptIsSoftest.png",true,false); 
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"cosh(etaGenEle[0])*ptGenEle[0]>>leadingGenEleenergy","leadingGenEleenergy","Energy leading Gen lepton","E (GeV)","g4",softLeadLept,"leadingGenEleenergy_leadingLeptIsSoftest.png",true,false);
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"(cosh(etaGenEle[1])*ptGenEle[1] - cosh(etaGenEle[0])*ptGenEle[0])>>GenEleenergyDifference","GenEleenergyDifference","#DeltaE(subleading - leading) Gen leptons","#DeltaE (GeV)","g5",softLeadLept,"GenEleenergyDifference_leadingLeptIsSoftest.png",true,false);
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"evtNumber>>GenEleevtNumber","GenEleevtNumber","Gen evt numbers no cuts","evt number","g6","","GenEleevtNumber_withPtEtaCuts.png",false,false);
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"ptGenEle[1]>>subleadingGenElept","subleadingGenElept","P_{T} subleading Gen lepton","PT (GeV)","g7",softLeadLept,"subleadingGenElept_leadingLeptIsSoftest.png",true,false); 
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"ptGenEle[0]>>leadingGenElept","leadingGenElept","P_{T} leading Gen lepton","PT (GeV)","g8",softLeadLept,"leadingGenElept_leadingLeptIsSoftest.png",true,false);
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"((ptGenEle[1]-ptGenEle[0])/ptGenEle[1])>>GenEleptDifferenceFraction","GenEleptDifferenceFraction","#DeltaP_{T}/highest P_{T} Gen leptons","#DeltaP_{T}/P_{T}","g9",softLeadLept,"GenEleptDifferenceFraction_leadingLeptIsSoftest.png",false,false);
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"((ptGenEle[1]*cosh(etaGenEle[1])-ptGenEle[0]*cosh(etaGenEle[0]))/ptGenEle[1]*cosh(etaGenEle[1]))>>GenEleenergyDifferenceFraction","GenEleenergyDifferenceFraction","#DeltaE/highest E Gen leptons","#DeltaE/E","g10",softLeadLept,"GenEleenergyDifferenceFraction_leadingLeptIsSoftest.png",false,false);
	//makeAndSaveSingleTreeHisto(matchedGenPtEtaCuts,"(ptGenEle[1]-ptGenEle[0])>>GenEleptDifference","GenEleptDifference","#DeltaP_{T} Gen leptons","#DeltaP_{T} (GeV)","g11",softLeadLept,"GenEleptDifference_leadingLeptIsSoftest.png",true,false);
	


	
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"etaGenEle[1]>>subleadingRecoEleeta","subleadingRecoEleeta","#eta subleading Reco lepton","#eta","r1",softLeadLept,"subleadingRecoEleeta_leadingLeptIsSoftest.png",false,false); 
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"etaGenEle[0]>>leadingRecoEleeta","leadingRecoEleeta","#eta leading Reco lepton","#eta","r2",softLeadLept,"leadingRecoEleeta_leadingLeptIsSoftest.png",false,false); 
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"cosh(etaGenEle[1])*ptGenEle[1]>>subleadingRecoEleenergy","subleadingRecoEleenergy","Energy subleading Reco lepton","E (GeV)","r3",softLeadLept,"subleadingRecoEleenergy_leadingLeptIsSoftest.png",true,false); 
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"cosh(etaGenEle[0])*ptGenEle[0]>>leadingRecoEleenergy","leadingRecoEleenergy","Energy leading Reco lepton","E (GeV)","r4",softLeadLept,"leadingRecoEleenergy_leadingLeptIsSoftest.png",true,false); 
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"(cosh(etaGenEle[1])*ptGenEle[1] - cosh(etaGenEle[0])*ptGenEle[0])>>RecoEleenergyDifference","RecoEleenergyDifference","#DeltaE(subleading - leading) Reco leptons","#DeltaE (GeV)","r5",softLeadLept,"RecoEleenergyDifference_leadingLeptIsSoftest.png",true,false);
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"evtNumber>>RecoEleevtNumber","RecoEleevtNumber","Reco evt numbers no cuts","evt number","r6","","RecoEleevtNumber_withPtEtaCuts.png",false,false);
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"ptGenEle[1]>>subleadingRecoElept","subleadingRecoElept","P_{T} subleading Reco lepton","PT (GeV)","r7",softLeadLept,"subleadingRecoElept_leadingLeptIsSoftest.png",true,false); 
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"ptGenEle[0]>>leadingRecoElept","leadingRecoElept","P_{T} leading Reco lepton","PT (GeV)","r8",softLeadLept,"leadingRecoElept_leadingLeptIsSoftest.png",true,false);
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"((ptGenEle[1]-ptGenEle[0])/ptGenEle[1])>>RecoEleptDifferenceFraction","RecoEleptDifferenceFraction","#DeltaP_{T}/highest P_{T} Reco leptons","#DeltaP_{T}/P_{T}","r9",softLeadLept,"RecoEleptDifferenceFraction_leadingLeptIsSoftest.png",false,false);
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"((ptGenEle[1]*cosh(etaGenEle[1])-ptGenEle[0]*cosh(etaGenEle[0]))/ptGenEle[1]*cosh(etaGenEle[1]))>>RecoEleenergyDifferenceFraction","RecoEleenergyDifferenceFraction","#DeltaE/highest E Reco leptons","#DeltaE/E","r10",softLeadLept,"RecoEleenergyDifferenceFraction_leadingLeptIsSoftest.png",false,false);
	//makeAndSaveSingleTreeHisto(matchedRecoPtEtaCuts,"(ptGenEle[1]-ptGenEle[0])>>RecoEleptDifference","RecoEleptDifference","#DeltaP_{T} Reco leptons","#DeltaP_{T} (GeV)","r11",softLeadLept,"RecoEleptDifference_leadingLeptIsSoftest.png",true,false);
	






}///end leadingIsSofterExploration()

