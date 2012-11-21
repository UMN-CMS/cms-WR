// -*- C++ -*-
//
// Package:    HeavyNu
// Class:      HeavyNuTop
// 
/**\class HeavyNuTop HeavyNuTop.cc HeavyNu/AnalyzerModules/src/HeavyNuTop.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: HeavyNuTop.cc,v 1.36 2012/10/26 23:25:53 pastika Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID 
// Order valid for 38X only.  Can be moved after Frameworkfwd.h in 39X
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
// Needed for 39X
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

// #include "RecoEgamma/EgammaElectronProducers/plugins/CorrectedGsfElectronProducer.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TVector3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
// #include "HeavyNu/AnalysisModules/src/HeavyNu_NNIF.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTopHist.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTree.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//////////////////////////////////////////////////////////////////

template <class T>
inline std::string int2str(T i)
{
    std::ostringstream ss;
    ss << i;
    return ss.str();
}

//============================================================

class comparePt
{
public:

    template <class T> bool operator() (const T& a, const T& b)
    {
        return a.pt() > b.pt() ;
    }
} ;

static std::string btagName;
static double      minBtagDiscVal; // for discriminating B-tagged jets.

class HeavyNuTop : public edm::EDFilter
{
public:
    explicit HeavyNuTop(const edm::ParameterSet&);
    ~HeavyNuTop();


private:

    virtual void beginJob          ();
    virtual bool filter            ( edm::Event&, const edm::EventSetup& );
    virtual void endJob            ();

    virtual TH1 *bookRunHisto      ( uint32_t runNumber );

    void fill(pat::MuonCollection muons,
              pat::ElectronCollection electrons,
              pat::JetCollection  jets,
              bool isMC,
              double wgt,
              HeavyNuHistSet* hnhs);

    // double GetCorrectedPt ( const pat::Electron& e );

    edm::InputTag muonTag_;
    edm::InputTag jetTag_;
    edm::InputTag metTag_;
    edm::InputTag elecTag_;

    edm::InputTag elecRhoTag_ ;
    double elecRho_ ;

    double applyMESfactor_;             // for Muon Energy Scale studies
    int    applyTrigEffsign_;           // for Trigger Efficiency studies
    bool   applyMuIDCorrections_ ;
    int applyMuIDEffsign_ ;
    bool correctEscale_ ; 
    // bool   applyEleScaleCorrections_ ; 
    // bool   applyEleIDWeightFactor_ ; 
    bool   studyScaleFactorEvolution_;  // for Top, Z+jets scale factors (by Mu1 pT) studies

    bool useDirtyElectrons_, useDirtyMuons_ ; 
    int  fakeFlavorToWeight_ ; 

    int    pileupEra_;

    int    heepVersion_ ;

    double EBscalefactor_, EEscalefactor_ ;
    double ebIDwgt_, eeIDwgt_ ;

    JetCorrectionUncertainty *jecuObj_;
    bool dolog_;
    bool firstEvent_;
    bool isPFJets_;

    // HeavyNu_NNIF *nnif_;
    HeavyNuTrigger *trig_;
    HeavyNuID *muid_ ;

    edm::LumiReWeighting MCweightByVertex_;

    std::map<uint32_t, TH1 *> m_runHistos_;

    // ----------member data ---------------------------

    bool init_;

    bool addSlopeTree_;
    HeavyNuTree *hnuTree_;

    // gf set of histo for all Z definitions in a stack

    struct HistStruct
    {
        TH1 *nelec, *njet, *nmet, *nmuAll, *nmuLoose, *nmuTight;
        TH1 *muPt, *muEta, *muPhi, *looseMuPt, *tightMuPt ;

        TH1* cutlevel;

        // Muon quality histos as a function of Pt
        TH2 *muNvalidHitsVsPt, *mudBvsPt, *muNormChi2vsPt, *muQualVsPt;
        TH2 *muNmatchesVsPt, *muNvalidMuonHitsVsPt, *muNvalidPixelHitsVsPt;
        TH2 *muTrckIsoVsPt, *muHcalIsoVsPt, *muEcalIsoVsPt, *muCaloIsoVsPt;

        TH1 *jetPt, *jetEta, *jetPhi, *jetID, *jecUncHi, *jecUncLo, *met ;
        TH2 *jetPtvsNum;
        TProfile2D *jecUncHiVsEtaPt, *jecUncLoVsEtaPt;

        TH1 *trkIsoStudy;

        TFileDirectory *rundir;
        HeavyNuHistSet* noCuts;
        // HeavyNuHistSet* LLptCuts;
        // HeavyNuHistSet* MuTightCuts;
        HeavyNuHistSet* LLJJptCuts;
        HeavyNuHistSet* TrigMatches;
        HeavyNuHistSet* VertexCuts;
        HeavyNuHistSet* Mu1HighPtCut;
        HeavyNuHistSet* Mu1HighPtCut_1bjet;
        HeavyNuHistSet* Mu1HighPtCut_2bjet;
        HeavyNuHistSet* Mu1Pt30GeVCut;
        HeavyNuHistSet* Mu1Pt40GeVCut;
        HeavyNuHistSet* Mu1Pt50GeVCut;
        HeavyNuHistSet* Mu1Pt60GeVCut;
        HeavyNuHistSet* Mu1Pt80GeVCut;
        HeavyNuHistSet* Mu1Pt100GeVCut;
        HeavyNuHistSet* Mu1HighPtCutVtxEq1;
        HeavyNuHistSet* Mu1HighPtCutVtx2to5;
        HeavyNuHistSet* Mu1HighPtCutVtxGt5;
        HeavyNuHistSet* diLmassCut;
        //HeavyNuHistSet* diLmassCutEChannel;
        HeavyNuHistSet* mWRmassCut;

    } hists;

    struct CutsStruct
    {
        double minimum_lep1_pt;
        double minimum_lep2_pt;
        double minimum_jet_pt;
        double maximum_mu_abseta;
        double maximum_elec_abseta;
        double maximum_jet_abseta;
        double minimum_mumu_mass;
        double minimum_mWR_mass;
        double minimum_lep_jet_dR;
        double muon_trackiso_limit;
        double maxVertexZsep;
        double maxJetVZsepCM;
    } cuts;

} ;


const int muonQualityFlags = 4 ;
const std::string muonQuality[] = {
    "All", "AllGlobalMuons", "AllStandAloneMuons", "AllTrackerMuons"
};

inline void labelMuonQualAxis(TAxis *ax)
{
    for (int i = 0; i < muonQualityFlags; i++)
    {
        ax->SetBinLabel(i + 1, muonQuality[i].c_str()) ;
        ax->SetBinLabel(i + 1, muonQuality[i].c_str()) ;
    }
}

void HeavyNuTop::fill(pat::MuonCollection muons,
                      pat::ElectronCollection electrons,
                      pat::JetCollection  jets,
                      bool isMC,
                      double wgt,
                      HeavyNuHistSet* hnhs)
{
    HeavyNuEvent hne(HeavyNuEvent::TOP);

    std::sort(muons.begin(), muons.end(), comparePt()) ;
    std::sort(electrons.begin(), electrons.end(), comparePt()) ;
    std::sort(jets.begin(), jets.end(), comparePt()) ;

    if(jets.size() >= 2)
    {
        hne.j1 = jets.at(0);
        hne.j2 = jets.at(1);
        hne.nJets = 2;
    }
    if(muons.size() >= 1)
    {
        hne.mu1 = muons.at(0);
        hne.nLeptons++;
    }
    if(electrons.size() >= 1)
    {
        hne.e1 = electrons.at(0);
        hne.nLeptons++;
    }

    hne.eventWgt = wgt;
    hne.isMC = isMC;

    //hne.regularize();
    //hne.calculateLL();
    hne.calculate();

    hnhs->fill(hne);

}// end of fill()

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HeavyNuTop::HeavyNuTop(const edm::ParameterSet& iConfig)
{
    // ==================== Get all parameters ====================
    //
    dolog_ = iConfig.getParameter<bool>("DoLog");

    muonTag_ = iConfig.getParameter< edm::InputTag > ( "muonTag" );
    jetTag_  = iConfig.getParameter< edm::InputTag > ( "jetTag"  );
    metTag_  = iConfig.getParameter< edm::InputTag > ( "metTag"  );
    elecTag_ = iConfig.getParameter< edm::InputTag > ( "electronTag" );

    elecRhoTag_ = iConfig.getParameter< edm::InputTag > ("electronRho");    

    btagName = iConfig.getParameter<std::string > ("BtagName");
    minBtagDiscVal = iConfig.getParameter<double>("minBtagDiscr");

    cuts.minimum_lep1_pt      = iConfig.getParameter<double>("minLep1pt");
    cuts.minimum_lep2_pt      = iConfig.getParameter<double>("minLep2pt");
    cuts.minimum_jet_pt       = iConfig.getParameter<double>("minJetPt");
    cuts.maximum_mu_abseta    = iConfig.getParameter<double>("maxMuAbsEta");
    cuts.maximum_elec_abseta  = iConfig.getParameter<double>("maxElecAbsEta");
    cuts.maximum_jet_abseta   = iConfig.getParameter<double>("maxJetAbsEta");
    cuts.minimum_mumu_mass    = iConfig.getParameter<double>("minMuMuMass");
    cuts.minimum_mWR_mass     = iConfig.getParameter<double>("min4objMass");
    cuts.minimum_lep_jet_dR   = iConfig.getParameter<double>("minLeptonJetdR");
    cuts.muon_trackiso_limit  = iConfig.getParameter<double>("muonTrackRelIsoLimit");
    cuts.maxVertexZsep        = iConfig.getParameter<double>("maxVertexZsepCM");
    cuts.maxJetVZsepCM        = iConfig.getParameter<double>("maxVertexJetVZsepCM");

    pileupEra_ = iConfig.getParameter<int>("pileupEra");

    // Default HEEP version is 4.0 (2012 selection)
    heepVersion_ = iConfig.getParameter<int>("heepVersion");
    if ( heepVersion_ < 30 || heepVersion_ > 41 ) 
    {
        std::cout << "!!!!!!!!INVALID HEEP VERSION: " << heepVersion_ << " (setting HEEP 4.1)!!!!!!!!" << std::endl;
        heepVersion_ = 41;
    }

    useDirtyElectrons_  = iConfig.getParameter<bool>("useDirtyElectrons") ; 
    useDirtyMuons_      = iConfig.getParameter<bool>("useDirtyMuons") ; 
    fakeFlavorToWeight_ = iConfig.getParameter<int>("fakeWeightFlavor") ; 

    EBscalefactor_ = iConfig.getParameter<double>("EBscalefactor") ;
    EEscalefactor_ = iConfig.getParameter<double>("EEscalefactor") ;
    ebIDwgt_ = iConfig.getParameter<double>("EBidWgt") ;
    eeIDwgt_ = iConfig.getParameter<double>("EEidWgt") ;

    applyMESfactor_ = 1.0 ; // Hardcoded for Top studies
    correctEscale_  = iConfig.getParameter<bool>("correctEscale") ; 
    applyTrigEffsign_ = iConfig.getParameter<int>("applyTrigEffsign");
    if(applyTrigEffsign_) applyTrigEffsign_ /= abs(applyTrigEffsign_); // ensure -1,0,+1

    isPFJets_ = iConfig.getParameter<bool>("isPFJets");
    addSlopeTree_ = iConfig.getUntrackedParameter<bool>("addSlopeTree");

    // applyEleScaleCorrections_ = iConfig.getParameter<bool>("applyEleEScale") ; ;
    // applyEleIDWeightFactor_   = iConfig.getParameter<bool>("applyEleIDweight") ;
    applyMuIDCorrections_ = iConfig.getParameter<bool>("applyMuIDEffcorr");
    applyMuIDEffsign_     = iConfig.getParameter<int>("applyMuIDEffsign");
    if (applyMuIDEffsign_) applyMuIDEffsign_ /= abs(applyMuIDEffsign_); // ensure -1,0,+1


    studyScaleFactorEvolution_ = iConfig.getParameter<bool>("studyScaleFactor");

    // ==================== Init other members ====================
    //

    // nnif_ = new HeavyNu_NNIF(iConfig);
    trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet > ("trigMatchPset"));
    muid_ = new HeavyNuID(iConfig.getParameter<edm::ParameterSet > ("muIDPset"));

    // ==================== Book the histos ====================
    //
    edm::Service<TFileService> fs;
    hists.nelec    = fs->make<TH1D > ("nelec",     "N(e^{#pm})", 10, -0.5, 9.5);

    //this must be after at least 1 call of fs->make<...>(...) in order for the directory to have been created.
    if(addSlopeTree_) hnuTree_ = new HeavyNuTree(*fs->getBareDirectory(), true);

    hists.nmuAll   = fs->make<TH1D > ("nmuAll",    "N(#mu^{#pm})", 10, -0.5, 9.5);
    hists.nmuLoose = fs->make<TH1D > ("nmuLoose",  "N(#mu^{#pm}) passes Loose", 10, -0.5, 9.5);
    hists.nmuTight = fs->make<TH1D > ("nmuTight",  "N(#mu^{#pm}) passes Tight", 10, -0.5, 9.5);
    hists.njet     = fs->make<TH1D > ("njet",      "N(Jet)", 50, -0.5, 49.5);
    hists.nmet     = fs->make<TH1D > ("nmet",      "N(MET)", 50, -0.5, 49.5);
    hists.met      = fs->make<TH1D > ("met",       "MET distribution", 100, 0, 2000);
    hists.muPt     = fs->make<TH1D > ("muPt",      "#mu p_{T} distribution", 100, 0, 2000);
    hists.muEta    = fs->make<TH1D > ("muEta",     "#mu #eta distribution", 50, -2.5, 2.5);
    hists.muPhi    = fs->make<TH1D > ("muPhi",     "#mu #phi distribution", 60, -3.14159, 3.14159);
    hists.jetPt    = fs->make<TH1D > ("jetPt",     "jet p_{T} distribution", 100, 0, 2000);
    hists.jetEta   = fs->make<TH1D > ("jetEta",    "jet #eta distribution", 50, -5, 5);
    hists.jetPhi   = fs->make<TH1D > ("jetPhi",    "jet #phi distribution", 60, -3.14159, 3.14159);
    hists.jetID    = fs->make<TH1I > ("jetID",     "Jet ID", 3, 0, 3);
    hists.jetPtvsNum = fs->make<TH2D > ("jetPtvsNum", "Jet P_{T} vs. Jet # ", 11, -0.5, 10.5, 200, 0., 2000.);

    hists.jetID->GetXaxis()->SetBinLabel(1, "Neither");
    hists.jetID->GetXaxis()->SetBinLabel(2, "PURE09 Loose");
    hists.jetID->GetXaxis()->SetBinLabel(3, "PURE09 Tight");

    hists.trkIsoStudy = fs->make<TH1D > ("trkIsoStudy", ";Tracker Relative Isolation", 100, 0., 1.);

    hists.cutlevel = fs->make<TH1D > ("cutlevel", "Cut Level", 11, -1.5, 9.5);
    hists.cutlevel->GetXaxis()->SetBinLabel(1, "Raw");
    hists.cutlevel->GetXaxis()->SetBinLabel(2, "No cuts");
    hists.cutlevel->GetXaxis()->SetBinLabel(3, "mmjj p_{T}");
    hists.cutlevel->GetXaxis()->SetBinLabel(4, "M1 (Trigger)");
    hists.cutlevel->GetXaxis()->SetBinLabel(5, "M2 (Vertex)");
    hists.cutlevel->GetXaxis()->SetBinLabel(6, "M3 (high p_{T})");
    hists.cutlevel->GetXaxis()->SetBinLabel(7, "M4 (M_{#mu#mu})");
    hists.cutlevel->GetXaxis()->SetBinLabel(8, "M5 (M(W_{R})");

    // Histos per cut:
    //
    hists.noCuts = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut0_none")),       "(no cuts)");
    hists.LLJJptCuts = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut1_LLJJpt")),     "(4 objects with ptcuts:1)");
    hists.TrigMatches = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:M1)");
    hists.VertexCuts = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut3_Vertex")),     "(vertex requirements:M2)");
    hists.Mu1HighPtCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut4_Mu1HighPt")),  "(Mu1 High pt cut:M3)");
    hists.Mu1HighPtCut_1bjet = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut4a_Mu1HighPt_1b")),  "(Mu1 High pt cut:M3a)");
    hists.Mu1HighPtCut_2bjet = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut4b_Mu1HighPt_2b")),  "(Mu1 High pt cut:M3b)");
    hists.diLmassCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut5_diLmass")),    "(mumu mass cut:M4)");
    hists.mWRmassCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("cut6_mWRmass")),    "(mumujj mass cut:M5)");

    if (studyScaleFactorEvolution_)
    {
        hists.Mu1Pt30GeVCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1Pt30GeV")),  "(Mu1 30 GeV pt cut)");
        hists.Mu1Pt40GeVCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1Pt40GeV")),  "(Mu1 40 GeV pt cut)");
        hists.Mu1Pt50GeVCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1Pt50GeV")),  "(Mu1 50 GeV pt cut)");
        hists.Mu1Pt60GeVCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1Pt60GeV")),  "(Mu1 60 GeV pt cut)");
        hists.Mu1Pt80GeVCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1Pt80GeV")),  "(Mu1 80 GeV pt cut)");
        hists.Mu1Pt100GeVCut = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1Pt100GeV")), "(Mu1 100 GeV pt cut)");

        hists.Mu1HighPtCutVtxEq1 = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1HighPtVtxEq1")), "(Mu1 60 GeV pt cut, 1 vtx)");
        hists.Mu1HighPtCutVtx2to5 = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1HighPtVtx2to5")), "(Mu1 60 GeV pt cut, 2-5 vtx)");
        hists.Mu1HighPtCutVtxGt5 = new HeavyNuTopHist( new TFileDirectory(fs->mkdir("Mu1HighPtVtxGt5")), "(Mu1 60 GeV pt cut, 6+ vtx)");
    }

    hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

    init_ = false;

    MCweightByVertex_ = edm::LumiReWeighting(hnu::generate_flat10_mc(pileupEra_),
                                             hnu::get_standard_pileup_data(pileupEra_));

    // For the record...
    std::cout << "Configurable cut values applied:" << std::endl;
    std::cout << "muonTag           = " << muonTag_                   << std::endl;
    std::cout << "jetTag            = " << jetTag_                    << std::endl;
    std::cout << "metTag            = " << metTag_                    << std::endl;
    std::cout << "electronTag       = " << elecTag_                   << std::endl;
    std::cout << "btagName          = " << btagName                   << std::endl;
    std::cout << "minBtagDiscr      = " << minBtagDiscVal             << std::endl;
    std::cout << "minLep1pt          = " << cuts.minimum_lep1_pt        << " GeV" << std::endl;
    std::cout << "minLep2pt          = " << cuts.minimum_lep2_pt        << " GeV" << std::endl;
    std::cout << "minJetPt          = " << cuts.minimum_jet_pt        << " GeV" << std::endl;
    std::cout << "maxMuAbsEta       = " << cuts.maximum_mu_abseta     << std::endl;
    std::cout << "maxElecAbsEta       = " << cuts.maximum_elec_abseta     << std::endl;
    std::cout << "maxJetAbsEta      = " << cuts.maximum_jet_abseta    << std::endl;
    std::cout << "minLeptonJetdR      = " << cuts.minimum_lep_jet_dR   << std::endl;
    std::cout << "muonTrackRelIso   = " << cuts.muon_trackiso_limit   << std::endl;
    std::cout << "minMuMuMass       = " << cuts.minimum_mumu_mass     << " GeV" << std::endl;
    std::cout << "min4objMass       = " << cuts.minimum_mWR_mass      << " GeV" << std::endl;
    std::cout << "applyMuIDEffcorr  = " << applyMuIDCorrections_      << std::endl;
    std::cout << "applyMuIDEffsign  = " << applyMuIDEffsign_ << std::endl;
    // std::cout << "applyEleEScale    = " << applyEleScaleCorrections_  << std::endl ;
    std::cout << "EB scale factor   = " << EBscalefactor_             << std::endl ;
    std::cout << "EE scale factor   = " << EEscalefactor_             << std::endl ;
    // std::cout << "applyEleIDweight  = " << applyEleIDWeightFactor_    << std::endl ;
    std::cout << "EB weight         = " << ebIDwgt_                   << std::endl ;
    std::cout << "EE weight         = " << eeIDwgt_                   << std::endl ;

    std::cout << "HEEP version      = " << heepVersion_               << std::endl ;

    std::cout << "pileup era        = " << pileupEra_ << std::endl;
    std::cout << "studyScaleFactor  = " << studyScaleFactorEvolution_ << std::endl;
    std::cout << "addSlopeTree      = " << addSlopeTree_ << std::endl;
}

HeavyNuTop::~HeavyNuTop(){

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
 }

TH1* HeavyNuTop::bookRunHisto(uint32_t runNumber)
{
    std::string runstr = int2str<uint32_t > (runNumber);
    return hists.rundir->make <TH1I > (runstr.c_str(), runstr.c_str(), 1, 1, 2);
}

// ------------ method called to for each event  ------------

bool HeavyNuTop::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    HeavyNuEvent hnuEvent(HeavyNuEvent::TOP);
    if(addSlopeTree_) hnuTree_->clear();

    hnuEvent.isMC = !iEvent.isRealData();
    hnuEvent.pfJets = isPFJets_;

    if (iEvent.isRealData())
    {
        if( (applyMESfactor_ != 1.0) )
            throw cms::Exception( "Energy scale studies not allowed on data currently");

        uint32_t runn = iEvent.id().run();
        std::map<uint32_t, TH1 *>::const_iterator it = m_runHistos_.find(runn);
        TH1 *runh;
        if (it == m_runHistos_.end())
        {if(pileupEra_ < 20100) hnuEvent.eventWgt = 1.0;
            runh = bookRunHisto(runn);
            m_runHistos_[runn] = runh;
        }
        else
            runh = it->second;
        runh->Fill(1);
    }

    edm::Handle<reco::JPTJetCollection> jptJets;
    iEvent.getByLabel("JetPlusTrackZSPCorJetAntiKt5", jptJets); // Some day we should make this readable from the cfg file

    edm::Handle<pat::MuonCollection> pMuons ;
    iEvent.getByLabel(muonTag_, pMuons) ;

    edm::Handle<pat::ElectronCollection> pElecs ;
    iEvent.getByLabel(elecTag_, pElecs) ;

    edm::Handle<double> electronRhoHandle ; 
    iEvent.getByLabel(elecRhoTag_, electronRhoHandle) ; 
    elecRho_ = ((electronRhoHandle.isValid()) ? (*(electronRhoHandle.product())) : 0.) ; 

    edm::Handle<pat::JetCollection> pJets ;
    iEvent.getByLabel(jetTag_, pJets) ;

    edm::Handle<pat::METCollection> pMET ;
    iEvent.getByLabel(metTag_, pMET) ;

    //Shirpa reweighting info
    edm::Handle<GenEventInfoProduct> geneventinfo;
    iEvent.getByLabel("generator", geneventinfo);

    if (hnuEvent.isMC)
    {
        edm::Handle<std::vector<PileupSummaryInfo> > pPU;
        iEvent.getByLabel("addPileupInfo", pPU);
        if(pileupEra_ < 20100) hnuEvent.eventWgt = 1.0;
        else
        {
            std::pair<float, double> pileup = hnu::pileupReweighting(pPU, MCweightByVertex_) ;
            hnuEvent.n_pue     = int(pileup.first) ;
            hnuEvent.eventWgt *= pileup.second ;
        }
        //Shirpa reweighting
        hnuEvent.eventWgt *= geneventinfo->weight();
    }
    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    hnuEvent.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle) ;

    if ( !pElecs.isValid() ||
        !pMuons.isValid() ||
        !pJets.isValid() )// ||
        //!(pMET.isValid() && pMET->size() > 0))
    {
        std::cout << "Exiting as valid PAT objects not found" << std::endl ;
        return false;
    }

    if (firstEvent_)
    {
        firstEvent_ = false;
        // Some sanity checks inserted
        // If running on Monte Carlo, expect
        //   - pileup configuration to make sense
        //   - to apply nominal trigger and ID corrections
        //   - to apply proper weighting to events for electron reco/ID
        // If running on data
        //   - no corrections are applied
        if ( studyScaleFactorEvolution_ ) std::cout << "Histograms for top scale factor cross checks will be created" << std::endl ;
        if ( hnuEvent.isMC )
        {
            if ( applyMuIDCorrections_ )
            {
                if (applyMuIDEffsign_) std::cout << "Studies will be used to estimate Mu ID uncertainty: " << applyMuIDEffsign_
                        << std::endl;
                else                   std::cout << "Nominal Mu ID corrections applied" << std::endl ;
            }
            else
            {
                std::cout << "WARNING: You have disabled Mu ID corrections.  Are you sure?" << std::endl ;
            }
            if ( ebIDwgt_ != 1.0 &&
                eeIDwgt_ != 1.0 )            std::cout << "Events will be weighted to account for data/MC reconstruction and ID differences" << std::endl ;
            else                              std::cout << "WARNING: You are not applying EB/EE weighting to MC.  Are you sure?" << std::endl ;

            int pileupYear = pileupEra_ ;
            int idYear     = muid_->idEra() ;

            bool allErasMatch = ( pileupYear == idYear ) ;
            if ( !allErasMatch )
            {
                std::cout << "WARNING: You do not appear to have consistent corrections applied!" << std::endl ;
                std::cout << "         pileup year is " << pileupEra_ << ", year for mu ID is " << idYear
                        << std::endl ;
            }
            else
            {
                std::cout << "Looking at corrections, I assume you are running with the " << pileupYear << " year settings" << std::endl ;
            }
        }
        else
        {
            if ( ebIDwgt_ != 1.0 && eeIDwgt_ != 1.0 )
                std::cout << "WARNING: You are applying EB/EE weighting to DATA.  Are you sure?" << std::endl ;
        }
    }

    hists.nelec ->Fill(pElecs->size()) ;
    hists.nmuAll->Fill(pMuons->size()) ;
    hists.njet  ->Fill(pJets->size()) ;
    //hists.nmet  ->Fill(pMET->size()) ;

    //if (pMET->size())
    //    hists.met->Fill(pMET->at(0).pt());
    //else
    //    hists.met->Fill(0);

    hists.cutlevel->Fill(-1.0, hnuEvent.eventWgt);

    // Basic selection requirements: Require at least two muons, two jets
    if ( pMuons->size() >= 1 && pElecs->size() >= 1 && pJets->size() >= 2 )
    {
        hists.cutlevel->Fill(0.0, hnuEvent.eventWgt);
        fill( *pMuons, *pElecs, *pJets, hnuEvent.isMC, hnuEvent.eventWgt, hists.noCuts);
    }
    else return false;

    if ( dolog_ ) std::cout << "Found an event with " << pMuons->size() << " muons and "
        << pElecs->size() << " electrons, and " << pJets->size() << " jets" << std::endl ;

    // Look for valid jets and put them in the event
    std::vector< std::pair<pat::Jet, float> > jetCands =
            hnu::getJetList(pJets, jecuObj_, cuts.minimum_jet_pt, cuts.maximum_jet_abseta, 0, 3, hnuEvent.isMC, 0) ;
    for (unsigned int i = 0; i < jetCands.size(); i++)
    {
        if ( hnuEvent.nJets == 2 ) break ;
        pat::Jet iJ = jetCands.at(i).first ;
        // double dRej = deltaR(iJ.eta(), iJ.phi(), hnuEvent.e1.eta(), hnuEvent.e1.phi()) ;
        // if (dRej > cuts.minimum_lep_jet_dR)
        // {
	hnuEvent.nJets++ ;
	if ( hnuEvent.nJets == 1 )
	  {
	    hnuEvent.j1 = iJ ;
	    hnuEvent.j1scale = 1.0 ; // No jet corrections are applied for Top!
	  }
	else if ( hnuEvent.nJets == 2 )
	  {
	    hnuEvent.j2 = iJ ;
	    hnuEvent.j2scale = 1.0 ;
	  }
	else
	  std::cout << "WARNING: Expected empty jet position" << std::endl ;
        // }
    }
    if ( hnuEvent.nJets < 2 ) return false ;
    hnuEvent.tjV1 = hnu::caloJetVertex(hnuEvent.j1, *jptJets);
    hnuEvent.tjV2 = hnu::caloJetVertex(hnuEvent.j2, *jptJets);

    // Look for valid electrons and put them in the event
    int eidVersion = abs(heepVersion_) ; 
    if ( useDirtyElectrons_ ) eidVersion = -1 * heepVersion_ ; 
    std::vector< std::pair<pat::Electron, float> > eCands =
      hnu::getElectronList(pElecs, cuts.maximum_elec_abseta,
			   cuts.minimum_lep2_pt, cuts.minimum_lep2_pt, eidVersion, pvHandle, elecRho_) ;

    for (unsigned int i=0; i<eCands.size(); i++) { 
      pat::Electron iE = eCands.at(i).first;
      double dRj1 = deltaR(iE.eta(), iE.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
      double dRj2 = deltaR(iE.eta(), iE.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
      if(dRj1 > cuts.minimum_lep_jet_dR && dRj2 > cuts.minimum_lep_jet_dR)
	{
	  hnuEvent.nLeptons++;
	  hnuEvent.e1 = iE;
	  hnuEvent.ElecScale = eCands.at(i).second ; 
	  break ; 
	}
    }
    if ( hnuEvent.nLeptons < 1 ) return false ;

    // Finally, look for valid muons and put them in the event
    std::vector<pat::Muon> muCands = hnu::getMuonList(pMuons, pvHandle, (int(muid_->idEra()/10)), 
						      cuts.minimum_lep2_pt, cuts.maximum_mu_abseta, applyMESfactor_) ;
    for (unsigned int i = 0; i < muCands.size(); i++)
    {
        pat::Muon iM = muCands.at(i) ;

	if ( useDirtyMuons_ ) { 
	  double caloIso = iM.ecalIso() + iM.hcalIso() ;
	  bool   isDirty = ( (iM.pt() < 100.) ? (caloIso > 10.) : ((caloIso/iM.pt()) > 0.10) ) ; 
	  if ( isDirty && deltaR(iM.eta(), iM.phi(), hnuEvent.e1.eta(), hnuEvent.e1.phi()) < 0.3 ) isDirty = false ; 

	  if ( isDirty ) { // Calo non-isolated muon 
	    double dRj1 = deltaR(iM.eta(), iM.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi()) ;
	    double dRj2 = deltaR(iM.eta(), iM.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi()) ;
            if (dRj1 > cuts.minimum_lep_jet_dR && dRj2 > cuts.minimum_lep_jet_dR)
            {
                hnuEvent.nLeptons++ ;
                hnuEvent.mu1 = iM ;
                break;
            }
	  }
	} else if ( hnu::muIsolation(iM) < cuts.muon_trackiso_limit ) {
            double dRj1 = deltaR(iM.eta(), iM.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi()) ;
            double dRj2 = deltaR(iM.eta(), iM.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi()) ;
            if (dRj1 > cuts.minimum_lep_jet_dR && dRj2 > cuts.minimum_lep_jet_dR)
            {
                hnuEvent.nLeptons++ ;
                hnuEvent.mu1 = iM ;
                break;
            }
        }
    }
    if ( hnuEvent.nLeptons < 2 ) return false ; // No muons

    // Determine weight factors for dirty-e + dirty-mu jj
    // NOTE: For emjj closure tests, only relative rates matter (absolute scale less important)
    if ( useDirtyElectrons_ && useDirtyMuons_ ) { 
      double fakeWeight = 1.0 ; 
      // Note: "fake" electrons must pass minimal cuts but FAIL full selection
      if      (fakeFlavorToWeight_ == 11) 
	fakeWeight = (hnu::fakeProbability(hnuEvent.e1)/(1.0-hnu::fakeProbability(hnuEvent.e1))) ; 
      else if (fakeFlavorToWeight_ == 13) fakeWeight = hnu::fakeProbability( hnuEvent.mu1 ) ; 
      // Special case: Consider trying to estimate QCD background to e/mu jj:
      else if (fakeFlavorToWeight_ == 1113 || fakeFlavorToWeight_ == 1311) 
	fakeWeight = (hnu::fakeProbability(hnuEvent.e1)/(1.0-hnu::fakeProbability(hnuEvent.e1))) * 
	  hnu::fakeProbability( hnuEvent.mu1 ) ; 
      hnuEvent.eventWgt *= fakeWeight ; 
    }

    if (applyMuIDCorrections_ && hnuEvent.isMC)
    {
        double mu1wgt = muid_->weightForMC((hnuEvent.mu1.pt()), applyMuIDEffsign_) ;
        double e1wgt  = muid_->weightElectronsForMC(hnu::getElectronSCEta(hnuEvent.e1), applyMuIDEffsign_) ;
        hnuEvent.eventWgt *= mu1wgt * e1wgt ;
    }
    
    hnuEvent.btagName = btagName;
    hnuEvent.isBJet1 = hnuEvent.j1.bDiscriminator(btagName) >= minBtagDiscVal;
    hnuEvent.isBJet2 = hnuEvent.j2.bDiscriminator(btagName) >= minBtagDiscVal;

    hists.cutlevel->Fill(1.0, hnuEvent.eventWgt); // Two highest pT muons that are isolated, separated from chosen jets

    //hnuEvent.met1 = pMET->at(0);

    //hnuEvent.regularize(); // assign internal primary lepton variables  //this is called in calculate now and no longer needs to be called here
    //hnuEvent.scaleMuE(applyMESfactor_,hnuEvent.ElecScale);
    hnuEvent.calculate(correctEscale_); // calculate various details

    // Fill slope fit tuple here
    if(addSlopeTree_)
    {
        hnuTree_->event_.mll = hnuEvent.mLL;
        hnuTree_->event_.mlljj = hnuEvent.mWR;
        hnuTree_->event_.l1pt = hnuEvent.l1pt;
        hnuTree_->event_.l2pt = hnuEvent.l2pt;
        hnuTree_->event_.weight = hnuEvent.eventWgt;
        hnuTree_->event_.flavor = hnuEvent.mode;
        hnuTree_->event_.n_pileup = hnuEvent.n_pue;
        hnuTree_->event_.n_primaryVertex = hnuEvent.n_primary_vertex;
        hnuTree_->event_.cutlevel = 1;
    }

    // Basic requirements on muon, electron, jets
    hists.LLJJptCuts->fill(hnuEvent);

    // require that one muon be BOTH tight and trigger-matched
    bool mu1trig = true ;
    if ( trig_->matchingEnabled() && iEvent.isRealData() )
    {
        mu1trig = mu1trig &&
                trig_->isTriggerMatched( hnuEvent.mu1, iEvent ) ;
    }
    else if (!iEvent.isRealData())
    {
        mu1trig = mu1trig &&
                trig_->simulateForMC( hnuEvent.mu1.pt(), hnuEvent.mu1.eta(), applyTrigEffsign_ );
    }

    if( !mu1trig ) 
    {
        if(addSlopeTree_) hnuTree_->fill();
        return false;
    }
    hists.TrigMatches->fill(hnuEvent);
    hists.cutlevel->Fill(2.0, hnuEvent.eventWgt); // Trigger
    if(addSlopeTree_) hnuTree_->event_.cutlevel = 2;

    //--- Impose vertex requirement here ---//
    float deltaVzJ1J2 = fabs(hnuEvent.tjV1 - hnuEvent.tjV2);
    float deltaVzJ1M1 = fabs(hnuEvent.tjV1 - hnuEvent.mu1.vertex().Z());
    float deltaVzJ2E1 = fabs(hnuEvent.tjV2 - hnuEvent.e1.vertex().Z());
    float deltaVzJ1E1 = fabs(hnuEvent.tjV1 - hnuEvent.e1.vertex().Z());
    float deltaVzJ2M1 = fabs(hnuEvent.tjV2 - hnuEvent.mu1.vertex().Z());
    float deltaVzM1E1 = fabs(hnuEvent.mu1.vertex().Z() - hnuEvent.e1.vertex().Z());
    if ((cuts.maxJetVZsepCM > 0 && cuts.maxVertexZsep > 0) &&
        ((deltaVzJ1J2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M1 >= cuts.maxJetVZsepCM) ||
        (deltaVzJ2E1 >= cuts.maxJetVZsepCM) || (deltaVzJ1E1 >= cuts.maxJetVZsepCM) ||
        (deltaVzJ2M1 >= cuts.maxJetVZsepCM) || (deltaVzM1E1 >= cuts.maxVertexZsep)) )
    {
        if(addSlopeTree_) hnuTree_->fill();
        return false ;
    }
    hists.VertexCuts->fill(hnuEvent);
    hists.cutlevel->Fill(3.0, hnuEvent.eventWgt); // Vertex
    if(addSlopeTree_) hnuTree_->event_.cutlevel = 3;

    double mu1pt = hnuEvent.mu1.pt() ;
    double e1pt  = hnu::getElectronEt(hnuEvent.e1,(heepVersion_ != 40 || heepVersion_ != 41)) ;
    double highestPt = std::max(mu1pt, e1pt);
    if ( studyScaleFactorEvolution_ && hnuEvent.mLL >= cuts.minimum_mumu_mass)
    {
        if ( highestPt > 30. )  hists.Mu1Pt30GeVCut->fill(hnuEvent);
        if ( highestPt > 40. )  hists.Mu1Pt40GeVCut->fill(hnuEvent);
        if ( highestPt > 50. )  hists.Mu1Pt50GeVCut->fill(hnuEvent);
        if ( highestPt > 60. )  hists.Mu1Pt60GeVCut->fill(hnuEvent);
        if ( highestPt > 80. )  hists.Mu1Pt80GeVCut->fill(hnuEvent);
        if ( highestPt > 100. ) hists.Mu1Pt100GeVCut->fill(hnuEvent);
    }

    if(highestPt < cuts.minimum_lep1_pt)
    {
        if(addSlopeTree_) hnuTree_->fill();
        return false;
    }

    hists.cutlevel->Fill(4.0, hnuEvent.eventWgt); // Primary lepton pT
    if(addSlopeTree_) hnuTree_->event_.cutlevel = 4;

    if ( studyScaleFactorEvolution_ )
    {
        if ( hnuEvent.n_primary_vertex == 1 )
            hists.Mu1HighPtCutVtxEq1->fill(hnuEvent);
        else if ( hnuEvent.n_primary_vertex <= 5 )
            hists.Mu1HighPtCutVtx2to5->fill(hnuEvent);
        else if ( hnuEvent.n_primary_vertex > 5 )
            hists.Mu1HighPtCutVtxGt5->fill(hnuEvent);
    }
    hists.Mu1HighPtCut->fill(hnuEvent);
    if(hnuEvent.isBJet1 || hnuEvent.isBJet2) hists.Mu1HighPtCut_1bjet->fill(hnuEvent);
    if(hnuEvent.isBJet1 && hnuEvent.isBJet2) hists.Mu1HighPtCut_2bjet->fill(hnuEvent);

    if(hnuEvent.mLL >= cuts.minimum_mumu_mass) { 
        hists.diLmassCut->fill(hnuEvent);
	hists.cutlevel->Fill(5.0, hnuEvent.eventWgt); // Dilepton mass

        if(addSlopeTree_ && !(hnuEvent.mWR >= cuts.minimum_mWR_mass))
        {
            hnuTree_->event_.cutlevel = 5;
            hnuTree_->fill();
        }
    }

    if ( iEvent.isRealData() )
    {
        std::cout << "\t" << iEvent.id() << std::endl;
        std::cout << "\tM(W_R)  = " << hnuEvent.mWR  << " GeV";
        std::cout << ", M(NuR1) = " << hnuEvent.mNuR1 << " GeV";
        std::cout << ", M(NuR2) = " << hnuEvent.mNuR2 << " GeV" << std::endl;
        std::cout << "\tM(mumu) = " << hnuEvent.mLL << " GeV";
        std::cout << ", M(JJ) = "  << hnuEvent.mJJ << " GeV" << std::endl;
        std::cout << "\tJets:   j1 "
                << "pt=" << hnuEvent.j1.pt() << " GeV, eta=" << hnuEvent.j1.eta() << ", phi=" << hnuEvent.j1.phi();
        std::cout <<        ", j2 "
                << "pt=" << hnuEvent.j2.pt() << " GeV, eta=" << hnuEvent.j2.eta() << ", phi=" << hnuEvent.j2.phi();
        std::cout << std::endl;
        std::cout << "\tLeptons: mu1 "
                << "pt=" << hnuEvent.mu1.pt() << " GeV, eta=" << hnuEvent.mu1.eta() << ", phi=" << hnuEvent.mu1.phi();
        std::cout <<          ", e1 "
                << "pt=" << hnuEvent.e1.pt() << " GeV, eta=" << hnuEvent.e1.eta() << ", phi=" << hnuEvent.e1.phi();
        std::cout << std::endl;
    }

    // Change the final logic of the filter for e/mu studies:
    // Keep anything passing the 60 GeV lepton 1 cut
    if ( hnuEvent.mLL >= cuts.minimum_mumu_mass && hnuEvent.mWR >= cuts.minimum_mWR_mass ) { 
      hists.cutlevel->Fill(6.0, hnuEvent.eventWgt); // Event meets W_R mass requirements
      hists.mWRmassCut->fill(hnuEvent);

      if(addSlopeTree_)
      {
          hnuTree_->event_.cutlevel = 5;
          hnuTree_->fill();
      }
    }

    return true;
}

// ------------ method called once each job just before starting event loop  ------------

void
HeavyNuTop::beginJob()
{
    // nnif_->beginJob();
    firstEvent_ = true;
}

// ------------ method called once each job just after ending the event loop  ------------

void
HeavyNuTop::endJob()
{
    // nnif_->endJob();
    trig_->endJob();
    muid_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuTop);
