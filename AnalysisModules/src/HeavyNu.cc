// -*- C++ -*-
//
// Package:    HeavyNu
// Class:      HeavyNu
// 
/**\class HeavyNu HeavyNu.cc HeavyNu/AnalyzerModules/src/HeavyNu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: HeavyNu.cc,v 1.94 2012/05/09 17:07:35 pastika Exp $
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
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TRandom.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuMuHist.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


template <class T>
inline std::string int2str(T i)
{
    std::ostringstream ss;
    ss << i;
    return ss.str();
}

const std::vector<hNuMassHypothesis> v_null;

template <class T> void outputCandidate(const T& p)
{
    std::cout << "pt=" << p.pt() << " GeV, eta=" << p.eta() << ", phi=" << p.phi();
}


static std::string btagName;
static double minBtagDiscVal; // for discriminating B-tagged jets.

const int muonQualityFlags = 4;
const std::string muonQuality[] = {
    "All", "AllGlobalMuons", "AllStandAloneMuons", "AllTrackerMuons"
};

inline void labelMuonQualAxis(TAxis *ax)
{
    //for(int i = 0; i < muonQualityFlags; i++)
    //{
    //    ax->SetBinLabel(i + 1, muonQuality[i].c_str());
    //    ax->SetBinLabel(i + 1, muonQuality[i].c_str());
    //}
}

class HeavyNu : public edm::EDFilter
{
public:
    explicit HeavyNu(const edm::ParameterSet&);
    ~HeavyNu();


private:

    virtual void respondToOpenInputFile(edm::FileBlock const& fb)
    {
        currentFile_ = fb.fileName();
    }

    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    virtual void fillBasicMuHistos(const pat::Muon& m);
    virtual void fillBasicJetHistos(const pat::Jet& j,
            int jetnum);
    virtual TH1 *bookRunHisto(uint32_t runNumber);

    virtual void studyMuonSelectionEff(const std::vector< std::pair<pat::Muon,pat::GenericParticle> >& trackTP, 
				       const std::vector< std::pair<pat::Muon,pat::Muon> >& tightTP, 
				       const std::vector< std::pair<pat::Muon,pat::Muon> >& cjTP, 
				       const std::vector< pat::Muon >& tightMuons, 
				       const edm::Handle<reco::TrackCollection>& gTracks,
				       const edm::Handle<reco::BeamSpot>& beamspot, 
				       double wgt, int nJets);

    virtual void studyJetVertex(edm::Handle<pat::JetCollection>& pJets,
                                edm::Handle<reco::JPTJetCollection>& jptJets,
                                edm::Handle<pat::MuonCollection>& pMuons, int npue);

    inline bool inZmassWindow(double mMuMu)
    {
        return(mMuMu <= ZwinMaxGeV_) && (mMuMu >= ZwinMinGeV_);
    }

    void fill(pat::MuonCollection muons,
              pat::ElectronCollection electrons,
              pat::JetCollection jets,
              pat::METCollection metc,
              bool isMC,
              double wgt,
              bool pfJets,
              int nPU, int nPV,
              HeavyNuHistSet* hnmh);

    bool isWrDaughter(const reco::Candidate* mother);

    edm::InputTag muonTag_;
    edm::InputTag trackTag_;
    edm::InputTag jetTag_;
    edm::InputTag metTag_;
    edm::InputTag elecTag_;

  edm::InputTag elecRhoTag_ ; 
  double elecRho_ ; 

    int evtCounter;

    double ZwinMinGeV_, ZwinMaxGeV_; // for trigger efficiency studies

    int jecVal_; // Jet correction eras: 0 = 2010(A+B), 1 = 2010A, 2 = 2010B, 3 = 2011A
    int applyJECUsign_; // for Jet Energy Correction Uncertainty studies
    int applyJERsign_; // for Jet Energy Resolution and Resolution Uncertainty studies
    double applyMESfactor_; // for Muon Energy Scale studies
    bool merUncertainty_ ; 
    int applyTrigEffsign_; // for Trigger Efficiency studies
    bool highestPtTriggerOnly_;
    bool applyMuIDCorrections_;
    int applyMuIDEffsign_; // for Loose+Iso efficiency correction
    bool studyMuonSelectionEff_; // for Muon Loose ID Efficiency studies
    bool oneTP_ ;
    bool studyScaleFactorEvolution_; // for Top, Z+jets scale factors (by Mu1 pT) studies
    int tpSeed_ ; 
    TRandom *tpRandom_;

    bool disableTriggerCorrection_; // turn off the trigger emulator
    int pileupEra_;
    double puShift_ ;
    reweight::PoissonMeanShifter poissonNvtxShifter_ ; 
    bool isPFJets_; // true if PFJets are used (turns off jet id requirement)
    bool studyRatePerRun_ ; // should only be true for data
    bool studyAlternativeSelection_ ; 
    
    std::string currentFile_;
    bool dolog_;
    bool firstEvent_;

    int heepVersion_;

    HeavyNuTrigger *trig_;
    HeavyNuID *muid_;
    JetCorrectionUncertainty *jecuObj_;

    edm::LumiReWeighting MCweightByVertex_;

    std::map<uint32_t, TH1 *> m_runHistos_;

    HeavyNuEvent::anal_type analysisMode_;

    bool doPDFreweight_;
    std::string pdfReweightBaseName, pdfReweightTargetName;
    int pdfReweightBaseId, pdfReweightTargetId;


    // ----------member data ---------------------------
    bool init_;

    // gf set of histo for all Z definitions in a stack

    struct HistStruct
    {
        TH1 *nelec, *njet, *nmet, *nmuAll, *nmuLoose, *nmuTight;
        TH1 *muPt, *muEta, *muPhi, *looseMuPt, *tightMuPt;
        TH1* mc_type;

        TH1* cutlevel;
        TH1* weights;

        TH1* z2jetPerRun; 
        
        // Muon quality histos as a function of Pt
        TH2 *muNvalidHitsVsPt, *mudBvsPt, *muNormChi2vsPt, *muQualVsPt;
        TH2 *muNmatchesVsPt, *muNvalidMuonHitsVsPt, *muNvalidPixelHitsVsPt;
        TH2 *muTrckIsoVsPt, *muHcalIsoVsPt, *muEcalIsoVsPt, *muCaloIsoVsPt;

        TH1 *jetPt, *jetEta, *jetPhi, *jetID, *jecUncHi, *jecUncLo, *met;
        TH2 *dVzMuJets, *dVzMuMus;
        TH2 *jetPtvsNum;
        TProfile2D *jecUncHiVsEtaPt, *jecUncLoVsEtaPt;

        TH1 *trkIsoStudy;
        TH1 *closejetMu2tagMu1probeInZwin, *closejetMu2tagMu1passInZwin;
        TH1 *closejetMu1tagMu2probeInZwin, *closejetMu1tagMu2passInZwin;

        TFileDirectory *rundir;
        // HistPerDef twoL;
        HeavyNuHistSet* noCuts;
        //HeavyNuHistSet* LLptCuts;
        HeavyNuHistSet* TrigMatches;
        HeavyNuHistSet* LLJJptCuts;
        //HeavyNuHistSet* VertexCuts;
        HeavyNuHistSet* Mu1HighPtCut;
        HeavyNuHistSet* Mu1Pt40GeVCut;
        HeavyNuHistSet* Mu1Pt50GeVCut;
        HeavyNuHistSet* Mu1Pt60GeVCut;
        HeavyNuHistSet* Mu1Pt80GeVCut;
        HeavyNuHistSet* Mu1Pt100GeVCut;
        HeavyNuHistSet* AlternativeElecChanPt;
        HeavyNuHistSet* AlternativeMu1Pt40;
        HeavyNuHistSet* AlternativeMu2Pt40;
        HeavyNuHistSet* AlternativeMu2Pt60;
        HeavyNuHistSet* AlternativeJetPt60;
        HeavyNuHistSet* AlternativeBarrelLoose;
        HeavyNuHistSet* AlternativeBarrelTight;
        HeavyNuHistSet* AlternativeAtLeastOneBjet;
        HeavyNuHistSet* AlternativeTwoBjets;
        HeavyNuHistSet* AlternativeDimuonMass120;
        HeavyNuHistSet* Mu1HighPtCutVtxEq1;
        HeavyNuHistSet* Mu1HighPtCutVtx2to5;
        HeavyNuHistSet* Mu1HighPtCutVtxGt5;
        HeavyNuHistSet* Mu1HighPtCutNoJets;
        HeavyNuHistSet* Mu1HighPtCut1Jet;
        HeavyNuHistSet* diLmassCut;
        //HeavyNuHistSet* loDiLmassCut;
        HeavyNuHistSet* mWRmassCut;
        HeavyNuHistSet* oneBtag;
        HeavyNuHistSet* twoBtag;
        // efficiency studies:
        // HeavyNuHistSet* Mu1tagInZwin;
        // HeavyNuHistSet* Mu2tagInZwin;
        // HeavyNuHistSet* Mu1tagMu2passesInZwin;
        // HeavyNuHistSet* Mu2tagMu1passesInZwin;
        HeavyNuHistSet* TightTagTrackProbeInZwin0jets;
        HeavyNuHistSet* TightTagTrackProbeInZwin1jet;
        HeavyNuHistSet* TightTagTrackProbeInZwin2jets;
        HeavyNuHistSet* TightTagTrackProbePassesInZwin0jets;
        HeavyNuHistSet* TightTagTrackProbePassesInZwin1jet;
        HeavyNuHistSet* TightTagTrackProbePassesInZwin2jets;
        HeavyNuHistSet* TightTagTightProbeInZwin0jets;
        HeavyNuHistSet* TightTagTightProbeInZwin1jet;
        HeavyNuHistSet* TightTagTightProbeInZwin2jets;
        HeavyNuHistSet* TightTagTightProbePassesInZwin0jets;
        HeavyNuHistSet* TightTagTightProbePassesInZwin1jet;
        HeavyNuHistSet* TightTagTightProbePassesInZwin2jets;
        HeavyNuHistSet* TightTagTightCJProbeInZwin0jets;
        HeavyNuHistSet* TightTagTightCJProbeInZwin1jet;
        HeavyNuHistSet* TightTagTightCJProbeInZwin2jets;
        HeavyNuHistSet* TightTagTightCJProbePassesInZwin0jets;
        HeavyNuHistSet* TightTagTightCJProbePassesInZwin1jet;
        HeavyNuHistSet* TightTagTightCJProbePassesInZwin2jets;
        // Trigger efficiency studies force two jets
        HeavyNuHistSet* TightTagTrigProbeInZwin;
        HeavyNuHistSet* TightTagTrigProbePassesInZwin;
        HeavyNuHistSet* Mu1TrigMatchesInZwin;
        HeavyNuHistSet* Mu2TrigMatchesInZwin;
        // HeavyNuHistSet* Mu1Mu2TrigMatchesInZwin;
    } hists;

    struct CutsStruct
    {
        double minimum_mu1_pt;
        double minimum_mu2_pt;
        double minimum_jet_pt;
        double maximum_mu_abseta;
        double maximum_jet_abseta;
        double minimum_mumu_mass;
        double minimum_mWR_mass;
        double minimum_muon_jet_dR;
        double muon_trackiso_limit;
        double maxVertexZsep;
        double maxJetVZsepCM;
    } cuts;

};

bool HeavyNu::isWrDaughter(const reco::Candidate* mother)
{
    for(size_t i = 0; i < mother->numberOfMothers(); i++)
    {
        if(mother->mother(i)->pdgId() == (analysisMode_ == HeavyNuEvent::HNUMU)?9900014:9900012 || isWrDaughter(mother->mother(i))) return true;
    }
    return false;
}

void HeavyNu::fill(pat::MuonCollection muons,
                   pat::ElectronCollection electrons,
                   pat::JetCollection jets,
                   pat::METCollection metc,
                   bool isMC,
                   double wgt,
                   bool pfJets,
                   int nPU, int nPV,
                   HeavyNuHistSet *hnmh)
{
    HeavyNuEvent hne(analysisMode_);
    hne.isMC = isMC;
    hne.eventWgt = wgt;
    hne.pfJets = pfJets;
    hne.nJets = jets.size();
    hne.n_pue = nPU ;
    hne.n_primary_vertex = nPV ;

    std::sort(muons.begin(), muons.end(), hnu::pTcompare());
    std::sort(electrons.begin(), electrons.end(), hnu::pTcompare());
    std::sort(jets.begin(), jets.end(), hnu::pTcompare());

    if((analysisMode_ == HeavyNuEvent::HNUMU || analysisMode_ == HeavyNuEvent::TOP) && muons.size() > 0)
    {
        hne.mu1 = muons[0];
        hne.nLeptons++;
    }
    if(analysisMode_ == HeavyNuEvent::HNUMU && muons.size() > 1)
    {
        hne.mu2 = muons[1];
        hne.nLeptons++;
    }

    if((analysisMode_ == HeavyNuEvent::HNUE || analysisMode_ == HeavyNuEvent::TOP) && electrons.size() > 0)
    {
        hne.e1 = electrons[0];
        hne.nLeptons++;
    }
    if(analysisMode_ == HeavyNuEvent::HNUE && electrons.size() > 1)
    {
        hne.e2 = electrons[1];
        hne.nLeptons++;
    }

    hne.regularize();

    hne.btagName = btagName;
    if(hne.nJets > 0)
    {
        hne.j1 = jets[0];
        double j1bdisc = hne.j1.bDiscriminator(btagName);
        hne.isBJet1 = j1bdisc >= minBtagDiscVal;
    }
    if(hne.nJets > 1)
    {
        hne.j2 = jets[1];
        double j2bdisc = hne.j2.bDiscriminator(btagName);
        hne.isBJet2 = j2bdisc >= minBtagDiscVal;
    }

    hne.met1 = metc[0];

    hne.calculateLL();
    hne.calculate();

    hnmh->fill(hne);
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HeavyNu::HeavyNu(const edm::ParameterSet& iConfig)
{
    // ==================== Get all parameters ====================
    //
    dolog_ = iConfig.getParameter<bool>("DoLog");

    muonTag_ = iConfig.getParameter< edm::InputTag > ("muonTag");
    jetTag_ = iConfig.getParameter< edm::InputTag > ("jetTag");
    metTag_ = iConfig.getParameter< edm::InputTag > ("metTag");
    elecTag_ = iConfig.getParameter< edm::InputTag > ("electronTag");
    trackTag_ = iConfig.getParameter< edm::InputTag > ("trackTag");

    elecRhoTag_ = iConfig.getParameter< edm::InputTag > ("electronRho");    

    btagName = iConfig.getParameter<std::string > ("BtagName");
    minBtagDiscVal = iConfig.getParameter<double>("minBtagDiscr");

    cuts.minimum_mu1_pt = iConfig.getParameter<double>("minMu1pt");
    cuts.minimum_mu2_pt = iConfig.getParameter<double>("minMu2pt");
    cuts.minimum_jet_pt = iConfig.getParameter<double>("minJetPt");
    cuts.maximum_mu_abseta = iConfig.getParameter<double>("maxMuAbsEta");
    cuts.maximum_jet_abseta = iConfig.getParameter<double>("maxJetAbsEta");
    cuts.minimum_mumu_mass = iConfig.getParameter<double>("minMuMuMass");
    cuts.minimum_mWR_mass = iConfig.getParameter<double>("min4objMass");
    cuts.minimum_muon_jet_dR = iConfig.getParameter<double>("minMuonJetdR");
    cuts.muon_trackiso_limit = iConfig.getParameter<double>("muonTrackRelIsoLimit");
    cuts.maxVertexZsep = iConfig.getParameter<double>("maxVertexZsepCM");
    cuts.maxJetVZsepCM = iConfig.getParameter<double>("maxJetVZsepCM");

    ZwinMinGeV_ = iConfig.getParameter<double>("ZmassWinMinGeV");
    ZwinMaxGeV_ = iConfig.getParameter<double>("ZmassWinMaxGeV");

    jecVal_ = iConfig.getParameter<int>("jecEra");

    pileupEra_ = iConfig.getParameter<int>("pileupEra");
    puShift_   = iConfig.getParameter<double>("systPileupShift") ;
    if ( fabs(puShift_) > 0.001 ) poissonNvtxShifter_ = reweight::PoissonMeanShifter( puShift_ ) ; 
    disableTriggerCorrection_ = iConfig.getParameter<bool>("DisableTriggerCorrection");

    applyJECUsign_ = iConfig.getParameter<int>("applyJECUsign");
    if(applyJECUsign_) applyJECUsign_ /= abs(applyJECUsign_); // ensure -1,0,+1

    applyJERsign_ = iConfig.getParameter<int>("applyJERsign");
    if(applyJERsign_) applyJERsign_ /= abs(applyJERsign_); // ensure -1,0,+1

    applyMESfactor_ = iConfig.getParameter<double>("applyMESfactor");
    merUncertainty_ = iConfig.getParameter<bool>("checkMERUnc");

    applyMuIDCorrections_ = iConfig.getParameter<bool>("applyMuIDEffcorr");
    applyMuIDEffsign_ = iConfig.getParameter<int>("applyMuIDEffsign");
    if(applyMuIDEffsign_) applyMuIDEffsign_ /= abs(applyMuIDEffsign_); // ensure -1,0,+1

    highestPtTriggerOnly_ = iConfig.getParameter<bool>("highestPtTriggerOnly");
    applyTrigEffsign_ = iConfig.getParameter<int>("applyTrigEffsign");
    if(applyTrigEffsign_) applyTrigEffsign_ /= abs(applyTrigEffsign_); // ensure -1,0,+1

    studyMuonSelectionEff_ = iConfig.getParameter<bool>("studyMuSelectEff");
    oneTP_  = iConfig.getParameter<bool>("oneTPcand");
    if (studyMuonSelectionEff_ && oneTP_) { 
      tpSeed_   = iConfig.getParameter<int>("tpRandomSeed");
      tpRandom_ = new TRandom(tpSeed_);
    }

    studyScaleFactorEvolution_ = iConfig.getParameter<bool>("studyScaleFactor");
    studyAlternativeSelection_ = iConfig.getParameter<bool>("alternativeSelections");
    
    isPFJets_ = iConfig.getParameter<bool>("isPFJets");
    studyRatePerRun_ = iConfig.getParameter<bool>("studyRatePerRun");

    // Default HEEP version is 4.0 (2012 selection)
    heepVersion_ = iConfig.getUntrackedParameter<int>("heepVersion");
    if ( heepVersion_ < 30 || heepVersion_ > 40 ) heepVersion_ = 40 ;     

    std::string am = iConfig.getUntrackedParameter<std::string>("analysisMode");
    if(!am.compare("HNUMU")) analysisMode_ = HeavyNuEvent::HNUMU;
    else if(!am.compare("HNUE")) analysisMode_ = HeavyNuEvent::HNUE;
    else if(!am.compare("TOP")) analysisMode_ = HeavyNuEvent::TOP;
    else if(!am.compare("QCD")) analysisMode_ = HeavyNuEvent::QCD;
    else if(!am.compare("CLO")) analysisMode_ = HeavyNuEvent::CLO;
    else std::cout << "!!!!!!!!INVALID ANALYSIS MODE : " << am << " !!!!!!!!\noptions are: HNUMU, HNUE, TOP, QCD, CLO" << std::endl;

    doPDFreweight_=iConfig.getUntrackedParameter<bool>("doPDFReweight",false);
    if (doPDFreweight_) {
      pdfReweightBaseName=iConfig.getUntrackedParameter<std::string>("pdfReweightBaseName");
      pdfReweightTargetName=iConfig.getUntrackedParameter<std::string>("pdfReweightTargetName");
      pdfReweightBaseId=iConfig.getUntrackedParameter<int>("pdfReweightBaseId",0);
      pdfReweightTargetId=iConfig.getUntrackedParameter<int>("pdfReweightTargetId",0);
      std::cout << "PDF Reweighting from " << pdfReweightBaseName << ":" << pdfReweightBaseId
		<< " to " << pdfReweightTargetName << ":" << pdfReweightTargetId
		<< std::endl;
    }
    
    // ==================== Init other members ====================
    //

    //    nnif_ = new HeavyNu_NNIF(iConfig);
    trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet > ("trigMatchPset"));
    muid_ = new HeavyNuID(iConfig.getParameter<edm::ParameterSet > ("muIDPset"));
    // ==================== Book the histos ====================
    //
    edm::Service<TFileService> fs;
    hists.mc_type = fs->make<TH1D > ("mc_type", "MC Type Code", 100, -0.5, 99.5);
    if ( studyRatePerRun_ )
        hists.z2jetPerRun = fs->make<TH1I > ("z2jetPerRun","M3 Z #to #mu#mu events per run",20000,160000,180000); 
    hists.nelec = fs->make<TH1D > ("nelec", "N(e^{#pm})", 10, -0.5, 9.5);
    hists.nmuAll = fs->make<TH1D > ("nmuAll", "N(#mu^{#pm})", 10, -0.5, 9.5);
    hists.nmuLoose = fs->make<TH1D > ("nmuLoose", "N(#mu^{#pm}) passes Loose", 10, -0.5, 9.5);
    hists.nmuTight = fs->make<TH1D > ("nmuTight", "N(#mu^{#pm}) passes Tight", 10, -0.5, 9.5);
    hists.cutlevel = fs->make<TH1D > ("cutlevel", "Cut Level", 11, -1.5, 9.5);
    hists.cutlevel->GetXaxis()->SetBinLabel(1, "Raw");
    hists.cutlevel->GetXaxis()->SetBinLabel(2, "No cuts");
    hists.cutlevel->GetXaxis()->SetBinLabel(3, "mmjj p_{T}");
    hists.cutlevel->GetXaxis()->SetBinLabel(4, "M1 (Trigger)");
    hists.cutlevel->GetXaxis()->SetBinLabel(5, "M2 (Vertex)");
    hists.cutlevel->GetXaxis()->SetBinLabel(6, "M3 (high p_{T})");
    hists.cutlevel->GetXaxis()->SetBinLabel(7, "M4 (M_{#mu#mu})");
    hists.cutlevel->GetXaxis()->SetBinLabel(8, "M5 (M(W_{R})");
    hists.weights = fs->make<TH1D > ("weights", "event weights", 1000, 0.0, 10.0);
    hists.njet = fs->make<TH1D > ("njet", "N(Jet)", 50, -0.5, 49.5);
    hists.nmet = fs->make<TH1D > ("nmet", "N(MET)", 50, -0.5, 49.5);
    hists.met = fs->make<TH1D > ("met", "MET distribution", 100, 0, 2000);
    hists.muPt = fs->make<TH1D > ("muPt", "#mu p_{T} distribution", 100, 0, 2000);
    hists.muEta = fs->make<TH1D > ("muEta", "#mu #eta distribution", 50, -2.5, 2.5);
    hists.muPhi = fs->make<TH1D > ("muPhi", "#mu #phi distribution", 60, -3.14159, 3.14159);
    hists.jetPt = fs->make<TH1D > ("jetPt", "jet p_{T} distribution", 100, 0, 2000);
    hists.jetEta = fs->make<TH1D > ("jetEta", "jet #eta distribution", 50, -5, 5);
    hists.jetPhi = fs->make<TH1D > ("jetPhi", "jet #phi distribution", 60, -3.14159, 3.14159);
    hists.jetID = fs->make<TH1I > ("jetID", "Jet ID", 3, 0, 3);
    hists.jetPtvsNum = fs->make<TH2D > ("jetPtvsNum", "Jet P_{T} vs. Jet # ", 11, -0.5, 10.5, 200, 0., 2000.);
    hists.dVzMuJets = fs->make<TH2D > ("dVzMuJets", "vertex dz Mu Jets vs npue", 500, -5.0, 5.0, 50, -0.5, 49.5);
    hists.dVzMuMus = fs->make<TH2D > ("dVzMuMus", "vertex dz Mu Mus vs npue", 500, -5.0, 5.0, 50, -0.5, 49.5);

    hists.jetID->GetXaxis()->SetBinLabel(1, "Neither");
    hists.jetID->GetXaxis()->SetBinLabel(2, "Loose");
    hists.jetID->GetXaxis()->SetBinLabel(3, "Tight");

    hists.trkIsoStudy = fs->make<TH1D > ("trkIsoStudy", "Tracker Relative Isolation", 100, 0., 1.);
    hists.closejetMu1tagMu2probeInZwin = fs->make<TH1D > ("closejetMu1tagMu2probeInZwin", "Probe (#mu_{2}) p_{T}", 100, 0., 1000.);
    hists.closejetMu1tagMu2passInZwin = fs->make<TH1D > ("closejetMu1tagMu2passInZwin", "Passing probe (#mu_{2}) p_{T}", 100, 0., 1000.);
    hists.closejetMu2tagMu1probeInZwin = fs->make<TH1D > ("closejetMu2tagMu1probeInZwin", "Probe (#mu_{1}) p_{T}", 100, 0., 1000.);
    hists.closejetMu2tagMu1passInZwin = fs->make<TH1D > ("closejetMu2tagMu1passInZwin", "Passing probe (#mu_{1}) p_{T}", 100, 0., 1000.);

    if(applyMESfactor_ == 1.0)
    { // otherwise don't bother

        // Loose/Tight vs Pt
        hists.looseMuPt = fs->make<TH1D > ("looseMuPt", "#mu p_{T}, passes VBTF Loose", 100, 0, 2000);
        hists.tightMuPt = fs->make<TH1D > ("tightMuPt", "#mu p_{T}, passes VBTF Tight", 100, 0, 2000);

        //  Muon quality variables vs muon p_T
        //
        hists.muNvalidHitsVsPt = fs->make<TH2D > ("muNvalidHitsVsPt",
                "#mu # Valid hits vs p_{T}; p_{T}(#mu) (GeV); # Tracker+Pixel Hits",
                150, 0, 3000, 50, 0, 50);
        hists.mudBvsPt = fs->make<TH2D > ("mudBvsPt",
                "#mu dXY vs p_{T}; p_{T}(#mu) (GeV); dXY(#mu)",
                150, 0, 3000, 40, 0., 20.);
        hists.muNormChi2vsPt = fs->make<TH2D > ("muNormChi2vsPt",
                "#mu Norm #chi^{2} vs p_{T}; p_{T}(#mu) (GeV); norm #chi^{2}",
                150, 0, 3000, 25, 0, 50);
        hists.muQualVsPt = fs->make<TH2D > ("muIDvsPt",
                "Qual(#mu) vs p_{T}(#mu); p_{T}(#mu) (GeV)",
                150, 0, 3000, 4, 0, 4);
        hists.muNmatchesVsPt = fs->make<TH2D > ("muNmatchesVsPt",
                "#mu # Matches vs p_{T}; p_{T}(#mu) (GeV); # Matches",
                150, 0, 3000, 10, 0, 10);
        hists.muNvalidMuonHitsVsPt = fs->make<TH2D > ("muNvalidMuonHitsVsPt",
                "# Valid Muon hits vs p_{T}; p_{T}(#mu) (GeV); # Valid Muon Hits",
                150, 0, 3000, 50, 0, 50);
        hists.muNvalidPixelHitsVsPt = fs->make<TH2D > ("muNvalidPixelHitsVsPt",
                "# Valid Pixel hits vs p_{T}; p_{T}(#mu) (GeV); # Valid Pixel Hits",
                150, 0, 3000, 10, 0, 10);

        labelMuonQualAxis(hists.muQualVsPt->GetXaxis());

        hists.muTrckIsoVsPt = fs->make<TH2D > ("muTrckIsoVsPt", "trackIso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);trackIso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muHcalIsoVsPt = fs->make<TH2D > ("muHcalIsoVsPt", "HCAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);HCAL Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muEcalIsoVsPt = fs->make<TH2D > ("muEcalIsoVsPt", "ECAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);ECAL Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muCaloIsoVsPt = fs->make<TH2D > ("muCaloIsoVsPt", "Calo Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);Calo Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
    }

    // Histos per cut:
    //
    // hists.twoL->book(new TFileDirectory(fs->mkdir("cutm1_LL")), "(two lepton:-1)");
    if(analysisMode_ == HeavyNuEvent::HNUMU)
    {
        hists.noCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)");
        hists.LLJJptCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)");
        hists.TrigMatches = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:2)");
        hists.Mu1HighPtCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4_L1HighPt")), "(Mu1 High pt cut:4)");
        hists.diLmassCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(mumu mass cut:5)");
        hists.mWRmassCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(mumujj mass cut:6)");

        if ( studyAlternativeSelection_ )
        {
            hists.AlternativeElecChanPt = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu1Pt80Mu2Pt40")), "(Alternative: Electron pT selection)");
            hists.AlternativeMu1Pt40 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu1Pt40")), "(Alternative: mu1 pT 40 GeV)");
            hists.AlternativeMu2Pt40 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu2Pt40")), "(Alternative: mu2 pT 40 GeV)");
            hists.AlternativeMu2Pt60 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu2Pt60")), "(Alternative: mu2 pT 60 GeV)");
            hists.AlternativeJetPt60 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltJetPt60")), "(Alternative: jet pT 60 GeV)");
            hists.AlternativeBarrelLoose = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltBarrelLoose")), "(Alternative: one barrel muon, loose)");
            hists.AlternativeBarrelTight = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltBarrelTight")), "(Alternative: one barrel muon, tight)");
            hists.AlternativeAtLeastOneBjet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltOneBjet")), "(Alternative: at least one b-jet)");
            hists.AlternativeTwoBjets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltTwoBjets")), "(Alternative: two b-jets)");
            hists.AlternativeDimuonMass120 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltDimuon120")), "(Alternative: dimuon mass 120 GeV)");
        }

        if(studyScaleFactorEvolution_)
        {
            hists.Mu1Pt40GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt40GeV")), "(Mu1 40 GeV pt cut)");
            hists.Mu1Pt50GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt50GeV")), "(Mu1 50 GeV pt cut)");
            hists.Mu1Pt60GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt60GeV")), "(Mu1 60 GeV pt cut)");
            hists.Mu1Pt80GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt80GeV")), "(Mu1 80 GeV pt cut)");
            hists.Mu1Pt100GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt100GeV")), "(Mu1 100 GeV pt cut)");

            hists.Mu1HighPtCutVtxEq1 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtVtxEq1")), "(Mu1 60 GeV pt cut, 1 vtx)");
            hists.Mu1HighPtCutVtx2to5 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtVtx2to5")), "(Mu1 60 GeV pt cut, 2-5 vtx)");
            hists.Mu1HighPtCutVtxGt5 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtVtxGt5")), "(Mu1 60 GeV pt cut, 6+ vtx)");

            hists.Mu1HighPtCutNoJets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtNoJets")), "(Mu1 60 GeV pt cut, no jets)");
            hists.Mu1HighPtCut1Jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPt1Jet")), "(Mu1 60 GeV pt cut, 1 jet)");
        }

        if(trig_->matchingEnabled())
        {
            TFileDirectory *tdirm1 = new TFileDirectory(fs->mkdir("Mu1TrigMatchesInZwin"));
            TFileDirectory *tdirm2 = new TFileDirectory(fs->mkdir("Mu2TrigMatchesInZwin"));
            hists.Mu1TrigMatchesInZwin = new HeavyNuMuHist(tdirm1, "(#mu1 trigger match in Z mass Window)");
            hists.Mu2TrigMatchesInZwin = new HeavyNuMuHist(tdirm2 , "(#mu2 Trigger match in Z mass Window)");
            // hists.Mu1Mu2TrigMatchesInZwin = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Mu2TrigMatchesInZwin")), "(#mu1,#mu2 Trigger match in Z mass Window)");
            trig_->book(*tdirm1, &(((HeavyNuMuHist*)hists.Mu1TrigMatchesInZwin)->trigHistos));
            trig_->book(*tdirm2, &(((HeavyNuMuHist*)hists.Mu2TrigMatchesInZwin)->trigHistos));
        }

        if(studyMuonSelectionEff_)
        {
            hists.TightTagTrackProbeInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbeInZwin0jets")), "(probe_{0}, ID in Z mass Window)", 2);
            hists.TightTagTrackProbeInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbeInZwin1jet")), "(probe_{1}, ID in Z mass Window)", 2);
            hists.TightTagTrackProbeInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbeInZwin2jets")), "(probe_{2}, ID in Z mass Window)", 2);
            hists.TightTagTrackProbePassesInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbePassesInZwin0jets")), "(probe_{0} passes, ID in Z mass Window)", 2);
            hists.TightTagTrackProbePassesInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbePassesInZwin1jet")), "(probe_{1} passes, ID in Z mass Window)", 2);
            hists.TightTagTrackProbePassesInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbePassesInZwin2jets")), "(probe_{2} passes, ID in Z mass Window)", 2);
            hists.TightTagTightProbeInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbeInZwin0jets")), "(probe_{0}, iso in Z mass Window)", 2);
            hists.TightTagTightProbeInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbeInZwin1jet")), "(probe_{1}, iso in Z mass Window)", 2);
            hists.TightTagTightProbeInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbeInZwin2jets")), "(probe_{2}, iso in Z mass Window)", 2);
            hists.TightTagTightProbePassesInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbePassesInZwin0jets")), "(probe_{0} passes, iso in Z mass Window)", 2);
            hists.TightTagTightProbePassesInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbePassesInZwin1jet")), "(probe_{1} passes, iso in Z mass Window)", 2);
            hists.TightTagTightProbePassesInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbePassesInZwin2jets")), "(probe_{2} passes, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbeInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbeInZwin0jets")), "(CJ probe_{0}, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbeInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbeInZwin1jet")), "(CJ probe_{1}, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbeInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbeInZwin2jets")), "(CJ probe_{2}, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbePassesInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbePassesInZwin0jets")), "(CJ probe_{0} passes, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbePassesInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbePassesInZwin1jet")), "(CJ probe_{1} passes, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbePassesInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbePassesInZwin2jets")), "(CJ probe_{2} passes, iso in Z mass Window)", 2);
            hists.TightTagTrigProbeInZwin = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrigProbeInZwin")), "(Trig probe, in Z mass Window)", 2);
            hists.TightTagTrigProbePassesInZwin = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrigProbePassesInZwin")), "(Trig probe passes, in Z mass Window)", 2);
        }

    }
    if(analysisMode_ == HeavyNuEvent::HNUE)
    {
        hists.noCuts =       new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)");
        hists.LLJJptCuts =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)");
        hists.TrigMatches =  new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:2)");
        hists.Mu1HighPtCut = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4_L1HighPt")), "(L1 High pt cut:4)");
        hists.diLmassCut =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(ee mass cut:5)");
        hists.mWRmassCut =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(eejj mass cut:6)");
    }

    hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

    init_ = false;

    if(applyJECUsign_ > 0)
    {
        hists.jecUncHi = fs->make<TH1F > ("jecUncHi", "JEC Uncertainty (high); Uncertainty (%)", 100, 0, 10);
        hists.jecUncHiVsEtaPt = fs->make<TProfile2D > ("jecUncHiVsEtaPt",
                "JEC Uncertainty (high)(%);Jet #eta;Jet p_{T} (GeV)",
                50, -2.5, 2.5, 100, 0, 1000);
    }
    else if(applyJECUsign_ < 0)
    {
        hists.jecUncLo = fs->make<TH1F > ("jecUncLo", "JEC Uncertainty (low);  Uncertainty (%)", 100, 0, 10);
        hists.jecUncLoVsEtaPt = fs->make<TProfile2D > ("jecUncLoVsEtaPt",
                "JEC Uncertainty (low)(%);Jet #eta;Jet p_{T} (GeV)",
                50, -2.5, 2.5, 100, 0, 1000);
    }

    MCweightByVertex_ = edm::LumiReWeighting(hnu::generate_flat10_mc(pileupEra_),
            hnu::get_standard_pileup_data(pileupEra_));
    
    // For the record...
    std::cout << "Configurable cut values applied:" << std::endl;
    std::cout << "Analysis Mode     = " << analysisMode_ << std::endl;
    std::cout << "muonTag           = " << muonTag_ << std::endl;
    std::cout << "jetTag            = " << jetTag_ << std::endl;
    std::cout << "metTag            = " << metTag_ << std::endl;
    std::cout << "electronTag       = " << elecTag_ << std::endl;
    std::cout << "heepVersion       = " << heepVersion_ << std::endl;
    std::cout << "trackTag          = " << trackTag_ << std::endl;
    std::cout << "btagName          = " << btagName << std::endl;
    std::cout << "minBtagDiscr      = " << minBtagDiscVal << std::endl;
    std::cout << "ZmassWinMinGeV    = " << ZwinMinGeV_ << " GeV" << std::endl;
    std::cout << "ZmassWinMaxGeV    = " << ZwinMaxGeV_ << " GeV" << std::endl;
    std::cout << "minMu1pt          = " << cuts.minimum_mu1_pt << " GeV" << std::endl;
    std::cout << "minMu2pt          = " << cuts.minimum_mu2_pt << " GeV" << std::endl;
    std::cout << "minJetPt          = " << cuts.minimum_jet_pt << " GeV" << std::endl;
    std::cout << "maxMuAbsEta       = " << cuts.maximum_mu_abseta << std::endl;
    std::cout << "maxJetAbsEta      = " << cuts.maximum_jet_abseta << std::endl;
    std::cout << "minMuonJetdR      = " << cuts.minimum_muon_jet_dR << std::endl;
    std::cout << "muonTrackRelIso   = " << cuts.muon_trackiso_limit << std::endl;
    std::cout << "minMuMuMass       = " << cuts.minimum_mumu_mass << " GeV" << std::endl;
    std::cout << "min4objMass       = " << cuts.minimum_mWR_mass << " GeV" << std::endl;
    std::cout << "jecEra            = " << jecVal_ << std::endl;
    std::cout << "applyJECUsign     = " << applyJECUsign_ << std::endl;
    std::cout << "applyMESfactor    = " << applyMESfactor_ << std::endl;
    std::cout << "applyMuIDEffcorr  = " << applyMuIDCorrections_ << std::endl;
    std::cout << "applyMuIDEffsign  = " << applyMuIDEffsign_ << std::endl;
    std::cout << "applyTrigEffsign  = " << applyTrigEffsign_ << std::endl;
    std::cout << "studyMuSelectEff  = " << studyMuonSelectionEff_ << std::endl;
    std::cout << "studyScaleFactor  = " << studyScaleFactorEvolution_ << std::endl;

    std::cout << "pileup era        = " << pileupEra_ << std::endl;
    std::cout << "DisableTriggerCorrection = " << disableTriggerCorrection_ << std::endl;
    std::cout << "isPFJets          = " << isPFJets_ << std::endl;

}

HeavyNu::~HeavyNu()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

void HeavyNu::fillBasicMuHistos(const pat::Muon& m)
{
    double mupt = m.pt();
    hists.muPt->Fill(applyMESfactor_ * mupt);
    hists.muEta->Fill(m.eta());
    hists.muPhi->Fill(m.phi());

    if(applyMESfactor_ == 1.0)
    {
        hists.mudBvsPt->Fill(mupt, m.dB());

        if(hnu::isVBTFloose(m))
        {
            hists.looseMuPt ->Fill(mupt);
            hists.muNvalidHitsVsPt ->Fill(mupt, m.numberOfValidHits());
            hists.muNormChi2vsPt ->Fill(mupt, m.normChi2());
            hists.muNmatchesVsPt ->Fill(mupt, m.numberOfMatches());

            reco::TrackRef gt = m.globalTrack();
            // gt.isNonnull() guaranteed at this point?
            hists.muNvalidMuonHitsVsPt ->Fill(mupt, gt->hitPattern().numberOfValidMuonHits());
            hists.muNvalidPixelHitsVsPt->Fill(mupt, gt->hitPattern().numberOfValidPixelHits());

            if(hnu::isVBTFtight(m)) hists.tightMuPt->Fill(mupt);
        }
        hists.muQualVsPt->Fill(mupt, 0);
        for(int i = 1; i < muonQualityFlags; i++)
            if(m.muonID(muonQuality[i]))
                hists.muQualVsPt->Fill(mupt, i);

        hists.muTrckIsoVsPt->Fill(mupt, m.trackIso());
        hists.muHcalIsoVsPt->Fill(mupt, m.hcalIso());
        hists.muEcalIsoVsPt->Fill(mupt, m.ecalIso());
        hists.muCaloIsoVsPt->Fill(mupt, m.caloIso());
    }
} // HeavyNu::fillBasicMuHistos

void HeavyNu::fillBasicJetHistos(const pat::Jet& j, int jetnum)
{
    double jpt = j.pt(), jeta = j.eta();
    float totalunc = 0.0f;

    int jpdgId = 0;
    if(j.genParton()) jpdgId = j.genParton()->pdgId();
    bool isBjet = (abs(jpdgId) == 5);

    float jecuHi = 0.;
    float jecuLo = 0.;

    if(applyJECUsign_)
    {
        jecuHi = hnu::jecTotalUncertainty(jpt, jeta, jecuObj_, jecVal_, isBjet, true);
        jecuLo = hnu::jecTotalUncertainty(jpt, jeta, jecuObj_, jecVal_, isBjet, false);
        totalunc = (applyJECUsign_ > 0)?jecuHi:jecuLo;
        jpt *= (1.0 + (applyJECUsign_ * totalunc));
    }

    hists.jetPt ->Fill(jpt);
    hists.jetEta->Fill(jeta);
    hists.jetPhi->Fill(j.phi());
    if(isPFJets_) hists.jetID ->Fill(hnu::jetID(j));
    hists.jetPtvsNum->Fill(jetnum, jpt);

    jecuHi *= 100.f;
    jecuLo *= 100.f;
    if(jecuHi != jecuLo) std::cout << jeta << "\t" << jpt << "\t" << jecuHi << "\t" << jecuLo << std::endl;

    if(applyJECUsign_ > 0)
    {
        hists.jecUncHi->Fill(jecuHi);
        hists.jecUncHiVsEtaPt->Fill(jeta, j.pt(), (double)jecuHi);
    }
    else if(applyJECUsign_ < 0)
    {
        hists.jecUncLo->Fill(jecuLo);
        hists.jecUncLoVsEtaPt->Fill(jeta, j.pt(), (double)jecuLo);
    }
} // HeavyNu::fillBasicJetHistos

TH1 * HeavyNu::bookRunHisto(uint32_t runNumber)
{
    std::string runstr = int2str<uint32_t > (runNumber);
    return hists.rundir->make <TH1I > (runstr.c_str(), runstr.c_str(), 1, 1, 2);
}

void HeavyNu::studyMuonSelectionEff(const std::vector< std::pair<pat::Muon,pat::GenericParticle> >& trackTP, 
				    const std::vector< std::pair<pat::Muon,pat::Muon> >& tightTP, 
				    const std::vector< std::pair<pat::Muon,pat::Muon> >& cjTP, 
				    const std::vector< pat::Muon >& tightMuons, 
                                    const edm::Handle<reco::TrackCollection>& gTracks,
                                    const edm::Handle<reco::BeamSpot>& beamspot, 
                                    double wgt, int nJets)
{
    
    reco::TrackCollection generalTracks = *(gTracks.product()); 

    // First case: Check muon ID
    // Probe is generic track meeting pT, eta requirements
    for (unsigned int i=0; i<trackTP.size(); i++) {
      pat::Muon            theTag   = trackTP.at(i).first ; 
      pat::GenericParticle theProbe = trackTP.at(i).second ; 
      // Calculate the isolation for the probe
      double trkSumPtIsoCone = 0. ; 
      for (unsigned int k=0; k<generalTracks.size(); k++) {
	reco::Track track = generalTracks.at(k) ; 
	double dR = deltaR(track.eta(), track.phi(), theProbe.eta(), theProbe.phi());
	if ( fabs(dR) > 0.3 || fabs(dR) < 0.01 ) continue ;
	if ( fabs(track.dxy(beamspot->position())) > 0.1 ) continue ;
	double dz = fabs( theProbe.vertex().Z() - track.vertex().Z() ) ;
	if ( dz > 0.2 ) continue ;
	trkSumPtIsoCone += track.pt() ;
      }
      if ( nJets >= 0 ) hists.TightTagTrackProbeInZwin0jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      if ( nJets >= 1 ) hists.TightTagTrackProbeInZwin1jet->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      if ( nJets >= 2 ) hists.TightTagTrackProbeInZwin2jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      // Does the probe pass tight ID selection?
      unsigned int muIdx = tightMuons.size() ;
      for (unsigned int k=0; k<tightMuons.size(); k++) {
	pat::Muon tightMuon = tightMuons.at(k) ; 
	if ( deltaR(tightMuon.eta(), tightMuon.phi(), theProbe.eta(), theProbe.phi()) < 0.02 &&
	     fabs(theProbe.pt()-tightMuon.pt())/theProbe.pt() < 0.05 ) {
	  muIdx = k ; break ; 
	}
      }
      if ( muIdx < tightMuons.size() ) {
          if ( nJets >= 0 ) hists.TightTagTrackProbePassesInZwin0jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
          if ( nJets >= 1 ) hists.TightTagTrackProbePassesInZwin1jet->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
          if ( nJets >= 2 ) hists.TightTagTrackProbePassesInZwin2jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      }
    }
        
    // Second case: Check isolation
    // Probe is tight muon separated from jets
    for (unsigned int i=0; i<tightTP.size(); i++) {
      pat::Muon theTag   = tightTP.at(i).first ; 
      pat::Muon theProbe = tightTP.at(i).second ; 
      // Confirmed: we have a valid probe
      if ( nJets >= 0 ) hists.TightTagTightProbeInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 1 ) hists.TightTagTightProbeInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 2 ) hists.TightTagTightProbeInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( hnu::muIsolation(theProbe,1.0) < cuts.muon_trackiso_limit ) { 
          if ( nJets >= 0 ) hists.TightTagTightProbePassesInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 1 ) hists.TightTagTightProbePassesInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 2 ) hists.TightTagTightProbePassesInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      }
    }

    // Third case:
    // Probe is tight muon with 0.5 < dR(mu,j) < 0.8
    for (unsigned int i=0; i<cjTP.size(); i++) {
      pat::Muon theTag   = cjTP.at(i).first ; 
      pat::Muon theProbe = cjTP.at(i).second ; 
      if ( nJets >= 0 ) hists.TightTagTightCJProbeInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 1 ) hists.TightTagTightCJProbeInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 2 ) hists.TightTagTightCJProbeInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( hnu::muIsolation(theProbe,1.0) < cuts.muon_trackiso_limit ) { 
          if ( nJets >= 0 ) hists.TightTagTightCJProbePassesInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 1 ) hists.TightTagTightCJProbePassesInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 2 ) hists.TightTagTightCJProbePassesInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      }
    }            
}


void HeavyNu::studyJetVertex(edm::Handle<pat::JetCollection>& pJets,
                             edm::Handle<reco::JPTJetCollection>& jptJets,
                             edm::Handle<pat::MuonCollection>& pMuons, int npue)
{
    if(pMuons->size() < 2) return;

    const pat::Muon& m1 = pMuons->at(0);
    const pat::Muon& m2 = pMuons->at(1);

    double mLL = (m1.p4() + m2.p4()).M();
    double dzll = m1.vz() - m2.vz();

    hists.dVzMuMus->Fill(dzll, npue);

    if(mLL < 71 || mLL > 111 || fabs(dzll) > 0.1) return;

    for(size_t iJet = 0; iJet < pJets->size(); iJet++)
    {
        pat::JetRef iJ = pat::JetRef(pJets, iJet);

        if(hnu::jetID(*iJ) > 1 && iJ->pt() > 40)
        {
            double jvz = hnu::caloJetVertex(*iJ, *jptJets);
            hists.dVzMuJets->Fill(m1.vz() - jvz, npue);
        }
    }
}

// ------------ method called to for each event  ------------

bool HeavyNu::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    HeavyNuEvent hnuEvent(analysisMode_);

    evtCounter++;

    hnuEvent.isMC = !iEvent.isRealData();
    hnuEvent.pfJets = isPFJets_;
    hnuEvent.scaleMuE(applyMESfactor_);

    if(iEvent.isRealData())
    {
        if((applyMESfactor_ != 1.0) || (applyJECUsign_ != 0.0))
            throw cms::Exception("Energy scale studies not allowed on data currently");

        uint32_t runn = iEvent.id().run();
        std::map<uint32_t, TH1 *>::const_iterator it = m_runHistos_.find(runn);
        TH1 *runh;
        if(it == m_runHistos_.end())
        {
            runh = bookRunHisto(runn);
            m_runHistos_[runn] = runh;
        }
        else
            runh = it->second;
        runh->Fill(1);
    }

    edm::Handle<reco::JPTJetCollection> jptJets;
    iEvent.getByLabel("JetPlusTrackZSPCorJetAntiKt5", jptJets); // Some day we should make this readable from the cfg file

    edm::Handle<pat::MuonCollection> pMuons;
    iEvent.getByLabel(muonTag_, pMuons);

    edm::Handle<pat::ElectronCollection> pElecs;
    iEvent.getByLabel(elecTag_, pElecs);

    edm::Handle<double> electronRhoHandle ; 
    iEvent.getByLabel(elecRhoTag_, electronRhoHandle) ; 
    elecRho_ = ((electronRhoHandle.isValid()) ? (*(electronRhoHandle.product())) : 0.) ; 

    edm::Handle<pat::JetCollection> pJets;
    iEvent.getByLabel(jetTag_, pJets);

    edm::Handle<reco::TrackCollection> gTracks ; 
    iEvent.getByLabel("generalTracks", gTracks) ; 
        
    edm::Handle<pat::GenericParticleCollection> pTracks;
    iEvent.getByLabel(trackTag_, pTracks);

    edm::Handle<pat::METCollection> pMET;
    iEvent.getByLabel(metTag_, pMET);

    edm::Handle<reco::MuonCollection> tevMuons;
    iEvent.getByLabel("refitMuons", tevMuons);

    //Shirpa reweighting info
    edm::Handle<GenEventInfoProduct> geneventinfo;
    iEvent.getByLabel("generator", geneventinfo);

    //Gen jets for gen matching
    edm::Handle<std::vector<reco::GenJet> > genjets;

    if(hnuEvent.isMC)
    {
        edm::Handle<std::vector<PileupSummaryInfo> > pPU;
        iEvent.getByLabel("addPileupInfo", pPU);
        std::pair<float, double> pileup = hnu::pileupReweighting(pPU, MCweightByVertex_);
        hnuEvent.n_pue = pileup.first ; // Will only be used for studies, thus no syst. correction necessary
        if(pileupEra_ < 20100) hnuEvent.eventWgt *= 1.0;
        else hnuEvent.eventWgt *= pileup.second;
        if ( fabs(puShift_) > 0.001 ) hnuEvent.eventWgt *= poissonNvtxShifter_.ShiftWeight( pileup.first ) ; 
        //Shirpa reweighting
        hnuEvent.eventWgt *= geneventinfo->weight();
	// PDF reweighting
#ifdef DO_LHAPDF
	if (doPDFreweight_) {
	  edm::Handle<GenEventInfoProduct> geip;
	  iEvent.getByLabel("generator",geip);
      
	  float Q=geip->pdf()->scalePDF;
	  int id1=geip->pdf()->id.first;
	  int id2=geip->pdf()->id.second;
	  float x1=geip->pdf()->x.first;
	  float x2=geip->pdf()->x.second;

	  hnuEvent.eventWgt *= hnu::getPDFWeight(Q,id1,x1,id2,x2,
                                                 doPDFreweight_,
                                                 pdfReweightBaseId,pdfReweightTargetId);
	}
#endif

        hists.weights->Fill(hnuEvent.eventWgt);
        // generator information
        edm::Handle<reco::GenParticleCollection> genInfo;
        if(iEvent.getByLabel("genParticles", genInfo))
        {
            hnuEvent.decayID(*genInfo);
        }
        iEvent.getByLabel("ak5GenJets", genjets);
    }
    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    hnuEvent.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle);

    if(!pElecs.isValid() || !pMuons.isValid() || !pJets.isValid() || !(pMET.isValid() && (pMET->size() > 0)))
    {
        std::cout << "Exiting as valid PAT objects not found" << std::endl;
        std::cout << "Electrons: " << pElecs.isValid() << std::endl;
        std::cout << "Muons:     " << pMuons.isValid() << std::endl;
        std::cout << "Jets:      " << pJets.isValid() << std::endl;
        std::cout << "MET:       " << pMET.isValid() << std::endl;
        return false;
    }

    if(firstEvent_)
    {
        // handle the jet corrector parameters collection,
        // get the jet corrector parameters collection from the global tag
        //
        edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
        iSetup.get<JetCorrectionsRecord > ().get("AK5Calo", JetCorParColl);

        // get the uncertainty parameters from the collection,
        // instantiate the jec uncertainty object
        if(applyJECUsign_)
        {
            JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
            jecuObj_ = new JetCorrectionUncertainty(JetCorPar);
        }

        // Some sanity checks inserted
        // If running on Monte Carlo, expect
        //   - pileup configuration to make sense
        //   - to apply trigger and ID corrections
        // If running on data
        //   - no corrections are applied
        // Of course, this changes completely when doing systematics checks
        if(studyMuonSelectionEff_) std::cout << "Histograms for studying muon reco/ID efficiency will be created" << std::endl;
        if(studyScaleFactorEvolution_) std::cout << "Histograms for Z scale factor cross checks will be created" << std::endl;
        if(applyJECUsign_)
        {
            std::cout << "Studies will be used to estimate JEC uncertainty" << std::endl;
            if(studyMuonSelectionEff_ || studyScaleFactorEvolution_)
                std::cout << "WARNING: You are performing studies with modified jets.  This is most likely wrong!" << std::endl;
        }
        else std::cout << "Nominal Jet corrections applied" << std::endl;
        if(applyMESfactor_ != 1.0)
        {
            std::cout << "Studies will be used to estimate MES uncertainty: " << applyMESfactor_ << std::endl;
            if(studyMuonSelectionEff_ || studyScaleFactorEvolution_)
                std::cout << "WARNING: You are performing studies with MES factor != 1.  This is most likely wrong!" << std::endl;
        }
        else std::cout << "No MES corrections applied" << std::endl;
        if(!disableTriggerCorrection_)
        {
            if(applyTrigEffsign_) std::cout << "Studies will be used to estimate trigger efficiency uncertainty" << std::endl;
            else std::cout << "Nominal trigger corrections applied" << std::endl;
        }

        if(hnuEvent.isMC)
        {
            int pileupYear = pileupEra_ ;
            int idYear = muid_->idEra();

            bool allErasMatch = true;
            if(applyMuIDCorrections_)
            {
                if(applyMuIDEffsign_) std::cout << "Studies will be used to estimate Mu ID uncertainty: " << applyMuIDEffsign_
                        << std::endl;
                else std::cout << "Nominal Mu ID corrections applied" << std::endl;
            }
            else std::cout << "You have disabled Mu ID corrections.  Are you sure?" << std::endl;
            if(disableTriggerCorrection_)
            {
                std::cout << "You have disabled the trigger correction.  Are you sure?" << std::endl;
                if(pileupYear != idYear) allErasMatch = false;
            }
            else
            {
                allErasMatch = (pileupYear == idYear);
            }
            if(!allErasMatch)
            {
                std::cout << "WARNING: You do not appear to have consistent corrections applied!" << std::endl;
                std::cout << "         pileup era is " << pileupEra_ << ", year for mu ID is " << idYear
                        << std::endl;
            }
            else
            {
                std::cout << "Looking at corrections, I assume you are running with the " << pileupYear << " year settings" << std::endl;
            }
        }
        else
        {
            if(disableTriggerCorrection_)
                std::cout << "WARNING: You have disabled the trigger correction in data.  What are you doing?!?" << std::endl;
        }
        std::cout << "===============================================" << std::endl;

        firstEvent_ = false;
    }

    hists.mc_type->Fill(hnuEvent.mc_class,hnuEvent.eventWgt);
    hists.nelec->Fill(pElecs->size());
    hists.nmuAll->Fill(pMuons->size());
    hists.njet->Fill(pJets->size());
    hists.nmet->Fill(pMET->size());

    studyJetVertex(pJets, jptJets, pMuons, hnuEvent.n_primary_vertex);

    if(pMET->size())
        hists.met->Fill(pMET->at(0).pt());
    else
        hists.met->Fill(0);

    for(size_t iJet = 0; iJet < pJets->size(); iJet++)
    {
        pat::JetRef iJ = pat::JetRef(pJets, iJet);
        fillBasicJetHistos(*iJ, iJet + 1);
    }

    int nloose = 0, ntight = 0;
    for(size_t iMuon = 0; iMuon < pMuons->size(); iMuon++)
    {
        pat::MuonRef iM = pat::MuonRef(pMuons, iMuon);
        if(!iM.isAvailable()) continue;
        fillBasicMuHistos(*iM);
        if(hnu::isVBTFloose(*iM))
        {
            nloose++;
            if(hnu::isVBTFtight(*iM))
                ntight++;
        }
    }
    hists.nmuLoose->Fill(nloose);
    hists.nmuTight->Fill(ntight);

    hists.cutlevel->Fill(-1.0, hnuEvent.eventWgt);

    // Basic selection requirements: Require at least two leptons, two jets
    switch(analysisMode_)
    {
        case HeavyNuEvent::HNUMU:
            if(pMuons->size() < 2) return false;
            break;
        case HeavyNuEvent::HNUE:
            if(pElecs->size() < 2) return false;
            break;
        case HeavyNuEvent::TOP:
            if(pElecs->size() < 1 || pMuons->size() < 1) return false;
            break;
        case HeavyNuEvent::QCD:
        case HeavyNuEvent::CLO:
            break;
    }
    if(pJets->size() >= 2)
    {
        hists.cutlevel->Fill(0.0, hnuEvent.eventWgt);
        fill(*pMuons, *pElecs, *pJets, *pMET, hnuEvent.isMC, hnuEvent.eventWgt, isPFJets_, hnuEvent.n_pue, hnuEvent.n_primary_vertex, hists.noCuts);
    }
    else return false;

    // Look for valid jets and put them in the event
    std::vector< std::pair<pat::Jet, float> > jetCands =
      hnu::getJetList(pJets, jecuObj_, cuts.minimum_jet_pt, cuts.maximum_jet_abseta, applyJECUsign_, jecVal_, hnuEvent.isMC, applyJERsign_);
    hnuEvent.nJets = jetCands.size();

    // Look for valid muons
    std::vector<pat::Muon> muCands =
      hnu::getMuonList(pMuons, tevMuons, pvHandle, (int(muid_->idEra()/10)), 
		       cuts.minimum_mu2_pt, cuts.maximum_mu_abseta, applyMESfactor_, merUncertainty_, false);

    // Look for valid electrons
    std::vector< std::pair<pat::Electron, float> > eCands =
      hnu::getElectronList(pElecs, cuts.maximum_mu_abseta, cuts.minimum_mu2_pt, cuts.minimum_mu2_pt, 
			   heepVersion_,elecRho_);

    if (analysisMode_ == HeavyNuEvent::HNUE)
      std::cout << "I have " << muCands.size() << " muons and " << eCands.size() << " electrons" << std::endl ; 

    // In order to avoid bias, it is necessary to perform muon studies
    // immediately after first creating the muon list
    //bool debuggingEvents = false ;
    if (studyMuonSelectionEff_ && analysisMode_ == HeavyNuEvent::HNUMU) {

        std::vector<pat::Muon> tagMuons ;
        for (unsigned int i=0; i<muCands.size(); i++) {
          // Muons already tight, in proper detector region, with sufficient pT
          pat::Muon tagCand = muCands.at(i) ;
          if (hnu::muIsolation(tagCand) < cuts.muon_trackiso_limit) {
            if ( (trig_->matchingEnabled() && iEvent.isRealData() && trig_->isTriggerMatched(tagCand, iEvent)) ||
                 (!iEvent.isRealData() && trig_->simulateForMC(tagCand.pt(), tagCand.eta())) ) {
              bool jetOverlap = false ;
              for (unsigned int j=0; j<jetCands.size(); j++) {
                pat::Jet jet = jetCands.at(j).first ;
                if ( hnu::jetID(jet) > 0 ) {
                  double dR = deltaR(tagCand.eta(), tagCand.phi(), jet.eta(), jet.phi()) ;
                  if (dR <= cuts.minimum_muon_jet_dR) { jetOverlap = true ; break ; }
                }
              }
              // If muon is separated from jets, is isolated, and trigger matched, it is a tag
              if ( !jetOverlap ) tagMuons.push_back( tagCand ) ;
            }
          }
        }

        edm::Handle<reco::BeamSpot> beamSpotHandle;
        if (!iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle))
            throw cms::Exception("Trying to do efficiency studies, cannot find beam spot");

        // Make sure that we have at least two jets
        // keeping in mind that we only study efficiency for the nominal case
        std::vector<pat::Jet> validJets ; 
        for (unsigned int i=0; i<jetCands.size(); i++) { 
            pat::Jet jet = jetCands.at(i).first ;
            if ( hnu::jetID(jet) > 0 && jet.pt() > cuts.minimum_jet_pt ) validJets.push_back( jet ) ;
        }
        // if ( validJets.size() >= 2 && tagMuons.size() > 0 ) {
        if ( tagMuons.size() > 0 ) { // Allow for checks of ID efficiency for different jet multiplicity
	  //--- Create the two lists of probe muons ---//
	  std::vector<pat::Muon> tightMuonProbes ; 
	  std::vector<pat::Muon> cjMuonProbes ; 
	  std::vector<pat::Muon> trigProbes ; 
	  for (unsigned int i=0; i<muCands.size(); i++) {
	    bool jetOverlap = false ; 
	    double minDeltaR = 0.8 ; 
	    for (unsigned int j=0; j<validJets.size(); j++) { 
	      pat::Jet jet = validJets.at(j) ; 
	      double dR = deltaR(muCands.at(i).eta(), muCands.at(i).phi(), jet.eta(), jet.phi()) ;
	      if (dR < minDeltaR) minDeltaR = dR ; 
	      if (dR <= cuts.minimum_muon_jet_dR) { jetOverlap = true ; break ; } 
	    }
	    if ( !jetOverlap ) {
	      tightMuonProbes.push_back( muCands.at(i) ) ; 
	      if ( minDeltaR < 0.8 ) cjMuonProbes.push_back( muCands.at(i) ) ; 
	      if ( fabs(muCands.at(i).eta()) < 2.1 && 
		   hnu::muIsolation(muCands.at(i)) < cuts.muon_trackiso_limit ) 
		trigProbes.push_back( muCands.at(i) ) ; 
	    }
	  }
	  std::vector< std::pair<pat::Muon,pat::Muon> > tagTightProbes = 
              hnu::getTagProbePair<pat::Muon>( tagMuons,tightMuonProbes,ZwinMinGeV_,ZwinMaxGeV_,
                                               ((oneTP_)?(tpRandom_->Uniform()):(-1.0)) ) ; 
	  std::vector< std::pair<pat::Muon,pat::Muon> > tagCJmuonProbes = 
              hnu::getTagProbePair<pat::Muon>( tagMuons,cjMuonProbes,ZwinMinGeV_,ZwinMaxGeV_,
                                               ((oneTP_)?(tpRandom_->Uniform()):(-1.0)),false ) ; 
	  std::vector< std::pair<pat::Muon,pat::Muon> > tagTrigProbes = 
              hnu::getTagProbePair<pat::Muon>( tagMuons,trigProbes,ZwinMinGeV_,ZwinMaxGeV_,
                                               ((oneTP_)?(tpRandom_->Uniform()):(-1.0)) ) ; 
	  
	  std::vector<pat::GenericParticle> trackProbes ; 
	  pat::GenericParticleCollection trackCands = *(pTracks.product());
	  std::sort(trackCands.begin(), trackCands.end(), hnu::pTcompare());
	  for (unsigned int i=0; i<trackCands.size(); i++) {
	    pat::GenericParticle trackCand = trackCands.at(i) ; 
            if ( trackCand.pt() <= cuts.minimum_mu2_pt ) break ; // Sorted collection: quit once below
            if ( fabs(trackCand.eta()) >= cuts.maximum_mu_abseta ) continue ;
	    trackProbes.push_back( trackCand ) ; 
	  } 
	  std::vector< std::pair<pat::Muon,pat::GenericParticle> > tagTrackProbes = 
              hnu::getTagProbePair<pat::GenericParticle>( tagMuons,trackProbes,ZwinMinGeV_,ZwinMaxGeV_,
                                                          ((oneTP_)?(tpRandom_->Uniform()):(-1.0)) ) ; 

	  //if ( tagTightProbes.size() > 0 || tagTrackProbes.size() > 0 ) debuggingEvents = true ;

	  // Defined both tag+probe collections
	  studyMuonSelectionEff(tagTrackProbes,tagTightProbes,tagCJmuonProbes, 
				muCands,gTracks,beamSpotHandle,hnuEvent.eventWgt,validJets.size());

	  // Study Trigger Matching efficiency for data only...here 2 or more jets is enforced
	  if ( trig_->matchingEnabled() && iEvent.isRealData() && validJets.size() >= 2 ) {
	    for (unsigned int i=0; i<tagTrigProbes.size(); i++) { 
	      pat::Muon theTag   = tagTrigProbes.at(i).first ; 
	      pat::Muon theProbe = tagTrigProbes.at(i).second ; 
	      hists.TightTagTrigProbeInZwin->tapfill( theTag,theProbe,theProbe.trackIso(),hnuEvent.eventWgt ) ;
	      if ( trig_->isTriggerMatched(theProbe, iEvent) )
		hists.TightTagTrigProbePassesInZwin->tapfill( theTag,theProbe,theProbe.trackIso(),hnuEvent.eventWgt ) ;
	    }
	  }
	}
    }

    if (analysisMode_ == HeavyNuEvent::HNUE)
      std::cout << "I have " << hnuEvent.nJets << " jets" << std::endl ; 

    if(hnuEvent.nJets < 2) return false;

    hnuEvent.j1 = jetCands.at(0).first;
    hnuEvent.j2 = jetCands.at(1).first;
    hnuEvent.j1scale = jetCands.at(0).second;
    hnuEvent.j2scale = jetCands.at(1).second;

    hnuEvent.btagName = btagName;
    hnuEvent.isBJet1 = hnuEvent.j1.bDiscriminator(btagName) >= minBtagDiscVal;
    hnuEvent.isBJet2 = hnuEvent.j2.bDiscriminator(btagName) >= minBtagDiscVal;

    //conduct gen study on jets
    if(hnuEvent.isMC)
    {
        //here we match the reco jets to gen jets
        double mindR1 = 100, mindR2 = 100;
        for(std::vector<reco::GenJet>::const_iterator igj = genjets->begin(); igj != genjets->end(); ++igj)
        {
            double dR1 = deltaR(hnuEvent.j1.p4(), igj->p4());
            double dR2 = deltaR(hnuEvent.j2.p4(), igj->p4());

            if(dR1 < mindR1)
            {
                mindR1 = dR1;
                hnuEvent.gj1 = *igj;
            }
            if(dR2 < mindR2)
            {
                mindR2 = dR2;
                hnuEvent.gj2 = *igj;
            }
        }

        //after finding the matching gen jets we track their parentage and try to match them to a Nu_mu
        bool gmj1 = false, gmj2 = false;
        std::vector<const reco::GenParticle*> mothers = hnuEvent.gj1.getGenConstituents();
        for(std::vector<const reco::GenParticle*>::const_iterator iM = mothers.begin(); iM != mothers.end(); ++iM)
        {
            if((gmj1 |= isWrDaughter(*iM))) break;
            //for(size_t i = 0; i < (*iM)->numberOfMothers(); i++)
            //{
            //    if(gmj1 |= isWrDaughter((*iM)->mother(i))) break;
            //}
            //if(gmj1) break;
        }
        mothers = hnuEvent.gj2.getGenConstituents();
        for(std::vector<const reco::GenParticle*>::const_iterator iM = mothers.begin(); iM != mothers.end(); ++iM)
        {
            if((gmj2 |= isWrDaughter(*iM))) break;
            //for(size_t i = 0; i < (*iM)->numberOfMothers(); i++)
            //{
            //    if(gmj2 |= isWrDaughter((*iM)->mother(i))) break;
            //}
            //if(gmj2) break;
        }
        hnuEvent.numNuLJetsMatched = (int)gmj1 + (int)gmj2;
    }

    hnuEvent.tjV1 = hnu::caloJetVertex(hnuEvent.j1, *jptJets);
    hnuEvent.tjV2 = hnu::caloJetVertex(hnuEvent.j2, *jptJets);

    bool l1trig = false;
    bool l2trig = false;

    if(analysisMode_ == HeavyNuEvent::HNUMU)
    {
        for(unsigned int i = 0; i < muCands.size(); i++)
        {
            if(hnuEvent.nLeptons == 2) break;
            pat::Muon iM = muCands.at(i);
            if(hnu::muIsolation(iM) < cuts.muon_trackiso_limit)
            {
                double dRj1 = deltaR(iM.eta(), iM.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
                double dRj2 = deltaR(iM.eta(), iM.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
                if(dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR)
                {
                    hnuEvent.nLeptons++;
                    if(hnuEvent.nLeptons == 1) hnuEvent.mu1 = iM;
                    else if(hnuEvent.nLeptons == 2) hnuEvent.mu2 = iM;
                    else std::cout << "WARNING: Expected empty muon position" << std::endl;
                }
            }
        }

        if(applyMuIDCorrections_ && hnuEvent.isMC)
        {
            double mu1wgt = (hnuEvent.nLeptons > 0)? (muid_->weightForMC((hnuEvent.mu1.pt()), applyMuIDEffsign_)):1.0;
            double mu2wgt = (hnuEvent.nLeptons > 1)? (muid_->weightForMC((hnuEvent.mu2.pt()), applyMuIDEffsign_)):1.0;

            hnuEvent.eventWgt *= (mu1wgt * mu2wgt);
        }

        //--- Trigger Matching needed for efficiency studies ---//
        if(trig_->matchingEnabled() && iEvent.isRealData())
        {
            l1trig = (hnuEvent.nLeptons > 0) &&
                    trig_->isTriggerMatched(hnuEvent.mu1, iEvent, &(((HeavyNuMuHist*)hists.Mu1TrigMatchesInZwin)->trigHistos));
            l2trig = (hnuEvent.nLeptons > 1) &&
                    trig_->isTriggerMatched(hnuEvent.mu2, iEvent, &(((HeavyNuMuHist*)hists.Mu2TrigMatchesInZwin)->trigHistos));
        }
        else if(!iEvent.isRealData())
        {
            if(disableTriggerCorrection_)
            {
                l1trig = true;
                l2trig = true;
            }
            else
            {
                l1trig = (hnuEvent.nLeptons > 0) && trig_->simulateForMC(hnuEvent.mu1.pt(), hnuEvent.mu1.eta(), applyTrigEffsign_) && hnuEvent.mu1.pt() > 40;
                l2trig = (hnuEvent.nLeptons > 1) && trig_->simulateForMC(hnuEvent.mu2.pt(), hnuEvent.mu2.eta(), applyTrigEffsign_) && hnuEvent.mu2.pt() > 40;
            }
        }
        // std::cout << "Trigger results: " << mu1trig << ", " << mu2trig << std::endl ;
    }
    else if(analysisMode_ == HeavyNuEvent::HNUE)
    {
      std::cout << "Examining electron trigger" << std::endl ; 
        for(unsigned int i = 0; i < eCands.size(); i++)
        {
            if(hnuEvent.nLeptons == 2) break;
            pat::Electron iE = eCands.at(i).first;
            double dRj1 = deltaR(iE.eta(), iE.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
            double dRj2 = deltaR(iE.eta(), iE.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
            if(dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR)
            {
                hnuEvent.nLeptons++;
                if(hnuEvent.nLeptons == 1) hnuEvent.e1 = iE;
                else if(hnuEvent.nLeptons == 2) hnuEvent.e2 = iE;
                else std::cout << "WARNING: Expected empty electron position" << std::endl;
            }
        }

	// bool l12trig = (hnuEvent.nLeptons > 1) &&
	//   trig_->isTriggerMatched(hnuEvent.e1, hnuEvent.e2, iEvent, &(((HeavyNuMuHist*)hists.Mu2TrigMatchesInZwin)->trigHistos));
	l1trig = l2trig = true ;
    }

    hnuEvent.regularize();

    if(hnuEvent.nLeptons < 2) return false;

    if(hnu::jetID(hnuEvent.j1) < 1 || hnu::jetID(hnuEvent.j2) < 1) return false;

    hists.cutlevel->Fill(1.0, hnuEvent.eventWgt); // Two highest pT muons that are isolated, separated from chosen jets
    hnuEvent.scaleMuE(applyMESfactor_);
    hnuEvent.calculate(); // calculate various details
    hists.LLJJptCuts->fill(hnuEvent);
    if(pMET->size()) hnuEvent.met1 = pMET->at(0);

    //--- Trigger code needs to be updated...placeholder for now ---//
    if(!l1trig && !l2trig) return false;
    hists.cutlevel->Fill(2.0, hnuEvent.eventWgt); // Event meets trigger requirements
    hists.TrigMatches->fill(hnuEvent);

    // hists.LLJJptCuts->fill(hnuEvent);

    //--- Impose vertex requirement here ---//
    float deltaVzJ1J2 = fabs(hnuEvent.tjV1 - hnuEvent.tjV2);
    float deltaVzJ1M1 = fabs(hnuEvent.tjV1 - hnuEvent.l1->vertex().Z());
    float deltaVzJ2M2 = fabs(hnuEvent.tjV2 - hnuEvent.l2->vertex().Z());
    float deltaVzJ1M2 = fabs(hnuEvent.tjV1 - hnuEvent.l2->vertex().Z());
    float deltaVzJ2M1 = fabs(hnuEvent.tjV2 - hnuEvent.l1->vertex().Z());
    if( (cuts.maxJetVZsepCM > 0) && 
        ((deltaVzJ1J2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M1 >= cuts.maxJetVZsepCM) ||
         (deltaVzJ2M2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M2 >= cuts.maxJetVZsepCM) ||
         (deltaVzJ2M1 >= cuts.maxJetVZsepCM)) )
        return false;
    float deltaVzM1M2 = fabs(hnuEvent.l1->vertex().Z() - hnuEvent.l2->vertex().Z());
    if(cuts.maxVertexZsep > 0 && deltaVzM1M2 >= cuts.maxVertexZsep) return false;

    hists.cutlevel->Fill(3.0, hnuEvent.eventWgt); // Event meets vertex requirements
    //hists.VertexCuts->fill(hnuEvent);

    //--- The "basic" object, trigger, and (possibly) vertex requirements should be done ---//
    //--- Consider alternative selection requirements ---//
    if ( studyAlternativeSelection_ ) {
        double l1pt = hnuEvent.l1->pt() ;
        double l2pt = hnuEvent.l2->pt() ;
        //double j1pt  = hnuEvent.j1scale * hnuEvent.j1.pt() ;
        double j2pt  = hnuEvent.j2.pt() ; 
        if ( hnuEvent.mLL >= cuts.minimum_mumu_mass ) { // Standard dimuon requirement
            if ( l1pt >= 40. ) hists.AlternativeMu1Pt40->fill(hnuEvent) ;
            if ( l1pt >=  cuts.minimum_mu1_pt ) { // Standard mu1 pT requirement
                if ( l2pt >= 40. ) {
                    hists.AlternativeMu2Pt40->fill(hnuEvent) ;
                    if ( l1pt >= 80. ) hists.AlternativeElecChanPt->fill(hnuEvent) ;
                }
                if ( l2pt >= 60. ) hists.AlternativeMu2Pt60->fill(hnuEvent) ;
                if ( j2pt >= 60. )  hists.AlternativeJetPt60->fill(hnuEvent) ;
                //--- Requirement for at least one muon to be barrel ---//
                bool atLeastOneBarrelMuonLoose = ( fabs(hnuEvent.l1->eta()) < 1.2 || fabs(hnuEvent.l2->eta()) < 1.2 ) ;
                bool atLeastOneBarrelMuonTight = ( fabs(hnuEvent.l1->eta()) < 0.8 || fabs(hnuEvent.l2->eta()) < 0.8 ) ;
                if ( atLeastOneBarrelMuonLoose )
                {
                    hists.AlternativeBarrelLoose->fill(hnuEvent) ;
                    if ( atLeastOneBarrelMuonTight ) hists.AlternativeBarrelTight->fill(hnuEvent) ;
                }
                //--- Checking for b-tagged jets using TCHE Loose ---//
                int nBtags = 0 ;
                if ( hnuEvent.j1.bDiscriminator(btagName) >= minBtagDiscVal ) nBtags++ ;  
                if ( hnuEvent.j2.bDiscriminator(btagName) >= minBtagDiscVal ) nBtags++ ;  
                if ( nBtags > 0 )  hists.AlternativeAtLeastOneBjet->fill(hnuEvent) ;
                if ( nBtags == 2 ) hists.AlternativeTwoBjets->fill(hnuEvent) ;
            }
        }
        if ( hnuEvent.mLL >= 120. ) {
            if ( l1pt >=  cuts.minimum_mu1_pt ) { // Standard mu1 pT requirement
                hists.AlternativeDimuonMass120->fill(hnuEvent) ;
            }
        }
    }
    
    if(studyScaleFactorEvolution_)
    {
        double l1pt = hnuEvent.l1->pt();
        if(l1pt > 40.) hists.Mu1Pt40GeVCut->fill(hnuEvent);
        if(l1pt > 50.) hists.Mu1Pt50GeVCut->fill(hnuEvent);
        if(l1pt > 60.) hists.Mu1Pt60GeVCut->fill(hnuEvent);
        if(l1pt > 80.) hists.Mu1Pt80GeVCut->fill(hnuEvent);
        if(l1pt > 100.) hists.Mu1Pt100GeVCut->fill(hnuEvent);
    }

    if(hnuEvent.l1->pt() < cuts.minimum_mu1_pt)
        return false;

    if(studyScaleFactorEvolution_)
    {
        if(hnuEvent.n_primary_vertex == 1)
            hists.Mu1HighPtCutVtxEq1->fill(hnuEvent);
        else if(hnuEvent.n_primary_vertex <= 5)
            hists.Mu1HighPtCutVtx2to5->fill(hnuEvent);
        else if(hnuEvent.n_primary_vertex > 5)
            hists.Mu1HighPtCutVtxGt5->fill(hnuEvent);
    }
    hists.cutlevel->Fill(4.0, hnuEvent.eventWgt); // Event meets high muon pT requirements
    hists.Mu1HighPtCut->fill(hnuEvent);

    if(hnuEvent.mLL < 40) return false; // Sanity check...remove low mass points
    //hists.loDiLmassCut->fill(hnuEvent);
    if ( studyRatePerRun_ && inZmassWindow(hnuEvent.mLL) )
        hists.z2jetPerRun->Fill( iEvent.id().run() ) ; 
    
    if(hnuEvent.mLL < cuts.minimum_mumu_mass) return false; // dimuon mass cut
    hists.cutlevel->Fill(5.0, hnuEvent.eventWgt); // Event meets dimuon mass requirements
    hists.diLmassCut->fill(hnuEvent);

    if(iEvent.isRealData())
    {
        bool mu1posChg = (hnuEvent.mu1.charge() > 0) ; 
        bool mu2posChg = (hnuEvent.mu2.charge() > 0) ;
        if ( hnuEvent.mu1.charge() == 0 || hnuEvent.mu2.charge() == 0 )
            std::cout << "WARNING: found muon with zero charge" << std::endl ; 

        std::cout << "\t" << iEvent.id() << std::endl;
        std::cout << "\tM(W_R)  = " << hnuEvent.mWR << " GeV";
        std::cout << ", M(NuR1) = " << hnuEvent.mNuR1 << " GeV";
        std::cout << ", M(NuR2) = " << hnuEvent.mNuR2 << " GeV" << std::endl;
        std::cout << "\tM(mumu) = " << hnuEvent.mLL << " GeV";
        std::cout << ", M(JJ) = " << hnuEvent.mJJ << " GeV" << std::endl;
        std::cout << "\tJets:   j1 ";
        outputCandidate(hnuEvent.j1);
        std::cout << ", j2 ";
        outputCandidate(hnuEvent.j2);
        std::cout << std::endl;
        std::cout << "\tMuons: mu1, mu" << (mu1posChg ? "+":"-");
        outputCandidate(*(hnuEvent.l1));
        std::cout << ", mu2, mu" << (mu2posChg ? "+":"-");
        outputCandidate(*(hnuEvent.l2));
        std::cout << std::endl;
    }

    // Change the final logic of the filter based on LQ meeting discussion:
    // Interest in seeing events that pass the dilepton mass requirement
    // if ( hnuEvent.mWR<cuts.minimum_mWR_mass ) return false;  // 4-object mass cut
    if(hnuEvent.mWR >= cuts.minimum_mWR_mass)
    {
        hists.cutlevel->Fill(6.0, hnuEvent.eventWgt); // Event meets W_R mass requirements
        hists.mWRmassCut->fill(hnuEvent);
    }
    return true;
}

// #include "LHAPDF/LHAPDF.h"

// double HeavyNu::getPDFWeight(float Q, int id1, float x1, int id2, float x2) {

  
//   if (!doPDFreweight_) return 1.0;
  
//   LHAPDF::usePDFMember(1,pdfReweightBaseId);
//   double pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
//   double pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
  
//   LHAPDF::usePDFMember(2,pdfReweightTargetId);
//   double newpdf1 = LHAPDF::xfx(2, x1, Q, id1)/x1;
//   double newpdf2 = LHAPDF::xfx(2, x2, Q, id2)/x2;
  
//   double w=(newpdf1/pdf1*newpdf2/pdf2);
  
//   //  printf("My weight is %f\n",w);
  
//   return w;
  
// }


// ------------ method called once each job just before starting event loop  ------------

void HeavyNu::beginJob()
{
  //    nnif_->beginJob();
    firstEvent_ = true;
    evtCounter = 0;

    if (doPDFreweight_) {
      hnu::initPDFSet(1,pdfReweightBaseName);
      hnu::initPDFSet(2,pdfReweightTargetName);
    }

}

// ------------ method called once each job just after ending the event loop  ------------

void HeavyNu::endJob()
{
  //    nnif_->endJob();
    trig_->endJob();
    muid_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNu);



