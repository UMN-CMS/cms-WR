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
// $Id: HeavyNu.cc,v 1.73 2011/10/14 23:01:07 bdahmes Exp $
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
#include "TProfile2D.h"
#include "TVector3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNu_NNIF.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


//////////////////////////////////////////////////////////////////
// generic maximum/minimum

template <class T> const T& max(const T& a, const T& b)
{
    return(b < a)?a:b;
}

template <class T> const T& min(const T& a, const T& b)
{
    return(b < a)?b:a;
}

template <class T>
inline std::string int2str(T i)
{
    std::ostringstream ss;
    ss << i;
    return ss.str();
}

//============================================================

// class compare
// {
// public:

//     template <class T> bool operator() (const T& a, const T& b)
//     {
//         return a.pt() > b.pt();
//     }
// };

//============================================================

template <class T> void outputCandidate(const T& p)
{
    std::cout << "pt=" << p.pt() << " GeV, eta=" << p.eta() << ", phi=" << p.phi();
}

//============================================================

// inline void dumpJetCorInfo(const pat::Jet& j)
// {
//     std::cout << "Available JEC sets and levels:\n";
//     const std::vector<std::string> jecsets = j.availableJECSets();
//     for(size_t i = 0; i < jecsets.size(); i++)
//     {
//         std::cout << jecsets[i] << ":";
//         const std::vector<std::string> jeclevs = j.availableJECLevels(i);
//         for(size_t j = 0; j < jeclevs.size(); j++)
//             std::cout << " " << jeclevs[i];
//         std::cout << std::endl;
//     }
//     std::cout << "current JEC set   : " << j.currentJECSet() << std::endl;
//     std::cout << "current JEC level : " << j.currentJECLevel() << std::endl;
//     std::cout << "current JEC flavor: " << j.currentJECFlavor() << std::endl;
// }

//============================================================

// Returns 0=neither, 1=loose or 2=tight, -1 if tight but not loose (!)

// int jetID(const pat::Jet& j)
// {
//     if(!j.isPFJet())
//     {
//         JetIDSelectionFunctor jetIDloose(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE);
//         JetIDSelectionFunctor jetIDtight(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::TIGHT);

//         pat::strbitset ret = jetIDloose.getBitTemplate();
//         ret.set(false);
//         bool loose = jetIDloose(j, ret);
//         ret.set(false);
//         bool tight = jetIDtight(j, ret);
//         return(tight?(loose?2:-1):(loose?1:0));
//     }
//     else return 2;
// }

inline void labelJetIDaxis(TAxis *ax)
{
    ax->SetBinLabel(1, "Neither");
    ax->SetBinLabel(2, "PURE09 Loose");
    ax->SetBinLabel(3, "PURE09 Tight");
}

//============================================================
// From the JetEnergyScale twiki:
//

inline float sqr(float x)
{
    return x*x;
}
//const float unc4bjetScale2 = sqr(0.05f);
//const float unc4pileup2    = sqr(0.01f);
const float unc4calib2 = sqr(0.015f);
const float unc4ecalCalib2 = sqr(0.025f);
const float pileup2010A = 1.2;
const float pileup2010B = 2.2;
const float pileup2011A = 5.0;
//const float totalOtherUnc2 = unc4bjetScale2+unc4pileup2+unc4calib2;

float jecTotalUncertainty(float jpt, float jeta,
        JetCorrectionUncertainty *jecuObj,
        int correctEra,
        bool isBjet,
        bool directionIsUp)
{
    float offunc; // the "official" eta/pt dependent uncertainty

    jecuObj->setJetPt((float)jpt);
    jecuObj->setJetEta((float)jeta);
    offunc = jecuObj->getUncertainty(directionIsUp);
    //return offunc ; // According to twiki, all uncertainties accounted for in GT

    float pileupCorrection = (2 * 0.75 * 0.8) / jpt;
    if(correctEra == 1) pileupCorrection *= pileup2010A;
    else if(correctEra == 2) pileupCorrection *= pileup2010B;
    else if(correctEra == 3) pileupCorrection *= pileup2011A;
    else
    { // Merge 2010A and 2010B corrections
        float pileup2010 = ((3.18 * pileup2010A) + (32.96 * pileup2010B)) / 36.14;
        pileupCorrection *= pileup2010;
    }
    //float unc4pileup2 = sqr( pileupCorrection ) ;

    //float uncBjetScale2 = 0. ;
    //if ( isBjet ) {
    //  if ( jpt > 50. && jpt < 200. && fabs(jeta) < 2.0 ) uncBjetScale2 = sqr( 0.04 ) ;
    //  else uncBjetScale2 = sqr( 0.06 ) ;
    //}

    // Calculations taken from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC, version 23
    // float totalUnc = sqrt( sqr(offunc) + unc4calib2 + unc4pileup2 + uncBjetScale2 ) ;
    float totalUnc = sqrt(sqr(offunc) + sqr(pileupCorrection) + 0.025 * 0.025);
    return totalUnc;
}

//============================================================


static std::string btagName;
static double minBtagDiscVal; // for discriminating B-tagged jets.

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
    // virtual void selectJets(edm::Handle<pat::JetCollection>& pJets,
    //         HeavyNuEvent& hne);
    // virtual bool muPassesSelection(const pat::Muon& m,
    //         const HeavyNuEvent& hne);
    // virtual void selectMuons(edm::Handle<pat::MuonCollection>& pMuons,
    //         HeavyNuEvent& hne);
    virtual TH1 *bookRunHisto(uint32_t runNumber);

    virtual void studyMuonSelectionEff(std::vector<pat::Muon>& tightMuons,
                                       std::vector<pat::Jet>& jets,
                                       edm::Handle<pat::GenericParticleCollection>& pTracks,
                                       edm::Handle<reco::TrackCollection>& gTracks,
                                       edm::Handle<reco::BeamSpot>& beamspot, 
                                       std::vector<bool> isTag,
                                       double wgt);
    virtual void studyIsolation(const std::vector<pat::Muon>& muons,
                                const std::vector< std::pair<pat::Jet, float> >& jets,
                                bool mu1tag, bool mu2tag, double weight);

    virtual void studyJetVertex(edm::Handle<pat::JetCollection>& pJets,
                                edm::Handle<reco::JPTJetCollection>& jptJets,
                                edm::Handle<pat::MuonCollection>& pMuons, int npue);

    bool passesTrigger(const double mu1pt, const double mu2pt,
                       const bool mu1trig, const bool mu2trig,
                       const uint32_t run);

    inline bool inZmassWindow(double mMuMu)
    {
        return(mMuMu <= ZwinMaxGeV_) && (mMuMu >= ZwinMinGeV_);
    }

    edm::InputTag muonTag_;
    edm::InputTag trackTag_;
    edm::InputTag jetTag_;
    edm::InputTag metTag_;
    edm::InputTag elecTag_;

    int evtCounter;

    double ZwinMinGeV_, ZwinMaxGeV_; // for trigger efficiency studies

    int jecVal_; // Jet correction eras: 0 = 2010(A+B), 1 = 2010A, 2 = 2010B, 3 = 2011A
    int applyJECUsign_; // for Jet Energy Correction Uncertainty studies
    double applyMESfactor_; // for Muon Energy Scale studies
    int applyTrigEffsign_; // for Trigger Efficiency studies
    bool highestPtTriggerOnly_;
    bool applyMuIDCorrections_;
    int applyMuIDEffsign_; // for Loose+Iso efficiency correction
    bool studyMuonSelectionEff_; // for Muon Loose ID Efficiency studies
    bool studyScaleFactorEvolution_; // for Top, Z+jets scale factors (by Mu1 pT) studies

    bool disableTriggerCorrection_; // turn off the trigger emulator
    int pileupEra_;
    bool isPFJets_; // true if PFJets are used (turns off jet id requirement)

    std::string currentFile_;
    bool dolog_;
    bool firstEvent_;

    HeavyNu_NNIF *nnif_;
    HeavyNuTrigger *trig_;
    HeavyNuID *muid_;
    JetCorrectionUncertainty *jecuObj_;

    edm::LumiReWeighting MCweightByVertex_;

    std::map<uint32_t, TH1 *> m_runHistos_;

    // ----------member data ---------------------------

    struct HistPerDef
    {
        //book histogram set w/ common suffix inside the provided TFileDirectory
        void book(TFileDirectory *, const std::string&, const std::vector<hNuMassHypothesis>&);
        void bookTagProbe(TFileDirectory *, const std::string&);
        // fill all histos of the set with the two electron candidates
        void fill(pat::MuonCollection muons,
                  pat::JetCollection jets,
                  pat::METCollection metc,
                  bool isMC,
                  double wgt,
                  bool pfJets);
        // fill all histos of the set with the two electron candidates
        void fill(const HeavyNuEvent& hne, const std::vector<hNuMassHypothesis>&);
        // Special fill for muon efficiency studies
        void fill(const pat::Muon& theTag, const pat::Muon& theProbe, const double wgt) ; 
        void fill(const pat::Muon& theTag, const pat::GenericParticle& theProbe, const double wgt) ; 
        
        TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2;
        TH2 *ptMu1VsPtMu2ss, *ptMu1VsPtMu2os;
        TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2;
        TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2;
        TH1 *dEtaMu, *dPhiMu, *dEtaJet, *dPhiJet;
        TH2 *dEtaPhiMu, *dEtaPhiJet;
        TH1 *dRminMu1jet, *dRminMu2jet, *dRminMuJet;
        TH1 *hptrelMu1, *hptrelMu2;
        TH2 *ptrelVsdRminMu1jet, *ptrelVsdRminMu2jet;
        TH2 *jetID2d;

        TH1 *ptTag, *etaTag, *phiTag ; 
        TH1 *ptProbe, *etaProbe, *phiProbe ; 
        TH1 *ptProbeRiso100, *etaProbeRiso100, *phiProbeRiso100 ; 
        TH1 *ptProbeRiso50, *etaProbeRiso50, *phiProbeRiso50 ; 
        TH1 *ptProbeRiso20, *etaProbeRiso20, *phiProbeRiso20 ; 
        TH1 *ptProbeRiso10, *etaProbeRiso10, *phiProbeRiso10 ; 
        TH1 *ptProbeRiso5, *etaProbeRiso5, *phiProbeRiso5 ; 
        TH1 *tagTrackIso, *tagTrackRelIso, *probeTrackIso, *probeTrackRelIso ; 
        TH1 *mMuMuTP, *mMuMuTPRiso100, *mMuMuTPRiso50, *mMuMuTPRiso20, *mMuMuTPRiso10, *mMuMuTPRiso5 ; 
        
        TH1 *dptMu1gen, *dptMu2gen;
        TH1 *dRMu1gen, *dRMu2gen;
        TH1 *qualMu1, *qualMu2;

        TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1caloIso, *mu1dB;
        TH1 *mu2trackIso, *mu2hcalIso, *mu2ecalIso, *mu2caloIso, *mu2dB;

        TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1caloRelIso;
        TH1 *mu2trackRelIso, *mu2hcalRelIso, *mu2ecalRelIso, *mu2caloRelIso;

        TH1 *mMuMu, *mMuMuOS, *mMuMuSS, *diMuCharge, *mMuMuZoom, *mMuMuGenZoom;
        TH1 *mWR, *mNuR1, *mNuR2, *mJJ;
        TH2 *mNuR2D, *jetPtvsNum;
        TH1 *mu1ptFracWRmass;

        TH1 *btagJet1, *btagJet2;
        TH1 *numBjets, *njets;

        TH1* met;

        TH1* czeta_mumu;
        TH1* czeta_mumu_zoom;

        TH1 *mu1jj_surfarea, *mu2jj_surfarea;

        // Jeremy's crazy angles...
        TH1* ctheta_mumu, *cthetaz_mumu;
        TH1* ctheta_jj, *cthetaz_jj;
        TH1* ctheta_mu1_jj, *cthetaz_mu1_jj;
        TH1* ctheta_mu2_jj, *cthetaz_mu2_jj;

        // Vertex plots
        TH1* vtx_mumu, *vtx_jj, *vtx_min_mu1j, *vtx_min_mu2j, *vtx_min_muj, *vtx_max_dist;

        TFileDirectory *mydir;
        TFileDirectory *nndir;

        HeavyNuTrigger::trigHistos_t trigHistos;

        // mc type
        TH1* mc_type;

    };

    bool init_;

    // gf set of histo for all Z definitions in a stack

    struct HistStruct
    {
        TH1 *nelec, *njet, *nmet, *nmuAll, *nmuLoose, *nmuTight;
        TH1 *muPt, *muEta, *muPhi, *looseMuPt, *tightMuPt;
        TH1* mc_type;

        TH1* cutlevel;
        TH1* weights;

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
        HistPerDef noCuts;
        HistPerDef LLptCuts;
        HistPerDef TrigMatches;
        HistPerDef LLJJptCuts;
        HistPerDef VertexCuts;
        HistPerDef Mu1HighPtCut;
        HistPerDef Mu1Pt30GeVCut;
        HistPerDef Mu1Pt40GeVCut;
        HistPerDef Mu1Pt50GeVCut;
        HistPerDef Mu1Pt60GeVCut;
        HistPerDef Mu1Pt80GeVCut;
        HistPerDef Mu1Pt100GeVCut;
        HistPerDef Mu1HighPtCutVtxEq1;
        HistPerDef Mu1HighPtCutVtx2to5;
        HistPerDef Mu1HighPtCutVtxGt5;
        HistPerDef Mu1HighPtCutNoJets;
        HistPerDef Mu1HighPtCut1Jet;
        HistPerDef diLmassCut;
        HistPerDef loDiLmassCut;
        HistPerDef mWRmassCut;
        HistPerDef oneBtag;
        HistPerDef twoBtag;
        // efficiency studies:
        // HistPerDef Mu1tagInZwin;
        // HistPerDef Mu2tagInZwin;
        // HistPerDef Mu1tagMu2passesInZwin;
        // HistPerDef Mu2tagMu1passesInZwin;
        HistPerDef TightTagTrackProbeInZwin;
        HistPerDef TightTagTrackProbePassesInZwin;
        HistPerDef TightTagTightProbeInZwin;
        HistPerDef TightTagTightProbePassesInZwin;
        HistPerDef Mu1TrigMatchesInZwin;
        HistPerDef Mu2TrigMatchesInZwin;
        HistPerDef Mu1Mu2TrigMatchesInZwin;
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



const int muonQualityFlags = 4;
const std::string muonQuality[] = {
    "All", "AllGlobalMuons", "AllStandAloneMuons", "AllTrackerMuons"
};

inline void labelMuonQualAxis(TAxis *ax)
{
    for(int i = 0; i < muonQualityFlags; i++)
    {
        ax->SetBinLabel(i + 1, muonQuality[i].c_str());
        ax->SetBinLabel(i + 1, muonQuality[i].c_str());
    }
}

inline std::string nnhistoname(int mwr, int mnu)
{
    return("WR" + int2str<int>(mwr) + "nuRmu" + int2str<int>(mnu));
}

void HeavyNu::HistPerDef::book(TFileDirectory *td, const std::string& post,
        const std::vector<hNuMassHypothesis>& v_masspts)
{
    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();

    mydir = td;

    // ----------  Muon histograms  ----------

    t = "p_{T}(#mu_{1}) " + post;
    ptMu1 = td->make<TH1D > ("ptMu1", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(#mu_{2}) " + post;
    ptMu2 = td->make<TH1D > ("ptMu2", t.c_str(), 100, 0., 1000.);
    t = "#eta(#mu_{1}) " + post;
    etaMu1 = td->make<TH1D > ("etaMu1", t.c_str(), 40, -2.5, 2.5);
    t = "#eta(#mu_{2}) " + post;
    etaMu2 = td->make<TH1D > ("etaMu2", t.c_str(), 40, -2.5, 2.5);
    t = "#phi(#mu_{1}) " + post;
    phiMu1 = td->make<TH1D > ("phiMu1", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(#mu_{2}) " + post;
    phiMu2 = td->make<TH1D > ("phiMu2", t.c_str(), 30, -3.14159, 3.14159);

    t = "p_{T}(#mu_{1}) vs. p_{T}(#mu_{2}) (SS) " + post + ";p_{T}(#mu_{1})(GeV);p_{T}(#mu_{2}(GeV))";

    ptMu1VsPtMu2ss = td->make<TH2D > ("ptMu1VsPtMu2ss", t.c_str(), 50, 0., 2000., 50, 0., 2000);

    t = "p_{T}(#mu_{1}) vs. p_{T}(#mu_{2}) (OS) " + post + ";p_{T}(#mu_{1})(GeV);p_{T}(#mu_{2}(GeV))";

    ptMu1VsPtMu2os = td->make<TH2D > ("ptMu1VsPtMu2os", t.c_str(), 50, 0., 2000., 50, 0., 2000);

    // delta angles

    t = "#Delta#eta(#mu_{1},#mu_{2}) " + post;
    dEtaMu = td->make<TH1D > ("dEtaMu", t.c_str(), 40, 0, 5);
    t = "#Delta#phi(#mu_{1},#mu_{2}) " + post;
    dPhiMu = td->make<TH1D > ("dPhiMu", t.c_str(), 30, 0, 3.14159);
    t = "#Delta p_{T}(#mu_{1},gen) " + post;
    dptMu1gen = td->make<TH1D > ("dptMu1gen", t.c_str(), 50, -0.50, 0.50);
    t = "#Delta p_{T}(#mu_{2},gen) " + post;
    dptMu2gen = td->make<TH1D > ("dptMu2gen", t.c_str(), 50, -0.50, 0.50);
    t = "#Delta R(#mu_{1},gen) " + post;
    dRMu1gen = td->make<TH1D > ("dRMu1gen", t.c_str(), 50, 0, 0.01);
    t = "#Delta R(#mu_{2},gen) " + post;
    dRMu2gen = td->make<TH1D > ("dRMu2gen", t.c_str(), 50, 0, 0.01);
    t = "#mu #Delta#eta vs. #Delta#phi " + post;
    t += ";#Delta#eta; #Delta#phi";
    dEtaPhiMu = td->make<TH2D > ("dEtaPhiMu", t.c_str(),
            50, 0, 5, 30, 0, 3.14159);

    t = "Quality (#mu_{1}) " + post;
    qualMu1 = td->make<TH1D > ("qualMu1", t.c_str(), muonQualityFlags, 0, muonQualityFlags);
    t = "Quality (#mu_{2}) " + post;
    qualMu2 = td->make<TH1D > ("qualMu2", t.c_str(), muonQualityFlags, 0, muonQualityFlags);
    labelMuonQualAxis(qualMu1->GetXaxis());
    labelMuonQualAxis(qualMu2->GetXaxis());

    // isolation

    t = "trackIso(#mu_{1}) " + post;
    mu1trackIso = td->make<TH1D > ("mu1trackIso", t.c_str(), 40, 0., 200.);
    t = "hcalIso(#mu_{1}) " + post;
    mu1hcalIso = td->make<TH1D > ("mu1hcalIso", t.c_str(), 40, 0., 200.);
    t = "ecalIso(#mu_{1}) " + post;
    mu1ecalIso = td->make<TH1D > ("mu1ecalIso", t.c_str(), 40, 0., 200.);
    t = "caloIso(#mu_{1}) " + post;
    mu1caloIso = td->make<TH1D > ("mu1caloIso", t.c_str(), 40, 0., 200.);
    t = "Dxy(#mu_{1}) " + post;
    mu1dB = td->make<TH1D > ("mu1dB", t.c_str(), 50, 0., 1.);

    t = "trackIso(#mu_{2}) " + post;
    mu2trackIso = td->make<TH1D > ("mu2trackIso", t.c_str(), 40, 0., 200.);
    t = "hcalIso(#mu_{2}) " + post;
    mu2hcalIso = td->make<TH1D > ("mu2hcalIso", t.c_str(), 40, 0., 200.);
    t = "ecalIso(#mu_{2}) " + post;
    mu2ecalIso = td->make<TH1D > ("mu2ecalIso", t.c_str(), 40, 0., 200.);
    t = "caloIso(#mu_{2}) " + post;
    mu2caloIso = td->make<TH1D > ("mu2caloIso", t.c_str(), 40, 0., 200.);
    t = "Dxy(#mu_{2}) " + post;
    mu2dB = td->make<TH1D > ("mu2dB", t.c_str(), 50, 0., 1.);

    t = "trackRelIso(#mu_{1}) " + post;
    mu1trackRelIso = td->make<TH1D > ("mu1trackRelIso", t.c_str(), 50, 0., 5.);
    t = "hcalRelIso(#mu_{1}) " + post;
    mu1hcalRelIso = td->make<TH1D > ("mu1hcalRelIso", t.c_str(), 50, 0., 5.);
    t = "ecalRelIso(#mu_{1}) " + post;
    mu1ecalRelIso = td->make<TH1D > ("mu1ecalRelIso", t.c_str(), 50, 0., 5.);
    t = "caloRelIso(#mu_{1}) " + post;
    mu1caloRelIso = td->make<TH1D > ("mu1caloRelIso", t.c_str(), 50, 0., 5.);

    t = "trackRelIso(#mu_{2}) " + post;
    mu2trackRelIso = td->make<TH1D > ("mu2trackRelIso", t.c_str(), 50, 0., 5.);
    t = "hcalRelIso(#mu_{2}) " + post;
    mu2hcalRelIso = td->make<TH1D > ("mu2hcalRelIso", t.c_str(), 50, 0., 5.);
    t = "ecalRelIso(#mu_{2}) " + post;
    mu2ecalRelIso = td->make<TH1D > ("mu2ecalRelIso", t.c_str(), 50, 0., 5.);
    t = "caloRelIso(#mu_{2}) " + post;
    mu2caloRelIso = td->make<TH1D > ("mu2caloRelIso", t.c_str(), 50, 0., 5.);

    // ----------  Jet histograms ----------

    t = "p_{T}(j_{1}) " + post;
    ptJet1 = td->make<TH1D > ("ptJet1", t.c_str(), 50, 0., 500.);
    t = "p_{T}(j_{2}) " + post;
    ptJet2 = td->make<TH1D > ("ptJet2", t.c_str(), 50, 0., 500.);
    t = "#eta(j_{1}) " + post;
    etaJet1 = td->make<TH1D > ("etaJet1", t.c_str(), 40, -5, 5);
    t = "#eta(j_{2}) " + post;
    etaJet2 = td->make<TH1D > ("etaJet2", t.c_str(), 40, -5, 5);
    t = "#phi(j_{1}) " + post;
    phiJet1 = td->make<TH1D > ("phiJet1", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(j_{2}) " + post;
    phiJet2 = td->make<TH1D > ("phiJet2", t.c_str(), 30, -3.14159, 3.14159);

    t = "#Delta#eta(j_{1},j_{2}) " + post;
    dEtaJet = td->make<TH1D > ("dEtaJet", t.c_str(), 40, 0, 5);
    t = "#Delta#phi(j_{1},j_{2}) " + post;
    dPhiJet = td->make<TH1D > ("dPhiJet", t.c_str(), 30, 0, 3.14159);

    t = "btag(j_{1}) " + post;
    btagJet1 = td->make<TH1D > ("btagJet1", t.c_str(), 40, 0, 5);
    t = "btag(j_{2}) " + post;
    btagJet2 = td->make<TH1D > ("btagJet2", t.c_str(), 40, 0, 5);

    t = "# B-tagged Jets in Event " + post;
    numBjets = td->make<TH1D > ("numBjets", t.c_str(), 3, -0.5, 2.5);
    t = "#  Jets in Event " + post;
    njets = td->make<TH1D > ("njets", t.c_str(), 100, -0.5, 99.5);

    t = "Jet #Delta#eta vs. #Delta#phi ";
    t += post + ";#Delta#eta; #Delta#phi";
    dEtaPhiJet = td->make<TH2D > ("dEtaPhiJet", t.c_str(),
            50, 0, 5, 30, 0, 3.14159);
    t = "Jet ID(j_{2}) vs. ID(j_{1}) ";
    t += post + "; ID(j_{1}); ID(j_{2})";
    jetID2d = td->make<TH2I > ("jetID2d", t.c_str(), 3, 0, 3, 3, 0, 3);
    labelJetIDaxis(jetID2d->GetXaxis());
    labelJetIDaxis(jetID2d->GetYaxis());

    // ----------  MET histograms     ----------

    t = "MET distribution " + post;
    met = td->make<TH1D > ("met", t.c_str(), 100, 0, 2000);


    t = "MC Type " + post;
    mc_type = td->make<TH1D > ("mc_type", "MC Type Code", 100, -0.5, 99.5);

    // ----------  Mu/Jet histograms  ----------

    t = "Minimum #Delta R(#mu,jet) " + post;
    dRminMuJet = td->make<TH1D > ("dRminMuJet", t.c_str(), 50, 0, 5.);
    t = "Minimum #Delta R(#mu_{1},jet) " + post;
    dRminMu1jet = td->make<TH1D > ("dRminMu1jet", t.c_str(), 50, 0, 5.);
    t = "Minimum #Delta R(#mu_{2},jet) " + post;
    dRminMu2jet = td->make<TH1D > ("dRminMu2jet", t.c_str(), 50, 0, 5.);

    t = "p_{T,rel}(#mu_{1},jet)" + post;
    hptrelMu1 = td->make<TH1D > ("ptrelMu1", t.c_str(), 50, 0, 1000.);
    t = "p_{T,rel}(#mu_{2},jet)" + post;
    hptrelMu2 = td->make<TH1D > ("ptrelMu2", t.c_str(), 50, 0, 1000.);

    t = "p_{T,rel}(#mu_{1},jet) vs #Delta R(#mu_{1},jet)" + post;
    t += "; #Delta R(#mu_{1},jet); p_{T,rel}(#mu_{1},jet)";
    ptrelVsdRminMu1jet = td->make<TH2D > ("ptrelVsdRminMu1jet",
            t.c_str(),
            50, 0, 5., 50, 0, 1000);
    t = "p_{T,rel}(#mu_{2},jet) vs #Delta R(#mu_{2},jet)" + post;
    t += "; #Delta R(#mu_{2},jet); p_{T,rel}(#mu_{2},jet)";
    ptrelVsdRminMu2jet = td->make<TH2D > ("ptrelVsdRminMu2jet",
            t.c_str(),
            50, 0, 5., 50, 0, 1000);

    // ----------  Composite histograms  ----------

    t = "M(W_{R}) " + post;
    mWR = td->make<TH1D > ("mWR", t.c_str(), 70, 0, 2800);
    t = "M(N_{R}) with #mu_{1} " + post;
    mNuR1 = td->make<TH1D > ("mNuR1", t.c_str(), 70, 0, 2800);
    t = "M(N_{R}) with #mu_{2} " + post;
    mNuR2 = td->make<TH1D > ("mNuR2", t.c_str(), 70, 0, 1400);
    t = "M(N_{R}) #mu_{1} vs. #mu_{2} " + post;
    mNuR2D = td->make<TH2D > ("mNuR2D", t.c_str(), 70, 0, 2800, 70, 0, 1400);

    t = "#mu_{1} p_{T} / M(W_{R}) " + post;
    mu1ptFracWRmass = td->make<TH1D > ("mu1ptFracWRmass", t.c_str(), 100, 0, 1.);

    t = "M(jj) " + post;
    mJJ = td->make<TH1D > ("mJJ", t.c_str(), 50, 0, 2000);
    t = "M(#mu #mu) " + post;
    mMuMu = td->make<TH1D > ("mMuMu", t.c_str(), 100, 0, 2000);
    t = "M(#mu #mu)(OS) " + post;
    mMuMuOS = td->make<TH1D > ("mMuMuOS", t.c_str(), 100, 0, 2000);
    t = "M(#mu #mu)(SS) " + post;
    mMuMuSS = td->make<TH1D > ("mMuMuSS", t.c_str(), 100, 0, 2000);
    t = "M(#mu #mu) " + post;
    mMuMuZoom = td->make<TH1D > ("mMuMuZoom", t.c_str(), 400, 0, 400);
    t = "M(#mu #mu)(generated) " + post;
    mMuMuGenZoom = td->make<TH1D > ("mMuMuGenZoom", t.c_str(), 400, 0, 400);
    t = "DiMuon Charge " + post;
    diMuCharge = td->make<TH1D > ("diMuCharge", t.c_str(), 2, -1, 1);

    diMuCharge->GetXaxis()->SetBinLabel(1, "OS");
    diMuCharge->GetXaxis()->SetBinLabel(2, "SS");

    t = "cZeta(mumu) " + post;
    czeta_mumu = td->make<TH1D > ("czMM", t.c_str(), 100, -1, 1);
    t = "cZeta(mumu) Zoom " + post;
    czeta_mumu_zoom = td->make<TH1D > ("czMMzoom", t.c_str(), 100, -1, -0.9);

    // crazy angles
    t = "cT(mumu) " + post;
    ctheta_mumu = td->make<TH1D > ("ctMM", t.c_str(), 50, 0, 1);
    t = "cT(jj) " + post;
    ctheta_jj = td->make<TH1D > ("ctJJ", t.c_str(), 50, 0, 1);
    t = "cT(mu1-jj) " + post;
    ctheta_mu1_jj = td->make<TH1D > ("ctM1JJ", t.c_str(), 50, 0, 1);
    t = "cT(mu2-jj) " + post;
    ctheta_mu2_jj = td->make<TH1D > ("ctM2JJ", t.c_str(), 50, 0, 1);

    t = "cTz(mumu) " + post;
    cthetaz_mumu = td->make<TH1D > ("ctzMM", t.c_str(), 50, 0, 1);
    t = "cTz(jj) " + post;
    cthetaz_jj = td->make<TH1D > ("ctzJJ", t.c_str(), 50, 0, 1);
    t = "cTz(mu1-jj) " + post;
    cthetaz_mu1_jj = td->make<TH1D > ("ctzM1JJ", t.c_str(), 50, 0, 1);
    t = "cTz(mu2-jj) " + post;
    cthetaz_mu2_jj = td->make<TH1D > ("ctzM2JJ", t.c_str(), 50, 0, 1);

    t = "surface area m1jj " + post;
    mu1jj_surfarea = td->make<TH1D > ("mu1jj_surfarea", t.c_str(), 80, 0., 6.2831853);
    t = "surface area m2jj " + post;
    mu2jj_surfarea = td->make<TH1D > ("mu2jj_surfarea", t.c_str(), 80, 0., 6.2831853);

    // vertex histograms
    t = "mumu vtx DeltaZ " + post;
    vtx_mumu = td->make<TH1D > ("vtxMM", t.c_str(), 200, 0., 2.);
    t = "jj vtx DeltaZ " + post;
    vtx_jj = td->make<TH1D > ("vtxJJ", t.c_str(), 200, 0., 2.);
    t = "m1j vtx min DeltaZ " + post;
    vtx_min_mu1j = td->make<TH1D > ("vtxM1Jmin", t.c_str(), 200, 0., 2.);
    t = "m2j vtx min DeltaZ " + post;
    vtx_min_mu2j = td->make<TH1D > ("vtxM2Jmin", t.c_str(), 200, 0., 2.);
    t = "mu-j vtx min DeltaZ " + post;
    vtx_min_muj = td->make<TH1D > ("vtxMuJmin", t.c_str(), 200, 0., 2.);
    t = "mmjj vtx max DeltaZ " + post;
    vtx_max_dist = td->make<TH1D > ("vtxDistmax", t.c_str(), 400, 0., 5.);


    // ----------  Neural Net histograms  ----------

    if(v_masspts.size())
    {
        nndir = new TFileDirectory(td->mkdir("_NNdata"));
        for(size_t i = 0; i < v_masspts.size(); i++)
        {
            int mwr = v_masspts[i].first;
            int mnu = v_masspts[i].second;
            std::string name = nnhistoname(mwr, mnu);
            nndir->make<TH1D > (name.c_str(), (name + post).c_str(), 51, -0.01, 1.01);
        }
    }
}// end of book()

void HeavyNu::HistPerDef::bookTagProbe(TFileDirectory *td, const std::string& post)
{
    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();

    mydir = td;

    // ----------  Muon histograms  ----------

    t = "p_{T}(tag) " + post;
    ptTag = td->make<TH1D > ("ptTag", t.c_str(), 100, 0., 1000.);
    t = "#eta(tag) " + post;
    etaTag = td->make<TH1D > ("etaTag", t.c_str(), 40, -2.5, 2.5);
    t = "#phi(tag) " + post;
    phiTag = td->make<TH1D > ("phiTag", t.c_str(), 30, -3.14159, 3.14159);

    t = "p_{T}(probe) " + post;
    ptProbe = td->make<TH1D > ("ptProbe", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 100% probe relIso " + post;
    ptProbeRiso100 = td->make<TH1D > ("ptProbeRiso100", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 50% probe relIso " + post;
    ptProbeRiso50 = td->make<TH1D > ("ptProbeRiso50", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 20% probe relIso " + post;
    ptProbeRiso20 = td->make<TH1D > ("ptProbeRiso20", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 10% probe relIso " + post;
    ptProbeRiso10 = td->make<TH1D > ("ptProbeRiso10", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 5% probe relIso " + post;
    ptProbeRiso5 = td->make<TH1D > ("ptProbeRiso5", t.c_str(), 100, 0., 1000.);

    t = "#eta(probe) " + post;
    etaProbe = td->make<TH1D > ("etaProbe", t.c_str(), 40, -2.5, 2.5);
    t = "#eta(probe) 100% probe relIso " + post;
    etaProbeRiso100 = td->make<TH1D > ("etaProbeRiso100", t.c_str(), 40, -2.5, 2.5);
    t = "#eta(probe) 50% probe relIso " + post;
    etaProbeRiso50 = td->make<TH1D > ("etaProbeRiso50", t.c_str(), 40, -2.5, 2.5);
    t = "#eta(probe) 20% probe relIso " + post;
    etaProbeRiso20 = td->make<TH1D > ("etaProbeRiso20", t.c_str(), 40, -2.5, 2.5);
    t = "#eta(probe) 10% probe relIso " + post;
    etaProbeRiso10 = td->make<TH1D > ("etaProbeRiso10", t.c_str(), 40, -2.5, 2.5);
    t = "#eta(probe) 5% probe relIso " + post;
    etaProbeRiso5 = td->make<TH1D > ("etaProbeRiso5", t.c_str(), 40, -2.5, 2.5);

    t = "#phi(probe) " + post;
    phiProbe = td->make<TH1D > ("phiProbe", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 100% probe relIso " + post;
    phiProbeRiso100 = td->make<TH1D > ("phiProbeRiso100", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 50% probe relIso " + post;
    phiProbeRiso50 = td->make<TH1D > ("phiProbeRiso50", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 20% probe relIso " + post;
    phiProbeRiso20 = td->make<TH1D > ("phiProbeRiso20", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 10% probe relIso " + post;
    phiProbeRiso10 = td->make<TH1D > ("phiProbeRiso10", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 5% probe relIso " + post;
    phiProbeRiso5 = td->make<TH1D > ("phiProbeRiso5", t.c_str(), 30, -3.14159, 3.14159);


    // isolation
    t = "trackIso(tag) " + post;
    tagTrackIso = td->make<TH1D > ("tagTrackIso", t.c_str(), 100, 0., 100.);
    t = "trackIso(tag) " + post;
    probeTrackIso = td->make<TH1D > ("probeTrackIso", t.c_str(), 100, 0., 100.);
    t = "trackRelIso(tag) " + post;
    tagTrackRelIso = td->make<TH1D > ("tagTrackRelIso", t.c_str(), 100, 0., 2.);
    t = "trackRelIso(probe) " + post;
    probeTrackRelIso = td->make<TH1D > ("probeTrackRelIso", t.c_str(), 100, 0., 2.);


    // ----------  Composite histograms  ----------
    t = "M(#mu #mu) " + post;
    mMuMuTP = td->make<TH1D > ("mMuMuTP", t.c_str(), 80, 80, 120);
    t = "M(#mu #mu) 5% probe relIso " + post;
    mMuMuTPRiso5 = td->make<TH1D > ("mMuMuTPRiso5", t.c_str(), 80, 80, 120);
    t = "M(#mu #mu) 10% probe relIso " + post;
    mMuMuTPRiso10 = td->make<TH1D > ("mMuMuTPRiso10", t.c_str(), 80, 80, 120);
    t = "M(#mu #mu) 20% probe relIso " + post;
    mMuMuTPRiso20 = td->make<TH1D > ("mMuMuTPRiso20", t.c_str(), 80, 80, 120);
    t = "M(#mu #mu) 50% probe relIso " + post;
    mMuMuTPRiso50 = td->make<TH1D > ("mMuMuTPRiso50", t.c_str(), 80, 80, 120);
    t = "M(#mu #mu) 100% probe relIso " + post;
    mMuMuTPRiso100 = td->make<TH1D > ("mMuMuTPRiso100", t.c_str(), 80, 80, 120);

}// end of book()

void HeavyNu::HistPerDef::fill(pat::MuonCollection muons,
                               pat::JetCollection jets,
                               pat::METCollection metc,
                               bool isMC,
                               double wgt,
                               bool pfJets)
{
    std::sort(muons.begin(), muons.end(), hnu::pTcompare());
    std::sort(jets.begin(), jets.end(), hnu::pTcompare());

    reco::Particle::LorentzVector vWR;

    const pat::Muon& m0 = muons.at(0);
    const pat::Muon& m1 = muons.at(1);

    // Muons
    // std::cout << "Muon collection with size " << muons.size() 
    // 	      << " and Jet collection with size " << jets.size() << std::endl ; 
    // std::cout << "No cuts, filling with weight: " << wgt << std::endl ; 
    ptMu1->Fill(m0.pt(), wgt);
    ptMu2->Fill(m1.pt(), wgt);

    etaMu1->Fill(m0.eta(), wgt);
    etaMu2->Fill(m1.eta(), wgt);

    phiMu1->Fill(m0.phi(), wgt);
    phiMu2->Fill(m1.phi(), wgt);

    dPhiMu->Fill(fabs(deltaPhi(m0.phi(), m1.phi())), wgt);
    dEtaMu->Fill(fabs(m0.eta() - m1.eta()), wgt);
    dEtaPhiMu->Fill(fabs(m0.eta() - m1.eta()),
            fabs(deltaPhi(m0.phi(), m1.phi())), wgt);

    mu1trackIso ->Fill(m0.trackIso(), wgt);
    mu1hcalIso ->Fill(m0.hcalIso(), wgt);
    mu1ecalIso ->Fill(m0.ecalIso(), wgt);
    mu1caloIso ->Fill(m0.caloIso(), wgt);
    mu1dB ->Fill(m0.dB(), wgt);

    mu2trackIso ->Fill(m1.trackIso(), wgt);
    mu2hcalIso ->Fill(m1.hcalIso(), wgt);
    mu2ecalIso ->Fill(m1.ecalIso(), wgt);
    mu2caloIso ->Fill(m1.caloIso(), wgt);
    mu2dB ->Fill(m1.dB(), wgt);

    mu1trackRelIso->Fill(m0.trackIso() / m0.pt(), wgt);
    mu1hcalRelIso ->Fill(m0.hcalIso() / m0.pt(), wgt);
    mu1ecalRelIso ->Fill(m0.ecalIso() / m0.pt(), wgt);
    mu1caloRelIso ->Fill(m0.caloIso() / m0.pt(), wgt);

    mu2trackRelIso->Fill(m1.trackIso() / m1.pt(), wgt);
    mu2hcalRelIso ->Fill(m1.hcalIso() / m1.pt(), wgt);
    mu2ecalRelIso ->Fill(m1.ecalIso() / m1.pt(), wgt);
    mu2caloRelIso ->Fill(m1.caloIso() / m1.pt(), wgt);

    if(isMC)
    {
        for(unsigned int i = 0; i < 2; i++)
        {
            if(muons.at(i).genLepton() != 0)
            {
                float dpt = muons.at(i).pt() - muons.at(i).genLepton()->pt();
                float dR = deltaR(muons.at(i).eta(), muons.at(i).phi(),
                        muons.at(i).genLepton()->eta(), muons.at(i).genLepton()->phi());
                if(i == 0)
                {
                    dptMu1gen->Fill(dpt / muons.at(i).genLepton()->pt());
                    dRMu1gen->Fill(dR);
                }
                else
                {
                    dptMu2gen->Fill(dpt / muons.at(i).genLepton()->pt());
                    dRMu2gen->Fill(dR);
                }
            }
        }
        if((m0.genLepton() != 0) &&
                (m1.genLepton() != 0))
        {
            reco::Particle::LorentzVector mu1gp4 = m0.genLepton()->p4();
            reco::Particle::LorentzVector mu2gp4 = m1.genLepton()->p4();
            mMuMuGenZoom->Fill((mu1gp4 + mu2gp4).M());
        }
    }
    for(int i = 0; i < muonQualityFlags; i++)
    {
        if(m0.muonID(muonQuality[i])) qualMu1->Fill(i, wgt);
        if(m1.muonID(muonQuality[i])) qualMu2->Fill(i, wgt);
    }

    // Jets
    const pat::Jet& j0 = jets.at(0);
    const pat::Jet& j1 = jets.at(1);

    ptJet1->Fill(j0.pt(), wgt);
    ptJet2->Fill(j1.pt(), wgt);

    etaJet1->Fill(j0.eta(), wgt);
    etaJet2->Fill(j1.eta(), wgt);

    phiJet1->Fill(j0.phi(), wgt);
    phiJet2->Fill(j1.phi(), wgt);

    double j0bdisc = j0.bDiscriminator(btagName);
    double j1bdisc = j1.bDiscriminator(btagName);

    btagJet1->Fill(j0bdisc, wgt);
    btagJet2->Fill(j1bdisc, wgt);

    njets->Fill(jets.size());

    if((j0bdisc >= minBtagDiscVal) &&
            (j1bdisc >= minBtagDiscVal)) numBjets->Fill(2., wgt);
    else if((j0bdisc >= minBtagDiscVal) ||
            (j1bdisc >= minBtagDiscVal)) numBjets->Fill(1., wgt);
    else numBjets->Fill(0., wgt);

    dPhiJet->Fill(fabs(deltaPhi(j0.phi(), j1.phi())), wgt);
    dEtaJet->Fill(fabs(j0.eta() - j1.eta()), wgt);
    dEtaPhiJet->Fill(fabs(j0.eta() - j1.eta()),
            fabs(deltaPhi(j0.phi(), j1.phi())), wgt);

    if(!pfJets) jetID2d->Fill(hnu::jetID(j0), hnu::jetID(j1), wgt);

    // met
    if(metc.size())
        met->Fill(metc.at(0).pt(), wgt);
    else
        met->Fill(0., wgt);

    // Muon-Jet plots
    float dRmu1jet1 = deltaR(m0.eta(), m0.phi(), j0.eta(), j0.phi());
    float dRmu1jet2 = deltaR(m0.eta(), m0.phi(), j1.eta(), j1.phi());
    float dRmu2jet1 = deltaR(m1.eta(), m1.phi(), j0.eta(), j0.phi());
    float dRmu2jet2 = deltaR(m1.eta(), m1.phi(), j1.eta(), j1.phi());

    const pat::Jet& j4mu1 = (dRmu1jet1 < dRmu1jet2)?j0:j1;
    const pat::Jet& j4mu2 = (dRmu2jet1 < dRmu2jet2)?j0:j1;

    TVector3 mu1vec(m0.momentum().X(), m0.momentum().Y(), m0.momentum().Z());
    TVector3 mu2vec(m1.momentum().X(), m1.momentum().Y(), m1.momentum().Z());

    TVector3 jt1vec(j4mu1.p4().Vect().X(), j4mu1.p4().Vect().Y(), j4mu1.p4().Vect().Z());
    TVector3 jt2vec(j4mu2.p4().Vect().X(), j4mu2.p4().Vect().Y(), j4mu2.p4().Vect().Z());

    double ptrelMu1 = mu1vec.Perp(jt1vec);
    double ptrelMu2 = mu2vec.Perp(jt2vec);

    dRminMu1jet->Fill(min(dRmu1jet1, dRmu1jet2), wgt);
    dRminMu2jet->Fill(min(dRmu2jet1, dRmu2jet2), wgt);

    dRminMuJet->Fill(min(min(dRmu1jet1, dRmu1jet2), min(dRmu2jet1, dRmu2jet2)), wgt);

    hptrelMu1->Fill(ptrelMu1, wgt);
    hptrelMu2->Fill(ptrelMu2, wgt);

    ptrelVsdRminMu1jet->Fill(min(dRmu1jet1, dRmu1jet2), ptrelMu1, wgt);
    ptrelVsdRminMu2jet->Fill(min(dRmu2jet1, dRmu1jet2), ptrelMu2, wgt);

    // Composite objects
    reco::Particle::LorentzVector vJJ = j0.p4() + j1.p4();
    vWR = vJJ + m0.p4() + m1.p4();

    mWR->Fill(vWR.M(), wgt);
    mNuR1 ->Fill((vJJ + m0.p4()).M(), wgt);
    mNuR2 ->Fill((vJJ + m1.p4()).M(), wgt);
    mNuR2D->Fill((vJJ + m0.p4()).M(), (vJJ + m1.p4()).M(), wgt);

    reco::Particle::LorentzVector mumu = m0.p4() + m1.p4();
    reco::Particle::LorentzVector jj = j0.p4() + j1.p4();

    mMuMu->Fill(mumu.M(), wgt);
    if(m0.charge() == m1.charge())
    {
        mMuMuSS->Fill(mumu.M(), wgt);
        ptMu1VsPtMu2ss->Fill(m0.pt(), m1.pt(), wgt);
    }
    else
    {
        mMuMuOS->Fill(mumu.M(), wgt);
        ptMu1VsPtMu2os->Fill(m0.pt(), m1.pt(), wgt);
    }

    diMuCharge->Fill(0.5 * m0.charge() * m1.charge(), wgt);

    mMuMuZoom->Fill(mumu.M(), wgt);
    mJJ->Fill(jj.M(), wgt);

    float deltaVzJ1J2 = fabs(hnu::avgVertex(j0, 1.0) - hnu::avgVertex(j1, 1.0));
    float deltaVzJ1M1 = fabs(hnu::avgVertex(j0, 1.0) - m0.vertex().Z());
    float deltaVzJ2M2 = fabs(hnu::avgVertex(j1, 1.0) - m1.vertex().Z());
    float deltaVzJ1M2 = fabs(hnu::avgVertex(j0, 1.0) - m1.vertex().Z());
    float deltaVzJ2M1 = fabs(hnu::avgVertex(j1, 1.0) - m0.vertex().Z());
    float deltaVzM1M2 = fabs(m0.vertex().Z() - m1.vertex().Z());

    vtx_jj->Fill(deltaVzJ1J2, wgt);
    float minDeltaVzMu1J = std::min(deltaVzJ1M1, deltaVzJ2M1);
    float minDeltaVzMu2J = std::min(deltaVzJ2M2, deltaVzJ2M2);
    vtx_min_mu1j->Fill(minDeltaVzMu1J, wgt);
    vtx_min_mu2j->Fill(minDeltaVzMu2J, wgt);
    vtx_min_muj->Fill(std::min(minDeltaVzMu1J, minDeltaVzMu2J), wgt);

    float maxDeltaVzMuJ1 = std::max(deltaVzJ1M1, deltaVzJ1M2);
    float maxDeltaVzMuJ2 = std::max(deltaVzJ2M1, deltaVzJ2M2);
    float maxDeltaVzMMJJ = std::max(deltaVzM1M2, deltaVzJ1J2);
    float maxDeltaVzMuJ = std::max(maxDeltaVzMuJ1, maxDeltaVzMuJ2);
    vtx_max_dist->Fill(std::max(maxDeltaVzMMJJ, maxDeltaVzMuJ), wgt);

}// end of fill()

void HeavyNu::HistPerDef::fill(const HeavyNuEvent& hne,
        const std::vector<hNuMassHypothesis>& v_masspts)
{
    double wgt = hne.eventWgt;

    mc_type->Fill(hne.mc_class, wgt);

    // Muons
    double mu1pt = hne.MuScale * hne.mu1.pt();
    double mu2pt = hne.MuScale * hne.mu2.pt();

    ptMu1->Fill(mu1pt, wgt);
    ptMu2->Fill(mu2pt, wgt);

    etaMu1->Fill(hne.mu1.eta(), wgt);
    etaMu2->Fill(hne.mu2.eta(), wgt);

    phiMu1->Fill(hne.mu1.phi(), wgt);
    phiMu2->Fill(hne.mu2.phi(), wgt);

    dPhiMu->Fill(fabs(deltaPhi(hne.mu1.phi(), hne.mu2.phi())), wgt);
    dEtaMu->Fill(fabs(hne.mu1.eta() - hne.mu2.eta()), wgt);
    dEtaPhiMu->Fill(fabs(hne.mu1.eta() - hne.mu2.eta()),
            fabs(deltaPhi(hne.mu1.phi(), hne.mu2.phi())), wgt);

    mu1trackIso->Fill(hne.mu1.trackIso(), wgt);
    mu1hcalIso ->Fill(hne.mu1.hcalIso(), wgt);
    mu1ecalIso ->Fill(hne.mu1.ecalIso(), wgt);
    mu1caloIso ->Fill(hne.mu1.caloIso(), wgt);
    mu1dB ->Fill(hne.mu1.dB(), wgt);
    mu2trackIso->Fill(hne.mu2.trackIso(), wgt);
    mu2hcalIso ->Fill(hne.mu2.hcalIso(), wgt);
    mu2ecalIso ->Fill(hne.mu2.ecalIso(), wgt);
    mu2caloIso ->Fill(hne.mu2.caloIso(), wgt);
    mu2dB ->Fill(hne.mu2.dB(), wgt);

    mu1trackRelIso->Fill(hne.mu1.trackIso() / mu1pt, wgt);
    mu1hcalRelIso ->Fill(hne.mu1.hcalIso() / mu1pt, wgt);
    mu1ecalRelIso ->Fill(hne.mu1.ecalIso() / mu1pt, wgt);
    mu1caloRelIso ->Fill(hne.mu1.caloIso() / mu1pt, wgt);
    mu2trackRelIso->Fill(hne.mu2.trackIso() / mu2pt, wgt);
    mu2hcalRelIso ->Fill(hne.mu2.hcalIso() / mu2pt, wgt);
    mu2ecalRelIso ->Fill(hne.mu2.ecalIso() / mu2pt, wgt);
    mu2caloRelIso ->Fill(hne.mu2.caloIso() / mu2pt, wgt);

    vtx_mumu->Fill(fabs(hne.mu1.vertex().Z() - hne.mu2.vertex().Z()), wgt);

    if(hne.isMC)
    {
        for(unsigned int i = 0; i < 2; i++)
        {
            if(hne.mu[i].genLepton() != 0)
            {
                float dpt = (hne.MuScale * hne.mu[i].pt()) - hne.mu[i].genLepton()->pt();
                float dR = deltaR(hne.mu[i].eta(), hne.mu[i].phi(),
                        hne.mu[i].genLepton()->eta(), hne.mu[i].genLepton()->phi());
                if(i == 0)
                {
                    dptMu1gen->Fill(dpt / hne.mu[i].genLepton()->pt());
                    dRMu1gen->Fill(dR);
                }
                else
                {
                    dptMu2gen->Fill(dpt / hne.mu[i].genLepton()->pt());
                    dRMu2gen->Fill(dR);
                }
            }
        }
        if((hne.mu1.genLepton() != 0) &&
                (hne.mu2.genLepton() != 0))
        {
            reco::Particle::LorentzVector mu1gp4 = hne.mu1.genLepton()->p4();
            reco::Particle::LorentzVector mu2gp4 = hne.mu2.genLepton()->p4();
            mMuMuGenZoom->Fill((mu1gp4 + mu2gp4).M());
        }
    }
    for(int i = 0; i < muonQualityFlags; i++)
    {
        if(hne.mu1.muonID(muonQuality[i])) qualMu1->Fill(i, wgt);
        if(hne.mu2.muonID(muonQuality[i])) qualMu2->Fill(i, wgt);
    }

    int jet1id = 0;
    int jet2id = 0;

    // Jets
    njets->Fill(hne.nJets);
    if(hne.nJets > 0)
    {
        if(!hne.pfJets) jet1id = hnu::jetID(hne.j1);

        double j1bdisc = hne.j1.bDiscriminator(btagName);

        ptJet1->Fill(hne.j1scale * hne.j1.pt(), wgt);
        etaJet1->Fill(hne.j1.eta(), wgt);
        phiJet1->Fill(hne.j1.phi(), wgt);
        btagJet1->Fill(j1bdisc, wgt);

        if(hne.nJets > 1)
        {
            double j2bdisc = hne.j2.bDiscriminator(btagName);

            ptJet2->Fill(hne.j2scale * hne.j2.pt(), wgt);
            etaJet2->Fill(hne.j2.eta(), wgt);
            phiJet2->Fill(hne.j2.phi(), wgt);
            btagJet2->Fill(j2bdisc, wgt);

            if((j1bdisc >= minBtagDiscVal) &&
                    (j2bdisc >= minBtagDiscVal)) numBjets->Fill(2., wgt);
            else if((j1bdisc >= minBtagDiscVal) ||
                    (j2bdisc >= minBtagDiscVal)) numBjets->Fill(1., wgt);
            else numBjets->Fill(0., wgt);

            if(!hne.pfJets) jet2id = hnu::jetID(hne.j2);

            dPhiJet->Fill(fabs(deltaPhi(hne.j1.phi(), hne.j2.phi())), wgt);
            dEtaJet->Fill(fabs(hne.j1.eta() - hne.j2.eta()), wgt);
            dEtaPhiJet->Fill(fabs(hne.j1.eta() - hne.j2.eta()),
                    fabs(deltaPhi(hne.j1.phi(), hne.j2.phi())), wgt);

            mWR->Fill(hne.mWR, wgt);
            mNuR1->Fill(hne.mNuR1, wgt);
            mNuR2->Fill(hne.mNuR2, wgt);
            mNuR2D->Fill(hne.mNuR1, hne.mNuR2, wgt);
            mJJ->Fill(hne.mJJ, wgt);

            mu1ptFracWRmass->Fill(mu1pt / hne.mWR, wgt);
            mu1jj_surfarea->Fill(hne.area_1jj, wgt);
            mu2jj_surfarea->Fill(hne.area_2jj, wgt);

            ctheta_jj->Fill(hne.ctheta_jj, wgt);
            ctheta_mu1_jj->Fill(hne.ctheta_mu1_jj, wgt);
            ctheta_mu2_jj->Fill(hne.ctheta_mu2_jj, wgt);
            cthetaz_jj->Fill(hne.cthetaz_jj, wgt);
            cthetaz_mu1_jj->Fill(hne.cthetaz_mu1_jj, wgt);
            cthetaz_mu2_jj->Fill(hne.cthetaz_mu2_jj, wgt);

            float deltaVzJ1J2 = fabs(hne.tjV1 - hne.tjV2);
            float deltaVzJ1M1 = fabs(hne.tjV1 - hne.mu1.vertex().Z());
            float deltaVzJ2M2 = fabs(hne.tjV2 - hne.mu2.vertex().Z());
            float deltaVzJ1M2 = fabs(hne.tjV1 - hne.mu2.vertex().Z());
            float deltaVzJ2M1 = fabs(hne.tjV2 - hne.mu1.vertex().Z());
            float deltaVzM1M2 = fabs(hne.mu1.vertex().Z() - hne.mu2.vertex().Z());

            vtx_jj->Fill(deltaVzJ1J2, wgt);
            float minDeltaVzMu1J = std::min(deltaVzJ1M1, deltaVzJ2M1);
            float minDeltaVzMu2J = std::min(deltaVzJ2M2, deltaVzJ2M2);
            vtx_min_mu1j->Fill(minDeltaVzMu1J, wgt);
            vtx_min_mu2j->Fill(minDeltaVzMu2J, wgt);
            vtx_min_muj->Fill(std::min(minDeltaVzMu1J, minDeltaVzMu2J), wgt);

            float maxDeltaVzMuJ1 = std::max(deltaVzJ1M1, deltaVzJ1M2);
            float maxDeltaVzMuJ2 = std::max(deltaVzJ2M1, deltaVzJ2M2);
            float maxDeltaVzMMJJ = std::max(deltaVzM1M2, deltaVzJ1J2);
            float maxDeltaVzMuJ = std::max(maxDeltaVzMuJ1, maxDeltaVzMuJ2);
            vtx_max_dist->Fill(std::max(maxDeltaVzMMJJ, maxDeltaVzMuJ), wgt);
        }

        dRminMu1jet->Fill(hne.dRminMu1jet, wgt);
        dRminMu2jet->Fill(hne.dRminMu2jet, wgt);

        dRminMuJet->Fill(min(hne.dRminMu1jet, hne.dRminMu2jet), wgt);

        hptrelMu1->Fill(hne.ptrelMu1, wgt);
        hptrelMu2->Fill(hne.ptrelMu2, wgt);

        ptrelVsdRminMu1jet->Fill(hne.dRminMu1jet, hne.ptrelMu1, wgt);
        ptrelVsdRminMu2jet->Fill(hne.dRminMu2jet, hne.ptrelMu2, wgt);
    }

    if(!hne.pfJets) jetID2d->Fill(jet1id, jet2id, wgt);

    // met
    met->Fill(hne.met1.pt(), wgt);

    mMuMu->Fill(hne.mMuMu, wgt);

    if(hne.mu1.charge() == hne.mu2.charge())
    {
        mMuMuSS->Fill(hne.mMuMu, wgt);
        ptMu1VsPtMu2ss->Fill(hne.mu1.pt(), hne.mu2.pt(), wgt);
    }
    else
    {
        mMuMuOS->Fill(hne.mMuMu, wgt);
        ptMu1VsPtMu2os->Fill(hne.mu1.pt(), hne.mu2.pt(), wgt);
    }

    mMuMuZoom->Fill(hne.mMuMu, wgt);

    diMuCharge->Fill(0.5 * hne.mu1.charge() * hne.mu2.charge(), wgt);

    czeta_mumu->Fill(hne.czeta_mumu, wgt);
    czeta_mumu_zoom->Fill(hne.czeta_mumu, wgt);
    ctheta_mumu->Fill(hne.ctheta_mumu, wgt);
    cthetaz_mumu->Fill(hne.cthetaz_mumu, wgt);

    // Neural net histos
    if(v_masspts.size())
    {
        // defines in HeavyNuCommon.h
#ifdef CMSSW_4XX
        TDirectory *nnrootdir = nndir->getBareDirectory("");
#endif
#ifdef CMSSW_3XX 
        TDirectory *nnrootdir = nndir->cd();
#endif
        for(size_t i = 0; i < v_masspts.size(); i++)
        {
            int mwr = v_masspts[i].first;
            int mnu = v_masspts[i].second;
            std::string name = nnhistoname(mwr, mnu);
            TH1D *nnh = (TH1D *)nnrootdir->Get(name.c_str());
            assert(nnh);
            nnh->Fill(hne.nnoutputs[i]);
        }
    }
}// end of fill()

void HeavyNu::HistPerDef::fill(const pat::Muon& theTag,
                               const pat::Muon& theProbe,
                               const double wgt)
{
    ptTag->Fill(theTag.pt(), wgt);
    etaTag->Fill(theTag.eta(), wgt);
    phiTag->Fill(theTag.phi(), wgt);

    ptProbe->Fill(theProbe.pt(), wgt);
    etaProbe->Fill(theProbe.eta(), wgt);
    phiProbe->Fill(theProbe.phi(), wgt);

    tagTrackIso ->Fill(theTag.trackIso(), wgt);
    probeTrackIso ->Fill(theProbe.trackIso(), wgt);

    tagTrackRelIso->Fill(theTag.trackIso() / theTag.pt(), wgt);
    probeTrackRelIso->Fill(theProbe.trackIso() / theProbe.pt(), wgt);

    reco::Particle::LorentzVector mumu = theTag.p4() + theProbe.p4();
    mMuMuTP->Fill(mumu.M(), wgt);

    double probe_relIso = theProbe.trackIso() / theProbe.pt() ; 
    if ( probe_relIso < 1.0 ) {
        ptProbeRiso100->Fill(theProbe.pt(), wgt);
        etaProbeRiso100->Fill(theProbe.eta(), wgt);
        phiProbeRiso100->Fill(theProbe.phi(), wgt);
        mMuMuTPRiso100->Fill(mumu.M(), wgt);
        if ( probe_relIso < 0.5 ) {
            ptProbeRiso50->Fill(theProbe.pt(), wgt);
            etaProbeRiso50->Fill(theProbe.eta(), wgt);
            phiProbeRiso50->Fill(theProbe.phi(), wgt);
            mMuMuTPRiso50->Fill(mumu.M(), wgt);
            if ( probe_relIso < 0.2 ) {
                ptProbeRiso20->Fill(theProbe.pt(), wgt);
                etaProbeRiso20->Fill(theProbe.eta(), wgt);
                phiProbeRiso20->Fill(theProbe.phi(), wgt);
                mMuMuTPRiso20->Fill(mumu.M(), wgt);
                if ( probe_relIso < 0.1 ) {
                    ptProbeRiso10->Fill(theProbe.pt(), wgt);
                    etaProbeRiso10->Fill(theProbe.eta(), wgt);
                    phiProbeRiso10->Fill(theProbe.phi(), wgt);
                    mMuMuTPRiso10->Fill(mumu.M(), wgt);
                    if ( probe_relIso < 0.05 ) {
                        ptProbeRiso5->Fill(theProbe.pt(), wgt);
                        etaProbeRiso5->Fill(theProbe.eta(), wgt);
                        phiProbeRiso5->Fill(theProbe.phi(), wgt);
                        mMuMuTPRiso5->Fill(mumu.M(), wgt);
                    }
                }
            }
        }
    }
    

}// end of fill()

void HeavyNu::HistPerDef::fill(const pat::Muon& theTag,
                               const pat::GenericParticle& theProbe,
                               double wgt)
{
    pat::Muon muonProbe ;
    muonProbe.setP4( theProbe.p4() ) ;
    muonProbe.setTrackIso( theProbe.trackIso() ) ;
    HistPerDef::fill( theTag,muonProbe,wgt ) ;
}// end of fill()


//
// constants, enums and typedefs
//
const std::vector<hNuMassHypothesis> v_null;

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
    disableTriggerCorrection_ = iConfig.getParameter<bool>("DisableTriggerCorrection");

    applyJECUsign_ = iConfig.getParameter<int>("applyJECUsign");
    if(applyJECUsign_) applyJECUsign_ /= abs(applyJECUsign_); // ensure -1,0,+1

    applyMESfactor_ = iConfig.getParameter<double>("applyMESfactor");

    applyMuIDCorrections_ = iConfig.getParameter<bool>("applyMuIDEffcorr");
    applyMuIDEffsign_ = iConfig.getParameter<int>("applyMuIDEffsign");
    if(applyMuIDEffsign_) applyMuIDEffsign_ /= abs(applyMuIDEffsign_); // ensure -1,0,+1

    highestPtTriggerOnly_ = iConfig.getParameter<bool>("highestPtTriggerOnly");
    applyTrigEffsign_ = iConfig.getParameter<int>("applyTrigEffsign");
    if(applyTrigEffsign_) applyTrigEffsign_ /= abs(applyTrigEffsign_); // ensure -1,0,+1

    studyMuonSelectionEff_ = iConfig.getParameter<bool>("studyMuSelectEff");

    studyScaleFactorEvolution_ = iConfig.getParameter<bool>("studyScaleFactor");

    isPFJets_ = iConfig.getParameter<bool>("isPFJets");

    // ==================== Init other members ====================
    //

    nnif_ = new HeavyNu_NNIF(iConfig);
    trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet > ("trigMatchPset"));
    muid_ = new HeavyNuID(iConfig.getParameter<edm::ParameterSet > ("muIDPset"));
    // ==================== Book the histos ====================
    //
    edm::Service<TFileService> fs;
    hists.mc_type = fs->make<TH1D > ("mc_type", "MC Type Code", 100, -0.5, 99.5);
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
    hists.weights = fs->make<TH1D > ("weights", "event weights", 100, 0.0, 1.0);
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

    labelJetIDaxis(hists.jetID->GetXaxis());

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

        labelMuonQualAxis(hists.muQualVsPt->GetYaxis());

        hists.muTrckIsoVsPt = fs->make<TH2D > ("muTrckIsoVsPt", "trackIso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);trackIso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muHcalIsoVsPt = fs->make<TH2D > ("muHcalIsoVsPt", "HCAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);HCAL Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muEcalIsoVsPt = fs->make<TH2D > ("muEcalIsoVsPt", "ECAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);ECAL Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muCaloIsoVsPt = fs->make<TH2D > ("muCaloIsoVsPt", "Calo Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);Calo Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
    }

    // Histos per cut:
    //
    // hists.twoL.book(new TFileDirectory(fs->mkdir("cutm1_LL")), "(two lepton:-1)", v_null);
    hists.noCuts.book(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)", v_null);
    hists.LLptCuts.book(new TFileDirectory(fs->mkdir("cutX_LLpt")), "(dileptons with ptcuts:X)", v_null);
    // hists.MuTightCuts.book(new TFileDirectory(fs->mkdir("cut2_MuTight")), "(Mu tight cuts:2)", v_null);
    hists.LLJJptCuts.book(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)", nnif_->masspts());
    hists.TrigMatches.book(new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:2)", v_null);
    hists.VertexCuts.book(new TFileDirectory(fs->mkdir("cut3_Vertex")), "(vertex requirements:3)", nnif_->masspts());
    hists.Mu1HighPtCut.book(new TFileDirectory(fs->mkdir("cut4_Mu1HighPt")), "(Mu1 High pt cut:4)", nnif_->masspts());
    hists.loDiLmassCut.book(new TFileDirectory(fs->mkdir("cut5a_loDiLmass")), "(mumu mass cut:5a)", nnif_->masspts());
    hists.diLmassCut.book(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(mumu mass cut:5)", nnif_->masspts());
    hists.mWRmassCut.book(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(mumujj mass cut:6)", nnif_->masspts());

    if(studyScaleFactorEvolution_)
    {
        hists.Mu1Pt30GeVCut.book(new TFileDirectory(fs->mkdir("Mu1Pt30GeV")), "(Mu1 30 GeV pt cut)", v_null);
        hists.Mu1Pt40GeVCut.book(new TFileDirectory(fs->mkdir("Mu1Pt40GeV")), "(Mu1 40 GeV pt cut)", v_null);
        hists.Mu1Pt50GeVCut.book(new TFileDirectory(fs->mkdir("Mu1Pt50GeV")), "(Mu1 50 GeV pt cut)", v_null);
        hists.Mu1Pt60GeVCut.book(new TFileDirectory(fs->mkdir("Mu1Pt60GeV")), "(Mu1 60 GeV pt cut)", v_null);
        hists.Mu1Pt80GeVCut.book(new TFileDirectory(fs->mkdir("Mu1Pt80GeV")), "(Mu1 80 GeV pt cut)", v_null);
        hists.Mu1Pt100GeVCut.book(new TFileDirectory(fs->mkdir("Mu1Pt100GeV")), "(Mu1 100 GeV pt cut)", v_null);

        hists.Mu1HighPtCutVtxEq1.book(new TFileDirectory(fs->mkdir("Mu1HighPtVtxEq1")), "(Mu1 60 GeV pt cut, 1 vtx)", v_null);
        hists.Mu1HighPtCutVtx2to5.book(new TFileDirectory(fs->mkdir("Mu1HighPtVtx2to5")), "(Mu1 60 GeV pt cut, 2-5 vtx)", v_null);
        hists.Mu1HighPtCutVtxGt5.book(new TFileDirectory(fs->mkdir("Mu1HighPtVtxGt5")), "(Mu1 60 GeV pt cut, 6+ vtx)", v_null);

        hists.Mu1HighPtCutNoJets.book(new TFileDirectory(fs->mkdir("Mu1HighPtNoJets")), "(Mu1 60 GeV pt cut, no jets)", v_null);
        hists.Mu1HighPtCut1Jet.book(new TFileDirectory(fs->mkdir("Mu1HighPt1Jet")), "(Mu1 60 GeV pt cut, 1 jet)", v_null);
    }

    if(trig_->matchingEnabled())
    {
        hists.Mu1TrigMatchesInZwin.book(new TFileDirectory(fs->mkdir("Mu1TrigMatchesInZwin")), "(#mu1 trigger match in Z mass Window)", v_null);
        hists.Mu2TrigMatchesInZwin.book(new TFileDirectory(fs->mkdir("Mu2TrigMatchesInZwin")), "(#mu2 Trigger match in Z mass Window)", v_null);
        hists.Mu1Mu2TrigMatchesInZwin.book(new TFileDirectory(fs->mkdir("Mu1Mu2TrigMatchesInZwin")), "(#mu1,#mu2 Trigger match in Z mass Window)", v_null);
        trig_->book(*(hists.Mu1TrigMatchesInZwin.mydir), &(hists.Mu1TrigMatchesInZwin.trigHistos));
        trig_->book(*(hists.Mu2TrigMatchesInZwin.mydir), &(hists.Mu2TrigMatchesInZwin.trigHistos));
    }

    if(studyMuonSelectionEff_)
    {
        hists.TightTagTrackProbeInZwin.bookTagProbe(new TFileDirectory(fs->mkdir("TightTagTrackProbeInZwin")), "(probe, ID in Z mass Window)");
        hists.TightTagTrackProbePassesInZwin.bookTagProbe(new TFileDirectory(fs->mkdir("TightTagTrackProbePassesInZwin")), "(probe passes, ID in Z mass Window)");
        hists.TightTagTightProbeInZwin.bookTagProbe(new TFileDirectory(fs->mkdir("TightTagTightProbeInZwin")), "(probe, iso in Z mass Window)");
        hists.TightTagTightProbePassesInZwin.bookTagProbe(new TFileDirectory(fs->mkdir("TightTagTightProbePassesInZwin")), "(probe passes, iso in Z mass Window)");
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

    MCweightByVertex_ = edm::LumiReWeighting(hnu::generate_flat10_mc(50),
            hnu::get_standard_pileup_data(pileupEra_, 50));

    // For the record...
    std::cout << "Configurable cut values applied:" << std::endl;
    std::cout << "muonTag           = " << muonTag_ << std::endl;
    std::cout << "jetTag            = " << jetTag_ << std::endl;
    std::cout << "metTag            = " << metTag_ << std::endl;
    std::cout << "electronTag       = " << elecTag_ << std::endl;
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
        jecuHi = jecTotalUncertainty(jpt, jeta, jecuObj_, jecVal_, isBjet, true);
        jecuLo = jecTotalUncertainty(jpt, jeta, jecuObj_, jecVal_, isBjet, false);
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

// Study the ID/isolation efficiency using generic tracks 
// Input HNE may not have two valid muons.  
// If one muon is tight then there is only one tag.  If both are tight, you have two tags
void HeavyNu::studyMuonSelectionEff(std::vector<pat::Muon>& tightMuons, 
                                    std::vector<pat::Jet>& jets,
                                    edm::Handle<pat::GenericParticleCollection>& pTracks,
                                    edm::Handle<reco::TrackCollection>& gTracks,
                                    edm::Handle<reco::BeamSpot>& beamspot, 
                                    std::vector<bool> isTag,
                                    double wgt) {
    
    if ( tightMuons.size() < 1 || jets.size() < 2 ) return ;

    pat::GenericParticleCollection trackCands = *(pTracks.product());
    std::sort(trackCands.begin(), trackCands.end(), hnu::pTcompare());

    reco::TrackCollection generalTracks = *(gTracks.product()); 

    int nTrackProbes = 0 ; int nMuonIsoProbes = 0 ; 
    for (unsigned int i=0; i<tightMuons.size(); i++) {

        // Tag is ALWAYS tight, trigger-matched, and isolated
        if ( !isTag.at(i) ) continue ;
        pat::Muon theTag = tightMuons.at(i) ;
        if ( hnu::muIsolation(theTag,1.0) >= cuts.muon_trackiso_limit ) continue ; 

        // First case:
        // Probe is isolated track
        for (unsigned int j=0; j<trackCands.size(); j++) {
            pat::GenericParticle theProbe = trackCands.at(j) ; 
            if ( theProbe.pt() <= cuts.minimum_mu2_pt ) break; // Sorted collection: quit once below
            if ( fabs(theProbe.eta()) >= cuts.maximum_mu_abseta ) continue;
            // Remove tag/probe overlaps
            if ( deltaR(theTag.eta(), theTag.phi(), theProbe.eta(), theProbe.phi()) < 0.02 &&
                 fabs(theTag.pt()-theProbe.pt())/theTag.pt() < 0.05 ) continue ;  
            // Only consider Z window
            double tpMass = (theTag.p4() + theProbe.p4()).M() ;
            if ( !inZmassWindow(tpMass) ) continue ; 
            // Remove probes that overlap with valid jets
            bool probeJetOverlap = false ; 
            for (unsigned int k=0; k<jets.size(); k++) 
                if ( deltaR(jets.at(k).eta(), jets.at(k).phi(),
                            theProbe.eta(), theProbe.phi()) < cuts.minimum_muon_jet_dR ) probeJetOverlap = true ;
            if ( probeJetOverlap ) continue ;

            // Confirmed: we have a valid probe
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
            theProbe.setTrackIso( trkSumPtIsoCone ) ; 
            hists.TightTagTrackProbeInZwin.fill( theTag,theProbe,wgt ) ;
            nTrackProbes++ ;
            // Does the probe pass tight ID selection?
            unsigned int muIdx = tightMuons.size() ;
            for (unsigned int k=0; k<tightMuons.size(); k++) {
                if ( muIdx < tightMuons.size() ) break ; 
                if ( k == i ) continue ;
                pat::Muon tightMuon = tightMuons.at(k) ; 
                if ( deltaR(tightMuon.eta(), tightMuon.phi(), theProbe.eta(), theProbe.phi()) < 0.02 &&
                     fabs(theProbe.pt()-tightMuon.pt())/theProbe.pt() < 0.05 ) {
                    muIdx = k ;
                    std::cout << "Isolation comparison for track and matched tight muon: "
                              << trkSumPtIsoCone << " vs. " << tightMuon.trackIso() << std::endl ;
                }
            }
            if ( muIdx < tightMuons.size() ) 
                hists.TightTagTrackProbePassesInZwin.fill( theTag,theProbe,wgt ) ;
        }
        
        // Second case:
        // Probe is tight muon
        for (unsigned int j=0; j<tightMuons.size(); j++) {
            // Remove tag/probe overlaps
            if ( j == i ) continue ;
            pat::Muon theProbe = tightMuons.at(j) ; 
            // Only consider Z window
            double tpMass = (theTag.p4() + theProbe.p4()).M() ;
            if ( !inZmassWindow(tpMass) ) continue ; 
            // Remove probes that overlap with valid jets
            bool probeJetOverlap = false ; 
            for (unsigned int k=0; k<jets.size(); k++) 
                if ( deltaR(jets.at(k).eta(), jets.at(k).phi(),
                            theProbe.eta(), theProbe.phi()) < cuts.minimum_muon_jet_dR ) probeJetOverlap = true ;
            if ( probeJetOverlap ) continue ;

            // Confirmed: we have a valid probe
            nMuonIsoProbes++ ; 
            hists.TightTagTightProbeInZwin.fill( theTag,theProbe,wgt ) ;

            if ( hnu::muIsolation(theProbe,1.0) < cuts.muon_trackiso_limit ) 
                hists.TightTagTightProbePassesInZwin.fill( theTag,theProbe,wgt ) ;
        }
    }
    if ( nTrackProbes > 2 )   std::cout << "WARNING: Found more track probes than expected --> " << nTrackProbes << std::endl ; 
    if ( nMuonIsoProbes > 2 ) std::cout << "WARNING: Found more muon probes than expected  --> " << nMuonIsoProbes << std::endl ; 
            
}

void HeavyNu::studyIsolation(const std::vector<pat::Muon>& muons,
                             const std::vector< std::pair<pat::Jet, float> >& jets,
                             bool mu1tag, bool mu2tag, double weight) {
    //
    // Now look for "close jet" effects
    //
    if(muons.size() < 2) return; // Need both muons for this
    pat::Muon mu1 = muons.at(0);
    pat::Muon mu2 = muons.at(1);
    double m12 = (mu1.p4() + mu2.p4()).M();
    if(!inZmassWindow(m12)) return;
    // With at least one tagged: same ordering in muon list, hne
    if(!mu1tag && !mu2tag) return;

    bool mu1close = false;
    bool mu2close = false;
    for(unsigned int i = 0; i < jets.size(); i++)
    {
        pat::Jet j = jets.at(i).first;
        double dRj1 = deltaR(j.eta(), j.phi(), mu1.eta(), mu1.phi());
        double dRj2 = deltaR(j.eta(), j.phi(), mu2.eta(), mu2.phi());

        if(mu2tag && !mu1close && dRj1 < 0.08) mu1close = true;
        if(mu1tag && !mu2close && dRj2 < 0.08) mu2close = true;
    }

    if(mu1close)
    {
        hists.closejetMu2tagMu1probeInZwin->Fill(mu1.pt(), weight);
        hists.trkIsoStudy->Fill(mu1.trackIso() / mu1.pt(), weight);
    }
    if(mu1close && hnu::muIsolation(mu1, 1.0) < cuts.muon_trackiso_limit)
        hists.closejetMu2tagMu1passInZwin->Fill(mu1.pt(), weight);
    if(mu2close)
    {
        hists.closejetMu1tagMu2probeInZwin->Fill(mu2.pt(), weight);
        hists.trkIsoStudy->Fill(mu2.trackIso() / mu2.pt(), weight);
    }
    if(mu2close && hnu::muIsolation(mu2, 1.0) < cuts.muon_trackiso_limit)
        hists.closejetMu1tagMu2passInZwin->Fill(mu2.pt(), weight);
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

bool HeavyNu::passesTrigger(const double mu1pt, const double mu2pt,
        const bool mu1trig, const bool mu2trig,
        const uint32_t run)
{

    // Low luminosity running: 24 GeV single muon trigger
    if(run < 165000) return( mu1trig || mu2trig);
    // 2011 bulk running: 40 GeV single muon trigger
    return( (mu1trig && mu1pt > 40.) || (mu2trig && mu2pt > 40.));
}

// ------------ method called to for each event  ------------

bool HeavyNu::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    HeavyNuEvent hnuEvent;

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

    // std::cout << "Init event weight is: " << hnuEvent.eventWgt << std::endl ; 
    // if ( hnuEvent.eventWgt < 0.0001 || fabs(hnuEvent.eventWgt) > 1000. ) std::cout << evtCounter << std::endl ; 

    if(hnuEvent.isMC)
    {
        edm::Handle<std::vector<PileupSummaryInfo> > pPU;
        iEvent.getByLabel("addPileupInfo", pPU);
        std::pair<float, double> pileup = hnu::pileupReweighting(pPU, MCweightByVertex_);
        hnuEvent.n_pue = pileup.first;
        hnuEvent.eventWgt *= pileup.second;
        //Shirpa reweighting
        hnuEvent.eventWgt *= geneventinfo->weight();
        hists.weights->Fill(hnuEvent.eventWgt);
        // generator information
        edm::Handle<reco::GenParticleCollection> genInfo;
        if(iEvent.getByLabel("genParticles", genInfo))
        {
            hnuEvent.decayID(*genInfo);
        }
    }
    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    hnuEvent.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle);

    // std::cout << "PU event weight is:   " << hnuEvent.eventWgt << std::endl ; 
    // if ( hnuEvent.eventWgt < 0.0001 || fabs(hnuEvent.eventWgt) > 1000. ) std::cout << evtCounter << std::endl ; 

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
        // For Fall10 corrections...not yet applied
        // edm::FileInPath corrFile("HeavyNu/AnalysisModules/test/Jec10V3_Uncertainty_AK5Calo.txt") ;
        // jecuObj_ = new JetCorrectionUncertainty(corrFile.fullPath()) ;

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
            int pileupYear = pileupEra_ / 10;
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

    hists.mc_type->Fill(hnuEvent.mc_class);
    hists.nelec ->Fill(pElecs->size());
    hists.nmuAll->Fill(pMuons->size());
    hists.njet ->Fill(pJets->size());
    hists.nmet ->Fill(pMET->size());

    if(hnuEvent.isMC) studyJetVertex(pJets, jptJets, pMuons, hnuEvent.n_pue);
    else studyJetVertex(pJets, jptJets, pMuons, hnuEvent.n_primary_vertex);

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

    // Basic selection requirements: Require at least two muons, two jets
    if(pMuons->size() >= 2 && pJets->size() >= 2)
    {
        hists.cutlevel->Fill(0.0, hnuEvent.eventWgt);
        hists.noCuts.fill(*pMuons, *pJets, *pMET, hnuEvent.isMC, hnuEvent.eventWgt, isPFJets_);
    }
    else return false;

    // Look for valid jets and put them in the event
    std::vector< std::pair<pat::Jet, float> > jetCands =
            hnu::getJetList(pJets, jecuObj_, cuts.minimum_jet_pt, cuts.maximum_jet_abseta, applyJECUsign_, jecVal_);
    hnuEvent.nJets = jetCands.size();

    // Look for valid muons and put them in the event
    std::vector<pat::Muon> muCands =
        hnu::getMuonList(pMuons, tevMuons, cuts.minimum_mu2_pt, cuts.maximum_mu_abseta, applyMESfactor_);
    if (muCands.size() >= 2) { 
        hists.LLptCuts.fill(muCands, *pJets, *pMET, hnuEvent.isMC, hnuEvent.eventWgt, isPFJets_);
    } else return false ; 
        
    if(hnuEvent.nJets < 2) return false;

    hnuEvent.j1 = jetCands.at(0).first;
    hnuEvent.j2 = jetCands.at(1).first;
    hnuEvent.j1scale = jetCands.at(0).second;
    hnuEvent.j2scale = jetCands.at(1).second;

    hnuEvent.tjV1 = hnu::caloJetVertex(hnuEvent.j1, *jptJets);
    hnuEvent.tjV2 = hnu::caloJetVertex(hnuEvent.j2, *jptJets);


    // for (unsigned int i=0; i<muCands.size(); i++) { 
    //   std::cout << "muon " << i+1 << " of " << muCands.size() << " with pT " << muCands.at(i).pt() << std::endl ; 
    // }

    for(unsigned int i = 0; i < muCands.size(); i++)
    {
        if(hnuEvent.nMuons == 2) break;
        pat::Muon iM = muCands.at(i);
        if(hnu::muIsolation(iM, applyMESfactor_) < cuts.muon_trackiso_limit)
        {
            double dRj1 = deltaR(iM.eta(), iM.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
            double dRj2 = deltaR(iM.eta(), iM.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
            if(dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR)
            {
                hnuEvent.nMuons++;
                if(hnuEvent.nMuons == 1) hnuEvent.mu1 = iM;
                else if(hnuEvent.nMuons == 2) hnuEvent.mu2 = iM;
                else std::cout << "WARNING: Expected empty muon position" << std::endl;
            }
        }
    }

    if(applyMuIDCorrections_ && hnuEvent.isMC)
    {
        double mu1wgt = (hnuEvent.nMuons > 0)?
                (muid_->weightForMC((hnuEvent.mu1.pt() * applyMESfactor_), applyMuIDEffsign_)):1.0;
        double mu2wgt = (hnuEvent.nMuons > 1)?
                (muid_->weightForMC((hnuEvent.mu2.pt() * applyMESfactor_), applyMuIDEffsign_)):1.0;

        hnuEvent.eventWgt *= (mu1wgt * mu2wgt);
    }

    // std::cout << "muid event weight is: " << hnuEvent.eventWgt << std::endl ; 
    // if ( hnuEvent.eventWgt < 0.0001 || fabs(hnuEvent.eventWgt) > 1000. ) std::cout << evtCounter << std::endl ; 

    //--- Trigger Matching needed for efficiency studies ---//
    bool mu1trig = false;
    bool mu2trig = false;
    if(trig_->matchingEnabled() && iEvent.isRealData())
    {
        mu1trig = (hnuEvent.nMuons > 0) &&
                trig_->isTriggerMatched(hnuEvent.mu1, iEvent,
                &(hists.Mu1TrigMatchesInZwin.trigHistos));
        mu2trig = (hnuEvent.nMuons > 1) &&
                trig_->isTriggerMatched(hnuEvent.mu2, iEvent,
                &(hists.Mu2TrigMatchesInZwin.trigHistos));
    }
    else if(!iEvent.isRealData())
    {
        if(disableTriggerCorrection_)
        {
            mu1trig = true;
            mu2trig = true;
        }
        else
        {
            mu1trig = (hnuEvent.nMuons > 0) &&
                    trig_->simulateForMC(applyMESfactor_ * hnuEvent.mu1.pt(), applyTrigEffsign_);
            mu2trig = (hnuEvent.nMuons > 1) &&
                    trig_->simulateForMC(applyMESfactor_ * hnuEvent.mu2.pt(), applyTrigEffsign_);
        }
    }
    // std::cout << "Trigger results: " << mu1trig << ", " << mu2trig << std::endl ; 

    if(studyMuonSelectionEff_)
    {

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
        if ( validJets.size() >= 2 && muCands.size() > 0 ) {
            std::vector<bool> isTrigMuon ; // Muon can only be a tag if it is both tight AND trigger matched
            int nTrigs = 0 ; 
            for (unsigned int i=0; i<muCands.size(); i++) {
                bool isTrig = true ; // Default true for MC, as trigger "matching" is random
                // Overwrite trigger decision for data with the REAL decision
                if ( trig_->matchingEnabled() && iEvent.isRealData() )
                    isTrig = trig_->isTriggerMatched(muCands.at(i), iEvent) ;
                if ( isTrig ) nTrigs++ ; 
                isTrigMuon.push_back( isTrig ) ;
            }
            if ( nTrigs > 0 ) {
                studyMuonSelectionEff( muCands,validJets,pTracks, gTracks, beamSpotHandle, isTrigMuon, hnuEvent.eventWgt);
            }
        }
        if (muCands.size() == 2)
            studyIsolation(muCands, jetCands, mu1trig, mu2trig, hnuEvent.eventWgt);
    }
    if(hnuEvent.nMuons < 2) return false;

    // std::cout << "2m2j event weight is: " << hnuEvent.eventWgt << std::endl ; 
    // if ( hnuEvent.eventWgt < 0.0001 || fabs(hnuEvent.eventWgt) > 1000. ) std::cout << evtCounter << std::endl ; 

    if(hnu::jetID(hnuEvent.j1) < 1 || hnu::jetID(hnuEvent.j2) < 1) return false;

    hists.cutlevel->Fill(1.0, hnuEvent.eventWgt); // Two highest pT muons that are isolated, separated from chosen jets
    hnuEvent.regularize(); // assign internal standards
    hnuEvent.scaleMuE(applyMESfactor_);
    hnuEvent.calculate(); // calculate various details
    hists.LLJJptCuts.fill(hnuEvent, v_null);
    if(pMET->size()) hnuEvent.met1 = pMET->at(0);

    //--- Trigger code needs to be updated...placeholder for now ---//
    if(!passesTrigger(hnuEvent.mu1.pt(), hnuEvent.mu2.pt(),
            mu1trig, mu2trig, iEvent.id().run())) return false;
    hists.cutlevel->Fill(2.0, hnuEvent.eventWgt); // Event meets trigger requirements
    hists.TrigMatches.fill(hnuEvent, v_null);
    //--- Study Trigger Matching efficiency ---//
    if(trig_->matchingEnabled() && iEvent.isRealData())
    {
        if(inZmassWindow(hnuEvent.mMuMu) &&
                (fabs(hnuEvent.mu1.eta()) < 2.1) &&
                (fabs(hnuEvent.mu2.eta()) < 2.1))
        {
            if(mu1trig)
            {
                hists.Mu1TrigMatchesInZwin.fill(hnuEvent, v_null);
                if(mu2trig)
                {
                    hists.Mu2TrigMatchesInZwin.fill(hnuEvent, v_null);
                    hists.Mu1Mu2TrigMatchesInZwin.fill(hnuEvent, v_null);
                }
            }
            else if(mu2trig) hists.Mu2TrigMatchesInZwin.fill(hnuEvent, v_null);
        }
    }

    nnif_->fillvector(hnuEvent);
    nnif_->output(hnuEvent.nnoutputs);

    // hists.LLJJptCuts.fill(hnuEvent, nnif_->masspts());

    //--- Impose vertex requirement here ---//
    float deltaVzJ1J2 = fabs(hnuEvent.tjV1 - hnuEvent.tjV2);
    float deltaVzJ1M1 = fabs(hnuEvent.tjV1 - hnuEvent.mu1.vertex().Z());
    float deltaVzJ2M2 = fabs(hnuEvent.tjV2 - hnuEvent.mu2.vertex().Z());
    float deltaVzJ1M2 = fabs(hnuEvent.tjV1 - hnuEvent.mu2.vertex().Z());
    float deltaVzJ2M1 = fabs(hnuEvent.tjV2 - hnuEvent.mu1.vertex().Z());
    if( (cuts.maxJetVZsepCM > 0) && 
        ((deltaVzJ1J2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M1 >= cuts.maxJetVZsepCM) ||
         (deltaVzJ2M2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M2 >= cuts.maxJetVZsepCM) ||
         (deltaVzJ2M1 >= cuts.maxJetVZsepCM)) )
        return false;
    float deltaVzM1M2 = fabs(hnuEvent.mu1.vertex().Z() - hnuEvent.mu2.vertex().Z());
    if(cuts.maxVertexZsep > 0 && deltaVzM1M2 >= cuts.maxVertexZsep) return false;

    hists.cutlevel->Fill(3.0, hnuEvent.eventWgt); // Event meets vertex requirements
    hists.VertexCuts.fill(hnuEvent, nnif_->masspts());

    if(studyScaleFactorEvolution_)
    {
        double mu1pt = applyMESfactor_ * hnuEvent.mu1.pt();
        if(mu1pt > 40.) hists.Mu1Pt40GeVCut.fill(hnuEvent, v_null);
        if(mu1pt > 50.) hists.Mu1Pt50GeVCut.fill(hnuEvent, v_null);
        if(mu1pt > 60.) hists.Mu1Pt60GeVCut.fill(hnuEvent, v_null);
        if(mu1pt > 80.) hists.Mu1Pt80GeVCut.fill(hnuEvent, v_null);
        if(mu1pt > 100.) hists.Mu1Pt100GeVCut.fill(hnuEvent, v_null);
    }

    if((applyMESfactor_ * hnuEvent.mu1.pt()) < cuts.minimum_mu1_pt)
        return false;

    if(studyScaleFactorEvolution_)
    {
        if(hnuEvent.n_primary_vertex == 1)
            hists.Mu1HighPtCutVtxEq1.fill(hnuEvent, v_null);
        else if(hnuEvent.n_primary_vertex <= 5)
            hists.Mu1HighPtCutVtx2to5.fill(hnuEvent, v_null);
        else if(hnuEvent.n_primary_vertex > 5)
            hists.Mu1HighPtCutVtxGt5.fill(hnuEvent, v_null);
    }
    hists.cutlevel->Fill(4.0, hnuEvent.eventWgt); // Event meets high muon pT requirements
    hists.Mu1HighPtCut.fill(hnuEvent, nnif_->masspts());

    if(hnuEvent.mMuMu < 40) return false; // Sanity check...remove low mass points
    hists.loDiLmassCut.fill(hnuEvent, nnif_->masspts());

    if(hnuEvent.mMuMu < cuts.minimum_mumu_mass) return false; // dimuon mass cut
    hists.cutlevel->Fill(5.0, hnuEvent.eventWgt); // Event meets dimuon mass requirements
    hists.diLmassCut.fill(hnuEvent, nnif_->masspts());

    if(iEvent.isRealData())
    {
        std::cout << "\t" << iEvent.id() << std::endl;
        std::cout << "\tM(W_R)  = " << hnuEvent.mWR << " GeV";
        std::cout << ", M(NuR1) = " << hnuEvent.mNuR1 << " GeV";
        std::cout << ", M(NuR2) = " << hnuEvent.mNuR2 << " GeV" << std::endl;
        std::cout << "\tM(mumu) = " << hnuEvent.mMuMu << " GeV";
        std::cout << ", M(JJ) = " << hnuEvent.mJJ << " GeV" << std::endl;
        std::cout << "\tJets:   j1 ";
        outputCandidate(hnuEvent.j1);
        std::cout << ", j2 ";
        outputCandidate(hnuEvent.j2);
        std::cout << std::endl;
        std::cout << "\tMuons: mu1 ";
        outputCandidate(hnuEvent.mu1);
        std::cout << ", mu2 ";
        outputCandidate(hnuEvent.mu2);
        std::cout << std::endl;
    }

    // Change the final logic of the filter based on LQ meeting discussion:
    // Interest in seeing events that pass the dilepton mass requirement
    // if ( hnuEvent.mWR<cuts.minimum_mWR_mass ) return false;  // 4-object mass cut
    if(hnuEvent.mWR >= cuts.minimum_mWR_mass)
    {
        hists.cutlevel->Fill(6.0, hnuEvent.eventWgt); // Event meets W_R mass requirements
        hists.mWRmassCut.fill(hnuEvent, nnif_->masspts());
    }
    return true;
}

// ------------ method called once each job just before starting event loop  ------------

void HeavyNu::beginJob()
{
    nnif_->beginJob();
    firstEvent_ = true;
    evtCounter = 0;
}

// ------------ method called once each job just after ending the event loop  ------------

void HeavyNu::endJob()
{
    nnif_->endJob();
    trig_->endJob();
    muid_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNu);


