// -*- C++ -*-
//
// Package:    HeavyNu
// Class:      HeavyNuFNAL
// 
/**\class HeavyNuFNAL HeavyNuFNAL.cc HeavyNu/AnalyzerModules/src/HeavyNuFNAL.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: HeavyNuFNAL.cc,v 1.2 2011/11/18 13:28:16 pastika Exp $
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
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


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

template <class T> void outputCandidate(const T& p) { 
// inline void outputCandidate(const reco::CandidateBaseRef& rc)
// {
    std::cout << "pt=" << p.pt() << " GeV, eta=" << p.eta() << ", phi=" << p.phi();
}

//============================================================
// From the JetEnergyScale twiki:
//

//============================================================


static std::string btagName;
static double minBtagDiscVal; // for discriminating B-tagged jets.

class HeavyNuFNAL : public edm::EDFilter
{
public:
    explicit HeavyNuFNAL(const edm::ParameterSet&);
    ~HeavyNuFNAL();


private:

    virtual void respondToOpenInputFile(edm::FileBlock const& fb)
    {
        currentFile_ = fb.fileName();
    }

    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    virtual TH1 *bookRunHisto(uint32_t runNumber);


  bool passesTrigger(const double mu1pt,const double mu2pt,
		     const bool mu1trig,const bool mu2trig, 
		     const uint32_t run) ; 


    edm::InputTag muonTag_;
    edm::InputTag jetTag_;
    edm::InputTag metTag_;

    int evtCounter ; 

    double ZwinMinGeV_, ZwinMaxGeV_; // for trigger efficiency studies

    int jecVal_; // Jet correction eras: 0 = 2010(A+B), 1 = 2010A, 2 = 2010B, 3 = 2011A
    bool highestPtTriggerOnly_;

    bool isPFJets_; // true if PFJets are used (turns off jet id requirement)
    bool useTrackerPt_ ; 

    std::string currentFile_;
    bool dolog_;

    HeavyNuTrigger *trig_;
    JetCorrectionUncertainty *jecuObj_;

    std::map<uint32_t, TH1 *> m_runHistos_;

    // ----------member data ---------------------------

    struct HistPerDef
    {
        //book histogram set w/ common suffix inside the provided TFileDirectory
        void book(TFileDirectory *, const std::string&);
        // fill all histos of the set with the two electron candidates
        void fill(pat::MuonCollection muons,
                pat::JetCollection jets,
                pat::METCollection metc,
                bool isMC,
                double wgt,
                bool pfJets);
        // fill all histos of the set with the two electron candidates
        void fill(const HeavyNuEvent& hne);

        TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2;
        TH2 *ptMu1VsPtMu2ss, *ptMu1VsPtMu2os;
        TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2;
        TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2;
        TH1 *dEtaMu, *dPhiMu, *dEtaJet, *dPhiJet;
        TH2 *dEtaPhiMu, *dEtaPhiJet;
        TH1 *dRminMu1jet, *dRminMu2jet, *dRminMuJet;
        TH1 *hptrelMu1, *hptrelMu2;
        TH2 *ptrelVsdRminMu1jet, *ptrelVsdRminMu2jet;

        TH1 *dptMu1gen, *dptMu2gen;
        TH1 *dRMu1gen, *dRMu2gen;

        TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1caloIso, *mu1dB;
        TH1 *mu2trackIso, *mu2hcalIso, *mu2ecalIso, *mu2caloIso, *mu2dB;

        TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1caloRelIso;
        TH1 *mu2trackRelIso, *mu2hcalRelIso, *mu2ecalRelIso, *mu2caloRelIso;

        TH1 *mMuMu, *mMuMuOS, *mMuMuSS, *diMuCharge, *mMuMuZoom, *mMuMuGenZoom;
        TH1 *mWR, *mWR_1b, *mWR_2b, *mNuR1, *mNuR1_1b, *mNuR1_2b, *mNuR2, *mNuR2_1b, *mNuR2_2b;
        TH1 *mJJ, *mJJ_1b, *mJJ_2b;
        TH2 *mNuR2D, *jetPtvsNum;
        TH1 *mu1ptFracWRmass; 

        TH1* btagJet1, *btagJet2;
        TH1* numBjets;

        TH1* met;

        TH1* czeta_mumu;
        TH1* czeta_mumu_zoom;

        TH1 *mu1jj_surfarea, *mu2jj_surfarea ; 

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
        TH1 *njet, *nmet, *nmuAll, *nmuLoose, *nmuTight;
        TH1 *muPt, *muEta, *muPhi, *looseMuPt, *tightMuPt;
        TH1* mc_type;

      TH1* cutlevel ; 


        TH1 *jetPt, *jetEta, *jetPhi, *jecUncHi, *jecUncLo, *met;
        TH2 *dVzMuJets, *dVzMuMus;
        TH2 *jetPtvsNum;
        TProfile2D *jecUncHiVsEtaPt, *jecUncLoVsEtaPt;

        TH1 *trkIsoStudy;
        TH1 *closejetMu2tagMu1probeInZwin, *closejetMu2tagMu1passInZwin ; 
        TH1 *closejetMu1tagMu2probeInZwin, *closejetMu1tagMu2passInZwin ; 

        TFileDirectory *rundir;
        HistPerDef Mu1HighPtCut;
        HistPerDef loDiLmassCut;
        HistPerDef diLmassNo70to100Cut;
        HistPerDef diLmassNo70to110Cut;
        HistPerDef diLmassNo60to120Cut;
        HistPerDef diLmassAbove100Cut;
        HistPerDef diLmassAbove120Cut;
        HistPerDef diLmassAbove150Cut;
        HistPerDef diLmassAbove200Cut;
        HistPerDef oneBtag;
        HistPerDef twoBtag;
    } hists;

    struct CutsStruct
    {
        double minimum_mu1_pt;
        double minimum_mu2_pt;
        double minimum_jet_pt;
        double minimum_jet1_pt;
        double maximum_mu_abseta;
        double maximum_jet_abseta;
        double minimum_muon_jet_dR;
        double muon_trackiso_limit;
        double maxVertexZsep;
        double maxJetVZsepCM;
    } cuts;

};



void HeavyNuFNAL::HistPerDef::book(TFileDirectory *td, const std::string& post)
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

    t = "Jet #Delta#eta vs. #Delta#phi ";
    t += post + ";#Delta#eta; #Delta#phi";
    dEtaPhiJet = td->make<TH2D > ("dEtaPhiJet", t.c_str(),
            50, 0, 5, 30, 0, 3.14159);

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
    t = "M(W_{R}), at least one b-jet " + post;
    mWR_1b = td->make<TH1D > ("mWR_1b", t.c_str(), 70, 0, 2800);
    t = "M(W_{R}), two b-jets " + post;
    mWR_2b = td->make<TH1D > ("mWR_2b", t.c_str(), 70, 0, 2800);
    t = "M(N_{R}) with #mu_{1} " + post;
    mNuR1 = td->make<TH1D > ("mNuR1", t.c_str(), 70, 0, 2800);
    t = "M(N_{R}) with #mu_{1} and at least on b-jet " + post;
    mNuR1_1b = td->make<TH1D > ("mNuR1_1b", t.c_str(), 70, 0, 2800);
    t = "M(N_{R}) with #mu_{1} and two b-jets " + post;
    mNuR1_2b = td->make<TH1D > ("mNuR1_2b", t.c_str(), 70, 0, 2800);
    t = "M(N_{R}) with #mu_{2} " + post;
    mNuR2 = td->make<TH1D > ("mNuR2", t.c_str(), 70, 0, 1400);
    t = "M(N_{R}) with #mu_{2} and at least on b-jet " + post;
    mNuR2_1b = td->make<TH1D > ("mNuR2_1b", t.c_str(), 70, 0, 1400);
    t = "M(N_{R}) with #mu_{2} and two b-jets " + post;
    mNuR2_2b = td->make<TH1D > ("mNuR2_2b", t.c_str(), 70, 0, 1400);
    t = "M(N_{R}) #mu_{1} vs. #mu_{2} " + post;
    mNuR2D = td->make<TH2D > ("mNuR2D", t.c_str(), 70, 0, 2800, 70, 0, 1400);

    t = "#mu_{1} p_{T} / M(W_{R}) " + post;
    mu1ptFracWRmass = td->make<TH1D > ("mu1ptFracWRmass", t.c_str(), 100, 0, 1.);

    t = "M(jj) " + post;
    mJJ = td->make<TH1D > ("mJJ", t.c_str(), 50, 0, 2000);
    t = "M(jj) at least on b-jet " + post;
    mJJ_1b = td->make<TH1D > ("mJJ_1b", t.c_str(), 50, 0, 2000);
    t = "M(jj) two b-jets " + post;
    mJJ_2b = td->make<TH1D > ("mJJ_2b", t.c_str(), 50, 0, 2000);
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
    mu1jj_surfarea = td->make<TH1D> ("mu1jj_surfarea", t.c_str(), 80,0.,6.2831853) ; 
    t = "surface area m2jj " + post;
    mu2jj_surfarea = td->make<TH1D> ("mu2jj_surfarea", t.c_str(), 80,0.,6.2831853) ; 

    // vertex histograms
    t = "mumu vtx DeltaZ " + post;
    vtx_mumu     = td->make<TH1D > ("vtxMM", t.c_str(), 200, 0., 2.);
    t = "jj vtx DeltaZ " + post;
    vtx_jj       = td->make<TH1D > ("vtxJJ", t.c_str(), 200, 0., 2.);
    t = "m1j vtx min DeltaZ " + post;
    vtx_min_mu1j = td->make<TH1D > ("vtxM1Jmin", t.c_str(), 200, 0., 2.);
    t = "m2j vtx min DeltaZ " + post;
    vtx_min_mu2j = td->make<TH1D > ("vtxM2Jmin", t.c_str(), 200, 0., 2.);
    t = "mu-j vtx min DeltaZ " + post;
    vtx_min_muj  = td->make<TH1D > ("vtxMuJmin", t.c_str(), 200, 0., 2.);
    t = "mmjj vtx max DeltaZ " + post;
    vtx_max_dist = td->make<TH1D > ("vtxDistmax", t.c_str(), 400, 0., 5.);


}// end of book()

void HeavyNuFNAL::HistPerDef::fill(pat::MuonCollection muons,
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

    if((j0bdisc >= minBtagDiscVal) &&
            (j1bdisc >= minBtagDiscVal)) numBjets->Fill(2., wgt);
    else if((j0bdisc >= minBtagDiscVal) ||
            (j1bdisc >= minBtagDiscVal)) numBjets->Fill(1., wgt);
    else numBjets->Fill(0., wgt);

    dPhiJet->Fill(fabs(deltaPhi(j0.phi(), j1.phi())), wgt);
    dEtaJet->Fill(fabs(j0.eta() - j1.eta()), wgt);
    dEtaPhiJet->Fill(fabs(j0.eta() - j1.eta()),
            fabs(deltaPhi(j0.phi(), j1.phi())), wgt);

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

    bool oneBjet  = (j0bdisc >= minBtagDiscVal) || (j1bdisc >= minBtagDiscVal) ; 
    bool twoBjets = (j0bdisc >= minBtagDiscVal) && (j1bdisc >= minBtagDiscVal) ; 
    if ( oneBjet ) {
      mWR_1b->Fill(vWR.M(), wgt);
      mJJ_1b->Fill(jj.M(), wgt);
      mNuR1_1b->Fill((vJJ + m0.p4()).M(), wgt);
      mNuR2_1b->Fill((vJJ + m1.p4()).M(), wgt);
      if ( twoBjets ) { 
	mWR_2b->Fill(vWR.M(), wgt);
	mJJ_2b->Fill(jj.M(), wgt);
	mNuR1_2b->Fill((vJJ + m0.p4()).M(), wgt);
	mNuR2_2b->Fill((vJJ + m1.p4()).M(), wgt);
      }
    }

}// end of fill()

void HeavyNuFNAL::HistPerDef::fill(const HeavyNuEvent& hne)
{
    double wgt = hne.eventWgt;

    mc_type->Fill(hne.mc_class, wgt);

    // Muons
    double mu1pt = hne.mu1.pt();
    double mu2pt = hne.mu2.pt();

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
                float dpt = (hne.mu[i].pt()) - hne.mu[i].genLepton()->pt();
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

    int jet1id = 0;
    int jet2id = 0;

    // Jets
    if(hne.nJets > 0)
    {
      if(!hne.pfJets) jet1id = hnu::jetID(hne.j1);

        double j1bdisc = hne.j1.bDiscriminator(btagName);

        ptJet1->Fill(hne.j1.pt(), wgt);
        etaJet1->Fill(hne.j1.eta(), wgt);
        phiJet1->Fill(hne.j1.phi(), wgt);
        btagJet1->Fill(j1bdisc, wgt);

        if(hne.nJets > 1)
        {
            double j2bdisc = hne.j2.bDiscriminator(btagName);

            ptJet2->Fill(hne.j2.pt(), wgt);
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

	    bool oneBjet  = (j1bdisc >= minBtagDiscVal) || (j2bdisc >= minBtagDiscVal) ; 
	    bool twoBjets = (j1bdisc >= minBtagDiscVal) && (j2bdisc >= minBtagDiscVal) ; 
	    if ( oneBjet ) {
	      mWR_1b->Fill(hne.mWR, wgt);
	      mJJ_1b->Fill(hne.mJJ, wgt);
	      mNuR1_1b->Fill(hne.mNuR1, wgt);
	      mNuR2_1b->Fill(hne.mNuR2, wgt);
	      if ( twoBjets ) { 
		mWR_2b->Fill(hne.mWR, wgt);
		mJJ_2b->Fill(hne.mJJ, wgt);
		mNuR1_2b->Fill(hne.mNuR1, wgt);
		mNuR2_2b->Fill(hne.mNuR2, wgt);
	      }
	    }

	    mu1ptFracWRmass->Fill( mu1pt/hne.mWR, wgt ) ; 
	    mu1jj_surfarea->Fill( hne.area_1jj, wgt ) ; 
	    mu2jj_surfarea->Fill( hne.area_2jj, wgt ) ; 

            ctheta_jj->Fill(hne.ctheta_jj, wgt);
            ctheta_mu1_jj->Fill(hne.ctheta_mu1_jj, wgt);
            ctheta_mu2_jj->Fill(hne.ctheta_mu2_jj, wgt);
            cthetaz_jj->Fill(hne.cthetaz_jj, wgt);
            cthetaz_mu1_jj->Fill(hne.cthetaz_mu1_jj, wgt);
            cthetaz_mu2_jj->Fill(hne.cthetaz_mu2_jj, wgt);

	    float deltaVzJ1J2 = fabs(hne.tjV1-hne.tjV2);
	    float deltaVzJ1M1 = fabs(hne.tjV1-hne.mu1.vertex().Z());
	    float deltaVzJ2M2 = fabs(hne.tjV2-hne.mu2.vertex().Z());
	    float deltaVzJ1M2 = fabs(hne.tjV1-hne.mu2.vertex().Z());
	    float deltaVzJ2M1 = fabs(hne.tjV2-hne.mu1.vertex().Z());
	    float deltaVzM1M2 = fabs(hne.mu1.vertex().Z()-hne.mu2.vertex().Z());

	    vtx_jj->Fill( deltaVzJ1J2, wgt) ; 
	    float minDeltaVzMu1J = std::min( deltaVzJ1M1,deltaVzJ2M1 ) ; 
	    float minDeltaVzMu2J = std::min( deltaVzJ2M2,deltaVzJ2M2 ) ; 
            vtx_min_mu1j->Fill( minDeltaVzMu1J, wgt );
            vtx_min_mu2j->Fill( minDeltaVzMu2J, wgt );
	    vtx_min_muj->Fill( std::min(minDeltaVzMu1J,minDeltaVzMu2J),wgt ) ; 

	    float maxDeltaVzMuJ1 = std::max( deltaVzJ1M1,deltaVzJ1M2 ) ; 
	    float maxDeltaVzMuJ2 = std::max( deltaVzJ2M1,deltaVzJ2M2 ) ; 
	    float maxDeltaVzMMJJ = std::max( deltaVzM1M2,deltaVzJ1J2 ) ;
	    float maxDeltaVzMuJ  = std::max( maxDeltaVzMuJ1,maxDeltaVzMuJ2 ) ; 
	    vtx_max_dist->Fill( std::max(maxDeltaVzMMJJ,maxDeltaVzMuJ),wgt ) ; 
        }

        dRminMu1jet->Fill(hne.dRminMu1jet, wgt);
        dRminMu2jet->Fill(hne.dRminMu2jet, wgt);

        dRminMuJet->Fill(min(hne.dRminMu1jet, hne.dRminMu2jet), wgt);

        hptrelMu1->Fill(hne.ptrelMu1, wgt);
        hptrelMu2->Fill(hne.ptrelMu2, wgt);

        ptrelVsdRminMu1jet->Fill(hne.dRminMu1jet, hne.ptrelMu1, wgt);
        ptrelVsdRminMu2jet->Fill(hne.dRminMu2jet, hne.ptrelMu2, wgt);
    }

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

HeavyNuFNAL::HeavyNuFNAL(const edm::ParameterSet& iConfig)
{
    // ==================== Get all parameters ====================
    //
    dolog_ = iConfig.getParameter<bool>("DoLog");

    muonTag_  = iConfig.getParameter< edm::InputTag > ("muonTag");
    jetTag_   = iConfig.getParameter< edm::InputTag > ("jetTag");
    metTag_   = iConfig.getParameter< edm::InputTag > ("metTag");

    btagName = iConfig.getParameter<std::string > ("BtagName");
    minBtagDiscVal = iConfig.getParameter<double>("minBtagDiscr");

    cuts.minimum_mu1_pt = iConfig.getParameter<double>("minMu1pt");
    cuts.minimum_mu2_pt = iConfig.getParameter<double>("minMu2pt");
    cuts.minimum_jet_pt  = iConfig.getParameter<double>("minJetPt");
    cuts.minimum_jet1_pt = iConfig.getParameter<double>("minJet1Pt");
    cuts.maximum_mu_abseta = iConfig.getParameter<double>("maxMuAbsEta");
    cuts.maximum_jet_abseta = iConfig.getParameter<double>("maxJetAbsEta");
    cuts.minimum_muon_jet_dR = iConfig.getParameter<double>("minMuonJetdR");
    cuts.muon_trackiso_limit = iConfig.getParameter<double>("muonTrackRelIsoLimit");
    cuts.maxVertexZsep = iConfig.getParameter<double>("maxVertexZsepCM");
    cuts.maxJetVZsepCM = iConfig.getParameter<double>("maxJetVZsepCM");

    jecVal_ = iConfig.getParameter<int>("jecEra");

    highestPtTriggerOnly_ = iConfig.getParameter<bool>("highestPtTriggerOnly");

    isPFJets_ = iConfig.getParameter<bool>("isPFJets");
    useTrackerPt_ = iConfig.getParameter<bool>("useTrackerPt") ; 

    // ==================== Init other members ====================
    //

    trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet > ("trigMatchPset"));
    // ==================== Book the histos ====================
    //
    edm::Service<TFileService> fs;
    hists.mc_type = fs->make<TH1D > ("mc_type", "MC Type Code", 100, -0.5, 99.5);
    hists.nmuAll = fs->make<TH1D > ("nmuAll", "N(#mu^{#pm})", 10, -0.5, 9.5);
    hists.nmuLoose = fs->make<TH1D > ("nmuLoose", "N(#mu^{#pm}) passes Loose", 10, -0.5, 9.5);
    hists.nmuTight = fs->make<TH1D > ("nmuTight", "N(#mu^{#pm}) passes Tight", 10, -0.5, 9.5);
    hists.cutlevel = fs->make<TH1D > ("cutlevel", "Cut Level", 10, -0.5, 9.5);
    hists.cutlevel->GetXaxis()->SetBinLabel(1,"No cuts") ; 
    hists.cutlevel->GetXaxis()->SetBinLabel(2,"mmjj p_{T}") ; 
    hists.cutlevel->GetXaxis()->SetBinLabel(3,"M1 (Trigger)") ; 
    hists.cutlevel->GetXaxis()->SetBinLabel(4,"M2 (Vertex)") ; 
    hists.cutlevel->GetXaxis()->SetBinLabel(5,"M3 (high p_{T})") ; 
    hists.cutlevel->GetXaxis()->SetBinLabel(6,"M4 (M_{#mu#mu})") ; 
    hists.cutlevel->GetXaxis()->SetBinLabel(7,"M5 (M(W_{R})") ; 
    hists.njet = fs->make<TH1D > ("njet", "N(Jet)", 50, -0.5, 49.5);
    hists.nmet = fs->make<TH1D > ("nmet", "N(MET)", 50, -0.5, 49.5);
    hists.met = fs->make<TH1D > ("met", "MET distribution", 100, 0, 2000);
    hists.muPt = fs->make<TH1D > ("muPt", "#mu p_{T} distribution", 100, 0, 2000);
    hists.muEta = fs->make<TH1D > ("muEta", "#mu #eta distribution", 50, -2.5, 2.5);
    hists.muPhi = fs->make<TH1D > ("muPhi", "#mu #phi distribution", 60, -3.14159, 3.14159);
    hists.jetPt = fs->make<TH1D > ("jetPt", "jet p_{T} distribution", 100, 0, 2000);
    hists.jetEta = fs->make<TH1D > ("jetEta", "jet #eta distribution", 50, -5, 5);
    hists.jetPhi = fs->make<TH1D > ("jetPhi", "jet #phi distribution", 60, -3.14159, 3.14159);
    hists.jetPtvsNum = fs->make<TH2D > ("jetPtvsNum", "Jet P_{T} vs. Jet # ", 11, -0.5, 10.5, 200, 0., 2000.);
    hists.dVzMuJets = fs->make<TH2D > ("dVzMuJets", "vertex dz Mu Jets vs npue", 500, -5.0, 5.0, 50, -0.5, 49.5);
    hists.dVzMuMus = fs->make<TH2D > ("dVzMuMus", "vertex dz Mu Mus vs npue", 500, -5.0, 5.0, 50, -0.5, 49.5);

    hists.trkIsoStudy = fs->make<TH1D > ("trkIsoStudy", "Tracker Relative Isolation", 100, 0., 1.);
    hists.closejetMu1tagMu2probeInZwin = fs->make<TH1D> ("closejetMu1tagMu2probeInZwin", "Probe (#mu_{2}) p_{T}", 100, 0., 1000.); 
    hists.closejetMu1tagMu2passInZwin  = fs->make<TH1D> ("closejetMu1tagMu2passInZwin",  "Passing probe (#mu_{2}) p_{T}", 100, 0., 1000.); 
    hists.closejetMu2tagMu1probeInZwin = fs->make<TH1D> ("closejetMu2tagMu1probeInZwin", "Probe (#mu_{1}) p_{T}", 100, 0., 1000.); 
    hists.closejetMu2tagMu1passInZwin  = fs->make<TH1D> ("closejetMu2tagMu1passInZwin",  "Passing probe (#mu_{1}) p_{T}", 100, 0., 1000.); 


    // Histos per cut:
    //
    // hists.LLptCuts.book(new TFileDirectory(fs->mkdir("cut1_LLpt")), "(dileptons with ptcuts:1)");
    // hists.MuTightCuts.book(new TFileDirectory(fs->mkdir("cut2_MuTight")), "(Mu tight cuts:2)");
    hists.Mu1HighPtCut.book(new TFileDirectory(fs->mkdir("cut4_Mu1HighPt")), "(Mu1 High pt cut:4)");
    hists.loDiLmassCut.book(new TFileDirectory(fs->mkdir("cut5a_loDiLmass")), "(mumu mass cut:5a)");
    hists.diLmassNo70to100Cut.book(new TFileDirectory(fs->mkdir("diLcut_no70to100")), "(mumu mass cut, no 70-100)");
    hists.diLmassNo70to110Cut.book(new TFileDirectory(fs->mkdir("diLcut_no70to110")), "(mumu mass cut, no 70-110)");
    hists.diLmassNo60to120Cut.book(new TFileDirectory(fs->mkdir("diLcut_no60to120")), "(mumu mass cut, no 60-120)");
    hists.diLmassAbove100Cut.book(new TFileDirectory(fs->mkdir("diLcut_above100")), "(mumu mass cut, 100 GeV)");
    hists.diLmassAbove120Cut.book(new TFileDirectory(fs->mkdir("diLcut_above120")), "(mumu mass cut, 120 GeV)");
    hists.diLmassAbove150Cut.book(new TFileDirectory(fs->mkdir("diLcut_above150")), "(mumu mass cut, 150 GeV)");
    hists.diLmassAbove200Cut.book(new TFileDirectory(fs->mkdir("diLcut_above200")), "(mumu mass cut, 200 GeV)");


    hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

    init_ = false;


    // For the record...
    std::cout << "Configurable cut values applied:" << std::endl;
    std::cout << "muonTag           = " << muonTag_ << std::endl;
    std::cout << "jetTag            = " << jetTag_ << std::endl;
    std::cout << "metTag            = " << metTag_ << std::endl;
    std::cout << "btagName          = " << btagName << std::endl;
    std::cout << "minBtagDiscr      = " << minBtagDiscVal << std::endl;
    std::cout << "minMu1pt          = " << cuts.minimum_mu1_pt << " GeV" << std::endl;
    std::cout << "minMu2pt          = " << cuts.minimum_mu2_pt << " GeV" << std::endl;
    std::cout << "minJetPt          = " << cuts.minimum_jet_pt << " GeV" << std::endl;
    std::cout << "minJet1Pt         = " << cuts.minimum_jet1_pt << " GeV" << std::endl;
    std::cout << "maxMuAbsEta       = " << cuts.maximum_mu_abseta << std::endl;
    std::cout << "maxJetAbsEta      = " << cuts.maximum_jet_abseta << std::endl;
    std::cout << "minMuonJetdR      = " << cuts.minimum_muon_jet_dR << std::endl;
    std::cout << "muonTrackRelIso   = " << cuts.muon_trackiso_limit << std::endl;
    std::cout << "jecEra            = " << jecVal_ << std::endl;

    std::cout << "isPFJets          = " << isPFJets_ << std::endl;
    std::cout << "useTrackerPt      = " << useTrackerPt_ << std::endl;

}

HeavyNuFNAL::~HeavyNuFNAL()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}



//
// member functions
//


TH1 * HeavyNuFNAL::bookRunHisto(uint32_t runNumber)
{
    std::string runstr = int2str<uint32_t > (runNumber);
    return hists.rundir->make <TH1I > (runstr.c_str(), runstr.c_str(), 1, 1, 2);
}


bool HeavyNuFNAL::passesTrigger(const double mu1pt,const double mu2pt,
			    const bool mu1trig,const bool mu2trig, 
			    const uint32_t run) { 
  
  // Changes made for the cross check
  // Low luminosity running: 15 GeV single isolated muon trigger
  if ( run <= 163261 ) return ( mu1trig || mu2trig ) ; 
  // 2011 bulk running: 24 GeV single isolated muon trigger
  return ( ( mu1trig && mu1pt > 25. ) || ( mu2trig && mu2pt > 25. ) ) ; 

  // // Low luminosity running: 24 GeV single muon trigger
  // if ( run <= 163500 ) return ( mu1trig || mu2trig ) ; 
  // // 2011 bulk running: 40 GeV single muon trigger
  // return ( ( mu1trig && mu1pt > 40. ) || ( mu2trig && mu2pt > 40. ) ) ; 
}

// ------------ method called to for each event  ------------

bool HeavyNuFNAL::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    HeavyNuEvent hnuEvent;
    bool isVerbose = false;

    evtCounter++ ; 

    if(isVerbose) std::cout << "Investigating event: " << iEvent.id() << std::endl ;

    hnuEvent.isMC = !iEvent.isRealData();
    hnuEvent.pfJets = isPFJets_;

    if(iEvent.isRealData())
    {
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

    edm::Handle<pat::JetCollection> pJets;
    iEvent.getByLabel(jetTag_, pJets);

//     edm::Handle<reco::PFJetCollection> pfJets;
//     iEvent.getByLabel(edm::InputTag("ak5PFJets","","RECO"),pfJets) ; 

    edm::Handle<pat::METCollection> pMET;
    iEvent.getByLabel(metTag_, pMET);

    edm::Handle<reco::MuonCollection> tevMuons;
    iEvent.getByLabel("refitMuons", tevMuons); 

    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    hnuEvent.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle) ;

    // std::cout << "PU event weight is:   " << hnuEvent.eventWgt << std::endl ; 
    // if ( hnuEvent.eventWgt < 0.0001 || fabs(hnuEvent.eventWgt) > 1000. ) std::cout << evtCounter << std::endl ; 

    if(!pMuons.isValid() || !pJets.isValid() || !(pMET.isValid() && (pMET->size() > 0)))
    {
        std::cout << "Exiting as valid PAT objects not found" << std::endl;
	std::cout << "Muons:     " << pMuons.isValid() << std::endl ; 
	std::cout << "Jets:      " << pJets.isValid() << std::endl ; 
	std::cout << "MET:       " << pMET.isValid() << std::endl ; 
        return false;
    }
    
    hists.mc_type->Fill(hnuEvent.mc_class);
    hists.nmuAll->Fill(pMuons->size());
    hists.njet ->Fill(pJets->size());
    hists.nmet ->Fill(pMET->size());

    if(pMET->size())
        hists.met->Fill(pMET->at(0).pt());
    else
        hists.met->Fill(0);

    if(isVerbose) std::cout << "Initial info before cuts: " << pMuons->size()
                      << " muons and " << pJets->size() << " jets " << std::endl ;

    // Basic selection requirements: Require at least two muons, two jets
    if (pMuons->size() < 2 || pJets->size() < 2) return false ; 

//     std::cout << "Running through the jet candidates (PF on-the-fly): " << std::endl ;
//     for (unsigned int i=0; i<pJets->size(); i++) {
//         std::cout << "Jet " << i+1 << " of " << pJets->size()
//                   << " with pT " << pJets->at(i).pt() << " eta " << pJets->at(i).eta() << " and phi " << pJets->at(i).phi() << std::endl ; 
//     }
    
//     std::cout << "Running through the jet candidates (AOD PF): " << std::endl ;
//     for (unsigned int i=0; i<pfJets->size(); i++) {
//         std::cout << "Jet " << i+1 << " of " << pfJets->size()
//                   << " with pT " << pfJets->at(i).pt() << " eta " << pfJets->at(i).eta() << " and phi " << pfJets->at(i).phi() << std::endl ; 
//     }
    
    // Look for valid jets and put them in the event
    std::vector< std::pair<pat::Jet,float> > jetCands = 
      hnu::getJetList(pJets,jecuObj_,cuts.minimum_jet_pt,cuts.maximum_jet_abseta,0,jecVal_) ; 
    if(isVerbose) std::cout << "Find " << jetCands.size() << " jets after initial pT, eta requirements" << std::endl ;
    for (unsigned int i=0; i<jetCands.size(); i++) {
        if(isVerbose) std::cout << "Jet " << i+1 << " of " << jetCands.size()
                         << " with pT " << jetCands.at(i).first.pt()
                         << " eta " << jetCands.at(i).first.eta()
                         << " phi " << jetCands.at(i).first.phi() << std::endl ;
    }
    if ( jetCands.size() < 2 ) return false ; 

    std::vector<pat::Muon> muCands = 
      hnu::getMuonList(pMuons,tevMuons,cuts.minimum_mu2_pt,cuts.maximum_mu_abseta,1.0,useTrackerPt_) ; 
    if(isVerbose) std::cout << "There are " << muCands.size() << " muons after initial pt, eta, tight id requirements" << std::endl ;
    if ( muCands.size() < 2 ) return false ; 
    for (unsigned int i=0; i<muCands.size(); i++) { 
      pat::Muon iM = muCands.at(i) ; 
      if ( hnu::muIsolation(iM,1.0) < cuts.muon_trackiso_limit ) {
	hnuEvent.nMuons++ ; 
	if      ( hnuEvent.nMuons == 1 ) hnuEvent.mu1 = iM ; 
	else if ( hnuEvent.nMuons == 2 ) hnuEvent.mu2 = iM ; 
	else    std::cout << "WARNING: Expected empty muon position" << std::endl ; 
      }
    }
    if ( hnuEvent.nMuons < 2 ) return false ; 
    if(isVerbose)
    {
        std::cout << "Muon candidates: " << std::endl ;
        std::cout << "muon 1 with pT " << hnuEvent.mu1.pt() << ", eta " << hnuEvent.mu1.eta()
                  << " and phi " << hnuEvent.mu1.phi() << std::endl ;
        std::cout << "muon 2 with pT " << hnuEvent.mu2.pt() << ", eta " << hnuEvent.mu2.eta()
                  << " and phi " << hnuEvent.mu2.phi() << std::endl ;
    }
    
    for (unsigned int i=0; i<jetCands.size(); i++) { 
      if ( hnuEvent.nJets == 2 ) break ; 
      pat::Jet iJ = jetCands.at(i).first ; 
      double dRj1 = deltaR(iJ.eta(), iJ.phi(), hnuEvent.mu1.eta(), hnuEvent.mu1.phi()) ; 
      double dRj2 = deltaR(iJ.eta(), iJ.phi(), hnuEvent.mu2.eta(), hnuEvent.mu2.phi()) ; 
      if(isVerbose) std::cout << "Jet separation for jet " << i+1 << " of " << jetCands.size()
                     << " with pT " << iJ.pt() << " GeV: dR = " << dRj1
                     << " for muon 1 and dR = " << dRj2 << " for muon 2"
                     << std::endl ;
      if (dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR) { 
	hnuEvent.nJets++ ; 
	if      ( hnuEvent.nJets == 1 ) { hnuEvent.j1 = iJ ; hnuEvent.j1scale = jetCands.at(i).second ; } 
	else if ( hnuEvent.nJets == 2 ) { hnuEvent.j2 = iJ ; hnuEvent.j2scale = jetCands.at(i).second ; } 
	else    std::cout << "WARNING: Expected empty jet position" << std::endl ; 
      }
    }
    if(isVerbose) std::cout << "At least " << hnuEvent.nJets << " jets after mu-jet separation requirements" << std::endl ;
    if ( hnuEvent.nJets < 2 ) return false ; 
    if(isVerbose)
    {
        std::cout << "Jet candidates: " << std::endl ;
        std::cout << "jet 1 with pT " << hnuEvent.j1.pt() << ", eta " << hnuEvent.j1.eta()
                  << " and phi " << hnuEvent.j1.phi() << std::endl ;
        std::cout << "jet 2 with pT " << hnuEvent.j2.pt() << ", eta " << hnuEvent.j2.eta()
                  << " and phi " << hnuEvent.j2.phi() << std::endl ;
    }

    hnuEvent.tjV1 = hnu::caloJetVertex(hnuEvent.j1, *jptJets);
    hnuEvent.tjV2 = hnu::caloJetVertex(hnuEvent.j2, *jptJets);

    //--- Trigger Matching needed for efficiency studies ---//
    bool mu1trig = false ; bool mu2trig = false ; 
    if (trig_->matchingEnabled() && iEvent.isRealData()) {
      mu1trig = (hnuEvent.nMuons > 0) && 
          trig_->isTriggerMatched(hnuEvent.mu1, iEvent) ; 
      mu2trig = (hnuEvent.nMuons > 1) && 
          trig_->isTriggerMatched(hnuEvent.mu2, iEvent) ; 
    }

    if(isVerbose)  std::cout << "Trigger results for muon 1: " << (mu1trig?"Y":"N")
                     << ", for muon 2: " << (mu2trig?"Y":"N") << std::endl ;
    //if ( !passesTrigger(hnuEvent.mu1.pt(),hnuEvent.mu2.pt(),
    //                    mu1trig,mu2trig,iEvent.id().run())) return false;
    
    if(isVerbose)  std::cout << "Jet ID values: " << hnu::jetID(hnuEvent.j1) << " (j1), "
                      << hnu::jetID(hnuEvent.j2) << " (j2)" << std::endl ;
    if ( hnu::jetID(hnuEvent.j1) < 1 || hnu::jetID(hnuEvent.j2) < 1 ) return false ; 

    if ( hnuEvent.j1.pt()< cuts.minimum_jet1_pt ) return false ; 

    hnuEvent.regularize(); // assign internal standards
    hnuEvent.scaleMuE();
    hnuEvent.calculate(); // calculate various details

    if(isVerbose)
    {
        std::cout << "*** Event survives basic requirements" << std::endl;
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

    if (pMET->size()) hnuEvent.met1 = pMET->at(0) ;

    if ( cuts.maxJetVZsepCM < 0 && cuts.maxVertexZsep < 0 ) {
        if(isVerbose) std::cout << "Ignoring vertex requirements" << std::endl ;
    } else {
        if(isVerbose) std::cout << "Checking vertex requirements" << std::endl ;

      //--- Impose vertex requirement here ---//
      float deltaVzJ1J2 = fabs(hnuEvent.tjV1-hnuEvent.tjV2);
      float deltaVzJ1M1 = fabs(hnuEvent.tjV1-hnuEvent.mu1.vertex().Z());
      float deltaVzJ2M2 = fabs(hnuEvent.tjV2-hnuEvent.mu2.vertex().Z());
      float deltaVzJ1M2 = fabs(hnuEvent.tjV1-hnuEvent.mu2.vertex().Z());
      float deltaVzJ2M1 = fabs(hnuEvent.tjV2-hnuEvent.mu1.vertex().Z());

      if(isVerbose)
        {
            std::cout << "j1j2 " << deltaVzJ1J2 << ", " << (deltaVzJ1J2 >= cuts.maxJetVZsepCM) << std::endl;
            std::cout << "j1m1 " << deltaVzJ1M1 << ", " << (deltaVzJ1M1 >= cuts.maxJetVZsepCM) << std::endl;
            std::cout << "j1m2 " << deltaVzJ1M2 << ", " << (deltaVzJ1M2 >= cuts.maxJetVZsepCM) << std::endl;
            std::cout << "j2m1 " << deltaVzJ2M1 << ", " << (deltaVzJ2M1 >= cuts.maxJetVZsepCM) << std::endl;
            std::cout << "j2m2 " << deltaVzJ2M2 << ", " << (deltaVzJ2M2 >= cuts.maxJetVZsepCM) << std::endl;
      }

      if((deltaVzJ1J2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M1 >= cuts.maxJetVZsepCM) ||
	 (deltaVzJ2M2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M2 >= cuts.maxJetVZsepCM) ||
	 (deltaVzJ2M1 >= cuts.maxJetVZsepCM))
	return false;
      float deltaVzM1M2 = fabs(hnuEvent.mu1.vertex().Z()-hnuEvent.mu2.vertex().Z());
      if(isVerbose) std::cout << "m1m2 " << deltaVzM1M2 << ", " << (deltaVzM1M2 >= cuts.maxVertexZsep) << std::endl ;
      if (deltaVzM1M2 >= cuts.maxVertexZsep) return false ; 
    
      if(isVerbose) std::cout << "Event passes vertex requirements" << std::endl ;
    }

    // Still want to keep events if WR candidate...just for later looks
    // if(iEvent.isRealData() && isWrCand)
    // {
    //     std::cout << "\t" << iEvent.id() << std::endl;
    //     std::cout << "\tM(W_R)  = " << hnuEvent.mWR << " GeV";
    //     std::cout << ", M(NuR1) = " << hnuEvent.mNuR1 << " GeV";
    //     std::cout << ", M(NuR2) = " << hnuEvent.mNuR2 << " GeV" << std::endl;
    //     std::cout << "\tM(mumu) = " << hnuEvent.mMuMu << " GeV";
    //     std::cout << ", M(JJ) = " << hnuEvent.mJJ << " GeV" << std::endl;
    //     std::cout << "\tJets:   j1 ";
    //     outputCandidate(hnuEvent.j1);
    //     std::cout << ", j2 ";
    //     outputCandidate(hnuEvent.j2);
    //     std::cout << std::endl;
    //     std::cout << "\tMuons: mu1 ";
    //     outputCandidate(hnuEvent.mu1);
    //     std::cout << ", mu2 ";
    //     outputCandidate(hnuEvent.mu2);
    //     std::cout << std::endl;
    // }

    if(hnuEvent.mu1.pt() < cuts.minimum_mu1_pt)
        return false;

    if(isVerbose)
    {
        std::cout << "---> Mu1 pT exceeds " << cuts.minimum_mu1_pt
	      << ": " << hnuEvent.mu1.pt() << std::endl ;
    }

    hists.cutlevel->Fill(4) ; // Event meets high muon pT requirements
    hists.Mu1HighPtCut.fill(hnuEvent);

    if(hnuEvent.mMuMu < 40) return false; // Sanity check...remove low mass points
    hists.loDiLmassCut.fill(hnuEvent);

    if ( hnuEvent.mMuMu < 60 || hnuEvent.mMuMu > 120 ) 
      hists.diLmassNo60to120Cut.fill(hnuEvent);
    if ( hnuEvent.mMuMu < 70 || hnuEvent.mMuMu > 100 )
    {
      if(isVerbose) std::cout << "Event outside the dimuon (70,100) veto region" << std::endl ;
      if(isVerbose) std::cout << "*** Event passes tight Z window veto, and mumujj mass is " << hnuEvent.mWR << " GeV" << std::endl ;

      double j1b = hnuEvent.j1.bDiscriminator(btagName);
      double j2b = hnuEvent.j2.bDiscriminator(btagName);
      int nBjets = 0 ;
      if (j1b >= minBtagDiscVal) nBjets++ ; 
      if (j2b >= minBtagDiscVal) nBjets++ ;      
      if ( nBjets > 0 ) if(isVerbose) std::cout << "*** Event passes tight Z window veto, and mumujj mass is " << hnuEvent.mWR << " GeV for one or more b-jets" << std::endl ;
      if ( nBjets > 1 ) if(isVerbose) std::cout << "*** Event passes tight Z window veto, and mumujj mass is " << hnuEvent.mWR << " GeV for two or more b-jets" << std::endl ;
      hists.diLmassNo70to100Cut.fill(hnuEvent);
    }
    if ( hnuEvent.mMuMu < 70 || hnuEvent.mMuMu > 110 ) {
      hists.diLmassNo70to110Cut.fill(hnuEvent);
    }
    if ( hnuEvent.mMuMu > 100 ) 
      hists.diLmassAbove100Cut.fill(hnuEvent);
    if ( hnuEvent.mMuMu > 120 ) { 
      if(isVerbose) std::cout << "Event has dimuon mass above 120 GeV" << std::endl ;
      if(isVerbose) std::cout << "*** Event away from Z, and mumujj mass is " << hnuEvent.mWR << " GeV" << std::endl ;
      hists.diLmassAbove120Cut.fill(hnuEvent);
    }
    if ( hnuEvent.mMuMu > 150 ) 
      hists.diLmassAbove150Cut.fill(hnuEvent);
    if ( hnuEvent.mMuMu > 200 ) 
      hists.diLmassAbove200Cut.fill(hnuEvent);

    // Change the final logic of the filter based on LQ meeting discussion:
    // Interest in seeing events that pass the dilepton mass requirement
    return true;
}

// ------------ method called once each job just before starting event loop  ------------

void HeavyNuFNAL::beginJob()
{
    evtCounter = 0 ; 
}

// ------------ method called once each job just after ending the event loop  ------------

void HeavyNuFNAL::endJob()
{
    trig_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuFNAL);


