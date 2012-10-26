#ifndef HEAVYNU_H
#define HEAVYNU_H

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TRandom.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTree.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuMuHist.h"

namespace hnu
{
   struct EffInfo
   {
      double ZwinMinGeV, ZwinMaxGeV;
      bool oneTP;
      TRandom *tpRandom;
   };
   
   struct CutsStruct
   {
      double minimum_mu1_pt;
      double minimum_mu2_pt;
      double minimum_jet_pt;
      double maximum_elec_abseta;
      double maximum_mu_abseta;
      double maximum_jet_abseta;
      double minimum_mumu_mass;
      double minimum_mWR_mass;
      double minimum_muon_jet_dR;
      double muon_trackiso_limit;
      double maxVertexZsep;
      double maxJetVZsepCM;
   };
}

class HeavyNu : public edm::EDFilter
{
public:
    explicit HeavyNu(const edm::ParameterSet&);
    ~HeavyNu();

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
      HeavyNuHistSet* LLptCuts;
      HeavyNuHistSet* JJptCuts;
      HeavyNuHistSet* TrigMatches;
      HeavyNuHistSet* LLJJptCuts;
      //HeavyNuHistSet* VertexCuts;
      HeavyNuHistSet* ZRegionCut;
      HeavyNuHistSet* Mu1HighPtCut;
      HeavyNuHistSet* Mu1HighPtCut_1bjet;
      HeavyNuHistSet* Mu1HighPtCut_2bjet;
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
      HeavyNuHistSet* HeepTagGsfProbeInZwin0jets;
      HeavyNuHistSet* HeepTagGsfProbeInZwin1jet;
      HeavyNuHistSet* HeepTagGsfProbeInZwin2jets;
      HeavyNuHistSet* HeepTagGsfProbePassesInZwin0jets;
      HeavyNuHistSet* HeepTagGsfProbePassesInZwin1jet;
      HeavyNuHistSet* HeepTagGsfProbePassesInZwin2jets;
      // HeavyNuHistSet* Mu1Mu2TrigMatchesInZwin;
   } hists;

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

    virtual void studyJetVertex(edm::Handle<pat::JetCollection>& pJets,
                                edm::Handle<reco::JPTJetCollection>& jptJets,
                                edm::Handle<pat::MuonCollection>& pMuons, int npue);

    inline bool inZmassWindow(double mMuMu)
    {
        return(mMuMu <= effinfo_.ZwinMaxGeV) && (mMuMu >= effinfo_.ZwinMinGeV);
    }

    void fill(pat::MuonCollection muons,
              pat::ElectronCollection electrons,
              pat::JetCollection jets,
              const HeavyNuEvent& hnuEvent,
              const bool goodJets,
              const bool goodLeps,
              HeavyNuHistSet *hnmh);

    bool isWrDaughter(const reco::Candidate* mother);
    
    void studyJetMatching(HeavyNuEvent& hnuEvent, edm::Handle<std::vector<reco::GenJet> > genjets);
    
    void selectMuons(std::vector<pat::Muon>& muCands, HeavyNuEvent& hnuEvent, int nDirtyCands);
    void selectElectrons(std::vector< std::pair<pat::Electron, float> >& eCands, edm::Handle<pat::ElectronCollection>& pElecs, 
                         edm::Handle<reco::VertexCollection> pvHandle, HeavyNuEvent& hnuEvent, int nDirtyCands);
    void alternativeSelection(HeavyNuEvent& hnuEvent);

    edm::InputTag muonTag_;
    edm::InputTag trackTag_;
    edm::InputTag jetTag_;
    edm::InputTag metTag_;
    edm::InputTag elecTag_;

    edm::InputTag elecRhoTag_ ; 
    double elecRho_ ; 

    int evtCounter;

    int jecVal_; // Jet correction eras: 0 = 2010(A+B), 1 = 2010A, 2 = 2010B, 3 = 2011A
    int applyJECUsign_; // for Jet Energy Correction Uncertainty studies
    int applyJERsign_; // for Jet Energy Resolution and Resolution Uncertainty studies
    double applyMESfactor_; // for Muon Energy Scale studies
    bool merUncertainty_ ; 
    bool correctEscale_ ; 
    int applyTrigEffsign_; // for Trigger Efficiency studies
    bool highestPtTriggerOnly_;
    bool applyMuIDCorrections_;
    int applyMuIDEffsign_; // for Loose+Iso efficiency correction
    bool studyMuonSelectionEff_; // for Muon Loose ID Efficiency studies
    bool studyScaleFactorEvolution_; // for Top, Z+jets scale factors (by Mu1 pT) studies
    int tpSeed_ ; 
    
    hnu::EffInfo effinfo_;

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
    bool addSlopeTree_;

    int nDirtyCands_ ; 

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

    HeavyNuTree *hnuTree_;


    // ----------member data ---------------------------
    bool init_;

    hnu::CutsStruct cuts;

};

#endif