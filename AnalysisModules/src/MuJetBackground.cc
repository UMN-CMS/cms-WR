// -*- C++ -*-
//
// Package:    MuJetBackground
// Class:      MuJetBackground
// 
/**\class MuJetBackground MuJetBackground.cc HeavyNu/AnalyzerModules/src/MuJetBackground.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: MuJetBackground.cc,v 1.13 2011/10/08 18:28:17 bdahmes Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// According to
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
// this must be included before Frameworkfwd.h
//
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNu_NNIF.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//////////////////////////////////////////////////////////////////
// generic maximum/minimum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;
}
template <class T> const T& min ( const T& a, const T& b ) {
  return (b<a)?b:a;
}

template <class T>
inline std::string int2str(T i) {
  std::ostringstream ss;
  ss << i;
  return ss.str();
}

//============================================================

class compare {
public:
  template <class T> bool operator() (const T& a, const T& b) { return a.pt() > b.pt() ; } 
};

//============================================================

inline void outputCandidate(const reco::CandidateBaseRef& rc) {
  std::cout<<"pt="<<rc->pt()<<" GeV, eta="<<rc->eta()<<", phi="<<rc->phi();
}

//============================================================

inline void dumpJetCorInfo(const pat::Jet& j) {
  std::cout << "Available JEC sets and levels:\n";
  const std::vector<std::string> jecsets = j.availableJECSets();
  for (size_t i=0; i<jecsets.size(); i++) {
    std::cout << jecsets[i]<<":";
    const std::vector<std::string> jeclevs = j.availableJECLevels(i);
    for (size_t j=0; j<jeclevs.size(); j++)
      std::cout << " " << jeclevs[i];
    std::cout << std::endl;
  }
  std::cout << "current JEC set   : " << j.currentJECSet()    << std::endl;
  std::cout << "current JEC level : " << j.currentJECLevel()  << std::endl;
  std::cout << "current JEC flavor: " << j.currentJECFlavor() << std::endl;
}

//============================================================


// Returns 0=neither, 1=loose or 2=tight, -1 if tight but not loose (!)
int qcdJetID(const pat::Jet& j)
{
  JetIDSelectionFunctor jetIDloose(JetIDSelectionFunctor::PURE09,JetIDSelectionFunctor::LOOSE);
  JetIDSelectionFunctor jetIDtight(JetIDSelectionFunctor::PURE09,JetIDSelectionFunctor::TIGHT);

  pat::strbitset ret = jetIDloose.getBitTemplate();
  ret.set(false);  bool loose = jetIDloose(j, ret);
  ret.set(false);  bool tight = jetIDtight(j, ret);
  return (tight ? (loose ? 2 : -1) : (loose ? 1 : 0));
}

inline void labelJetIDaxis(TAxis *ax)
{
  ax->SetBinLabel(1,"Neither");
  ax->SetBinLabel(2,"PURE09 Loose");
  ax->SetBinLabel(3,"PURE09 Tight");
}


static std::string btagName;

class MuJetBackground : public edm::EDFilter {
public:
  explicit MuJetBackground(const edm::ParameterSet&);
  ~MuJetBackground();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }
  
  virtual void beginJob          ();
  virtual bool filter            ( edm::Event&, const edm::EventSetup& );
  virtual void endJob            ();
  // virtual bool isVBTFloose       ( const pat::Muon& m );
  // virtual bool isVBTFtight       ( const pat::Muon& m );
  // virtual void fillBasicMuHistos ( const pat::Muon& m );
  // virtual void fillBasicJetHistos( const pat::Jet& j,
  // 				   int jetnum );

  void initializeHNE(HeavyNuEvent& hne, edm::Handle< std::vector<PileupSummaryInfo> >& pPU, 
		     edm::Handle<reco::VertexCollection>& pvHandle, bool isMC, bool isPF) ; 

  // virtual bool selectJetsStd     ( edm::Handle<pat::JetCollection>& pJets,
  // 				   HeavyNuEvent& hne );
  virtual bool selectJets        ( std::vector< std::pair<pat::Jet,float> >& jets,
				   HeavyNuEvent& hne );
  // virtual bool selectMuons       ( edm::Handle<pat::MuonCollection>& pMuons,
  // 				   HeavyNuEvent& hne );
  virtual void selectMuonsInJets ( std::vector<pat::Muon>& muons,
				   std::vector< std::pair<pat::Jet,float> >& jets,
				   HeavyNuEvent& hne );
  bool findQCDmuon               ( const std::vector<pat::Muon>& muons, 
				   edm::Handle<pat::MuonCollection>& pMuons,
				   HeavyNuEvent& hne );
  bool secondQualityMuon         (const pat::Muon& mu1, const pat::Muon& mu2);
  virtual bool findQCDjet        ( const std::vector< std::pair<pat::Jet,float> >& jets,
				   edm::Handle<pat::MuonCollection>& pMuons,
				   HeavyNuEvent& hne );
  virtual bool isDijetCandidate  ( HeavyNuEvent& hne,pat::MET& theMET ); 
  virtual TH1 *bookRunHisto      ( uint32_t runNumber );
  
  inline bool inZmassWindow( double mMuMu ) {
    return (mMuMu <= ZwinMaxGeV_) && (mMuMu >= ZwinMinGeV_);
  }

  edm::InputTag muonTag_;
  edm::InputTag jetTag_;
  edm::InputTag elecTag_;
  edm::InputTag photTag_;
  edm::InputTag hybridSClabel_ ; 
  edm::InputTag multiSClabel_ ; 

  double ZwinMinGeV_, ZwinMaxGeV_; // for trigger efficiency studies
  double applyMESfactor_;             // for Muon Energy Scale studies

  std::string currentFile_;
  bool dolog_;
  bool firstEvent_;
  int theMETtype ; 
  std::vector<double> rwLowPtbin, rwHighPtbin ; 
  std::vector<double> rwLoose, rwTight ; 

  bool calcSurvival_ ; 
  bool doClosure_, doQuadJet_ ; 

  JetCorrectionUncertainty *jecuObj_;

  bool   applyMuIDCorrections_ ;
    
  // HeavyNu_NNIF *nnif_;
  HeavyNuTrigger *trig_;
  HeavyNuID *muid_ ; 

  int    pileupEra_;
  bool isPFJets_; 
  edm::LumiReWeighting MCweightByVertex_;
  std::map<uint32_t,TH1 *> m_runHistos_;

  // ----------member data ---------------------------

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory *, const std::string&, const std::vector<hNuMassHypothesis>&) ;
    // fill all histos of the set with the mu-jet candidate
    void fill(pat::Muon& theMuon, pat::Jet& theJet, pat::MET& theMET, bool isMC, double trkIso) ;
    // void fill(reco::SuperCluster theSC) ; 
    // void fill(pat::MuonCollection muons, pat::JetCollection jets,bool isMC) ;
    // fill all histos of the set with the two muon candidates
    void fill(const HeavyNuEvent& hne, double w1=1.0, double w2=1.0) ;

    TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2;
    TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2 ;
    TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2 ;
    TH1 *ptMET, *phiMET ; 
    TH1 *ptMu1trkIso, *etaMu1trkIso, *phiMu1trkIso ; 

    TH1* muQual ; 
    TH1 *ptEB, *ptEE ; 

    // TH1 *dEtaMu, *dPhiMu, *dEtaJet, *dPhiJet ;
    // TH1 *dEtaPhiMu, *dEtaPhiJet ; 
    // TH1 *dRminMu1jet, *dRminMu2jet ; 
    // TH1 *hptrelMu1, *hptrelMu2 ; 
    // TH2 *ptrelVsdRminMu1jet, *ptrelVsdRminMu2jet ;
    // TH2 *jetID2d;

    // TH1 *dptMu1gen, *dptMu2gen ; 
    // TH1 *dRMu1gen, *dRMu2gen ; 
    // TH1 *qualMu1, *qualMu2 ; 

    // TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1caloIso, *mu1dB;
    // TH1 *mu2trackIso, *mu2hcalIso, *mu2ecalIso, *mu2caloIso, *mu2dB;

    // TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1caloRelIso;
    // TH1 *mu2trackRelIso, *mu2hcalRelIso, *mu2ecalRelIso, *mu2caloRelIso;

    TH1 *mWR, *mNuR1, *mNuR2, *mMuMu, *mMuMuZoom, *mJJ ; 
    TH2 *mNuR2D, *mNuR2D_raw, *jetPtvsNum ; 

    TH1 *mWR_raw, *mNuR1_raw, *mNuR2_raw, *mMuMu_raw, *mMuMuZoom_raw, *mJJ_raw ; 

    // TH1* btagJet1, *btagJet2;

    TH1* czeta_mumu;
    TH1* czeta_mumu_zoom;


    // Jeremy's crazy angles...
    // TH1* ctheta_mumu, *cthetaz_mumu;
    // TH1* ctheta_jj, *cthetaz_jj;
    // TH1* ctheta_mu1_jj, *cthetaz_mu1_jj;
    // TH1* ctheta_mu2_jj, *cthetaz_mu2_jj;

    TFileDirectory *mydir;
    // TFileDirectory *nndir;

    // HeavyNuTrigger::trigHistos_t trigHistos;
  };

  bool init_;

  // gf set of histo for all Z definitions in a stack
  struct HistStruct {
    TH1 *nelec, *nmu, *njet ;
    TH1 *muPt, *muEta, *muPhi ; 

    // Muon quality histos as a function of Pt
    TH2 *muNvalidHitsVsPt, *mudBvsPt, *muNormChi2vsPt, *muQualVsPt;
    TH2 *muNmatchesVsPt, *muNvalidMuonHitsVsPt, *muNvalidPixelHitsVsPt;
    TH2 *muTrckIsoVsPt, *muHcalIsoVsPt, *muEcalIsoVsPt, *muCaloIsoVsPt;

    TH1 *jetPt, *jetEta, *jetPhi, *jetID ; 
    TH2 *jetPtvsNum;

    TFileDirectory *rundir;

    HistPerDef NoCuts;
    HistPerDef dPhiCuts;
    HistPerDef dPhi10pctCuts;

    HistPerDef LLJJpTCuts;
    HistPerDef TrigMatches;
    HistPerDef VertexCuts;
    HistPerDef Mu1HighPtCut;
    HistPerDef diLmassCut;

    HistPerDef LJJJClosure; 
    HistPerDef L2JClosure;
  } hists;

  struct CutsStruct {
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
    double minimum_dijet_dPhi; 
    double minimum_dijet_pt; 
    double minimum_SCEt; 
    double minimum_extraJet_dR; 
  } cuts;
  
};

//======================================================================

const int muonQualityFlags = 4 ;
const std::string muonQuality[] = {
  "All","AllGlobalMuons","AllStandAloneMuons","AllTrackerMuons"
};

// bool
// MuJetBackground::isVBTFloose(const pat::Muon& m)
// {
//   // std::cout << "Global muon: " << m.muonID("AllGlobalMuons") ; 
//   // if ( m.muonID("AllGlobalMuons") ) std::cout << " " << m.numberOfValidHits() << " hits" ; 
//   // std::cout << std::endl ; 
//   return (m.muonID("AllGlobalMuons") &&
// 	  (m.numberOfValidHits() > 10));
// }

// bool
// MuJetBackground::isVBTFtight(const pat::Muon& m)
// {
//   if (!isVBTFloose(m)) return false;

//   reco::TrackRef gt = m.globalTrack();
//   if (gt.isNull()) {
//     std::cerr << "Mu global track reference is NULL" << std::endl;
//     return false;
//   }
//   return (m.muonID("AllTrackerMuons") &&
// 	  (m.dB() < 0.2) &&
// 	  (m.normChi2() < 10) &&
// 	  (m.numberOfMatches() > 1) &&
// 	  (gt->hitPattern().numberOfValidMuonHits()>0) &&
// 	  (gt->hitPattern().numberOfValidPixelHits()>0) );

// }                                                // MuJetBackground::isVBTFtight

void MuJetBackground::initializeHNE(HeavyNuEvent& hne, 
				    edm::Handle< std::vector<PileupSummaryInfo> >& pPU, 
				    edm::Handle<reco::VertexCollection>& pvHandle, 
				    bool isMC, bool isPF) {  
  hne.isMC   = isMC ; 
  hne.pfJets = isPF ; 

  if (hne.isMC) { 
    std::pair<float,double> pileup = hnu::pileupReweighting(pPU,MCweightByVertex_) ; 
    hne.n_pue     = int(pileup.first) ; 
    hne.eventWgt *= pileup.second ; 
  }
  hne.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle) ;
}


inline void labelMuonQualAxis(TAxis *ax)
{
  for (int i=0; i<muonQualityFlags; i++) {
    ax->SetBinLabel(i+1,muonQuality[i].c_str()) ;
    ax->SetBinLabel(i+1,muonQuality[i].c_str()) ;
  }
}

inline std::string nnhistoname(int mwr,int mnu) {
  return ("WR"+int2str<int>(mwr)+"nuRmu"+int2str<int>(mnu));
}

double qcdScaleFactor(double pt, 
		   std::vector<double> ptLow, 
		   std::vector<double> ptHigh, 
		   std::vector<double> corr) {
  double correction = 1.0 ; 
  for (unsigned int i=0; i<ptLow.size(); i++) 
    if ( (pt >= ptLow.at(i)) && (pt < ptHigh.at(i)) ) return corr.at(i) ; 
  return correction ; 
}

//======================================================================

void
MuJetBackground::HistPerDef::book(TFileDirectory *td,
				  const std::string& post,
				  const std::vector<hNuMassHypothesis>& v_masspts )
{
  // std::cout << "Beginning booking: " << post << std::endl ; 
  std::string t; // histogram title string;
  
  TH1::SetDefaultSumw2();

  mydir = td;

  // ----------  Muon histograms  ----------

  t="p_{T}(#mu_{1}) "+post;   ptMu1=td->make<TH1D>("ptMu1",t.c_str(),200,0.,1000.);
  t="p_{T}(#mu_{2}) "+post;   ptMu2=td->make<TH1D>("ptMu2",t.c_str(),200,0.,1000.);
  t="#eta(#mu_{1}) " +post;  etaMu1=td->make<TH1D>("etaMu1",t.c_str(),40,-2.5,2.5);
  t="#eta(#mu_{2}) " +post;  etaMu2=td->make<TH1D>("etaMu2",t.c_str(),40,-2.5,2.5);
  t="#phi(#mu_{1}) " +post;  phiMu1=td->make<TH1D>("phiMu1",t.c_str(),30,-3.14159,3.14159);
  t="#phi(#mu_{2}) " +post;  phiMu2=td->make<TH1D>("phiMu2",t.c_str(),30,-3.14159,3.14159);

  t="p_{T}(#mu_{1}) "+post;   ptMu1trkIso=td->make<TH1D>("ptMu1trkIso",t.c_str(),200,0.,1000.);
  t="#eta(#mu_{1}) " +post;  etaMu1trkIso=td->make<TH1D>("etaMu1trkIso",t.c_str(),40,-2.5,2.5);
  t="#phi(#mu_{1}) " +post;  phiMu1trkIso=td->make<TH1D>("phiMu1trkIso",t.c_str(),30,-3.14159,3.14159);

  t="mu quality " + post; muQual=td->make<TH1D>("muQual",t.c_str(),3,0.,3.);

  // ----------  Electron histograms  ----------
  t="p_{T}(SC_{EB}) "+post;   ptEB=td->make<TH1D>("ptEB",t.c_str(),200,0.,1000.);
  t="p_{T}(SC_{EE}) "+post;   ptEE=td->make<TH1D>("ptEE",t.c_str(),200,0.,1000.);

  // ----------  Jet histograms ----------

  t="p_{T}(j_{1}) "            +post;     ptJet1=td->make<TH1D>("ptJet1",  t.c_str(),50,0.,500.);
  t="p_{T}(j_{2}) "            +post;     ptJet2=td->make<TH1D>("ptJet2",  t.c_str(),50,0.,500.);
  t= "#eta(j_{1}) "            +post;    etaJet1=td->make<TH1D>("etaJet1", t.c_str(),40,-5,5);
  t= "#eta(j_{2}) "            +post;    etaJet2=td->make<TH1D>("etaJet2", t.c_str(),40,-5,5);
  t= "#phi(j_{1}) "            +post;    phiJet1=td->make<TH1D>("phiJet1", t.c_str(),30,-3.14159,3.14159);
  t= "#phi(j_{2}) "            +post;    phiJet2=td->make<TH1D>("phiJet2", t.c_str(),30,-3.14159,3.14159);

  t="Missing E_{T} "+post;        ptMET=td->make<TH1D>("ptMET",t.c_str(),50,0.,100.);
  t="#phi (Missing E_{T}) "+post; phiMET=td->make<TH1D>("phiMET",t.c_str(),30,-3.14159,3.14159);
  
  // ----------  Composite histograms  ----------
  t="M(W_{R}) "                    +post;        mWR=td->make<TH1D>("mWR",   t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{1} "       +post;      mNuR1=td->make<TH1D>("mNuR1", t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{2} "       +post;      mNuR2=td->make<TH1D>("mNuR2", t.c_str(),70,0,1400);
  t="M(W_{R}) "                    +post;    mWR_raw=td->make<TH1D>("mWR_raw",   t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{1} "       +post;  mNuR1_raw=td->make<TH1D>("mNuR1_raw", t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{2} "       +post;  mNuR2_raw=td->make<TH1D>("mNuR2_raw", t.c_str(),70,0,1400);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post;     mNuR2D=td->make<TH2D>("mNuR2D",t.c_str(),70,0,2800,70,0,1400);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post; mNuR2D_raw=td->make<TH2D>("mNuR2D_raw",t.c_str(),70,0,2800,70,0,1400);

  mNuR2D->Sumw2() ; 

  t="M(#mu #mu)"                   +post;         mMuMu=td->make<TH1D>("mMuMu",    t.c_str(),50,0,2000);
  t="M(#mu #mu)"                   +post;     mMuMuZoom=td->make<TH1D>("mMuMuZoom",t.c_str(),50,0,200);
  t="M(jj)"                        +post;           mJJ=td->make<TH1D>("mJJ",      t.c_str(),50,0,2000);
  t="M(#mu #mu)"                   +post;     mMuMu_raw=td->make<TH1D>("mMuMu_raw",    t.c_str(),50,0,2000);
  t="M(#mu #mu)"                   +post; mMuMuZoom_raw=td->make<TH1D>("mMuMuZoom_raw",t.c_str(),50,0,200);
  t="M(jj)"                        +post;       mJJ_raw=td->make<TH1D>("mJJ_raw",      t.c_str(),50,0,2000);

  t="cZeta(mumu)"                 +post; czeta_mumu=td->make<TH1D>("czMM",    t.c_str(),100,-1,1);
  t="cZeta(mumu) Zoom"       +post; czeta_mumu_zoom=td->make<TH1D>("czMMzoom",t.c_str(),100,-1,-0.9);

}// end of book()

//======================================================================

void MuJetBackground::HistPerDef::fill(pat::Muon& theMuon,
				       pat::Jet& theJet, 
				       pat::MET& theMET, 
				       bool isMC,
				       double trkIsoLimit)
{  
  // std::cout << "Begin fill" << std::endl ; 

  // Muons 
  ptMu1->Fill(theMuon.pt()) ; 
  etaMu1->Fill(theMuon.eta()) ; 
  phiMu1->Fill(theMuon.phi()) ; 

  // Jets 
  ptJet1->Fill(theJet.pt()) ; 
  etaJet1->Fill(theJet.eta()) ; 
  phiJet1->Fill(theJet.phi()) ; 

  // MET 
  ptMET->Fill(theMET.pt()) ; 
  phiMET->Fill(theMET.phi()) ; 

  if ( (theMuon.trackIso()/theMuon.pt()) < trkIsoLimit ) { // Relative track isolation
    ptMu1trkIso->Fill(theMuon.pt()) ; 
    etaMu1trkIso->Fill(theMuon.eta()) ; 
    phiMu1trkIso->Fill(theMuon.phi()) ; 
  }
  // std::cout << "end fill" << std::endl ; 
}// end of fill()

// void MuJetBackground::HistPerDef::fill(reco::SuperCluster theSC)
// {  
//   float scEt = theSC.energy()/cosh(theSC.eta()) ; 

//   if ( fabs(theSC.eta()) < 1.442 ) ptEB->Fill(scEt) ; 
//   else if ( fabs(theSC.eta()) > 1.56 && fabs(theSC.eta()) < 2.5 ) ptEE->Fill(scEt) ; 
// }// end of fill()


//======================================================================

void
MuJetBackground::HistPerDef::fill(const HeavyNuEvent& hne, double w1, double w2)
{
  double weight = w1 * w2 ; 

  // Muons 
  if (hne.nMuons > 0) { 
    ptMu1->Fill(hne.mu1.pt(),hne.eventWgt) ;     
    etaMu1->Fill(hne.mu1.eta(),hne.eventWgt) ; 
    phiMu1->Fill(hne.mu1.phi(),hne.eventWgt) ; 

    if (hne.nMuons > 1 ) { 
      ptMu2->Fill(hne.mu2.pt(),hne.eventWgt) ; 
      etaMu2->Fill(hne.mu2.eta(),hne.eventWgt) ; 
      phiMu2->Fill(hne.mu2.phi(),hne.eventWgt) ; 
    }
  }

  // Jets 
  if (hne.nJets > 0) {

    ptJet1->Fill(hne.j1.pt(),hne.eventWgt) ; 
    etaJet1->Fill(hne.j1.eta(),hne.eventWgt) ; 
    phiJet1->Fill(hne.j1.phi(),hne.eventWgt) ; 

    if (hne.nJets > 1) {
      ptJet2->Fill(hne.j2.pt(),hne.eventWgt) ; 
      etaJet2->Fill(hne.j2.eta(),hne.eventWgt) ; 
      phiJet2->Fill(hne.j2.phi(),hne.eventWgt) ; 

      mWR->Fill   ( hne.mWR,weight*hne.eventWgt   ) ; 
      mNuR1->Fill ( hne.mNuR1,w1*hne.eventWgt ) ; 
      mNuR2->Fill ( hne.mNuR2,w2*hne.eventWgt ) ; 
      mNuR2D->Fill( hne.mNuR1, hne.mNuR2,weight*hne.eventWgt );
      mJJ->Fill   ( hne.mJJ,weight*hne.eventWgt   );

      mWR_raw->Fill   ( hne.mWR,hne.eventWgt   ) ; 
      mNuR1_raw->Fill ( hne.mNuR1,hne.eventWgt ) ; 
      mNuR2_raw->Fill ( hne.mNuR2,hne.eventWgt ) ; 
      mNuR2D_raw->Fill( hne.mNuR1, hne.mNuR2,hne.eventWgt );
      mJJ_raw->Fill   ( hne.mJJ,hne.eventWgt   );

    }

  }

  mMuMu->Fill( hne.mMuMu,weight*hne.eventWgt );
  mMuMuZoom->Fill( hne.mMuMu,weight*hne.eventWgt );

  mMuMu_raw->Fill( hne.mMuMu,hne.eventWgt );
  mMuMuZoom_raw->Fill( hne.mMuMu,hne.eventWgt );

}// end of fill()

//======================================================================

//
// constants, enums and typedefs
//
const std::vector<hNuMassHypothesis> v_null;

//
// static data member definitions
//

//======================================================================

//
// constructors and destructor
//
MuJetBackground::MuJetBackground(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  dolog_=iConfig.getParameter<bool>("DoLog");

  muonTag_ = iConfig.getParameter< edm::InputTag >( "muonTag" );
  jetTag_  = iConfig.getParameter< edm::InputTag >( "jetTag"  );
  elecTag_ = iConfig.getParameter< edm::InputTag >( "electronTag" );
  photTag_ = iConfig.getParameter< edm::InputTag >( "photonTag" );
  
  hybridSClabel_ = iConfig.getParameter< edm::InputTag >( "hybridSCs" );
  multiSClabel_  = iConfig.getParameter< edm::InputTag >( "multi5x5SCs" );

  btagName=iConfig.getParameter<std::string>("BtagName");

  // nnif_ = new HeavyNu_NNIF(iConfig);

  trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet>("trigMatchPset"));
  muid_ = new HeavyNuID(iConfig.getParameter<edm::ParameterSet>("muIDPset"));

  edm::Service<TFileService> fs;
  hists.nelec    = fs->make<TH1D>("nelec", "N(e^{#pm})",10,-0.5,9.5);
  hists.nmu      = fs->make<TH1D>("nmu",   "N(#mu^{#pm})",10,-0.5,9.5);
  hists.njet     = fs->make<TH1D>("njet",  "N(Jet)",50,-0.5,49.5);
  hists.muPt     = fs->make<TH1D>("muPt",  "#mu p_{T} distribution",100,0,2000);
  hists.muEta    = fs->make<TH1D>("muEta", "#mu #eta distribution",50,-2.5,2.5);
  hists.muPhi    = fs->make<TH1D>("muPhi", "#mu #phi distribution",60,-3.14159,3.14159);
  hists.jetPt    = fs->make<TH1D>("jetPt", "jet p_{T} distribution",100,0,2000);
  hists.jetEta   = fs->make<TH1D>("jetEta","jet #eta distribution",50,-5,5);
  hists.jetPhi   = fs->make<TH1D>("jetPhi","jet #phi distribution",60,-3.14159,3.14159);
  hists.jetID    = fs->make<TH1I>("jetID", "Jet ID",3,0,3);
  hists.jetPtvsNum=fs->make<TH2D>("jetPtvsNum","Jet P_{T} vs. Jet # ",11,-0.5,10.5,200,0.,2000.);

  labelJetIDaxis(hists.jetID->GetXaxis());

  applyMESfactor_ = iConfig.getParameter<double>("applyMESfactor");

  if (applyMESfactor_==1.0) { // otherwise don't bother
    //  Muon quality variables vs muon p_T
    //
    hists.muNvalidHitsVsPt      = fs->make<TH2D>("muNvalidHitsVsPt",
						 "#mu # Valid hits vs p_{T}; p_{T}(#mu) (GeV); # Tracker+Pixel Hits",
						 150,0,3000,50,0,50);
    hists.mudBvsPt              = fs->make<TH2D>("mudBvsPt",
						 "#mu dXY vs p_{T}; p_{T}(#mu) (GeV); dXY(#mu)",
						 150,0,3000,40,0.,20.);
    hists.muNormChi2vsPt        = fs->make<TH2D>("muNormChi2vsPt",
						 "#mu Norm #chi^{2} vs p_{T}; p_{T}(#mu) (GeV); norm #chi^{2}",
						 150,0,3000,25,0,50);
    hists.muQualVsPt            = fs->make<TH2D>("muIDvsPt",
						 "Qual(#mu) vs p_{T}(#mu); p_{T}(#mu) (GeV)",
						 150,0,3000,4,0,4);
    hists.muNmatchesVsPt        = fs->make<TH2D>("muNmatchesVsPt",
						 "#mu # Matches vs p_{T}; p_{T}(#mu) (GeV); # Matches",
						 150,0,3000,10,0,10);
    hists.muNvalidMuonHitsVsPt  = fs->make<TH2D>("muNvalidMuonHitsVsPt",
						 "# Valid Muon hits vs p_{T}; p_{T}(#mu) (GeV); # Valid Muon Hits",
						 150,0,3000,50,0,50);
    hists.muNvalidPixelHitsVsPt = fs->make<TH2D>("muNvalidPixelHitsVsPt",
						 "# Valid Pixel hits vs p_{T}; p_{T}(#mu) (GeV); # Valid Pixel Hits",
						 150,0,3000,10,0,10);
    
    labelMuonQualAxis(hists.muQualVsPt->GetYaxis());
    
    hists.muTrckIsoVsPt=fs->make<TH2D>("muTrckIsoVsPt","trackIso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);trackIso (GeV)",150,0.,3000.,30,0.,300.);
    hists.muHcalIsoVsPt=fs->make<TH2D>("muHcalIsoVsPt","HCAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);HCAL Iso (GeV)",150,0.,3000.,30,0.,300.);
    hists.muEcalIsoVsPt=fs->make<TH2D>("muEcalIsoVsPt","ECAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);ECAL Iso (GeV)",150,0.,3000.,30,0.,300.);
    hists.muCaloIsoVsPt=fs->make<TH2D>("muCaloIsoVsPt","Calo Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);Calo Iso (GeV)",150,0.,3000.,30,0.,300.);
  }

  // Histos per cut:
  //
  // hists.noCuts.book                ( new TFileDirectory(fs->mkdir("NoCuts")),       
  //  				     "(no cuts)", v_null );
  // // hists.dPhiCuts.book              ( new TFileDirectory(fs->mkdir("dPhiCuts")), 
  // // 				     "(VBTF loose)", v_null );
  // // hists.dPhiTightCuts.book         ( new TFileDirectory(fs->mkdir("dPhiTightCuts")), 
  // // 				     "(VBTF tight)", v_null );

  calcSurvival_ = iConfig.getParameter<bool>("getSurvivalRate") ;
  doClosure_    = iConfig.getParameter<bool>("doClosureTest") ;
  doQuadJet_    = iConfig.getParameter<bool>("doQuadJetTest") ;

  if ( calcSurvival_ ) { 
    hists.NoCuts.book        ( new TFileDirectory(fs->mkdir("NoCuts")), 
                               "(VBTF tight, mu+jet)", v_null );
    hists.dPhiCuts.book      ( new TFileDirectory(fs->mkdir("dPhiCuts")), 
                               "(VBTF tight, NI Calo)", v_null );
    hists.dPhi10pctCuts.book ( new TFileDirectory(fs->mkdir("dPhi10pctCuts")), 
                               "(VBTF tight, NI 10% Calo)", v_null );
  } else { 
    if ( doQuadJet_ ) { 
      hists.LLJJpTCuts.book   ( new TFileDirectory(fs->mkdir("LLJJpTCuts")), 
				"(two muons, two jets)", v_null );
      hists.TrigMatches.book  ( new TFileDirectory(fs->mkdir("TrigMatches")), 
				"(at least one trigger matched muon)", v_null );
      hists.VertexCuts.book   ( new TFileDirectory(fs->mkdir("VertexCuts")), 
				"(all objects share common vtx)", v_null );
      hists.Mu1HighPtCut.book ( new TFileDirectory(fs->mkdir("Mu1HighPtCut")), 
				"(passes mu1 pT rqmt)", v_null );
      hists.diLmassCut.book   ( new TFileDirectory(fs->mkdir("diLmassCuts")), 
				"(passes mumu mass cut)", v_null );
    }
    if ( doClosure_ ) { 
      hists.LJJJClosure.book ( new TFileDirectory(fs->mkdir("LJJJClosure")), 
			       "(3 jets, 1 muon in jet)", v_null );
      hists.L2JClosure.book  ( new TFileDirectory(fs->mkdir("L2JClosure")), 
			       "(2 jets, tight muon)", v_null );
    }
  }

  hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

  init_=false;

  cuts.minimum_mu1_pt       = iConfig.getParameter<double>("minMu1pt");
  cuts.minimum_mu2_pt       = iConfig.getParameter<double>("minMu2pt");
  cuts.minimum_jet_pt       = iConfig.getParameter<double>("minJetPt");
  cuts.maximum_mu_abseta    = iConfig.getParameter<double>("maxMuAbsEta");
  cuts.maximum_jet_abseta   = iConfig.getParameter<double>("maxJetAbsEta");
  cuts.minimum_mumu_mass    = iConfig.getParameter<double>("minMuMuMass");
  cuts.minimum_mWR_mass     = iConfig.getParameter<double>("min4objMass");
  cuts.minimum_muon_jet_dR  = iConfig.getParameter<double>("minMuonJetdR");
  cuts.muon_trackiso_limit  = iConfig.getParameter<double>("muonTrackIsoLimitGeV");
  cuts.maxVertexZsep        = iConfig.getParameter<double>("dimuonMaxVertexZsepCM");
  cuts.maxJetVZsepCM        = iConfig.getParameter<double>("maxJetVZsepCM");
  cuts.minimum_dijet_dPhi   = iConfig.getParameter<double>("minimumMuJetdPhi");
  cuts.minimum_dijet_pt     = iConfig.getParameter<double>("minimumJetPtForDijets");
  cuts.minimum_extraJet_dR  = iConfig.getParameter<double>("minimumDeltaRforExtraJets");
  cuts.minimum_SCEt         = iConfig.getParameter<double>("minimumSuperClusterEt");

  theMETtype  = iConfig.getParameter<int>("METvariety") ; 
  rwLowPtbin  = iConfig.getParameter< std::vector<double> >("reweightPtLow") ; 
  rwHighPtbin = iConfig.getParameter< std::vector<double> >("reweightPtHigh") ; 
  // rwLoose     = iConfig.getParameter< std::vector<double> >("reweightLoose") ; 
  rwTight     = iConfig.getParameter< std::vector<double> >("reweightTight") ; 

  // Special check: Make sure all vectors are of the same size
  unsigned int vecsize = rwLowPtbin.size() ; 
  if ( ( doClosure_ || doQuadJet_ ) && 
       ( rwHighPtbin.size() != vecsize ||
 	 // rwLoose.size() != vecsize ||
 	 rwTight.size() != vecsize ) )
    throw cms::Exception( "Please ensure that all QCD reweighting vectors are equal size");

  ZwinMinGeV_ = iConfig.getParameter<double>("ZmassWinMinGeV");
  ZwinMaxGeV_ = iConfig.getParameter<double>("ZmassWinMaxGeV");

  pileupEra_ = iConfig.getParameter<int>("pileupEra");
  isPFJets_ = iConfig.getParameter<bool>("isPFJets") ; 
  MCweightByVertex_ = edm::LumiReWeighting(hnu::generate_flat10_mc(50),
					   hnu::get_standard_pileup_data(pileupEra_,50));

  // For the record...
  std::cout << "Configurable cut values applied:" << std::endl;
  std::cout << "muonTag          = " << muonTag_                 << std::endl;
  std::cout << "jetTag           = " << jetTag_                  << std::endl;
  std::cout << "electronTag      = " << elecTag_                 << std::endl;
  std::cout << "photonTag        = " << photTag_                 << std::endl;
  std::cout << "ZmassWinMinGeV   = " << ZwinMinGeV_              << " GeV" << std::endl;
  std::cout << "ZmassWinMaxGeV   = " << ZwinMaxGeV_              << " GeV" << std::endl;
  std::cout << "minMu1pt         = " << cuts.minimum_mu1_pt      << " GeV" << std::endl;
  std::cout << "minMu2pt         = " << cuts.minimum_mu2_pt      << " GeV" << std::endl;
  std::cout << "minJetPt         = " << cuts.minimum_jet_pt      << " GeV" << std::endl;
  std::cout << "maxMuAbsEta      = " << cuts.maximum_mu_abseta   << std::endl;
  std::cout << "maxJetAbsEta     = " << cuts.maximum_jet_abseta  << std::endl;
  std::cout << "minMuonJetdR     = " << cuts.minimum_muon_jet_dR << std::endl;
  std::cout << "muonTrackIso     = " << cuts.muon_trackiso_limit << " GeV" << std::endl;
  std::cout << "minMuMuMass      = " << cuts.minimum_mumu_mass   << " GeV" << std::endl;
  std::cout << "min4objMass      = " << cuts.minimum_mWR_mass    << " GeV" << std::endl;
  std::cout << "applyMESfactor   = " << applyMESfactor_          << std::endl;
  std::cout << "minimumMuJetdPhi = " << cuts.minimum_dijet_dPhi  << std::endl; 
  std::cout << "minimumQCDjetPt  = " << cuts.minimum_dijet_pt   << " GeV " << std::endl; 
  std::cout << "minExtraJetdR    = " << cuts.minimum_extraJet_dR << std::endl; 
  std::cout << "pileup era       = " << pileupEra_ << std::endl;

}
  
MuJetBackground::~MuJetBackground()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//======================================================================

//
// member functions
//
//======================================================================

//======================================================================

TH1 *
MuJetBackground::bookRunHisto(uint32_t runNumber)
{
  // std::cout << "Begin run book" << std::endl ; 

  std::string runstr = int2str<uint32_t>(runNumber);
  return hists.rundir->make <TH1I> (runstr.c_str(), runstr.c_str(),1,1,2);
}

//======================================================================

bool 
MuJetBackground::selectJets(std::vector< std::pair<pat::Jet,float> >& jets,
			    HeavyNuEvent& hne) {
  for (unsigned int i=0; i<jets.size(); i++) { 
    if ( hne.nJets == 2 ) break ; 
    pat::Jet iJ = jets.at(i).first ;

    // Must have jets of sufficient transverse momentum
    if ( iJ.pt() < cuts.minimum_jet_pt ) continue ; 

    // Jets must be separated from candidate muon(s)
    double dRm1j = deltaR(iJ.eta(), iJ.phi(), hne.mu1.eta(), hne.mu1.phi()) ; 
    double dRm2j = ( hne.nMuons > 1 ) ? 
      ( deltaR(iJ.eta(), iJ.phi(), hne.mu2.eta(), hne.mu2.phi()) ) : ( 100.0 ) ; 
    if (std::min(dRm1j,dRm2j) > cuts.minimum_muon_jet_dR) { 
      hne.nJets++ ; 
      if      ( hne.nJets == 1 ) hne.j1 = iJ ; 
      else if ( hne.nJets == 2 ) hne.j2 = iJ ; 
      else    std::cout << "WARNING: Expected empty jet position" << std::endl ; 
    }
  }
  return ( hne.nJets >= 2 ) ; 
}


bool
MuJetBackground::findQCDjet(const std::vector< std::pair<pat::Jet,float> >& jets,
			    edm::Handle<pat::MuonCollection>& pMuons,
			    HeavyNuEvent& hne) {
  pat::Jet iJqcd ; 
  bool foundJet = false ; unsigned int jetLoc = jets.size() ;  
  float maxdPhi_jet_mu = cuts.minimum_dijet_dPhi ; 
  for (unsigned int i=0; i<jets.size(); i++) {
    pat::Jet iJ = jets.at(i).first ; 
    float dPhi = fabs( deltaPhi(iJ.phi(),hne.mu1.phi()) ) ; 
    if ( dPhi > maxdPhi_jet_mu ) { 
      foundJet = true ; 
      maxdPhi_jet_mu = dPhi ; 
      iJqcd = iJ ; 
      jetLoc = i ; 
    } 
  }
  if ( !foundJet ) return false ; 

  // Reject the event if extra high-energy jets are present outside mu/jet axis
  for (unsigned int i=0; i<jets.size(); i++) {
    if ( i == jetLoc ) continue ; 
    pat::Jet iJ = jets.at(i).first ; 
    if ( iJ.pt() > cuts.minimum_jet_pt ) {  
      float dRm = deltaR(iJ.eta(),iJ.phi(),hne.mu1.eta(),hne.mu1.phi()) ; 
      float dRj = deltaR(iJ.eta(),iJ.phi(),iJqcd.eta(),iJqcd.phi()) ; 
      if ( std::min(dRm,dRj) > cuts.minimum_extraJet_dR ) return false ; 
    }
  }

  // Reject event if muon found with most of jet energy 
  for (unsigned int i=0; i<pMuons->size(); i++) {
    pat::MuonRef iM=pat::MuonRef( pMuons,i ) ;
    float dR  = deltaR(iM->eta(),iM->phi(),hne.mu1.eta(),hne.mu1.phi()) ; 
    if ( dR < 0.001 ) continue ; // Matched to original muon
    dR = deltaR(iM->eta(),iM->phi(),iJqcd.eta(),iJqcd.phi()) ; 
    if ( dR < 0.5 ) { // muon is in the jet 
      double pTratio = iM->pt() / iJqcd.pt() ; 
      if ( pTratio > 0.75 ) return false ; 
    }
  }

  hne.j1 = iJqcd ; 
  return true ; 
} //MuJetBackground::findQCDjet

//======================================================================

void
MuJetBackground::selectMuonsInJets(std::vector<pat::Muon>& muons,
				   std::vector< std::pair<pat::Jet,float> >& jets,
				   HeavyNuEvent& hne) {
  double mu1wgt = 1.0 ; 
  double mu2wgt = 1.0 ; 
  for (unsigned int i=0; i<muons.size(); i++) {
    pat::Muon iM = muons.at(i) ; 
    bool muInJet = false ; 
    for (unsigned int j=0; j<jets.size(); j++) {
      if ( muInJet ) break ; 
      pat::Jet iJ = jets.at(j).first;
      double dR = deltaR(iJ.eta(),iJ.phi(),iM.eta(),iM.phi()) ;
      if ( dR < cuts.minimum_muon_jet_dR ) muInJet = true ; 
    }

    if ( muInJet ) { 
      hne.nMuons++ ; 
      if ( hne.nMuons == 1 ) { 
	hne.mu1 = iM ; 
	if ( hne.isMC && applyMuIDCorrections_ ) mu1wgt = muid_->weightForMC( iM.pt(),0 ) ; 
      } 
      if ( hne.nMuons == 2 ) { 
	hne.mu2 = iM ; 
	if ( hne.isMC && applyMuIDCorrections_ ) mu2wgt = muid_->weightForMC( iM.pt(),0 ) ; 
      } 
    }
  }

  hne.eventWgt *= mu1wgt * mu2wgt ;
}

bool MuJetBackground::secondQualityMuon(const pat::Muon& mu1, const pat::Muon& mu2) { 

  if ( !hnu::isVBTFloose(mu2) ) return false ; 

  float relIso = ( mu2.trackIso() + mu2.hcalIso() + mu2.ecalIso() ) / mu2.pt() ;
  if ( relIso < 0.15 ) return true ;

  reco::Particle::LorentzVector vMuMu = mu1.p4() + mu2.p4() ;
  if ( vMuMu.M() > 70. && vMuMu.M() < 110. ) return true ; 

  return false ; 
}

bool
MuJetBackground::findQCDmuon(const std::vector<pat::Muon>& muons, 
			     edm::Handle<pat::MuonCollection>& pMuons,
			     HeavyNuEvent& hne) {

  pat::Muon qcdCand = muons.at(0) ; 
  for (unsigned int i=1; i<muons.size(); i++) { // Skip first candidate
    if ( secondQualityMuon(qcdCand,muons.at(i)) ) return false ;  
  }
  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    pat::MuonRef iM=pat::MuonRef(pMuons,iMuon);
    if ( deltaR(qcdCand.eta(),qcdCand.phi(),iM->eta(),iM->phi()) < 0.001 ) continue ; // same muon 
    if ( secondQualityMuon(qcdCand,(*iM)) ) return false ;  
  }
  hne.mu1 = qcdCand ; 
  return true ; 
} // MuJetBackground::findQCDmuon

bool MuJetBackground::isDijetCandidate(HeavyNuEvent& hne,pat::MET& theMET) {

  if ( theMET.et() < 20. ) return false ; // Absolute MET requirement
  if ( hnu::jetID(hne.j1) < 1 ) return false ; // Require at least PURE09 Loose
  return true ; 

} // MuJetBackground::isDijetCandidate

//======================================================================

// ------------ method called to for each event  ------------
bool
MuJetBackground::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  HeavyNuEvent hnuEvent;

  bool keepThisEvent = false ; 

  applyMuIDCorrections_ = !iEvent.isRealData();

  // std::cout << "Test 0" << std::endl ; 

  if (iEvent.isRealData())
  {
    if( (applyMESfactor_ != 1.0) ) 
      throw cms::Exception( "Energy scale studies not allowed on data currently");

    uint32_t runn = iEvent.id().run();
    std::map<uint32_t,TH1 *>::const_iterator it = m_runHistos_.find(runn);
    TH1 *runh;
    if (it == m_runHistos_.end()) {
      runh = bookRunHisto(runn);
      m_runHistos_[runn] = runh;
    } else
      runh = it->second;
    runh->Fill(1);
  }
  // std::cout << "Test 0" << std::endl ; 

  edm::Handle<pat::MuonCollection> pMuons ; 
  iEvent.getByLabel(muonTag_,pMuons) ; 

  edm::Handle<pat::ElectronCollection> pElecs ;
  iEvent.getByLabel(elecTag_, pElecs) ;

  edm::Handle<pat::PhotonCollection> pGammas ;
  iEvent.getByLabel(photTag_, pGammas) ;

  edm::Handle<pat::JetCollection> pJets ;
  iEvent.getByLabel(jetTag_, pJets) ;

  edm::Handle<reco::SuperClusterCollection> hybridClusters ; 
  iEvent.getByLabel(hybridSClabel_, hybridClusters) ; 

  edm::Handle<reco::SuperClusterCollection> multi5x5Clusters ; 
  iEvent.getByLabel(multiSClabel_, multi5x5Clusters) ; 

  edm::Handle<reco::MuonCollection> tevMuons;
  iEvent.getByLabel("refitMuons", tevMuons); 

  edm::Handle<std::vector<PileupSummaryInfo> > pPU;
  iEvent.getByLabel("addPileupInfo", pPU);    

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);

  // if (hnuEvent.isMC) { 
  //   std::pair<int,double> pileup = hnu::pileupReweighting(MCweightByVertex_) ; 
  //   hnuEvent.n_pue     = pileup.first ; 
  //   hnuEvent.eventWgt *= pileup.second ; 
  // }
  // hnuEvent.n_primary_vertex = hnu::numberOfPrimaryVertices() ;

  if ( !pElecs.isValid() || 
       !pGammas.isValid() || 
       !pMuons.isValid() || 
       !pJets.isValid() ) {
    std::cout << "Exiting as valid PAT objects not found" << std::endl ;
    return false; 
  }

  edm::Handle<pat::METCollection> patMetCollection ; 
  iEvent.getByLabel("patMETs", patMetCollection) ; 

  edm::Handle<reco::PFMETCollection> pfMetCollection ; 
  iEvent.getByLabel("pfMet", pfMetCollection) ; 

  if ( ( theMETtype == 1 && !patMetCollection.isValid() ) ||
       ( theMETtype == 2 && !pfMetCollection.isValid() ) ) { 
    std::cout << "Exiting as valid MET collection not found" << std::endl ; 
    return false ; 
  }
  std::auto_ptr<pat::METCollection> pMET(new pat::METCollection) ; 
  if ( theMETtype == 1 ) pMET->push_back(patMetCollection.product()->at(0)) ;
  else if ( theMETtype == 2 ) {
    reco::MET thePFMET(pfMetCollection.product()->at(0).p4(),pfMetCollection.product()->at(0).vertex()) ; 
    pMET->push_back(pat::MET(thePFMET)) ;
  }  

  if (firstEvent_) {
    // handle the jet corrector parameters collection,
    // get the jet corrector parameters collection from the global tag
    //
    
    // get the uncertainty parameters from the collection,
    // instantiate the jec uncertainty object
    //
    if ( !iEvent.isRealData() ) {
        int pileupYear = pileupEra_ / 10 ;
        int idYear     = muid_->idEra() ;
        
        bool allErasMatch = ( pileupYear == idYear ) ;
        if ( !allErasMatch ) {
            std::cout << "WARNING: You do not appear to have consistent corrections applied!" << std::endl ;
            std::cout << "         pileup year is " << pileupEra_ << ", year for mu ID is " << idYear
		      << std::endl ; 
        } else {
            std::cout << "Looking at corrections, I assume you are running with the " << pileupYear << " year settings" << std::endl ; 
        }
        std::cout << "==================================" << std::endl ; 
    }  
    firstEvent_ = false;
  }

  hists.nelec->Fill(pElecs->size()) ;
  hists.nmu  ->Fill(pMuons->size()) ;
  hists.njet ->Fill(pJets->size()) ;

  std::vector<pat::Muon> muCands = 
    hnu::getMuonList(pMuons,tevMuons,cuts.minimum_mu2_pt,cuts.maximum_mu_abseta,applyMESfactor_) ; 
  std::vector< std::pair<pat::Jet,float> > jetCands = 
    hnu::getJetList(pJets,jecuObj_,cuts.minimum_dijet_pt,cuts.maximum_jet_abseta,0) ; 
  if ( muCands.size() < 1 || jetCands.size() < 1 ) return false ;
  
  if ( calcSurvival_ ) { 
    HeavyNuEvent hnuDijet(HeavyNuEvent::QCD) ;  
    initializeHNE(hnuDijet,pPU,pvHandle,!iEvent.isRealData(),isPFJets_) ; 
    hists.NoCuts.fill( muCands.at(0),jetCands.at(0).first,pMET->at(0),!iEvent.isRealData(),
                       cuts.muon_trackiso_limit ) ; 

    if ( findQCDmuon( muCands, pMuons, hnuDijet ) ) { 
      if ( findQCDjet( jetCands, pMuons, hnuDijet ) ) { 
	if ( isDijetCandidate( hnuDijet,pMET->at(0) ) ) {	  
	  if ( (hnuDijet.mu1.ecalIso()+hnuDijet.mu1.hcalIso()) > 10. ) {
	    hists.dPhiCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
				 cuts.muon_trackiso_limit ) ; 
	  }
	  if ( ((hnuDijet.mu1.ecalIso()+hnuDijet.mu1.hcalIso())/hnuDijet.mu1.pt()) > 0.10 ) {
	    hists.dPhi10pctCuts.fill( hnuDijet.mu1,hnuDijet.j1,pMET->at(0),!iEvent.isRealData(),
                                      cuts.muon_trackiso_limit ) ; 
	  }
	}
      }
    }
  } else { 
    if ( doQuadJet_ ) { 
      if ( muCands.size() >= 2 && jetCands.size() >= 4 ) { 
	HeavyNuEvent hnuQuadjet(HeavyNuEvent::QCD) ;  
	initializeHNE(hnuQuadjet,pPU,pvHandle,!iEvent.isRealData(),isPFJets_) ; 
	
	selectMuonsInJets( muCands,jetCands,hnuQuadjet ) ; 
	if ( hnuQuadjet.nMuons >= 2 ) {  
	  if ( selectJets( jetCands,hnuQuadjet ) ) { 
	    keepThisEvent = true ; 
	    hnuQuadjet.regularize();
	    hnuQuadjet.calculate() ; 

	    double mu1scaleFactor, mu2scaleFactor ; 
	    mu1scaleFactor = qcdScaleFactor(hnuQuadjet.mu1.pt(),
					    rwLowPtbin,rwHighPtbin,rwTight) ; 
	    mu2scaleFactor = qcdScaleFactor(hnuQuadjet.mu2.pt(),
					    rwLowPtbin,rwHighPtbin,rwTight) ; 
	  	  
	    hists.LLJJpTCuts.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 
	    
	    bool mu1trig = false ; bool mu2trig = false ;
	    if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
	      mu1trig = trig_->isTriggerMatched( hnuQuadjet.mu1, iEvent) ; 
	      mu2trig = trig_->isTriggerMatched( hnuQuadjet.mu2, iEvent) ; 
	    } else if ( !iEvent.isRealData() ) {
	      mu1trig = trig_->simulateForMC( hnuQuadjet.mu1.pt(),hnuQuadjet.mu1.eta(),0 );
	      mu2trig = trig_->simulateForMC( hnuQuadjet.mu2.pt(),hnuQuadjet.mu2.eta(),0 );
	    }

	    if ( mu1trig || mu2trig ) { 
	      hists.TrigMatches.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 
		  
	      //--- Impose vertex requirement here ---//
	      float deltaVzJ1J2 = fabs(hnuQuadjet.tjV1-hnuQuadjet.tjV2);
	      float deltaVzJ1M1 = fabs(hnuQuadjet.tjV1-hnuQuadjet.mu1.vertex().Z());
	      float deltaVzJ2M2 = fabs(hnuQuadjet.tjV2-hnuQuadjet.mu2.vertex().Z());
	      float deltaVzJ1M2 = fabs(hnuQuadjet.tjV1-hnuQuadjet.mu2.vertex().Z());
	      float deltaVzJ2M1 = fabs(hnuQuadjet.tjV2-hnuQuadjet.mu1.vertex().Z());
	      float deltaVzM1M2 = fabs(hnuQuadjet.mu1.vertex().Z()-hnuQuadjet.mu2.vertex().Z());
	      if ((deltaVzJ1J2 < cuts.maxJetVZsepCM) && (deltaVzJ1M1 < cuts.maxJetVZsepCM) &&
		  (deltaVzJ2M2 < cuts.maxJetVZsepCM) && (deltaVzJ1M2 < cuts.maxJetVZsepCM) &&
		  (deltaVzJ2M1 < cuts.maxJetVZsepCM) && (deltaVzM1M2 < cuts.maxVertexZsep) ) {

		hists.VertexCuts.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 

		if ( hnuQuadjet.mu1.pt() > cuts.minimum_mu1_pt ) {
		  hists.Mu1HighPtCut.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ;
		  if ( hnuQuadjet.mMuMu > cuts.minimum_mumu_mass ) { // dimuon mass rqmt
		    hists.diLmassCut.fill( hnuQuadjet,mu1scaleFactor,mu2scaleFactor ) ; 
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // 
    // Closure test
    // 
    if ( doClosure_ ) { 
      if ( pMuons->size() >= 1 && pJets->size() >= 3
	   && pMET->at(0).pt() < 20. ) { 

	HeavyNuEvent hnuClosure3jet1muon(HeavyNuEvent::CLO) ;  
	initializeHNE(hnuClosure3jet1muon,pPU,pvHandle,!iEvent.isRealData(),isPFJets_) ; 
      
	selectMuonsInJets( muCands,jetCands,hnuClosure3jet1muon ) ; 
	if ( hnuClosure3jet1muon.nMuons == 1 ) {  
	  if ( selectJets( jetCands,hnuClosure3jet1muon ) ) { 
	    // hnuClosure3jet1muon.mu2 = hnuClosure3jet1muon.mu1 ; // HACK

	    hnuClosure3jet1muon.regularize(); 
	    hnuClosure3jet1muon.scaleMuE(applyMESfactor_);
	    hnuClosure3jet1muon.calculate() ; 

	    bool mu1trig = false ; 
	    if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
	      mu1trig = trig_->isTriggerMatched( hnuClosure3jet1muon.mu1, iEvent) ; 
	    } else if (!iEvent.isRealData()) {
	      mu1trig = trig_->simulateForMC( hnuClosure3jet1muon.mu1.pt(),
					      hnuClosure3jet1muon.mu1.eta(),0 );
	    }

	    if ( mu1trig ) { 
	      double mu1scaleFactor = qcdScaleFactor(hnuClosure3jet1muon.mu1.pt(),
						     rwLowPtbin,rwHighPtbin,rwTight) ; 

	      //--- Impose vertex requirement here ---//
	      float deltaVzJ1J2 = fabs(hnuClosure3jet1muon.tjV1-hnuClosure3jet1muon.tjV2);
	      float deltaVzJ1M1 = fabs(hnuClosure3jet1muon.tjV1-hnuClosure3jet1muon.mu1.vertex().Z());
	      float deltaVzJ2M1 = fabs(hnuClosure3jet1muon.tjV2-hnuClosure3jet1muon.mu1.vertex().Z());
	      if ( (deltaVzJ1J2 < cuts.maxJetVZsepCM) && (deltaVzJ1M1 < cuts.maxJetVZsepCM) &&
		   (deltaVzJ2M1 < cuts.maxJetVZsepCM) ) { 
		hists.LJJJClosure.fill( hnuClosure3jet1muon, mu1scaleFactor ) ; 
	      }
	    }
	  }
	}
      }
    
      if ( muCands.size() >= 1 && jetCands.size() >= 2 
	   && pMET->at(0).pt() < 20. ) { // MET cut to remove W+2 jets

	HeavyNuEvent hnuClosure2jet1muon(HeavyNuEvent::CLO) ;
	initializeHNE(hnuClosure2jet1muon,pPU,pvHandle,!iEvent.isRealData(),isPFJets_);
      
	hnuClosure2jet1muon.j1 = jetCands.at(0).first ; 
	hnuClosure2jet1muon.j2 = jetCands.at(1).first ; 

	for (unsigned int i=0; i<muCands.size(); i++) { 
	  if ( hnuEvent.nMuons == 1 ) break ; 
	  pat::Muon iM = muCands.at(i) ; 
	  if ( hnu::muIsolation(iM,applyMESfactor_) < cuts.muon_trackiso_limit ) {
	    double dRj1 = deltaR(iM.eta(), iM.phi(), hnuClosure2jet1muon.j1.eta(), hnuClosure2jet1muon.j1.phi()) ; 
	    double dRj2 = deltaR(iM.eta(), iM.phi(), hnuClosure2jet1muon.j2.eta(), hnuClosure2jet1muon.j2.phi()) ; 
	    if (dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR) { 
	      hnuClosure2jet1muon.nMuons++ ; 
	      if   ( hnuClosure2jet1muon.nMuons == 1 ) hnuClosure2jet1muon.mu1 = iM ; 
	      else std::cout << "WARNING: Expected empty muon position" << std::endl ; 
	    }
	  }
	}
	if ( hnuClosure2jet1muon.nMuons > 0 ) { 
	  hnuClosure2jet1muon.eventWgt *= muid_->weightForMC((hnuClosure2jet1muon.mu1.pt()*applyMESfactor_),0) ;

	  // hnuClosure2jet1muon.mu2 = hnuClosure2jet1muon.mu1 ; // HACK
	  hnuClosure2jet1muon.regularize(); 
	  hnuClosure2jet1muon.scaleMuE(applyMESfactor_);
	  hnuClosure2jet1muon.calculate() ; 

	  bool mu1trig = false ; 
	  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
	    mu1trig = trig_->isTriggerMatched( hnuClosure2jet1muon.mu1, iEvent) ; 
	  } else if (!iEvent.isRealData()) {
	    mu1trig = trig_->simulateForMC( hnuClosure2jet1muon.mu1.pt(),hnuClosure2jet1muon.mu1.eta(),0 ) ; 
	  }

	  if ( mu1trig ) { 
	    //--- Impose vertex requirement here ---//
	    float deltaVzJ1J2 = fabs(hnuClosure2jet1muon.tjV1-hnuClosure2jet1muon.tjV2);
	    float deltaVzJ1M1 = fabs(hnuClosure2jet1muon.tjV1-hnuClosure2jet1muon.mu1.vertex().Z());
	    float deltaVzJ2M1 = fabs(hnuClosure2jet1muon.tjV2-hnuClosure2jet1muon.mu1.vertex().Z());
	    if ( (deltaVzJ1J2 < cuts.maxJetVZsepCM) && (deltaVzJ1M1 < cuts.maxJetVZsepCM) &&
		 (deltaVzJ2M1 < cuts.maxJetVZsepCM) ) { 
	      hists.L2JClosure.fill( hnuClosure2jet1muon ) ; 
	    }
	  }
	}
      }
    }
  }
  
  // Keep 4-jet events only...drop the rest
  if ( !calcSurvival_ && !doClosure_ ) return keepThisEvent ; 
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuJetBackground::beginJob() {
  // std::cout << "Beginning the job" << std::endl ; 
  // nnif_->beginJob();
  firstEvent_ = true;
  // std::cout << "leaving beginJob" << std::endl ; 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuJetBackground::endJob() {
  // nnif_->endJob();
  trig_->endJob();
  muid_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuJetBackground);
