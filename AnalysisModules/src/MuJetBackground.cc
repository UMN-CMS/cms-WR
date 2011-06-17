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
// $Id: MuJetBackground.cc,v 1.5 2011/06/09 19:16:32 bdahmes Exp $
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

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNu_NNIF.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"


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
  virtual bool isVBTFloose       ( const pat::Muon& m );
  virtual bool isVBTFtight       ( const pat::Muon& m );
  virtual void fillBasicMuHistos ( const pat::Muon& m );
  virtual void fillBasicJetHistos( const pat::Jet& j,
				   int jetnum );
  virtual bool selectJetsStd     ( edm::Handle<pat::JetCollection>& pJets,
				   HeavyNuEvent& hne );
  virtual bool selectJets        ( edm::Handle<pat::JetCollection>& pJets,
				   HeavyNuEvent& hne );
  virtual bool selectMuons       ( edm::Handle<pat::MuonCollection>& pMuons,
				   HeavyNuEvent& hne );
  virtual bool selectMuonsInJets ( edm::Handle<pat::MuonCollection>& pMuons,
				   edm::Handle<pat::JetCollection>& pJets,
				   HeavyNuEvent& hne );
  virtual bool selectQCDmuon     ( edm::Handle<pat::MuonCollection>& pMuons,
				   HeavyNuEvent& hne );
  virtual bool selectQCDjet      ( edm::Handle<pat::JetCollection>& pJets,
				   edm::Handle<pat::MuonCollection>& pMuons,
				   HeavyNuEvent& hne );
  virtual bool isDijetCandidate  ( HeavyNuEvent& hne,pat::MET theMET ); 
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
  
  HeavyNu_NNIF *nnif_;
  HeavyNuTrigger *trig_;

  std::map<uint32_t,TH1 *> m_runHistos_;

  // ----------member data ---------------------------

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory *, const std::string&, const std::vector<hNuMassHypothesis>&) ;
    // fill all histos of the set with the mu-jet candidate
    void fill(pat::Muon theMuon, pat::Jet theJet, pat::MET theMET, bool isMC, double trkIso) ;
    void fill(reco::SuperCluster theSC) ; 
    // void fill(pat::MuonCollection muons, pat::JetCollection jets,bool isMC) ;
    // fill all histos of the set with the two muon candidates
    void fill(const HeavyNuEvent& hne, double w1=1.0, double w2=1.0, int nTight=-1) ;

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

    HistPerDef dPhiCutsNICalo10; 
    HistPerDef dPhiTightCutsNICalo10; 

    HistPerDef dPhiHeepSC ; 
    HistPerDef dPhiHeepEle ; 

    HistPerDef FourJetLooseCuts;
    HistPerDef FourJetTightTrigCuts;
    HistPerDef FourJetVertexCuts;
    HistPerDef FourJetMu1PtCuts;
    HistPerDef FourJetZMassCuts;

    HistPerDef LJJJLooseCuts; 
    HistPerDef LJJJTightCuts; 
    HistPerDef LJJLooseCuts; 
    HistPerDef LJJTightCuts; 
    HistPerDef LJJZrejLooseCuts; 
    HistPerDef LJJZrejTightCuts; 

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
    double minimum_dijet_dPhi; 
    double minimum_QCDjet_pt; 
    double minimum_SCEt; 
    double minimum_extraJet_dR; 
  } cuts;
  
};

//======================================================================

const int muonQualityFlags = 4 ;
const std::string muonQuality[] = {
  "All","AllGlobalMuons","AllStandAloneMuons","AllTrackerMuons"
};

bool
MuJetBackground::isVBTFloose(const pat::Muon& m)
{
  // std::cout << "Global muon: " << m.muonID("AllGlobalMuons") ; 
  // if ( m.muonID("AllGlobalMuons") ) std::cout << " " << m.numberOfValidHits() << " hits" ; 
  // std::cout << std::endl ; 
  return (m.muonID("AllGlobalMuons") &&
	  (m.numberOfValidHits() > 10));
}

bool
MuJetBackground::isVBTFtight(const pat::Muon& m)
{
  if (!isVBTFloose(m)) return false;

  reco::TrackRef gt = m.globalTrack();
  if (gt.isNull()) {
    std::cerr << "Mu global track reference is NULL" << std::endl;
    return false;
  }
  return (m.muonID("AllTrackerMuons") &&
	  (m.dB() < 0.2) &&
	  (m.normChi2() < 10) &&
	  (m.numberOfMatches() > 1) &&
	  (gt->hitPattern().numberOfValidMuonHits()>0) &&
	  (gt->hitPattern().numberOfValidPixelHits()>0) );

}                                                // MuJetBackground::isVBTFtight

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
  t="M(W_{R}) "                    +post;        mWR=td->make<TH1D>("mWR",   t.c_str(),40,0,2000);
  t="M(N_{R}) with #mu_{1} "       +post;      mNuR1=td->make<TH1D>("mNuR1", t.c_str(),40,0,2000);
  t="M(N_{R}) with #mu_{2} "       +post;      mNuR2=td->make<TH1D>("mNuR2", t.c_str(),20,0,1000);
  t="M(W_{R}) "                    +post;    mWR_raw=td->make<TH1D>("mWR_raw",   t.c_str(),40,0,2000);
  t="M(N_{R}) with #mu_{1} "       +post;  mNuR1_raw=td->make<TH1D>("mNuR1_raw", t.c_str(),40,0,2000);
  t="M(N_{R}) with #mu_{2} "       +post;  mNuR2_raw=td->make<TH1D>("mNuR2_raw", t.c_str(),20,0,1000);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post;     mNuR2D=td->make<TH2D>("mNuR2D",t.c_str(),40,0,2000,20,0,1000);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post; mNuR2D_raw=td->make<TH2D>("mNuR2D_raw",t.c_str(),40,0,2000,20,0,1000);

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

void MuJetBackground::HistPerDef::fill(pat::Muon theMuon,
				       pat::Jet theJet, 
				       pat::MET theMET, 
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

void MuJetBackground::HistPerDef::fill(reco::SuperCluster theSC)
{  
  float scEt = theSC.energy()/cosh(theSC.eta()) ; 

  if ( fabs(theSC.eta()) < 1.442 ) ptEB->Fill(scEt) ; 
  else if ( fabs(theSC.eta()) > 1.56 && fabs(theSC.eta()) < 2.5 ) ptEE->Fill(scEt) ; 
}// end of fill()


//======================================================================

void
MuJetBackground::HistPerDef::fill(const HeavyNuEvent& hne, double w1, double w2, int nTight)
{
  double weight = w1 * w2 ; 

  const pat::Muon theMu = *(hne.mu1) ; 
  muQual->Fill( double(nTight) ) ; 

  // Muons 
  if (hne.mu1.isAvailable()) { 
    ptMu1->Fill(hne.mu1->pt()) ; 
    ptMu2->Fill(hne.mu2->pt()) ; 
    
    etaMu1->Fill(hne.mu1->eta()) ; 
    etaMu2->Fill(hne.mu2->eta()) ; 

    phiMu1->Fill(hne.mu1->phi()) ; 
    phiMu2->Fill(hne.mu2->phi()) ; 
  }

  // Jets 
  if (hne.j1.isAvailable()) {

    ptJet1->Fill(hne.j1->pt()) ; 
    etaJet1->Fill(hne.j1->eta()) ; 
    phiJet1->Fill(hne.j1->phi()) ; 

    if (hne.j2.isAvailable()) {
      ptJet2->Fill(hne.j2->pt()) ; 
      etaJet2->Fill(hne.j2->eta()) ; 
      phiJet2->Fill(hne.j2->phi()) ; 

      mWR->Fill   ( hne.mWR,weight   ) ; 
      mNuR1->Fill ( hne.mNuR1,w1 ) ; 
      mNuR2->Fill ( hne.mNuR2,w2 ) ; 
      mNuR2D->Fill( hne.mNuR1, hne.mNuR2,weight );
      mJJ->Fill   ( hne.mJJ,weight   );

      mWR_raw->Fill   ( hne.mWR   ) ; 
      mNuR1_raw->Fill ( hne.mNuR1 ) ; 
      mNuR2_raw->Fill ( hne.mNuR2 ) ; 
      mNuR2D_raw->Fill( hne.mNuR1, hne.mNuR2 );
      mJJ_raw->Fill   ( hne.mJJ   );

    }

  }

  mMuMu->Fill( hne.mMuMu,weight );
  mMuMuZoom->Fill( hne.mMuMu,weight );

  mMuMu_raw->Fill( hne.mMuMu );
  mMuMuZoom_raw->Fill( hne.mMuMu );

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

  nnif_ = new HeavyNu_NNIF(iConfig);

  trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet>("trigMatchPset"));

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
    hists.dPhiCutsNICalo10.book      ( new TFileDirectory(fs->mkdir("dPhiCutsNICalo10")), 
				       "(VBTF loose, NI Calo)", v_null );
    hists.dPhiTightCutsNICalo10.book ( new TFileDirectory(fs->mkdir("dPhiTightCutsNICalo10")), 
				       "(VBTF tight, NI Calo)", v_null );
    hists.dPhiHeepSC.book            ( new TFileDirectory(fs->mkdir("dPhiHeepSC")), 
				       "(Dijet Supercluster)", v_null );
    hists.dPhiHeepEle.book           ( new TFileDirectory(fs->mkdir("dPhiHeepEle")), 
				       "(Dijet HEEP electron)", v_null );
  }
  if ( !calcSurvival_ && doQuadJet_ ) { 
    hists.FourJetLooseCuts.book      ( new TFileDirectory(fs->mkdir("FourJetLooseCuts")), 
				       "(two loose muons, two jets)", v_null );
    hists.FourJetTightTrigCuts.book  ( new TFileDirectory(fs->mkdir("FourJetTightTrigCuts")), 
				       "(at least one tight + trigger matched muon)", v_null );
    hists.FourJetVertexCuts.book     ( new TFileDirectory(fs->mkdir("FourJetVertexCuts")), 
				       "(all objects share common vtx)", v_null );
    hists.FourJetMu1PtCuts.book      ( new TFileDirectory(fs->mkdir("FourJetMu1PtCuts")), 
				       "(passes mu1 pT rqmt)", v_null );
    hists.FourJetZMassCuts.book      ( new TFileDirectory(fs->mkdir("FourJetZMassCuts")), 
				       "(passes mumu mass cut)", v_null );
  }

  if ( !calcSurvival_ && doClosure_ ) { 
    hists.LJJJLooseCuts.book         ( new TFileDirectory(fs->mkdir("LJJJLooseCuts")), 
				       "(3 jets, loose muon in jet)", v_null );
    hists.LJJJTightCuts.book         ( new TFileDirectory(fs->mkdir("LJJJTightCuts")), 
				       "(3 jets, tight muon in jet)", v_null );
    hists.LJJLooseCuts.book          ( new TFileDirectory(fs->mkdir("LJJLooseCuts")), 
				       "(2 jets, loose muon)", v_null );
    hists.LJJTightCuts.book          ( new TFileDirectory(fs->mkdir("LJJTightCuts")), 
				       "(2 jets, tight muon)", v_null );
    hists.LJJZrejLooseCuts.book      ( new TFileDirectory(fs->mkdir("LJJZrejLooseCuts")), 
				       "(2 jets, loose muon, Z rejection)", v_null );
    hists.LJJZrejTightCuts.book      ( new TFileDirectory(fs->mkdir("LJJZrejTightCuts")), 
				       "(2 jets, tight muon, Z rejection)", v_null );
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
  cuts.minimum_dijet_dPhi   = iConfig.getParameter<double>("minimumMuJetdPhi");
  cuts.minimum_QCDjet_pt    = iConfig.getParameter<double>("minimumJetPtForDijets");
  cuts.minimum_extraJet_dR  = iConfig.getParameter<double>("minimumDeltaRforExtraJets");
  cuts.minimum_SCEt         = iConfig.getParameter<double>("minimumSuperClusterEt");

  theMETtype  = iConfig.getParameter<int>("METvariety") ; 
  rwLowPtbin  = iConfig.getParameter< std::vector<double> >("reweightPtLow") ; 
  rwHighPtbin = iConfig.getParameter< std::vector<double> >("reweightPtHigh") ; 
  rwLoose     = iConfig.getParameter< std::vector<double> >("reweightLoose") ; 
  rwTight     = iConfig.getParameter< std::vector<double> >("reweightTight") ; 

  // Special check: Make sure all vectors are of the same size
  unsigned int vecsize = rwLowPtbin.size() ; 
  if ( ( doClosure_ || doQuadJet_ ) && 
       ( rwHighPtbin.size() != vecsize ||
 	 rwLoose.size() != vecsize ||
 	 rwTight.size() != vecsize ) )
    throw cms::Exception( "Please ensure that all QCD reweighting vectors are equal size");

  ZwinMinGeV_ = iConfig.getParameter<double>("ZmassWinMinGeV");
  ZwinMaxGeV_ = iConfig.getParameter<double>("ZmassWinMaxGeV");

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
  std::cout << "minimumQCDjetPt  = " << cuts.minimum_QCDjet_pt   << " GeV " << std::endl; 
  std::cout << "minExtraJetdR    = " << cuts.minimum_extraJet_dR << std::endl; 

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

void
MuJetBackground::fillBasicMuHistos(const pat::Muon& m)
{

  //std::cout << "Begin fill basic" << std::endl ; 

  double mupt = m.pt();
  hists.muPt->Fill( applyMESfactor_*mupt ) ; 
  hists.muEta->Fill( m.eta() ) ; 
  hists.muPhi->Fill( m.phi() ) ; 

  if (applyMESfactor_==1.0) {
    hists.mudBvsPt->Fill( mupt, m.dB() );

    if (isVBTFloose(m)) {
      hists.muNvalidHitsVsPt->Fill     ( mupt, m.numberOfValidHits() );
      hists.muNormChi2vsPt->Fill       ( mupt, m.normChi2() );
      hists.muNmatchesVsPt->Fill       ( mupt, m.numberOfMatches() );
      
      reco::TrackRef gt = m.globalTrack();
      // gt.isNonnull() guaranteed at this point?
      hists.muNvalidMuonHitsVsPt->Fill ( mupt, gt->hitPattern().numberOfValidMuonHits() );
      hists.muNvalidPixelHitsVsPt->Fill( mupt, gt->hitPattern().numberOfValidPixelHits() );
    }
    hists.muQualVsPt->Fill( mupt, 0 );
    for (int i=1; i<muonQualityFlags; i++)
      if (m.muonID(muonQuality[i]))
	hists.muQualVsPt->Fill( mupt, i ) ; 
    
    hists.muTrckIsoVsPt->Fill( mupt, m.trackIso() );
    hists.muHcalIsoVsPt->Fill( mupt, m.hcalIso()  );
    hists.muEcalIsoVsPt->Fill( mupt, m.ecalIso()  );
    hists.muCaloIsoVsPt->Fill( mupt, m.caloIso()  );
  }
  //std::cout << "end basic fill" << std::endl ; 
}                                          // MuJetBackground::fillBasicMuHistos

//======================================================================

void
MuJetBackground::fillBasicJetHistos( const pat::Jet& j,
			     int jetnum )
{
  //std::cout << "Begin basic fill" << std::endl ; 

  double jpt=j.pt(),jeta=j.eta();
  float totalunc = 0.0f;

  hists.jetPt ->Fill( jpt  ) ; 
  hists.jetEta->Fill( jeta ) ; 
  hists.jetPhi->Fill( j.phi() ) ; 
  hists.jetID ->Fill( qcdJetID( j ) );
  hists.jetPtvsNum->Fill( jetnum, jpt ) ;
}                                         // MuJetBackground::fillBasicJetHistos

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
MuJetBackground::selectJetsStd(edm::Handle<pat::JetCollection>& pJets,
			       HeavyNuEvent& hne)
{
  for (size_t iJet=0; iJet<pJets->size(); iJet++) {
    pat::JetRef iJ=pat::JetRef( pJets,iJet );
    float jpt       = (*iJ).pt();
    float jeta      = (*iJ).eta();
    float jecuscale = 1.0f;
    // std::cout << "JES: " << jecuscale << std::endl ; 
    // std::cout << "JES: " << jecuscale << std::endl ; 
    if( (jpt       > cuts.minimum_jet_pt)   && // more later!
	(fabs(jeta)<=cuts.maximum_jet_abseta) ) {

      if(hne.j1.isNull()) {
	hne.j2=hne.j1;
	hne.j1=iJ;
	hne.j2scale=hne.j1scale;
	hne.j1scale=jecuscale;
      } else {
	float j1pt = hne.j1->pt() * hne.j1scale;
	if (j1pt < jpt) {
	  hne.j2=hne.j1;
	  hne.j1=iJ;
	  hne.j2scale=hne.j1scale;
	  hne.j1scale=jecuscale;
	}
	else if( hne.j2.isNull() ) {
	  hne.j2=iJ;
	  hne.j2scale=jecuscale;
	} else {
	  float j2pt = hne.j2->pt() * hne.j2scale;
	  if( j2pt < jpt ) {
	    hne.j2=iJ;
	    hne.j2scale=jecuscale;
	  }
	}
      } // if jet supplants one of the jets selected so far
    } // if jet passes pt/eta cuts
  } // jet loop
  if ( hne.j2.isNull() ) return false ; 
  return true ; 
} //HeavyNu::selectJets

bool
MuJetBackground::selectMuons(edm::Handle<pat::MuonCollection>& pMuons,
			     HeavyNuEvent& hne)
{
  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    pat::MuonRef iM=pat::MuonRef(pMuons,iMuon);
    double dr1=(hne.j1.isNull())?(10.0):(deltaR((*iM).eta(),(*iM).phi(),hne.j1->eta(),hne.j1->phi()));
    double dr2=(hne.j2.isNull())?(10.0):(deltaR((*iM).eta(),(*iM).phi(),hne.j2->eta(),hne.j2->phi()));

    double mupt = (*iM).pt();
    if( (mupt > cuts.minimum_mu2_pt)
	&& isVBTFloose(*iM)
	&& (fabs((*iM).eta()) < cuts.maximum_mu_abseta)
	&& (std::min(dr1,dr2) > cuts.minimum_muon_jet_dR)
	// && ((*iM).trackIso()  < cuts.muon_trackiso_limit)
	&& (((*iM).trackIso()/mupt)  < cuts.muon_trackiso_limit) // relative track isolation
	) {
      if( (hne.mu1.isNull()) ||
	  (hne.mu1->pt()<(*iM).pt()) ) { 
	hne.mu1=iM;
	hne.mu2=iM;
      }
    }
  }
  if ( hne.mu1.isNull() ) return false ; 
  return true ; 
}                                                // HeavyNu::selectMuons

bool 
MuJetBackground::selectJets(edm::Handle<pat::JetCollection>& pJets,
			    HeavyNuEvent& hne)
{
  for (size_t iJet=0; iJet<pJets->size(); iJet++) {
    pat::JetRef iJ=pat::JetRef( pJets,iJet );
    float jpt       = (*iJ).pt();
    float jeta      = (*iJ).eta();
    float jecuscale = 1.0f;
    // std::cout << "JES: " << jecuscale << std::endl ; 
    // std::cout << "JES: " << jecuscale << std::endl ; 

    // std::cout << "BMD: Looking at jet with pT " << jpt << " and eta " << jeta << std::endl ; 

    if ( (jpt       > cuts.minimum_jet_pt)   && // more later!
	 (fabs(jeta)<=cuts.maximum_jet_abseta) ) {
      double dr1=deltaR(hne.mu1->eta(),hne.mu1->phi(),(*iJ).eta(),(*iJ).phi()) ;
      double dr2=(hne.mu2.isNull()) ? (10.0) : 
	(deltaR(hne.mu2->eta(),hne.mu2->phi(),(*iJ).eta(),(*iJ).phi())) ;

      // std::cout << "BMD: dr values are " << dr1 << ", " << dr2 << std::endl ; 
      if ( std::min(dr1,dr2) > cuts.minimum_muon_jet_dR ) {
	// std::cout << "BMD: Jet is far enough away...keeping" << std::endl ; 
	if ( (hne.j1.isNull()) ||
	     (hne.j1->pt()<(*iJ).pt()) ) { 
	  // if ( hne.j1.isNull() ) std::cout << "BMD: Setting j1" << std::endl ; 
	  // else std::cout << "BMD: Resetting j1" << std::endl ; 
	  hne.j2=hne.j1;
	  hne.j1=iJ;
	  hne.j2scale=hne.j1scale;
	  hne.j1scale=jecuscale;
	} else if (hne.j2.isNull() ||
		   hne.j2->pt()<(*iJ).pt()) { 
	  hne.j2=iJ;
	  hne.j2scale=jecuscale;
	  // std::cout << "BMD: Setting j2" << std::endl ; 
	}
      }
    }
  }
  if (hne.j2.isNull()) return false ; 
  return true ; 
} //MuJetBackground::selectJets

bool
MuJetBackground::selectQCDjet(edm::Handle<pat::JetCollection>& pJets,
			      edm::Handle<pat::MuonCollection>& pMuons,
			      HeavyNuEvent& hne)
{
  pat::JetRef iJqcd ;
  if ( !iJqcd.isNull() ) {
    std::cout << "Warning: Non-null initial value" << std::endl ; 
    return false ; 
  } 
  float maxdPhi = cuts.minimum_dijet_dPhi ; 
  for (size_t iJet=0; iJet<pJets->size(); iJet++) {
    pat::JetRef iJ=pat::JetRef( pJets,iJet ) ;
    float jpt       = (*iJ).pt();
    float jeta      = (*iJ).eta();
    if( (jpt       > cuts.minimum_QCDjet_pt) && 
	(fabs(jeta)<=cuts.maximum_jet_abseta) ) {
      float dPhi = fabs(deltaPhi((*iJ).phi(),hne.mu1->phi())) ; 
      if ( dPhi > maxdPhi ) { maxdPhi = dPhi ; iJqcd = iJ ; } 
    }
  }

  if ( iJqcd.isNull() ) return false ; 
  
  for (size_t iJet=0; iJet<pJets->size(); iJet++) {
    pat::JetRef iJ=pat::JetRef( pJets,iJet ) ;
    float jpt       = (*iJ).pt();
    float jeta      = (*iJ).eta();
    if( (jpt       > 2*cuts.minimum_QCDjet_pt) && 
	(fabs(jeta)<=cuts.maximum_jet_abseta) ) {
      float dRm = deltaR((*iJ).eta(),(*iJ).phi(),hne.mu1->eta(),hne.mu1->phi()) ; 
      float dRj = deltaR((*iJ).eta(),(*iJ).phi(),(*iJqcd).eta(),(*iJqcd).phi()) ; 
      if ( std::min(dRm,dRj) > cuts.minimum_extraJet_dR ) return false ; 
    }
  }

  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    pat::MuonRef iM=pat::MuonRef( pMuons,iMuon ) ;
    double dR  = deltaR((*iM).eta(),(*iM).phi(),hne.mu1->eta(),hne.mu1->phi()) ; 
    double dpT = fabs( (*iM).pt() - hne.mu1->pt() ) ; 
    if ( dR < 0.001 && dpT < 0.001 ) continue ; // Matched to original muon
    dR = deltaR((*iM).eta(),(*iM).phi(),(*iJqcd).eta(),(*iJqcd).phi()) ; 
    if ( dR < 0.5 ) { // muon is in the jet 
      double pTratio = (*iM).pt() / (*iJqcd).pt() ; 
      if ( pTratio > 0.75 ) return false ; 
    }
  }

  hne.j1 = iJqcd ; 
  return true ; 
} //MuJetBackground::selectQCDjet

//======================================================================

bool 
MuJetBackground::selectMuonsInJets(edm::Handle<pat::MuonCollection>& pMuons,
				   edm::Handle<pat::JetCollection>& pJets,
				   HeavyNuEvent& hne)
{
  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    pat::MuonRef iM=pat::MuonRef(pMuons,iMuon);
    if ( (*iM).pt() < cuts.minimum_mu2_pt ) continue ; 
    // std::cout << "BMD: Found muon candidate with pT: " << (*iM).pt() 
    // 	      << ", eta " << (*iM).eta() 
    // 	      << std::endl ; 
    if ( isVBTFloose(*iM) && 
	 (*iM).pt() > cuts.minimum_mu2_pt &&
	 fabs((*iM).eta()) < cuts.maximum_mu_abseta ) { 
      // std::cout << "BMD: Muon is LOOSE" << std::endl ; 
      bool muInJet = false ; 
      for (size_t iJet=0; iJet<pJets->size(); iJet++) {
	pat::JetRef iJ=pat::JetRef(pJets,iJet);
	if ( qcdJetID(*iJ) < 1 || (*iJ).pt() < cuts.minimum_QCDjet_pt ) continue ;
	// std::cout << "BMD: Jet with pT " << (*iJ).pt() << std::endl ; 
	double dR = deltaR((*iJ).eta(),(*iJ).phi(),(*iM).eta(),(*iM).phi()) ;
	if ( dR < cuts.minimum_muon_jet_dR ) { 
	  // if ( ((*iM).pt()/(*iJ).pt()) > 0.75 ) 
	  // std::cout << "BMD: Muon/jet energy: " << ((*iM).pt()/(*iJ).pt()) << std::endl ; 
	  // if ( ((*iM).pt()/(*iJ).pt()) > 0.75 ) continue ; // muon reco'd as jet
	  muInJet = true ; 
	}
	if ( muInJet ) break ; 
      }
      if ( muInJet ) { 
	// std::cout << "BMD: Dirty muon candidate with pT: " << (*iM).pt() << std::endl ; 
	if ( hne.mu1.isNull() || hne.mu1->pt()<(*iM).pt() ) {
	  // if ( !hne.mu1.isNull() ) std::cout << "BMD: Resetting mu1" << std::endl ; 
	  // else std::cout << "BMD: Setting mu1" << std::endl ; 
	  hne.mu2=hne.mu1;
	  hne.mu1=iM;
	} else if ( hne.mu2.isNull() || hne.mu2->pt()<(*iM).pt() ) {
	  // std::cout << "BMD: Setting mu2 only" << std::endl ; 
	  hne.mu2=iM;
	}
      }
    }
  }
  if ( hne.mu1.isNull() ) return false ; 
  return true ; 
}

bool
MuJetBackground::selectQCDmuon(edm::Handle<pat::MuonCollection>& pMuons,
			       HeavyNuEvent& hne)
{
  pat::MuonRef iMqcd ; 
  if ( !iMqcd.isNull() ) {
    std::cout << "Warning: Non-null initial value" << std::endl ; 
    return false ; 
  } 
  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    pat::MuonRef iM=pat::MuonRef(pMuons,iMuon);
    if( isVBTFloose(*iM) ) {
      if( (iMqcd.isNull()) ||
	  (iMqcd->pt()<(*iM).pt()) ) { 
	iMqcd = iM ;
      }
    }
  }

  if (iMqcd.isNull()) return false ; 

  // Reject non-dijet events 
  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
    pat::MuonRef iM=pat::MuonRef(pMuons,iMuon);
    if( isVBTFloose(*iM) && (iM != iMqcd) ) {
      float relIso = ( (*iM).trackIso() + (*iM).hcalIso() + (*iM).ecalIso() ) / (*iM).pt() ;
      if ( relIso < 0.15 ) return false ; // 2nd quality muon found
      reco::Particle::LorentzVector vMuMu = (*iMqcd).p4() + (*iM).p4() ;
      if ( vMuMu.M() > 70. && vMuMu.M() < 110. ) return false ; // Z->mumu 
    }
  }
  
  hne.mu1 = iMqcd ; 
  return true ; 

} // MuJetBackground::selectQCDmuon

bool MuJetBackground::isDijetCandidate(HeavyNuEvent& hne,pat::MET theMET) {

  if ( theMET.et() < 20. ) return false ; // Absolute MET requirement
  if ( qcdJetID(*(hne.j1)) < 1 ) return false ; // Require at least PURE09 Loose
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

  // float maxdPhi = 1.570796 ; float minPtQCD = 10. ; 

  hnuEvent.isMC = !iEvent.isRealData();
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
    firstEvent_ = false;
  }

  // std::cout << "Test 1" << std::endl ; 

  hists.nelec->Fill(pElecs->size()) ;
  hists.nmu  ->Fill(pMuons->size()) ;
  hists.njet ->Fill(pJets->size()) ;

  if ( pMuons->size() < 1 || pJets->size() < 1 ) return false ; 

  if ( calcSurvival_ ) { 
    if ( selectQCDmuon( pMuons, hnuEvent ) ) {
      if ( selectQCDjet( pJets,pMuons,hnuEvent ) ) { 
	if ( isDijetCandidate( hnuEvent,pMET->at(0) ) ) {
	  
	  bool nonIsoCalo10 = false ; 
	  if ( hnuEvent.mu1->ecalIso()+hnuEvent.mu1->hcalIso() > 10. ) 
	    nonIsoCalo10 = true ; 
	  
	  if ( isVBTFtight(*(hnuEvent.mu1)) ) { 
	    if ( nonIsoCalo10 ) 
	      hists.dPhiTightCutsNICalo10.fill( *(hnuEvent.mu1),*(hnuEvent.j1),
						pMET->at(0),hnuEvent.isMC,
						cuts.muon_trackiso_limit ) ; 
	  } else { 
	    if ( nonIsoCalo10 ) 
	      hists.dPhiCutsNICalo10.fill( *(hnuEvent.mu1),*(hnuEvent.j1),
	   				   pMET->at(0),hnuEvent.isMC,
					   cuts.muon_trackiso_limit) ; 
	  }
	}
      }
    }

    // Search for jet containing muon
    HeavyNuEvent hnuHeepEvent;
    hnuHeepEvent.isMC = !iEvent.isRealData();
    if ( selectMuonsInJets(pMuons,pJets,hnuHeepEvent) ) { 
      //std::cout << "Found a muon in jet" << std::endl ; 
      if ( hnuHeepEvent.mu2.isNull() ) { 
	// Find a jet (with a muon) and a supercluster back-to-back
	pat::JetRef iJheep ; 
	for (size_t iJet=0; iJet<pJets->size(); iJet++) {
	  pat::JetRef iJ=pat::JetRef( pJets,iJet ) ;
	  if ( qcdJetID(*iJ) < 1 ) continue ; // Skip bad jets
	  float dRmin = 0.5 ; float pTratio = 0.75 ; 
	  float dR = deltaR((*iJ).eta(),(*iJ).phi(),hnuHeepEvent.mu1->eta(),hnuHeepEvent.mu1->phi()) ; 
	  if ( (dR < dRmin) && 
	       ((hnuHeepEvent.mu1->pt()/(*iJ).pt()) < pTratio) ) { 
	    iJheep = iJ ; dRmin = dR ;
	    //std::cout << "Found matching jet: " << dR << std::endl ; 
	  } // Get the closest "valid" jet 
	}

	// Search for supercluster back-to-back
	if ( !iJheep.isNull() ) {

	  float maxSCJetdPhi = cuts.minimum_dijet_dPhi ; 
	  reco::SuperClusterCollection::const_iterator scIter ; 
	  reco::SuperCluster theSC ; 
	  bool foundSC = false ; 
	  for (scIter=hybridClusters->begin(); scIter!=hybridClusters->end(); scIter++) {
	    float scEt = scIter->energy()/cosh(scIter->eta()) ; 
	    if ( scEt < cuts.minimum_SCEt ) continue ;  
	    if ( fabs(scIter->eta()) < 1.442 || 
		 ( fabs(scIter->eta()) > 1.56 && fabs(scIter->eta()) < 2.5 ) ) { 
	      float dPhi = fabs(deltaPhi((*iJheep).phi(),scIter->phi())) ;
	      if ( dPhi > maxSCJetdPhi ) { foundSC = true ; theSC = (*scIter) ; maxSCJetdPhi = dPhi ; }
	    }
	  }
	  for (scIter=multi5x5Clusters->begin(); scIter!=multi5x5Clusters->end(); scIter++) {
	    float scEt = scIter->energy()/cosh(scIter->eta()) ; 
	    if ( scEt < cuts.minimum_SCEt ) continue ;  
	    if ( fabs(scIter->eta()) < 1.442 || 
		 ( fabs(scIter->eta()) > 1.56 && fabs(scIter->eta()) < 2.5 ) ) { 
	      float dPhi = fabs(deltaPhi((*iJheep).phi(),scIter->phi())) ;
	      if ( dPhi > maxSCJetdPhi ) { foundSC = true ; theSC = (*scIter) ; maxSCJetdPhi = dPhi ; }
	    }
	  }

	  if ( foundSC ) { // Supercluster found...now impose additional requirements
	    pat::PhotonRef thePhoton ;
	    bool matched = false ; 
	    float scEt  = theSC.energy()/cosh(theSC.eta()) ; 
	    float dRmin = 0.01 ; float dpTmaxRatio = 0.01 ;
	    for (size_t iGam=0; iGam<pGammas->size(); iGam++) {
	      pat::PhotonRef iP=pat::PhotonRef( pGammas,iGam ) ;
	      float dR  = deltaR((*iP).superCluster()->eta(),(*iP).superCluster()->phi(),theSC.eta(),theSC.phi()) ; 
	      float gamEt = iP->superCluster()->energy()/cosh(iP->superCluster()->eta()) ; 
	      float dpTratio = fabs(scEt-gamEt)/gamEt ;
	      if ( dR < dRmin && dpTratio < dpTmaxRatio ) { 
		dRmin = dR ; dpTmaxRatio = dpTratio ; matched = true ; 
		thePhoton = iP ; 
	      }
	    }
	    foundSC = foundSC && matched ; 
	    if ( matched ) { 
	      bool passHoE     = (thePhoton->hadronicOverEm() < 0.05 ) ; 
	      float threshold = 6. + 0.01 * scEt ; 
	      if ( fabs(theSC.eta()) > 1.56 && fabs(theSC.eta()) < 2.5 ) { 
		float scale = scEt - 50. ; 
		if ( scale < 0. ) scale = 0. ; 
		threshold = 6. + 0.01 * scale ; 
	      }
	      bool passEcalIso = (thePhoton->ecalRecHitSumEtConeDR03() < threshold ) ;
	      foundSC = foundSC && ( passHoE && passEcalIso ) ; 
	    // } else { 
	    //   std::cout << "Did not find a photon match to the SC" << std::endl ; 
	    }
	  }
	  
	  if ( foundSC ) { // Reject event if additional jets are found
	    bool noExtraJets = true ; 
	    for (size_t iJet=0; iJet<pJets->size(); iJet++) {
	      pat::JetRef iJ=pat::JetRef( pJets,iJet ) ;
	      float jpt       = (*iJ).pt();
	      float jeta      = (*iJ).eta();
	      if( (jpt       > 2*cuts.minimum_QCDjet_pt) && 
		  (fabs(jeta)<=cuts.maximum_jet_abseta) ) {
		float dRsc = deltaR((*iJ).eta(),(*iJ).phi(),theSC.eta(),theSC.phi()) ; 
		float dRj  = deltaR((*iJ).eta(),(*iJ).phi(),(*iJheep).eta(),(*iJheep).phi()) ; 
		if ( std::min(dRsc,dRj) > cuts.minimum_extraJet_dR ) noExtraJets = false ; 
	      }
	    }
	    foundSC = foundSC && noExtraJets ; 
	  }

	  if ( foundSC ) { // Found a supercluster opposite a mu-in-jet
	    float scEt  =  theSC.energy()/cosh(theSC.eta()) ; 
	    // std::cout << "Supercluster found with ET " << scEt
	    // 	      << " and dPhi " << maxSCJetdPhi << std::endl ; 
	    // std::cout << "Checking MET: " << pMET->at(0).et() << std::endl ; 
	    if ( pMET->at(0).et() < 20. ) { // Decent di-jet candidate
	      //std::cout << "Filling" << std::endl ; 
	      hists.dPhiHeepSC.fill( theSC ) ; 
	      // Apply HEEP selection
	      pat::ElectronRef theEle ; 
	      //std::cout << "Searching for the electron: " << pElecs->size() << std::endl ; 
	      float dRmin = 0.01 ; float dpTmaxRatio = 0.01 ; 
	      for (size_t iEle=0; iEle<pElecs->size(); iEle++) {
		pat::ElectronRef iE=pat::ElectronRef( pElecs,iEle ) ;
		float dR  = deltaR((*iE).superCluster()->eta(),(*iE).superCluster()->phi(),theSC.eta(),theSC.phi()) ; 
		//std::cout << "Computing dR: " << dR << std::endl ; 
		float eleEt =  iE->superCluster()->energy()/cosh(iE->superCluster()->eta()) ; 
		//std::cout << "Computing eleEt: " << eleEt << std::endl ; 
		float dpTratio = fabs(scEt-eleEt)/eleEt ; 
		if ( (dR < dRmin) && (dpTratio < dpTmaxRatio) ) { 
		  dRmin = dR ; dpTmaxRatio = dpTratio ; 
		  theEle = iE ; 
		  //std::cout << "Electron match with " << dR << std::endl ; 
		}
	      }

	      bool validHEEP = false ; 
	      if ( !theEle.isNull() && theEle->ecalDrivenSeed() ) { 
		float eleEt =  theEle->superCluster()->energy()/cosh(theEle->superCluster()->eta()) ; 
		//std::cout << "Et: " << eleEt << std::endl ; 
		if ( fabs(theSC.eta()) < 1.442 ) { // barrel 
		  //std::cout << "EB" << std::endl ; 
		  if ( theEle->deltaEtaSuperClusterTrackAtVtx() < 0.005 ) { 
		    if ( theEle->deltaPhiSuperClusterTrackAtVtx() < 0.09 ) { 
		      if ( theEle->hadronicOverEm() < 0.05 ) { 
			if ( (theEle->e2x5Max()/theEle->e5x5() > 0.94 ) ||
			     (theEle->e1x5()/theEle->e5x5() > 0.83 ) ) { 
			  float caloIso = theEle->dr03EcalRecHitSumEt() + theEle->dr03HcalDepth1TowerSumEt() ; 
			  float caloIsoLimit = 2.0 + (0.03 * eleEt) ;  
			  if ( caloIso < caloIsoLimit ) {
			    if ( theEle->dr03TkSumPt() < 7.5 ) validHEEP = true ; 
			  }
			}
		      }
		    }
		  }
		} else if ( fabs(theSC.eta()) > 1.56 && fabs(theSC.eta()) < 2.5 ) { // endcap
		  //std::cout << "EE" << std::endl ; 
		  if ( theEle->deltaEtaSuperClusterTrackAtVtx() < 0.007 ) { 
		    if ( theEle->deltaPhiSuperClusterTrackAtVtx() < 0.09 ) { 
		      if ( theEle->hadronicOverEm() < 0.05 ) { 
			if ( theEle->sigmaIetaIeta() < 0.03 ) { 
			  float caloIso = theEle->dr03EcalRecHitSumEt() + theEle->dr03HcalDepth1TowerSumEt() ; 
			  float scale = eleEt - 50. ; if ( scale < 0. ) scale = 0. ; 
			  float caloIsoLimit = 2.5 + (0.03 * scale) ;  
			  if ( (caloIso < caloIsoLimit) && 
			       (theEle->dr03HcalDepth2TowerSumEt() < 0.5) ) {
			    if ( theEle->dr03TkSumPt() < 15. ) validHEEP = true ; 
			  }
			}
		      }
		    }
		  }
		}
	      } 
	      //if ( validHEEP ) std::cout << "Filling electron" << std::endl ; 
	      if ( validHEEP ) hists.dPhiHeepEle.fill( *(theEle->superCluster()) ) ; 
	    }
	  }
	}
      }
    }
  }

  // Multi-jet analysis  
  if ( !calcSurvival_ && doQuadJet_ ) { 
    if ( pMuons->size() >= 2 && pJets->size() >= 4 ) { 
      
      HeavyNuEvent hnuEvent4jet;
      hnuEvent4jet.isMC = !iEvent.isRealData();
      int nTightMuons = 0 ; 

      if ( selectMuonsInJets( pMuons,pJets,hnuEvent4jet ) && !hnuEvent4jet.mu2.isNull() ) { 
	if ( selectJets( pJets,hnuEvent4jet ) ) { 
	  keepThisEvent = true ; 
	  
	  hnuEvent4jet.regularize(); 
	  hnuEvent4jet.calculateMuMu(applyMESfactor_);
	  hnuEvent4jet.calculate() ; 
	  bool mu1isTight = isVBTFtight(*(hnuEvent4jet.mu1));
	  bool mu2isTight = isVBTFtight(*(hnuEvent4jet.mu2));
	  double mu1scaleFactor, mu2scaleFactor ; 
	  if ( mu1isTight ) { 
	    nTightMuons++ ; 
	    mu1scaleFactor = qcdScaleFactor(hnuEvent4jet.mu1->pt(),
					    rwLowPtbin,rwHighPtbin,rwTight) ; 
	  } else mu1scaleFactor = qcdScaleFactor(hnuEvent4jet.mu1->pt(),
						 rwLowPtbin,rwHighPtbin,rwLoose) ; 
	  if ( mu2isTight ) { 
	    nTightMuons++ ; 
	    mu2scaleFactor = qcdScaleFactor(hnuEvent4jet.mu2->pt(),
					    rwLowPtbin,rwHighPtbin,rwTight) ; 
	  } else mu2scaleFactor = qcdScaleFactor(hnuEvent4jet.mu2->pt(),
						 rwLowPtbin,rwHighPtbin,rwLoose) ; 	    
	  // std::cout << "BMD: 4 jj   mass = " << hnuEvent4jet.mJJ << std::endl ; 
	  // std::cout << "BMD: 4 hnu1 mass = " << hnuEvent4jet.mNuR1 << std::endl ; 
	  // std::cout << "BMD: 4 wr   mass = " << hnuEvent4jet.mWR << std::endl ; 
	  	  
	  if ( fabs(hnuEvent4jet.mu1->vertex().Z()-hnuEvent4jet.mu2->vertex().Z()) 
	       < cuts.maxVertexZsep ) {
	    
	    if ( (qcdJetID(*(hnuEvent4jet.j1)) >= 1) &&   
		 (qcdJetID(*(hnuEvent4jet.j2)) >= 1) ) { 

	      //std::cout << "Candidate passes 2 mu + 2 jet requirement" << std::endl ; 
	      hists.FourJetLooseCuts.fill( hnuEvent4jet,mu1scaleFactor,mu2scaleFactor,nTightMuons ) ; 
	    
	      if ( mu1isTight || mu2isTight ) { 
		
		//std::cout << "One of the muons passes tight requirements" << std::endl ; 
	      
		bool mu1trig = mu1isTight ; bool mu2trig = mu2isTight ; 
		if ( trig_->matchingEnabled() &&
		     iEvent.isRealData() ) {
		  // require that one muon be BOTH tight and trigger-matched
		  mu1trig = mu1trig && trig_->isTriggerMatched( hnuEvent4jet.mu1, iEvent) ; 
		  mu2trig = mu2trig && trig_->isTriggerMatched( hnuEvent4jet.mu2, iEvent) ; 
		} else if ( !iEvent.isRealData() ) {
		  mu1trig = mu1trig && trig_->simulateForMC( hnuEvent4jet.mu1->pt(),hnuEvent4jet.mu1->eta(),0 );
		  mu2trig = mu2trig && trig_->simulateForMC( hnuEvent4jet.mu2->pt(),hnuEvent4jet.mu2->eta(),0 );
		}

		if ( mu1trig || mu2trig ) { 		  
		  hists.FourJetTightTrigCuts.fill( hnuEvent4jet,mu1scaleFactor,mu2scaleFactor,nTightMuons ) ; 
		  
		  // bool commonVertex = ( fabs(hnuEvent4jet.mu1->vertex().Z()-hnuEvent4jet.mu2->vertex().Z()) 
		  // 			< cuts.maxVertexZsep ) ; 
		  // commonVertex = commonVertex && ( fabs(hnuEvent4jet.mu1->vertex().Z()-hnuEvent4jet.j1->vertex().Z()) 
		  // 				   < cuts.maxVertexZsep ) ; 
		  // commonVertex = commonVertex && ( fabs(hnuEvent4jet.mu1->vertex().Z()-hnuEvent4jet.j2->vertex().Z()) 
		  // 				   < cuts.maxVertexZsep ) ; 
		  // commonVertex = commonVertex && ( fabs(hnuEvent4jet.mu2->vertex().Z()-hnuEvent4jet.j1->vertex().Z()) 
		  // 				   < cuts.maxVertexZsep ) ; 
		  // commonVertex = commonVertex && ( fabs(hnuEvent4jet.mu2->vertex().Z()-hnuEvent4jet.j2->vertex().Z()) 
		  // 				   < cuts.maxVertexZsep ) ; 
		  // commonVertex = commonVertex && ( fabs(hnuEvent4jet.j1->vertex().Z()-hnuEvent4jet.j2->vertex().Z()) 
		  // 				   < cuts.maxVertexZsep ) ; 		    
		  // if ( commonVertex ) { 
		  //   hists.FourJetVertexCuts.fill( hnuEvent4jet,mu1scaleFactor,mu2scaleFactor ) ; 

		  if ( hnuEvent4jet.mu1->pt() > cuts.minimum_mu1_pt ) {
		    hists.FourJetMu1PtCuts.fill( hnuEvent4jet,mu1scaleFactor,mu2scaleFactor,nTightMuons ) ; 

		    if ( hnuEvent4jet.mMuMu > cuts.minimum_mumu_mass ) { // dimuon mass rqmt
		      hists.FourJetZMassCuts.fill( hnuEvent4jet,mu1scaleFactor,mu2scaleFactor,nTightMuons ) ; 
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // Closure test, using three jets
  if ( !calcSurvival_ && doClosure_ ) { 
    if ( pMuons->size() >= 1 && pJets->size() >= 3
	 && pMET->at(0).pt() < 20. ) { 

      HeavyNuEvent hnuEvent3jet;
      hnuEvent3jet.isMC = !iEvent.isRealData();
      
      if ( selectMuonsInJets( pMuons,pJets,hnuEvent3jet ) ) { 
	if ( selectJets( pJets,hnuEvent3jet ) ) { 
	  
	  hnuEvent3jet.mu2 = hnuEvent3jet.mu1 ; // HACK
	  hnuEvent3jet.regularize(); // assign internal standards
	  hnuEvent3jet.calculateMuMu(applyMESfactor_);
	  hnuEvent3jet.calculate() ; 
	  // std::cout << "BMD: 3 jj =   " << hnuEvent3jet.mJJ << std::endl ; 
	  // std::cout << "BMD: 3 hnu1 = " << hnuEvent3jet.mNuR1 << std::endl ; 

	  bool mu1isTight = isVBTFtight(*(hnuEvent3jet.mu1));
	  bool mu1trig = mu1isTight ; 
	  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
	    mu1trig = mu1trig && trig_->isTriggerMatched( hnuEvent3jet.mu1, iEvent) ; 
	  } else if (!iEvent.isRealData()) {
	    mu1trig = mu1trig &&  trig_->simulateForMC( hnuEvent3jet.mu1->pt(),hnuEvent3jet.mu1->eta(),0 );
	  }
	  if ( (qcdJetID(*(hnuEvent3jet.j1)) >= 1) &&   
	       (qcdJetID(*(hnuEvent3jet.j2)) >= 1) ) {
	    double mu1scaleFactor ; 
	    if ( mu1isTight ) mu1scaleFactor = qcdScaleFactor(hnuEvent3jet.mu1->pt(),
							      rwLowPtbin,rwHighPtbin,rwTight) ; 
	    else mu1scaleFactor = qcdScaleFactor(hnuEvent3jet.mu1->pt(),
						 rwLowPtbin,rwHighPtbin,rwLoose) ; 

	    // Establish common vertex for the events
	    // bool commonVertex = ( fabs(hnuEvent3jet.mu1->vertex().Z()-hnuEvent3jet.j1->vertex().Z()) 
	    // 			  < cuts.maxVertexZsep ) ; 
	    // commonVertex = commonVertex && ( fabs(hnuEvent3jet.mu1->vertex().Z()-hnuEvent3jet.j2->vertex().Z()) 
	    // 				     < cuts.maxVertexZsep ) ; 
	    // commonVertex = commonVertex && ( fabs(hnuEvent3jet.j1->vertex().Z()-hnuEvent3jet.j2->vertex().Z()) 
	    // 				     < cuts.maxVertexZsep ) ; 

	    // if ( commonVertex ) { 
	    if ( mu1trig ) hists.LJJJTightCuts.fill( hnuEvent3jet, mu1scaleFactor ) ; 
	    if ( !mu1isTight ) hists.LJJJLooseCuts.fill( hnuEvent3jet, mu1scaleFactor ) ; 
	    // }
	  }
	}
      }
    }

    // Closure test continued
    if ( pMuons->size() >= 1 && pJets->size() >= 2 
	 && pMET->at(0).pt() < 20. ) { // MET cut to remove W+2 jets

      HeavyNuEvent hnuEvent2jet;
      hnuEvent2jet.isMC = !iEvent.isRealData();
      
      if ( selectJetsStd( pJets,hnuEvent2jet ) ) { 
	if ( selectMuons( pMuons,hnuEvent2jet ) ) { 
	  
	  hnuEvent2jet.mu2 = hnuEvent2jet.mu1 ; // HACK
	  hnuEvent2jet.regularize(); // assign internal standards
	  hnuEvent2jet.calculateMuMu(applyMESfactor_);
	  hnuEvent2jet.calculate() ; 
	  // std::cout << "BMD: 2 jj =   " << hnuEvent2jet.mJJ << std::endl ; 
	  // std::cout << "BMD: 2 hnu1 = " << hnuEvent2jet.mNuR1 << std::endl ; 

	  // Special: Check for a Z boson
	  bool zMuMu = false ; 
	  for (size_t iMuon=0; iMuon<pMuons->size(); iMuon++) {
	    pat::MuonRef iM=pat::MuonRef(pMuons,iMuon);
	    double dR  = deltaR((*iM).eta(),(*iM).phi(),hnuEvent2jet.mu1->eta(),hnuEvent2jet.mu1->phi()) ; 
	    double dpT = fabs( (*iM).pt() - hnuEvent2jet.mu1->pt() ) ; 
	    if ( dR > 0.001 || dpT > 0.001 ) { // Not the original muon
	      if( isVBTFloose(*iM) ) {
		float relIso = ( (*iM).trackIso() + (*iM).hcalIso() + (*iM).ecalIso() ) / (*iM).pt() ;
		if ( relIso < 0.15 ) { 
		  reco::Particle::LorentzVector vMuMu = hnuEvent2jet.mu1->p4() + (*iM).p4() ;
		  if ( vMuMu.M() > 70. && vMuMu.M() < 110. ) zMuMu = true ; 
		}
	      }
	    }
	  }

	  bool mu1isTight = isVBTFtight(*(hnuEvent2jet.mu1));
	  bool mu1trig = mu1isTight ; 
	  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
	    mu1trig = mu1trig && trig_->isTriggerMatched( hnuEvent2jet.mu1, iEvent) ; 
	  } else if (!iEvent.isRealData()) {
	    mu1trig = mu1trig && trig_->simulateForMC( hnuEvent2jet.mu1->pt(),hnuEvent2jet.mu1->eta(),0 ) ; 
	  }
	  if ( (qcdJetID(*(hnuEvent2jet.j1)) >= 1) &&   
	       (qcdJetID(*(hnuEvent2jet.j2)) >= 1) ) {
	    double mu1scaleFactor ; 
	    if ( mu1isTight ) mu1scaleFactor = qcdScaleFactor(hnuEvent2jet.mu1->pt(),
							      rwLowPtbin,rwHighPtbin,rwTight) ; 
	    else mu1scaleFactor = qcdScaleFactor(hnuEvent2jet.mu1->pt(),
						 rwLowPtbin,rwHighPtbin,rwLoose) ; 

	    // Establish common vertex for the events
	    // bool commonVertex = ( fabs(hnuEvent2jet.mu1->vertex().Z()-hnuEvent2jet.j1->vertex().Z()) 
	    // 			  < cuts.maxVertexZsep ) ; 
	    // commonVertex = commonVertex && ( fabs(hnuEvent2jet.mu1->vertex().Z()-hnuEvent2jet.j2->vertex().Z()) 
	    // 				     < cuts.maxVertexZsep ) ; 
	    // commonVertex = commonVertex && ( fabs(hnuEvent2jet.j1->vertex().Z()-hnuEvent2jet.j2->vertex().Z()) 
	    // 				     < cuts.maxVertexZsep ) ; 
	    // if ( commonVertex ) { 
	    if ( mu1trig ) hists.LJJTightCuts.fill( hnuEvent2jet,mu1scaleFactor ) ; 
	    if ( !mu1isTight ) hists.LJJLooseCuts.fill( hnuEvent2jet,mu1scaleFactor ) ; 

	    if ( !zMuMu ) { 
	      if ( mu1trig ) hists.LJJZrejTightCuts.fill( hnuEvent2jet,mu1scaleFactor ) ; 
	      if ( !mu1isTight ) hists.LJJZrejLooseCuts.fill( hnuEvent2jet,mu1scaleFactor ) ; 
	    }
	    // }
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
  // trig_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuJetBackground);
