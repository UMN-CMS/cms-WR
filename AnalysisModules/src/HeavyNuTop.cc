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
// $Id: HeavyNuTop.cc,v 1.18 2011/12/24 01:54:14 pastika Exp $
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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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

class comparePt {
public:
  template <class T> bool operator() (const T& a, const T& b) { return a.pt() > b.pt() ; } 
};

static std::string btagName;
static double      minBtagDiscVal; // for discriminating B-tagged jets.

class HeavyNuTop : public edm::EDFilter {
public:
  explicit HeavyNuTop(const edm::ParameterSet&);
  ~HeavyNuTop();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }
  
  virtual void beginJob          ();
  virtual bool filter            ( edm::Event&, const edm::EventSetup& );
  virtual void endJob            ();

  virtual TH1 *bookRunHisto      ( uint32_t runNumber );
  
  // double GetCorrectedPt ( const pat::Electron& e );

  edm::InputTag muonTag_;
  edm::InputTag jetTag_;
  edm::InputTag metTag_;
  edm::InputTag elecTag_;

  double applyMESfactor_;             // for Muon Energy Scale studies
  int    applyTrigEffsign_;           // for Trigger Efficiency studies
  bool   applyMuIDCorrections_ ; 
    // bool   applyEleScaleCorrections_ ; 
    // bool   applyEleIDWeightFactor_ ; 
  bool   studyScaleFactorEvolution_;  // for Top, Z+jets scale factors (by Mu1 pT) studies

  int    pileupEra_;

  double EBscalefactor_, EEscalefactor_ ; 
  double ebIDwgt_, eeIDwgt_ ; 

  std::string currentFile_;
  JetCorrectionUncertainty *jecuObj_;
  bool dolog_;
  bool firstEvent_;
  bool isPFJets_;

  // HeavyNu_NNIF *nnif_;
  HeavyNuTrigger *trig_;
  HeavyNuID *muid_ ; 

  edm::LumiReWeighting MCweightByVertex_;

  std::map<uint32_t,TH1 *> m_runHistos_;

  // ----------member data ---------------------------

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory *, const std::string&, const std::vector<hNuMassHypothesis>&) ;
    // fill all histos of the set with the two electron candidates
    void fill(pat::MuonCollection muons, 
	      pat::ElectronCollection electrons,
	      pat::JetCollection jets,
	      pat::METCollection metc,
	      bool isMC,
	      double wgt) ;
    // fill all histos of the set with the two electron candidates
    void fill(const HeavyNuEvent& hne, const std::vector<hNuMassHypothesis>&) ;

    TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2;
    TH2 *ptMu1VsPtMu2ss,  *ptMu1VsPtMu2os;
    TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2 ;
    TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2 ;
    TH1 *dEtaMu, *dPhiMu, *dEtaJet, *dPhiJet ;
    TH2 *dEtaPhiMu, *dEtaPhiJet ; 
    TH1 *dRminMu1jet, *dRminMu2jet ; 
    TH1 *hptrelMu1, *hptrelMu2 ; 
    TH2 *ptrelVsdRminMu1jet, *ptrelVsdRminMu2jet ;
    TH2 *jetID2d;

    TH1 *dptMu1gen, *dptMu2gen ; 
    TH1 *dRMu1gen, *dRMu2gen ; 
    TH1 *qualMu1, *qualMu2 ; 

    TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1caloIso, *mu1dB;
    TH1 *mu2trackIso, *mu2hcalIso, *mu2ecalIso, *mu2caloIso, *mu2dB;

    TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1caloRelIso;
    TH1 *mu2trackRelIso, *mu2hcalRelIso, *mu2ecalRelIso, *mu2caloRelIso;

    TH1 *mMuMu, *mMuMuOS, *mMuMuSS, *diMuCharge, *mMuMuZoom, *mMuMuGenZoom, *mMuMuBarrel, *mMuMuBB;
    TH1 *mWR, *mNuR1, *mNuR2, *mJJ, *mWRBarrel, *mWRBB ;
    TH2 *mNuR2D, *jetPtvsNum, *mWRvsminLPt, *mWRvsNPV;

    TH1* btagJet1, *btagJet2;
    TH1* numBjets;

    TH1* met;

    TH1* czeta_mumu;
    TH1* czeta_mumu_zoom;

    // Jeremy's crazy angles...
    TH1* ctheta_mumu, *cthetaz_mumu;
    TH1* ctheta_jj, *cthetaz_jj;
    TH1* ctheta_mu1_jj, *cthetaz_mu1_jj;
    TH1* ctheta_mu2_jj, *cthetaz_mu2_jj;

    // Vertex plots
    TH1* vtx_mumu, *vtx_jj, *vtx_mu1jmin, *vtx_mu2jmin ; 

    TFileDirectory *mydir;
    TFileDirectory *nndir;

  };

  bool init_;

  // gf set of histo for all Z definitions in a stack
  struct HistStruct {
    TH1 *nelec, *njet, *nmet, *nmuAll, *nmuLoose, *nmuTight;
    TH1 *muPt, *muEta, *muPhi, *looseMuPt, *tightMuPt ; 

    // Muon quality histos as a function of Pt
    TH2 *muNvalidHitsVsPt, *mudBvsPt, *muNormChi2vsPt, *muQualVsPt;
    TH2 *muNmatchesVsPt, *muNvalidMuonHitsVsPt, *muNvalidPixelHitsVsPt;
    TH2 *muTrckIsoVsPt, *muHcalIsoVsPt, *muEcalIsoVsPt, *muCaloIsoVsPt;

    TH1 *jetPt, *jetEta, *jetPhi, *jetID, *jecUncHi, *jecUncLo, *met ; 
    TH2 *jetPtvsNum;
    TProfile2D *jecUncHiVsEtaPt,*jecUncLoVsEtaPt;

    TH1 *trkIsoStudy;

    TFileDirectory *rundir;
    HistPerDef noCuts; 
    // HistPerDef LLptCuts;
    // HistPerDef MuTightCuts;
    HistPerDef LLJJptCuts;
    HistPerDef TrigMatches;
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
    HistPerDef diLmassCut;
    //HistPerDef diLmassCutEChannel;
    HistPerDef mWRmassCut;

  } hists;

  struct CutsStruct {
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
  
};

//======================================================================

const int muonQualityFlags = 4 ;
const std::string muonQuality[] = {
  "All","AllGlobalMuons","AllStandAloneMuons","AllTrackerMuons"
};

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

//======================================================================

void
HeavyNuTop::HistPerDef::book(TFileDirectory *td,
			  const std::string& post,
			  const std::vector<hNuMassHypothesis>& v_masspts )
{
  std::string t; // histogram title string;
  
  TH1::SetDefaultSumw2();

  mydir = td;

  // ----------  Muon histograms  ----------

  t="p_{T}(#mu_{1}) "+post;   ptMu1=td->make<TH1D>("ptMu1",t.c_str(),100,0.,1000.);
  t="p_{T}(#mu_{2}) "+post;   ptMu2=td->make<TH1D>("ptMu2",t.c_str(),100,0.,1000.);
  t="#eta(#mu_{1}) " +post;  etaMu1=td->make<TH1D>("etaMu1",t.c_str(),40,-2.5,2.5);
  t="#eta(#mu_{2}) " +post;  etaMu2=td->make<TH1D>("etaMu2",t.c_str(),40,-2.5,2.5);
  t="#phi(#mu_{1}) " +post;  phiMu1=td->make<TH1D>("phiMu1",t.c_str(),30,-3.14159,3.14159);
  t="#phi(#mu_{2}) " +post;  phiMu2=td->make<TH1D>("phiMu2",t.c_str(),30,-3.14159,3.14159);

  t="p_{T}(#mu_{1}) vs. p_{T}(#mu_{2}) (SS) "+post+";p_{T}(#mu_{1})(GeV);p_{T}(#mu_{2}(GeV))";

  ptMu1VsPtMu2ss=td->make<TH2D>("ptMu1VsPtMu2ss",t.c_str(),50,0.,2000.,50,0.,2000);

  t="p_{T}(#mu_{1}) vs. p_{T}(#mu_{2}) (OS) "+post+";p_{T}(#mu_{1})(GeV);p_{T}(#mu_{2}(GeV))";

  ptMu1VsPtMu2os=td->make<TH2D>("ptMu1VsPtMu2os",t.c_str(),50,0.,2000.,50,0.,2000);

  // delta angles

  t="#Delta#eta(#mu_{1},#mu_{2}) "  +post;    dEtaMu=td->make<TH1D>("dEtaMu",t.c_str(),40,0,5);
  t="#Delta#phi(#mu_{1},#mu_{2}) "  +post;    dPhiMu=td->make<TH1D>("dPhiMu",t.c_str(),30,0,3.14159);
  t="#Delta p_{T}(#mu_{1},gen) "    +post; dptMu1gen=td->make<TH1D>("dptMu1gen",t.c_str(),50,-0.50,0.50);
  t="#Delta p_{T}(#mu_{2},gen) "    +post; dptMu2gen=td->make<TH1D>("dptMu2gen",t.c_str(),50,-0.50,0.50);
  t="#Delta R(#mu_{1},gen) "        +post;  dRMu1gen=td->make<TH1D>("dRMu1gen", t.c_str(),50,0,0.01);
  t="#Delta R(#mu_{2},gen) "        +post;  dRMu2gen=td->make<TH1D>("dRMu2gen", t.c_str(),50,0,0.01);
  t="#mu #Delta#eta vs. #Delta#phi "+post;
  t+=  ";#Delta#eta; #Delta#phi";          dEtaPhiMu=td->make<TH2D>("dEtaPhiMu",t.c_str(),
								   50,0,5,30,0,3.14159);

  t="Quality (#mu_{1}) "+post; qualMu1=td->make<TH1D>("qualMu1",t.c_str(),muonQualityFlags,0,muonQualityFlags);
  t="Quality (#mu_{2}) "+post; qualMu2=td->make<TH1D>("qualMu2",t.c_str(),muonQualityFlags,0,muonQualityFlags);
  labelMuonQualAxis(qualMu1->GetXaxis());
  labelMuonQualAxis(qualMu2->GetXaxis());

  // isolation

  t=   "trackIso(#mu_{1}) "+post;    mu1trackIso=td->make<TH1D>("mu1trackIso",   t.c_str(),40,0.,200.);
  t=    "hcalIso(#mu_{1}) "+post;     mu1hcalIso=td->make<TH1D>("mu1hcalIso",    t.c_str(),40,0.,200.);
  t=    "ecalIso(#mu_{1}) "+post;     mu1ecalIso=td->make<TH1D>("mu1ecalIso",    t.c_str(),40,0.,200.);
  t=    "caloIso(#mu_{1}) "+post;     mu1caloIso=td->make<TH1D>("mu1caloIso",    t.c_str(),40,0.,200.);
  t=        "Dxy(#mu_{1}) "+post;          mu1dB=td->make<TH1D>("mu1dB",         t.c_str(),50,0.,1.);

  t=   "trackIso(#mu_{2}) "+post;    mu2trackIso=td->make<TH1D>("mu2trackIso",   t.c_str(),40,0.,200.);
  t=    "hcalIso(#mu_{2}) "+post;     mu2hcalIso=td->make<TH1D>("mu2hcalIso",    t.c_str(),40,0.,200.);
  t=    "ecalIso(#mu_{2}) "+post;     mu2ecalIso=td->make<TH1D>("mu2ecalIso",    t.c_str(),40,0.,200.);
  t=    "caloIso(#mu_{2}) "+post;     mu2caloIso=td->make<TH1D>("mu2caloIso",    t.c_str(),40,0.,200.);
  t=        "Dxy(#mu_{2}) "+post;          mu2dB=td->make<TH1D>("mu2dB",         t.c_str(),50,0.,1.);

  t="trackRelIso(#mu_{1}) "+post; mu1trackRelIso=td->make<TH1D>("mu1trackRelIso",t.c_str(),50,0.,5.);
  t= "hcalRelIso(#mu_{1}) "+post;  mu1hcalRelIso=td->make<TH1D>("mu1hcalRelIso", t.c_str(),50,0.,5.);
  t= "ecalRelIso(#mu_{1}) "+post;  mu1ecalRelIso=td->make<TH1D>("mu1ecalRelIso", t.c_str(),50,0.,5.);
  t= "caloRelIso(#mu_{1}) "+post;  mu1caloRelIso=td->make<TH1D>("mu1caloRelIso", t.c_str(),50,0.,5.);

  t="trackRelIso(#mu_{2}) "+post; mu2trackRelIso=td->make<TH1D>("mu2trackRelIso",t.c_str(),50,0.,5.);
  t= "hcalRelIso(#mu_{2}) "+post;  mu2hcalRelIso=td->make<TH1D>("mu2hcalRelIso", t.c_str(),50,0.,5.);
  t= "ecalRelIso(#mu_{2}) "+post;  mu2ecalRelIso=td->make<TH1D>("mu2ecalRelIso", t.c_str(),50,0.,5.);
  t= "caloRelIso(#mu_{2}) "+post;  mu2caloRelIso=td->make<TH1D>("mu2caloRelIso", t.c_str(),50,0.,5.);

  // ----------  Jet histograms ----------

  t="p_{T}(j_{1}) "            +post;     ptJet1=td->make<TH1D>("ptJet1",  t.c_str(),50,0.,500.);
  t="p_{T}(j_{2}) "            +post;     ptJet2=td->make<TH1D>("ptJet2",  t.c_str(),50,0.,500.);
  t= "#eta(j_{1}) "            +post;    etaJet1=td->make<TH1D>("etaJet1", t.c_str(),40,-5,5);
  t= "#eta(j_{2}) "            +post;    etaJet2=td->make<TH1D>("etaJet2", t.c_str(),40,-5,5);
  t= "#phi(j_{1}) "            +post;    phiJet1=td->make<TH1D>("phiJet1", t.c_str(),30,-3.14159,3.14159);
  t= "#phi(j_{2}) "            +post;    phiJet2=td->make<TH1D>("phiJet2", t.c_str(),30,-3.14159,3.14159);

  t="#Delta#eta(j_{1},j_{2}) " +post;    dEtaJet=td->make<TH1D>("dEtaJet", t.c_str(),40,0,5);
  t="#Delta#phi(j_{1},j_{2}) " +post;    dPhiJet=td->make<TH1D>("dPhiJet", t.c_str(),30,0,3.14159);

  t= "btag(j_{1}) "            +post;   btagJet1=td->make<TH1D>("btagJet1",t.c_str(),40,0,5);
  t= "btag(j_{2}) "            +post;   btagJet2=td->make<TH1D>("btagJet2",t.c_str(),40,0,5);

  t="# B-tagged Jets in Event "+post;   numBjets=td->make<TH1D>("numBjets",t.c_str(),3,-0.5,2.5);

  t ="Jet #Delta#eta vs. #Delta#phi ";
  t+= post+";#Delta#eta; #Delta#phi";  dEtaPhiJet=td->make<TH2D>("dEtaPhiJet",t.c_str(),
								50,0,5,30,0,3.14159);
  t ="Jet ID(j_{2}) vs. ID(j_{1}) ";
  t+=post+"; ID(j_{1}); ID(j_{2})";         jetID2d=td->make<TH2I>("jetID2d",t.c_str(),3,0,3,3,0,3);
  jetID2d->GetXaxis()->SetBinLabel(1,"Neither");
  jetID2d->GetXaxis()->SetBinLabel(2,"PURE09 Loose");
  jetID2d->GetXaxis()->SetBinLabel(3,"PURE09 Tight");
  jetID2d->GetYaxis()->SetBinLabel(1,"Neither");
  jetID2d->GetYaxis()->SetBinLabel(2,"PURE09 Loose");
  jetID2d->GetYaxis()->SetBinLabel(3,"PURE09 Tight");

  // ----------  MET histograms     ----------

  t="MET distribution "             +post;          met=td->make<TH1D>("met",t.c_str(),100,0,2000);

  // ----------  Mu/Jet histograms  ----------

  t="Minimum #Delta R(#mu_{1},jet) "+post;  dRminMu1jet=td->make<TH1D>("dRminMu1jet",t.c_str(),50,0,5.);
  t="Minimum #Delta R(#mu_{2},jet) "+post;  dRminMu2jet=td->make<TH1D>("dRminMu2jet",t.c_str(),50,0,5.);

  t="p_{T,rel}(#mu_{1},jet)"+post;            hptrelMu1=td->make<TH1D>("ptrelMu1",t.c_str(),50,0,1000.);
  t="p_{T,rel}(#mu_{2},jet)"+post;            hptrelMu2=td->make<TH1D>("ptrelMu2",t.c_str(),50,0,1000.);

  t="p_{T,rel}(#mu_{1},jet) vs #Delta R(#mu_{1},jet)"+post;
  t+="; #Delta R(#mu_{1},jet); p_{T,rel}(#mu_{1},jet)"; ptrelVsdRminMu1jet=td->make<TH2D>("ptrelVsdRminMu1jet",
											 t.c_str(),
											 50,0,5.,50,0,1000);
  t="p_{T,rel}(#mu_{2},jet) vs #Delta R(#mu_{2},jet)"+post;
  t+="; #Delta R(#mu_{2},jet); p_{T,rel}(#mu_{2},jet)"; ptrelVsdRminMu2jet=td->make<TH2D>("ptrelVsdRminMu2jet",
											 t.c_str(),
											 50,0,5.,50,0,1000);
  
  // ----------  Composite histograms  ----------

  t="M(W_{R}) "                    +post;       mWR=td->make<TH1D>("mWR",   t.c_str(),70,0,2800);
  t="M(W_{R}) e in Barrel"         +post; mWRBarrel=td->make<TH1D>("mWReBarrel", t.c_str(),70,0,2800);
  t="M(W_{R}) e and #mu in Barrel" +post;     mWRBB=td->make<TH1D>("mWRemuBarrel", t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{1} "       +post;     mNuR1=td->make<TH1D>("mNuR1", t.c_str(),70,0,2800);
  t="M(N_{R}) with #mu_{2} "       +post;     mNuR2=td->make<TH1D>("mNuR2", t.c_str(),70,0,1400);
  t="M(N_{R}) #mu_{1} vs. #mu_{2} "+post;    mNuR2D=td->make<TH2D>("mNuR2D",t.c_str(),70,0,2800,70,0,1400);

  t="M(jj) "                       +post;          mJJ=td->make<TH1D>("mJJ",         t.c_str(),50,0,2000);
  t="M(#mu #mu) "                  +post;        mMuMu=td->make<TH1D>("mMuMu",       t.c_str(),50,0,2000);
  t="M(#mu #mu) e in Barrel"       +post;  mMuMuBarrel=td->make<TH1D>("mMuMueBarrel", t.c_str(),50,0,2000);
  t="M(#mu #mu) e and #mu in Barrel" +post;    mMuMuBB=td->make<TH1D>("mMuMuemuBarrel", t.c_str(),50,0,2000);
  t="M(#mu #mu)(OS) "              +post;      mMuMuOS=td->make<TH1D>("mMuMuOS",     t.c_str(),50,0,2000);
  t="M(#mu #mu)(SS) "              +post;      mMuMuSS=td->make<TH1D>("mMuMuSS",     t.c_str(),50,0,2000);
  t="M(#mu #mu) "                  +post;    mMuMuZoom=td->make<TH1D>("mMuMuZoom",   t.c_str(),75,0,300);
  t="M(#mu #mu)(generated) "       +post; mMuMuGenZoom=td->make<TH1D>("mMuMuGenZoom",t.c_str(),75,0,300);
  t="DiMuon Charge "               +post;   diMuCharge=td->make<TH1D>("diMuCharge",  t.c_str(),2,-1,1);

  double edges[] = {0,20,30,40,50,60,1500};
  t="M(W_{R}) vs. min lepton p_{t} " +post;  mWRvsminLPt=td->make<TH2D>("mWRvsminLPt",  t.c_str(),70,0,2800,6,edges);
  t="M(W_{R}) vs. N primary vertices " +post;  mWRvsNPV=td->make<TH2D>("mWRvsNPV",  t.c_str(),70,0,2800,30, 0, 30);

  diMuCharge->GetXaxis()->SetBinLabel(1,"OS");
  diMuCharge->GetXaxis()->SetBinLabel(2,"SS");

  t="cZeta(mumu) "                +post; czeta_mumu=td->make<TH1D>("czMM",    t.c_str(),100,-1,1);
  t="cZeta(mumu) Zoom "      +post; czeta_mumu_zoom=td->make<TH1D>("czMMzoom",t.c_str(),100,-1,-0.9);

  // crazy angles
  t="cT(mumu) "   +post;     ctheta_mumu=td->make<TH1D>("ctMM",   t.c_str(),50,0,1);
  t="cT(jj) "     +post;       ctheta_jj=td->make<TH1D>("ctJJ",   t.c_str(),50,0,1);
  t="cT(mu1-jj) " +post;   ctheta_mu1_jj=td->make<TH1D>("ctM1JJ", t.c_str(),50,0,1);
  t="cT(mu2-jj) " +post;   ctheta_mu2_jj=td->make<TH1D>("ctM2JJ", t.c_str(),50,0,1);

  t="cTz(mumu) "  +post;    cthetaz_mumu=td->make<TH1D>("ctzMM",  t.c_str(),50,0,1);
  t="cTz(jj) "    +post;      cthetaz_jj=td->make<TH1D>("ctzJJ",  t.c_str(),50,0,1);
  t="cTz(mu1-jj) "+post;  cthetaz_mu1_jj=td->make<TH1D>("ctzM1JJ",t.c_str(),50,0,1);
  t="cTz(mu2-jj) "+post;  cthetaz_mu2_jj=td->make<TH1D>("ctzM2JJ",t.c_str(),50,0,1);

  // vertex histograms 
  t="mumu vtx DeltaZ "    +post; vtx_mumu=td->make<TH1D>("vtxMM", t.c_str(),100,0.,1.);
  t="jj vtx DeltaZ "      +post; vtx_jj=td->make<TH1D>("vtxJJ", t.c_str(),100,0.,1.);
  t="m1j vtx min DeltaZ " +post; vtx_mu1jmin=td->make<TH1D>("vtxM1Jmin", t.c_str(),100,0.,1.);
  t="m2j vtx min DeltaZ " +post; vtx_mu2jmin=td->make<TH1D>("vtxM2Jmin", t.c_str(),100,0.,1.);

  // ----------  Neural Net histograms  ----------

  if (v_masspts.size()) {
    nndir= new TFileDirectory(td->mkdir("_NNdata"));
    for (size_t i=0; i<v_masspts.size(); i++) {
      int mwr = v_masspts[i].first;
      int mnu = v_masspts[i].second;
      std::string name = nnhistoname(mwr,mnu);
      nndir->make<TH1D>(name.c_str(),(name+post).c_str(),51,-0.01,1.01);
    }
  }
}// end of book()

//======================================================================

void HeavyNuTop::HistPerDef::fill(pat::MuonCollection muons,
				  pat::ElectronCollection electrons,
				  pat::JetCollection  jets,
				  pat::METCollection  metc,
				  bool isMC,
				  double wgt)
{  
  std::sort(muons.begin(),muons.end(),comparePt()) ; 
  std::sort(electrons.begin(),electrons.end(),comparePt()) ; 
  std::sort(jets.begin(),jets.end(),comparePt()) ; 

  reco::Particle::LorentzVector vWR; 

  const pat::Muon& m0=muons.at(0);
  const pat::Electron& m1=electrons.at(0); //kluge

  // Muons 
  ptMu1->Fill(m0.pt(),wgt) ; 
  ptMu2->Fill(m1.pt(),wgt) ; 

  etaMu1->Fill(m0.eta(),wgt) ; 
  etaMu2->Fill(m1.eta(),wgt) ; 

  phiMu1->Fill(m0.phi(),wgt) ; 
  phiMu2->Fill(m1.phi(),wgt) ; 

  dPhiMu->Fill( fabs(deltaPhi(m0.phi(),m1.phi())),wgt ) ; 
  dEtaMu->Fill( fabs(m0.eta() - m1.eta()),wgt ) ; 
  dEtaPhiMu->Fill(fabs(m0.eta()-m1.eta()),
		  fabs(deltaPhi(m0.phi(),m1.phi())),wgt) ; 
  
  mu1trackIso   ->Fill(m0.trackIso(),wgt);
  mu1hcalIso    ->Fill(m0.hcalIso(),wgt);
  mu1ecalIso    ->Fill(m0.ecalIso(),wgt);
  mu1caloIso    ->Fill(m0.caloIso(),wgt);
  mu1dB         ->Fill(m0.dB(),wgt);

  mu2trackIso   ->Fill(m1.trackIso(),wgt);
  mu2hcalIso    ->Fill(m1.hcalIso(),wgt);
  mu2ecalIso    ->Fill(m1.ecalIso(),wgt);
  mu2caloIso    ->Fill(m1.caloIso(),wgt);
  mu2dB         ->Fill(m1.dB(),wgt);
  
  mu1trackRelIso->Fill(m0.trackIso()/m0.pt(),wgt);
  mu1hcalRelIso ->Fill(m0.hcalIso() /m0.pt(),wgt);
  mu1ecalRelIso ->Fill(m0.ecalIso() /m0.pt(),wgt);
  mu1caloRelIso ->Fill(m0.caloIso() /m0.pt(),wgt);

  mu2trackRelIso->Fill(m1.trackIso()/m1.pt(),wgt);
  mu2hcalRelIso ->Fill(m1.hcalIso() /m1.pt(),wgt);
  mu2ecalRelIso ->Fill(m1.ecalIso() /m1.pt(),wgt);
  mu2caloRelIso ->Fill(m1.caloIso() /m1.pt(),wgt);

  for (int i=0; i<muonQualityFlags; i++) { 
    if (m0.muonID(muonQuality[i])) qualMu1->Fill( i,wgt ) ; 
  }

  // Jets 
  const pat::Jet& j0 = jets.at(0);
  const pat::Jet& j1 = jets.at(1);

  ptJet1->Fill(j0.pt(),wgt) ; 
  ptJet2->Fill(j1.pt(),wgt) ; 

  etaJet1->Fill(j0.eta(),wgt) ; 
  etaJet2->Fill(j1.eta(),wgt) ; 

  phiJet1->Fill(j0.phi(),wgt) ; 
  phiJet2->Fill(j1.phi(),wgt) ; 

  double j0bdisc = j0.bDiscriminator(btagName);
  double j1bdisc = j1.bDiscriminator(btagName);

  btagJet1->Fill(j0bdisc,wgt);
  btagJet2->Fill(j1bdisc,wgt);

  if( (j0bdisc >= minBtagDiscVal) && (j1bdisc >= minBtagDiscVal)   ) numBjets->Fill(2.,wgt);
  else if( (j0bdisc >= minBtagDiscVal) || (j1bdisc >= minBtagDiscVal)   ) numBjets->Fill(1.,wgt);
  else numBjets->Fill(0.,wgt);

  dPhiJet->Fill( fabs(deltaPhi(j0.phi(),j1.phi())),wgt ) ; 
  dEtaJet->Fill( fabs(j0.eta() - j1.eta()),wgt ) ; 
  dEtaPhiJet->Fill( fabs(j0.eta()-j1.eta()),
		    fabs(deltaPhi(j0.phi(),j1.phi())),wgt ) ;

  jetID2d->Fill(hnu::jetID(j0),hnu::jetID(j1),wgt);

  // met
  if( metc.size() )
    met->Fill(metc.at(0).pt(),wgt);
  else
    met->Fill(0.,wgt);

  // Muon-Jet plots
  float dRmu1jet1 = deltaR(m0.eta(),m0.phi(),j0.eta(),j0.phi()) ; 
  float dRmu1jet2 = deltaR(m0.eta(),m0.phi(),j1.eta(),j1.phi()) ; 
  float dRmu2jet1 = deltaR(m1.eta(),m1.phi(),j0.eta(),j0.phi()) ; 
  float dRmu2jet2 = deltaR(m1.eta(),m1.phi(),j1.eta(),j1.phi()) ; 

  const pat::Jet&  j4mu1 = (dRmu1jet1 < dRmu1jet2) ? j0 : j1;
  const pat::Jet&  j4mu2 = (dRmu2jet1 < dRmu2jet2) ? j0 : j1;

  TVector3 mu1vec(m0.momentum().X(), m0.momentum().Y(), m0.momentum().Z());
  TVector3 mu2vec(m1.momentum().X(), m1.momentum().Y(), m1.momentum().Z());

  TVector3 jt1vec(j4mu1.p4().Vect().X(), j4mu1.p4().Vect().Y(), j4mu1.p4().Vect().Z() );
  TVector3 jt2vec(j4mu2.p4().Vect().X(), j4mu2.p4().Vect().Y(), j4mu2.p4().Vect().Z() );

  double ptrelMu1 = mu1vec.Perp(jt1vec);
  double ptrelMu2 = mu2vec.Perp(jt2vec);

  dRminMu1jet->Fill( min(dRmu1jet1,dRmu1jet2),wgt ) ; 
  dRminMu2jet->Fill( min(dRmu2jet1,dRmu2jet2),wgt ) ; 

  hptrelMu1->Fill( ptrelMu1,wgt );
  hptrelMu2->Fill( ptrelMu2,wgt );

  ptrelVsdRminMu1jet->Fill(min(dRmu1jet1,dRmu1jet2),ptrelMu1,wgt);
  ptrelVsdRminMu2jet->Fill(min(dRmu2jet1,dRmu1jet2),ptrelMu2,wgt);

  // Composite objects
  reco::Particle::LorentzVector vJJ = j0.p4() + j1.p4() ; 
  vWR = vJJ + m0.p4() + m1.p4() ; 

  mWR->Fill( vWR.M(),wgt ) ;
  if(m1.isEB()) 
  {
      mWRBarrel->Fill(vWR.M(),wgt );
      if(m0.eta() < 1.44) mWRBB->Fill(vWR.M(),wgt );
  }
  mNuR1 ->Fill( (vJJ + m0.p4()).M(),wgt ) ; 
  mNuR2 ->Fill( (vJJ + m1.p4()).M(),wgt ) ; 
  mNuR2D->Fill( (vJJ + m0.p4()).M(),(vJJ + m1.p4()).M(),wgt ) ;

  mWRvsminLPt->Fill(vWR.M(), std::min(m0.pt(), m1.pt()), wgt);

  reco::Particle::LorentzVector mumu=m0.p4()+m1.p4();
  reco::Particle::LorentzVector jj=j0.p4()+j1.p4();

  mMuMu->Fill(mumu.M(),wgt );
  if(m1.isEB())
  {
      mMuMuBarrel->Fill(mumu.M(),wgt);
      if(m0.eta() < 1.44) mMuMuBB->Fill(mumu.M(),wgt);
  }
  if (m0.charge() == m1.charge()) {
    mMuMuSS->Fill( mumu.M(),wgt );
    ptMu1VsPtMu2ss->Fill( m0.pt(), m1.pt(),wgt );
  }
  else {
    mMuMuOS->Fill( mumu.M(),wgt );
    ptMu1VsPtMu2os->Fill( m0.pt(), m1.pt(),wgt );
  }

  diMuCharge->Fill(0.5*m0.charge()*m1.charge(),wgt);

  mMuMuZoom->Fill(mumu.M(),wgt );
  mJJ->Fill(jj.M(),wgt );

}// end of fill()

//======================================================================

void
HeavyNuTop::HistPerDef::fill(const HeavyNuEvent& hne,
			  const std::vector<hNuMassHypothesis>& v_masspts)
{
  double wgt = hne.eventWgt ;

  // Muons 
  double mu1pt = hne.mu1.pt();
  double e1pt  = hne.ElecScale*hne.e1.pt();

  ptMu1->Fill(mu1pt,wgt) ; 
  ptMu2->Fill(e1pt,wgt) ; 

  etaMu1->Fill(hne.mu1.eta(),wgt) ; 
  etaMu2->Fill(hne.e1.eta(),wgt) ; 

  phiMu1->Fill(hne.mu1.phi(),wgt) ; 
  phiMu2->Fill(hne.e1.phi(),wgt) ; 

  dPhiMu->Fill( fabs(deltaPhi(hne.mu1.phi(),hne.e1.phi())),wgt ) ; 
  dEtaMu->Fill( fabs(hne.mu1.eta() - hne.e1.eta()),wgt ) ; 
  dEtaPhiMu->Fill(fabs(hne.mu1.eta()-hne.e1.eta()),
		  fabs(deltaPhi(hne.mu1.phi(),hne.e1.phi())),wgt) ; 

  mu1trackIso->Fill(hne.mu1.trackIso(),wgt);
  mu1hcalIso ->Fill(hne.mu1.hcalIso(),wgt);
  mu1ecalIso ->Fill(hne.mu1.ecalIso(),wgt);
  mu1caloIso ->Fill(hne.mu1.caloIso(),wgt);
  mu1dB      ->Fill(hne.mu1.dB(),wgt);
  mu2trackIso->Fill(hne.e1.trackIso(),wgt);
  mu2hcalIso ->Fill(hne.e1.hcalIso(),wgt);
  mu2ecalIso ->Fill(hne.e1.ecalIso(),wgt);
  mu2caloIso ->Fill(hne.e1.caloIso(),wgt);
  mu2dB      ->Fill(hne.e1.dB(),wgt);
  
  mu1trackRelIso->Fill(hne.mu1.trackIso()/mu1pt,wgt);
  mu1hcalRelIso ->Fill(hne.mu1.hcalIso() /mu1pt,wgt);
  mu1ecalRelIso ->Fill(hne.mu1.ecalIso() /mu1pt,wgt);
  mu1caloRelIso ->Fill(hne.mu1.caloIso() /mu1pt,wgt);
  mu2trackRelIso->Fill(hne.e1.trackIso()/e1pt,wgt);
  mu2hcalRelIso ->Fill(hne.e1.hcalIso() /e1pt,wgt);
  mu2ecalRelIso ->Fill(hne.e1.ecalIso() /e1pt,wgt);
  mu2caloRelIso ->Fill(hne.e1.caloIso() /e1pt,wgt);

  vtx_mumu->Fill(fabs(hne.mu1.vertex().Z()-hne.e1.vertex().Z()),wgt);

  for (int i=0; i<muonQualityFlags; i++) { 
    if (hne.mu1.muonID(muonQuality[i])) qualMu1->Fill( i,wgt ) ; 
  }

  int jet1id = 0;
  int jet2id = 0;

  // Jets 
  if (hne.nJets > 0) {
    jet1id = hnu::jetID(hne.j1);

    double j1bdisc = hne.j1.bDiscriminator(btagName);

    ptJet1->Fill(hne.j1.pt(),wgt) ; 
    etaJet1->Fill(hne.j1.eta(),wgt) ; 
    phiJet1->Fill(hne.j1.phi(),wgt) ; 
    btagJet1->Fill(j1bdisc,wgt);

    if (hne.nJets > 1) {
      double j2bdisc = hne.j2.bDiscriminator(btagName);

      ptJet2->Fill(hne.j2.pt(),wgt) ; 
      etaJet2->Fill(hne.j2.eta(),wgt) ; 
      phiJet2->Fill(hne.j2.phi(),wgt) ; 
      btagJet2->Fill(j2bdisc,wgt);

      if((j1bdisc >= minBtagDiscVal) && (j2bdisc >= minBtagDiscVal)) numBjets->Fill(2.,wgt);
      else if((j1bdisc >= minBtagDiscVal) || (j2bdisc >= minBtagDiscVal)) numBjets->Fill(1.,wgt);
      else numBjets->Fill(0.,wgt);

      jet2id = hnu::jetID(hne.j2);

      dPhiJet->Fill( fabs(deltaPhi(hne.j1.phi(),hne.j2.phi())),wgt ) ; 
      dEtaJet->Fill( fabs(hne.j1.eta() - hne.j2.eta()),wgt ) ; 
      dEtaPhiJet->Fill(fabs(hne.j1.eta()-hne.j2.eta()),
		       fabs(deltaPhi(hne.j1.phi(),hne.j2.phi())),wgt) ;

      mWR->Fill(hne.mWR,wgt); 
      if(hne.e1.isEB()) 
      {
          mWRBarrel->Fill(hne.mWR,wgt);
          if(hne.mu1.eta() < 1.44) mWRBB->Fill(hne.mWR, wgt);
      }
      mNuR1->Fill ( hne.mNuR1,wgt ) ; 
      mNuR2->Fill ( hne.mNuR2,wgt ) ; 
      mNuR2D->Fill( hne.mNuR1, hne.mNuR2,wgt );
      mJJ->Fill   ( hne.mJJ,wgt   );

      mWRvsminLPt->Fill(hne.mWR, std::min(hne.mu1.pt(), hne.e1.pt()), wgt);
      mWRvsNPV->Fill(hne.mWR, hne.n_primary_vertex, wgt);
      
      ctheta_jj->Fill(hne.ctheta_jj,wgt);
      ctheta_mu1_jj->Fill(hne.ctheta_mu1_jj,wgt);
      ctheta_mu2_jj->Fill(hne.ctheta_mu2_jj,wgt);
      cthetaz_jj->Fill(hne.cthetaz_jj,wgt);
      cthetaz_mu1_jj->Fill(hne.cthetaz_mu1_jj,wgt);
      cthetaz_mu2_jj->Fill(hne.cthetaz_mu2_jj,wgt);

      vtx_jj->Fill(fabs(hne.tjV1-hne.tjV2),wgt);
      vtx_mu1jmin->Fill(min(fabs(hne.mu1.vertex().Z()-hne.tjV1),
			    fabs(hne.mu1.vertex().Z()-hne.tjV2)),wgt);
      vtx_mu2jmin->Fill(min(fabs(hne.e1.vertex().Z()-hne.tjV1),
			    fabs(hne.e1.vertex().Z()-hne.tjV2)),wgt);
    }

    dRminMu1jet->Fill(hne.dRminMu1jet,wgt);
    dRminMu2jet->Fill(hne.dRminMu2jet,wgt);

    hptrelMu1->Fill(hne.ptrelMu1,wgt);
    hptrelMu2->Fill(hne.ptrelMu2,wgt);

    ptrelVsdRminMu1jet->Fill(hne.dRminMu1jet,hne.ptrelMu1,wgt);
    ptrelVsdRminMu2jet->Fill(hne.dRminMu2jet,hne.ptrelMu2,wgt);
  }

  jetID2d->Fill(jet1id,jet2id,wgt);

  // met
  met->Fill( hne.met1.pt(),wgt );

  mMuMu->Fill( hne.mMuMu,wgt );
  if(hne.e1.isEB())
  {
      mMuMuBarrel->Fill(hne.mMuMu,wgt );
      if(hne.mu1.eta() < 1.44) mMuMuBB->Fill(hne.mMuMu,wgt );
  }
  if (hne.mu1.charge() == hne.e1.charge()) {
    mMuMuSS->Fill( hne.mMuMu,wgt );
    ptMu1VsPtMu2ss->Fill( hne.mu1.pt(), hne.e1.pt(),wgt );
  }  else {
    mMuMuOS->Fill( hne.mMuMu,wgt );
    ptMu1VsPtMu2os->Fill( hne.mu1.pt(), hne.e1.pt(),wgt );
  }

  mMuMuZoom->Fill( hne.mMuMu,wgt );

  diMuCharge->Fill(0.5*hne.mu1.charge()*hne.e1.charge(),wgt);

  czeta_mumu->Fill(hne.czeta_mumu,wgt);
  czeta_mumu_zoom->Fill(hne.czeta_mumu,wgt);
  ctheta_mumu->Fill(hne.ctheta_mumu,wgt);
  cthetaz_mumu->Fill(hne.cthetaz_mumu,wgt);
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
HeavyNuTop::HeavyNuTop(const edm::ParameterSet& iConfig)
{
  // ==================== Get all parameters ====================
  //
  dolog_=iConfig.getParameter<bool>("DoLog");

  muonTag_ = iConfig.getParameter< edm::InputTag >( "muonTag" );
  jetTag_  = iConfig.getParameter< edm::InputTag >( "jetTag"  );
  metTag_  = iConfig.getParameter< edm::InputTag >( "metTag"  );
  elecTag_ = iConfig.getParameter< edm::InputTag >( "electronTag" );

  btagName=iConfig.getParameter<std::string>("BtagName");
  minBtagDiscVal=iConfig.getParameter<double>("minBtagDiscr");

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

  EBscalefactor_ = iConfig.getParameter<double>("EBscalefactor") ; 
  EEscalefactor_ = iConfig.getParameter<double>("EEscalefactor") ; 
  ebIDwgt_ = iConfig.getParameter<double>("EBidWgt") ; 
  eeIDwgt_ = iConfig.getParameter<double>("EEidWgt") ; 

  applyMESfactor_ = 1.0 ; // Hardcoded for Top studies
  applyTrigEffsign_ = 0. ; // Hardcoded for Top studies 

  isPFJets_ = iConfig.getParameter<bool>("isPFJets");

  // applyEleScaleCorrections_ = iConfig.getParameter<bool>("applyEleEScale") ; ;
  // applyEleIDWeightFactor_   = iConfig.getParameter<bool>("applyEleIDweight") ; 
  applyMuIDCorrections_ = iConfig.getParameter<bool>("applyMuIDEffcorr");


  studyScaleFactorEvolution_ = iConfig.getParameter<bool>("studyScaleFactor");

  // ==================== Init other members ====================
  //

  // nnif_ = new HeavyNu_NNIF(iConfig);
  trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet>("trigMatchPset"));
  muid_ = new HeavyNuID(iConfig.getParameter<edm::ParameterSet>("muIDPset"));

  // ==================== Book the histos ====================
  //
  edm::Service<TFileService> fs;
  hists.nelec    = fs->make<TH1D>("nelec",     "N(e^{#pm})",10,-0.5,9.5);
  hists.nmuAll   = fs->make<TH1D>("nmuAll",    "N(#mu^{#pm})",10,-0.5,9.5);
  hists.nmuLoose = fs->make<TH1D>("nmuLoose",  "N(#mu^{#pm}) passes Loose",10,-0.5,9.5);
  hists.nmuTight = fs->make<TH1D>("nmuTight",  "N(#mu^{#pm}) passes Tight",10,-0.5,9.5);
  hists.njet     = fs->make<TH1D>("njet",      "N(Jet)",50,-0.5,49.5);
  hists.nmet     = fs->make<TH1D>("nmet",      "N(MET)",50,-0.5,49.5);
  hists.met      = fs->make<TH1D>("met",       "MET distribution",100,0,2000);
  hists.muPt     = fs->make<TH1D>("muPt",      "#mu p_{T} distribution",100,0,2000);
  hists.muEta    = fs->make<TH1D>("muEta",     "#mu #eta distribution",50,-2.5,2.5);
  hists.muPhi    = fs->make<TH1D>("muPhi",     "#mu #phi distribution",60,-3.14159,3.14159);
  hists.jetPt    = fs->make<TH1D>("jetPt",     "jet p_{T} distribution",100,0,2000);
  hists.jetEta   = fs->make<TH1D>("jetEta",    "jet #eta distribution",50,-5,5);
  hists.jetPhi   = fs->make<TH1D>("jetPhi",    "jet #phi distribution",60,-3.14159,3.14159);
  hists.jetID    = fs->make<TH1I>("jetID",     "Jet ID",3,0,3);
  hists.jetPtvsNum=fs->make<TH2D>("jetPtvsNum","Jet P_{T} vs. Jet # ",11,-0.5,10.5,200,0.,2000.);

  hists.jetID->GetXaxis()->SetBinLabel(1,"Neither");
  hists.jetID->GetXaxis()->SetBinLabel(2,"PURE09 Loose");
  hists.jetID->GetXaxis()->SetBinLabel(3,"PURE09 Tight");

  hists.trkIsoStudy = fs->make<TH1D>("trkIsoStudy",";Tracker Relative Isolation", 100,0.,1.);


  // Histos per cut:
  //
  hists.noCuts.book      ( new TFileDirectory(fs->mkdir("cut0_none")),       "(no cuts)",                 v_null );
  hists.LLJJptCuts.book  ( new TFileDirectory(fs->mkdir("cut1_LLJJpt")),     "(4 objects with ptcuts:1)", v_null );
  hists.TrigMatches.book ( new TFileDirectory(fs->mkdir("cut2_TrigMatches")),"(Trigger match:M1)",        v_null );
  hists.VertexCuts.book  ( new TFileDirectory(fs->mkdir("cut3_Vertex")),     "(vertex requirements:M2)",  v_null );
  hists.Mu1HighPtCut.book( new TFileDirectory(fs->mkdir("cut4_Mu1HighPt")),  "(Mu1 High pt cut:M3)",      v_null );
  hists.diLmassCut.book  ( new TFileDirectory(fs->mkdir("cut5_diLmass")),    "(mumu mass cut:M4)",        v_null );
  hists.mWRmassCut.book  ( new TFileDirectory(fs->mkdir("cut6_mWRmass")),    "(mumujj mass cut:M5)",      v_null );

  if (studyScaleFactorEvolution_) { 
    hists.Mu1Pt30GeVCut.book ( new TFileDirectory(fs->mkdir("Mu1Pt30GeV")),  "(Mu1 30 GeV pt cut)",  v_null );
    hists.Mu1Pt40GeVCut.book ( new TFileDirectory(fs->mkdir("Mu1Pt40GeV")),  "(Mu1 40 GeV pt cut)",  v_null );
    hists.Mu1Pt50GeVCut.book ( new TFileDirectory(fs->mkdir("Mu1Pt50GeV")),  "(Mu1 50 GeV pt cut)",  v_null );
    hists.Mu1Pt60GeVCut.book ( new TFileDirectory(fs->mkdir("Mu1Pt60GeV")),  "(Mu1 60 GeV pt cut)",  v_null );
    hists.Mu1Pt80GeVCut.book ( new TFileDirectory(fs->mkdir("Mu1Pt80GeV")),  "(Mu1 80 GeV pt cut)",  v_null );
    hists.Mu1Pt100GeVCut.book( new TFileDirectory(fs->mkdir("Mu1Pt100GeV")), "(Mu1 100 GeV pt cut)", v_null );

    hists.Mu1HighPtCutVtxEq1.book( new TFileDirectory(fs->mkdir("Mu1HighPtVtxEq1")), "(Mu1 60 GeV pt cut, 1 vtx)", v_null );
    hists.Mu1HighPtCutVtx2to5.book( new TFileDirectory(fs->mkdir("Mu1HighPtVtx2to5")), "(Mu1 60 GeV pt cut, 2-5 vtx)", v_null );
    hists.Mu1HighPtCutVtxGt5.book( new TFileDirectory(fs->mkdir("Mu1HighPtVtxGt5")), "(Mu1 60 GeV pt cut, 6+ vtx)", v_null );
  }

  hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

  init_=false;

  MCweightByVertex_ = edm::LumiReWeighting(hnu::generate_flat10_mc(),
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
  // std::cout << "applyEleEScale    = " << applyEleScaleCorrections_  << std::endl ; 
  std::cout << "EB scale factor   = " << EBscalefactor_             << std::endl ; 
  std::cout << "EE scale factor   = " << EEscalefactor_             << std::endl ; 
  // std::cout << "applyEleIDweight  = " << applyEleIDWeightFactor_    << std::endl ; 
  std::cout << "EB weight         = " << ebIDwgt_                   << std::endl ; 
  std::cout << "EE weight         = " << eeIDwgt_                   << std::endl ; 

  std::cout << "pileup era        = " << pileupEra_ << std::endl;
  std::cout << "studyScaleFactor  = " << studyScaleFactorEvolution_ << std::endl;
}
  



HeavyNuTop::~HeavyNuTop()
{
  
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}

TH1 *
HeavyNuTop::bookRunHisto(uint32_t runNumber)
{
  std::string runstr = int2str<uint32_t>(runNumber);
  return hists.rundir->make <TH1I> (runstr.c_str(), runstr.c_str(),1,1,2);
}

// ------------ method called to for each event  ------------
bool
HeavyNuTop::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  HeavyNuEvent hnuEvent(HeavyNuEvent::TOP);

  hnuEvent.isMC = !iEvent.isRealData();
  hnuEvent.pfJets = isPFJets_;

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

  edm::Handle<reco::JPTJetCollection> jptJets;
  iEvent.getByLabel("JetPlusTrackZSPCorJetAntiKt5", jptJets); // Some day we should make this readable from the cfg file

  edm::Handle<pat::MuonCollection> pMuons ; 
  iEvent.getByLabel(muonTag_,pMuons) ; 

  edm::Handle<pat::ElectronCollection> pElecs ;
  iEvent.getByLabel(elecTag_, pElecs) ;

  edm::Handle<pat::JetCollection> pJets ;
  iEvent.getByLabel(jetTag_, pJets) ;

  edm::Handle<pat::METCollection> pMET ;
  iEvent.getByLabel(metTag_, pMET) ;

  edm::Handle<reco::MuonCollection> tevMuons;
  iEvent.getByLabel("refitMuons", tevMuons);

  //Shirpa reweighting info
  edm::Handle<GenEventInfoProduct> geneventinfo;
  iEvent.getByLabel("generator", geneventinfo);

  if (hnuEvent.isMC) { 
    edm::Handle<std::vector<PileupSummaryInfo> > pPU;
    iEvent.getByLabel("addPileupInfo", pPU); 
    std::pair<float,double> pileup = hnu::pileupReweighting(pPU,MCweightByVertex_) ; 
    hnuEvent.n_pue     = int(pileup.first) ; 
    hnuEvent.eventWgt *= pileup.second ;
    //Shirpa reweighting
    hnuEvent.eventWgt *= geneventinfo->weight();
  }
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
  hnuEvent.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle) ;

  if ( !pElecs.isValid() || 
       !pMuons.isValid() || 
       !pJets.isValid()  ||
       !(pMET.isValid() && pMET->size() > 0)) {
    std::cout << "Exiting as valid PAT objects not found" << std::endl ;
    return false; 
  }

  if (firstEvent_) {      
    firstEvent_ = false;
    // Some sanity checks inserted 
    // If running on Monte Carlo, expect
    //   - pileup configuration to make sense
    //   - to apply nominal trigger and ID corrections
    //   - to apply proper weighting to events for electron reco/ID
    // If running on data
    //   - no corrections are applied
    if ( studyScaleFactorEvolution_ ) std::cout << "Histograms for top scale factor cross checks will be created" << std::endl ; 
    if ( hnuEvent.isMC ) {
        if ( applyMuIDCorrections_ )      std::cout << "Nominal Mu ID corrections applied" << std::endl ;
        else                              std::cout << "WARNING: You have disabled Mu ID corrections.  Are you sure?" << std::endl ;
        if ( ebIDwgt_ != 1.0 &&
             eeIDwgt_ != 1.0 )            std::cout << "Events will be weighted to account for data/MC reconstruction and ID differences" << std::endl ;
        else                              std::cout << "WARNING: You are not applying EB/EE weighting to MC.  Are you sure?" << std::endl ; 

        int pileupYear = pileupEra_ ;
        int idYear     = muid_->idEra() ;
        
        bool allErasMatch = ( pileupYear == idYear ) ;
        if ( !allErasMatch ) {
            std::cout << "WARNING: You do not appear to have consistent corrections applied!" << std::endl ;
            std::cout << "         pileup year is " << pileupEra_ << ", year for mu ID is " << idYear
		      << std::endl ; 
	} else {
            std::cout << "Looking at corrections, I assume you are running with the " << pileupYear << " year settings" << std::endl ; 
        }
    } else {
        if ( ebIDwgt_ != 1.0 && eeIDwgt_ != 1.0 )
            std::cout << "WARNING: You are applying EB/EE weighting to DATA.  Are you sure?" << std::endl ; 
    }
  }

  hists.nelec ->Fill(pElecs->size()) ;
  hists.nmuAll->Fill(pMuons->size()) ;
  hists.njet  ->Fill(pJets->size()) ;
  hists.nmet  ->Fill(pMET->size()) ;

  if (pMET->size())
    hists.met->Fill(pMET->at(0).pt());
  else
    hists.met->Fill(0);

  // Basic selection requirements: Require at least two muons, two jets
  if ( pMuons->size() >= 1 && pElecs->size() >= 1 && pJets->size() >= 2 ) {
    hists.noCuts.fill( *pMuons,*pElecs,*pJets,*pMET, hnuEvent.isMC, hnuEvent.eventWgt ) ; 
  } else return false;

  if ( dolog_ ) std::cout << "Found an event with " << pMuons->size() << " muons and " 
			  << pElecs->size() << " electrons, and " << pJets->size() << " jets" << std::endl ; 

  // Look for valid electrons and put them in the event
  std::vector< std::pair<pat::Electron,float> > eCands = 
    hnu::getElectronList(pElecs,cuts.maximum_elec_abseta, 
			 cuts.minimum_lep2_pt,cuts.minimum_lep2_pt) ; 
  hnuEvent.nElectrons = eCands.size() ; 
  if ( hnuEvent.nElectrons < 1 ) return false ; 
  hnuEvent.e1        = eCands.at(0).first ; 
  hnuEvent.ElecScale = eCands.at(0).second ; 

  // Look for valid jets and put them in the event
  std::vector< std::pair<pat::Jet,float> > jetCands = 
    hnu::getJetList(pJets,jecuObj_,cuts.minimum_jet_pt,cuts.maximum_jet_abseta,0) ; 
  for (unsigned int i=0; i<jetCands.size(); i++) { 
    if ( hnuEvent.nJets == 2 ) break ; 
    pat::Jet iJ = jetCands.at(i).first ; 
    double dRej = deltaR(iJ.eta(), iJ.phi(), hnuEvent.e1.eta(), hnuEvent.e1.phi()) ; 
    if (dRej > cuts.minimum_lep_jet_dR) { 
      hnuEvent.nJets++ ; 
      if ( hnuEvent.nJets == 1 ) {
        hnuEvent.j1 = iJ ;
        hnuEvent.j1scale = 1.0 ; // No jet corrections are applied for Top!
      } else if ( hnuEvent.nJets == 2 ) {
        hnuEvent.j2 = iJ ;
        hnuEvent.j2scale = 1.0 ; 
      } else
        std::cout << "WARNING: Expected empty jet position" << std::endl ; 
    }
  }
  if ( hnuEvent.nJets < 2 ) return false ; 
  hnuEvent.tjV1 = hnu::caloJetVertex(hnuEvent.j1, *jptJets);
  hnuEvent.tjV2 = hnu::caloJetVertex(hnuEvent.j2, *jptJets);

  // Finally, look for valid muons and put them in the event
  std::vector<pat::Muon> muCands = 
    hnu::getMuonList(pMuons,tevMuons,cuts.minimum_lep2_pt,cuts.maximum_mu_abseta,applyMESfactor_) ; 
  for (unsigned int i=0; i<muCands.size(); i++) { 
    if ( hnuEvent.nMuons > 0 ) break ; 
    pat::Muon iM = muCands.at(i) ; 
    if ( hnu::muIsolation(iM) < cuts.muon_trackiso_limit ) {
      double dRj1 = deltaR(iM.eta(), iM.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi()) ; 
      double dRj2 = deltaR(iM.eta(), iM.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi()) ; 
      if (dRj1 > cuts.minimum_lep_jet_dR && dRj2 > cuts.minimum_lep_jet_dR) { 
	hnuEvent.nMuons++ ; 
	if      ( hnuEvent.nMuons == 1 ) hnuEvent.mu1 = iM ; 
	else    std::cout << "WARNING: Expected empty muon position" << std::endl ; 
      }
    }
  }
  if ( hnuEvent.nMuons < 1 ) return false ; // No muons
  if (applyMuIDCorrections_ && hnuEvent.isMC) {
    double mu1wgt = muid_->weightForMC((hnuEvent.mu1.pt()),0) ;
    hnuEvent.eventWgt *= mu1wgt ; 
  }

  hnuEvent.met1 = pMET->at(0);

  hnuEvent.regularize(); // assign internal standards
  //hnuEvent.scaleMuE(applyMESfactor_,hnuEvent.ElecScale);
  hnuEvent.calculate(); // calculate various details
  
  // Require mu1 meets tight requirements
  bool mu1isTight = hnu::isVBTFtight(hnuEvent.mu1);
  if ( !mu1isTight ) return false;

  // Basic requirements on muon, electron, jets
  hists.LLJJptCuts.fill( hnuEvent,v_null );

  // require that one muon be BOTH tight and trigger-matched
  bool mu1trig=mu1isTight;
  if ( trig_->matchingEnabled() && iEvent.isRealData() ) {
    mu1trig = mu1trig &&
      trig_->isTriggerMatched( hnuEvent.mu1, iEvent ) ; 
  } else if (!iEvent.isRealData()) {
    mu1trig = mu1trig &&  
      trig_->simulateForMC( hnuEvent.mu1.pt(),hnuEvent.mu1.eta(),0 );
  }

  if( !mu1trig ) return false;
  hists.TrigMatches.fill( hnuEvent,v_null );

  //--- Impose vertex requirement here ---//
  float deltaVzJ1J2 = fabs(hnuEvent.tjV1-hnuEvent.tjV2);
  float deltaVzJ1M1 = fabs(hnuEvent.tjV1-hnuEvent.mu1.vertex().Z());
  float deltaVzJ2E1 = fabs(hnuEvent.tjV2-hnuEvent.e1.vertex().Z());
  float deltaVzJ1E1 = fabs(hnuEvent.tjV1-hnuEvent.e1.vertex().Z());
  float deltaVzJ2M1 = fabs(hnuEvent.tjV2-hnuEvent.mu1.vertex().Z());
  float deltaVzM1E1 = fabs(hnuEvent.mu1.vertex().Z()-hnuEvent.e1.vertex().Z());
  if ((cuts.maxJetVZsepCM > 0 && cuts.maxVertexZsep > 0) && 
      ((deltaVzJ1J2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M1 >= cuts.maxJetVZsepCM) ||
       (deltaVzJ2E1 >= cuts.maxJetVZsepCM) || (deltaVzJ1E1 >= cuts.maxJetVZsepCM) ||
       (deltaVzJ2M1 >= cuts.maxJetVZsepCM) || (deltaVzM1E1 >= cuts.maxVertexZsep)) ) 
    return false ; 
  hists.VertexCuts.fill( hnuEvent,v_null );

  double mu1pt = hnuEvent.mu1.pt() ;
  double e1pt  = hnu::getElectronEt(hnuEvent.e1) ; 
  double highestPt = std::max(mu1pt, e1pt);
  if ( studyScaleFactorEvolution_ && hnuEvent.mMuMu>=cuts.minimum_mumu_mass) {
    if ( highestPt > 30. )  hists.Mu1Pt30GeVCut.fill( hnuEvent,v_null );
    if ( highestPt > 40. )  hists.Mu1Pt40GeVCut.fill( hnuEvent,v_null );
    if ( highestPt > 50. )  hists.Mu1Pt50GeVCut.fill( hnuEvent,v_null );
    if ( highestPt > 60. )  hists.Mu1Pt60GeVCut.fill( hnuEvent,v_null );
    if ( highestPt > 80. )  hists.Mu1Pt80GeVCut.fill( hnuEvent,v_null );
    if ( highestPt > 100. ) hists.Mu1Pt100GeVCut.fill( hnuEvent,v_null );
  }

  if( highestPt < cuts.minimum_lep1_pt )
    return false;

  if ( studyScaleFactorEvolution_ ) { 
    if ( hnuEvent.n_primary_vertex == 1 ) 
      hists.Mu1HighPtCutVtxEq1.fill( hnuEvent,v_null );
    else if ( hnuEvent.n_primary_vertex <= 5 ) 
      hists.Mu1HighPtCutVtx2to5.fill( hnuEvent,v_null );
    else if ( hnuEvent.n_primary_vertex > 5 ) 
      hists.Mu1HighPtCutVtxGt5.fill( hnuEvent,v_null );
  }
  hists.Mu1HighPtCut.fill( hnuEvent,v_null );

  if ( hnuEvent.mMuMu>=cuts.minimum_mumu_mass ) 
    hists.diLmassCut.fill( hnuEvent,v_null );

  if ( iEvent.isRealData() ) {
    std::cout<<"\t"<<iEvent.id() << std::endl;
    std::cout<<"\tM(W_R)  = "<<hnuEvent.mWR  <<" GeV";
    std::cout<<", M(NuR1) = "<<hnuEvent.mNuR1<<" GeV";
    std::cout<<", M(NuR2) = "<<hnuEvent.mNuR2<<" GeV"<<std::endl;
    std::cout<<"\tM(mumu) = "<<hnuEvent.mMuMu<<" GeV"; 
    std::cout<<", M(JJ) = "  <<hnuEvent.mJJ<<" GeV"<<std::endl; 
    std::cout<<"\tJets:   j1 "
	     <<"pt="<<hnuEvent.j1.pt()<<" GeV, eta="<<hnuEvent.j1.eta()<<", phi="<<hnuEvent.j1.phi();
    std::cout<<        ", j2 " 
	     <<"pt="<<hnuEvent.j2.pt()<<" GeV, eta="<<hnuEvent.j2.eta()<<", phi="<<hnuEvent.j2.phi();
    std::cout<<std::endl;
    std::cout<<"\tLeptons: mu1 "
	     <<"pt="<<hnuEvent.mu1.pt()<<" GeV, eta="<<hnuEvent.mu1.eta()<<", phi="<<hnuEvent.mu1.phi();
    std::cout<<          ", e1 " 
	     <<"pt="<<hnuEvent.e1.pt()<<" GeV, eta="<<hnuEvent.e1.eta()<<", phi="<<hnuEvent.e1.phi();
    std::cout<<std::endl;
  }

  // Change the final logic of the filter for e/mu studies:
  // Keep anything passing the 60 GeV lepton 1 cut
  if ( hnuEvent.mMuMu>=cuts.minimum_mumu_mass && hnuEvent.mWR>=cuts.minimum_mWR_mass )
    hists.mWRmassCut.fill( hnuEvent,v_null );

  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNuTop::beginJob() {
  // nnif_->beginJob();
  firstEvent_ = true;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNuTop::endJob() {
  // nnif_->endJob();
  trig_->endJob();
  muid_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuTop);
