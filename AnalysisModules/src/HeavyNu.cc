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
// $Id: HeavyNu.cc,v 1.10 2010/10/13 16:50:59 dudero Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNu_NNIF.h"


//////////////////////////////////////////////////////////////////
// generic maximum/minimum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     
}
template <class T> const T& min ( const T& a, const T& b ) {
  return (b<a)?b:a;     
}

class compare {
public:
  template <class T> bool operator() (const T& a, const T& b) { return a.pt() > b.pt() ; } 
};

class HeavyNu : public edm::EDFilter {
public:
  explicit HeavyNu(const edm::ParameterSet&);
  ~HeavyNu();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }
  
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual bool isVBTFloose(const pat::Muon& m);
  virtual bool isVBTFtight(const pat::Muon& m);

  std::string currentFile_;
  bool dolog_;
  HeavyNu_NNIF *nnif_;

  // ----------member data ---------------------------
  static const std::string muonQuality[] ; 
  static const int muonQualityFlags ; 

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory td, const std::string&) ;
    // fill all histos of the set with the two electron candidates
    void fill(pat::MuonCollection muons, pat::JetCollection jets,bool isMC) ;
    // fill all histos of the set with the two electron candidates
    void fill(const HeavyNuEvent& hne) ;

    TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2 ;  
    TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2 ;  
    TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2 ;  
    TH1 *dEtaMu, *dPhiMu, *dEtaJet, *dPhiJet ;
    TH1 *dEtaPhiMu, *dEtaPhiJet ; 
    TH1 *dRminMu1jet, *dRminMu2jet ; 
    TH1 *hptrelMu1, *hptrelMu2 ; 
    TH2 *ptrelVsdRminMu1jet, *ptrelVsdRminMu2jet ;

    TH1 *dptMu1gen, *dptMu2gen ; 
    TH1 *dRMu1gen, *dRMu2gen ; 
    TH1 *qualMu1, *qualMu2 ; 

    TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1caloIso, *mu1dB;
    TH1 *mu2trackIso, *mu2hcalIso, *mu2ecalIso, *mu2caloIso, *mu2dB;

    TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1caloRelIso;
    TH1 *mu2trackRelIso, *mu2hcalRelIso, *mu2ecalRelIso, *mu2caloRelIso;

    TH1 *mWR, *mNuR1, *mNuR2, *mMuMu, *mMuMuZoom, *mJJ ; 
    TH2 *mNuR2D, *jetPtvsNum ; 

    // Jeremy's crazy angles...
    TH1* ctheta_mumu, *cthetaz_mumu;
    TH1* ctheta_jj, *cthetaz_jj;
    TH1* ctheta_mu1_jj, *cthetaz_mu1_jj;
    TH1* ctheta_mu2_jj, *cthetaz_mu2_jj;

  };

  bool init_;

  // gf set of histo for all Z definitions in a stack
  struct HistStruct {
    TH1 *nelec, *nmu, *njet ;
    TH1 *muPt, *muEta, *muPhi ; 
    TH1 *jetPt, *jetEta, *jetPhi ; 
    TH2 *jetPtvsNum;
    HistPerDef noCuts ; 
    HistPerDef LLptCuts;
    HistPerDef LLJJptCuts;
    HistPerDef massCut;
  } hists;

  struct CutsStruct {
    double minimum_mu1_pt;
    double minimum_mu2_pt;
    double minimum_jet_pt;
    double maximum_jet_abseta;
    double minimum_mumu_mass;
    double minimum_muon_jet_dR;
  } cuts;
  
};

const std::string HeavyNu::muonQuality[] = {"AllGlobalMuons","AllStandAloneMuons","AllTrackerMuons"} ; 
const int HeavyNu::muonQualityFlags = 3 ;

void HeavyNu::HistPerDef::book(TFileDirectory td, const std::string& post) {
  std::string title;

  // Muon histograms 
  title=std::string("p_{T}(#mu_{1}) ")+post;
  ptMu1=td.make<TH1D>("ptMu1",title.c_str(),50,0.,1000.);  
  title=std::string("p_{T}(#mu_{2}) ")+post;
  ptMu2=td.make<TH1D>("ptMu2",title.c_str(),50,0.,500.);  
  title=std::string("#eta(#mu_{1}) ")+post;
  etaMu1=td.make<TH1D>("etaMu1",title.c_str(),40,-2.5,2.5);  
  title=std::string("#eta(#mu_{2}) ")+post;
  etaMu2=td.make<TH1D>("etaMu2",title.c_str(),40,-2.5,2.5);  
  title=std::string("#phi(#mu_{1}) ")+post;
  phiMu1=td.make<TH1D>("phiMu1",title.c_str(),30,-3.14159,3.14159);  
  title=std::string("#phi(#mu_{2}) ")+post;
  phiMu2=td.make<TH1D>("phiMu2",title.c_str(),30,-3.14159,3.14159);  
  title=std::string("#Delta#eta(#mu_{1},#mu_{2}) ")+post;
  dEtaMu=td.make<TH1D>("dEtaMu",title.c_str(),40,0,5);  
  title=std::string("#Delta#phi(#mu_{1},#mu_{2}) ")+post;
  dPhiMu=td.make<TH1D>("dPhiMu",title.c_str(),30,0,3.14159);  
  title=std::string("#mu #Delta#eta vs. #Delta#phi ")+post;
  dEtaPhiMu=td.make<TH2D>("dEtaPhiMu",title.c_str(),50,0,5.,30,0,3.14159);  

  title=std::string("#Delta p_{T}(#mu_{1},gen) ")+post;
  dptMu1gen=td.make<TH1D>("dptMu1gen",title.c_str(),50,-0.50,0.50);  
  title=std::string("#Delta p_{T}(#mu_{2},gen) ")+post;
  dptMu2gen=td.make<TH1D>("dptMu2gen",title.c_str(),50,-0.50,0.50);  
  title=std::string("#Delta R(#mu_{1},gen) ")+post;
  dRMu1gen=td.make<TH1D>("dRMu1gen",title.c_str(),50,0,0.01);  
  title=std::string("#Delta R(#mu_{2},gen) ")+post;
  dRMu2gen=td.make<TH1D>("dRMu2gen",title.c_str(),50,0,0.01);  

  title=std::string("Quality (#mu_{1}) ")+post;
  qualMu1=td.make<TH1D>("qualMu1",title.c_str(),muonQualityFlags,0,muonQualityFlags) ;
  title=std::string("Quality (#mu_{2}) ")+post;
  qualMu2=td.make<TH1D>("qualMu2",title.c_str(),muonQualityFlags,0,muonQualityFlags) ;
  for (int i=0; i<muonQualityFlags; i++) {
    qualMu1->GetXaxis()->SetBinLabel(i+1,muonQuality[i].c_str()) ;
    qualMu2->GetXaxis()->SetBinLabel(i+1,muonQuality[i].c_str()) ;
  }

  title=std::string("trackIso(#mu_{1}) ")+post;
  mu1trackIso=td.make<TH1D>("mu1trackIso",title.c_str(),40,0.,200.);  
  title=std::string("hcalIso(#mu_{1}) ")+post;
  mu1hcalIso=td.make<TH1D>("mu1hcalIso",title.c_str(),40,0.,200.);  
  title=std::string("ecalIso(#mu_{1}) ")+post;
  mu1ecalIso=td.make<TH1D>("mu1ecalIso",title.c_str(),40,0.,200.);  
  title=std::string("caloIso(#mu_{1}) ")+post;
  mu1caloIso=td.make<TH1D>("mu1caloIso",title.c_str(),40,0.,200.);  
  title=std::string("Dxy(#mu_{1}) ")+post;
  mu1dB=td.make<TH1D>("mu1dB",title.c_str(),50,-5.,5.);  

  title=std::string("trackIso(#mu_{2}) ")+post;
  mu2trackIso=td.make<TH1D>("mu2trackIso",title.c_str(),40,0.,200.);  
  title=std::string("hcalIso(#mu_{2}) ")+post;
  mu2hcalIso=td.make<TH1D>("mu2hcalIso",title.c_str(),40,0.,200.);  
  title=std::string("ecalIso(#mu_{2}) ")+post;
  mu2ecalIso=td.make<TH1D>("mu2ecalIso",title.c_str(),40,0.,200.);  
  title=std::string("caloIso(#mu_{2}) ")+post;
  mu2caloIso=td.make<TH1D>("mu2caloIso",title.c_str(),40,0.,200.);  
  title=std::string("Dxy(#mu_{2}) ")+post;
  mu2dB=td.make<TH1D>("mu2dB",title.c_str(),50,-5.,5.);  

  title=std::string("trackRelIso(#mu_{1}) ")+post;
  mu1trackRelIso=td.make<TH1D>("mu1trackRelIso",title.c_str(),50,0.,5.);  
  title=std::string("hcalRelIso(#mu_{1}) ")+post;
  mu1hcalRelIso=td.make<TH1D>("mu1hcalRelIso",title.c_str(),50,0.,5.);  
  title=std::string("ecalRelIso(#mu_{1}) ")+post;
  mu1ecalRelIso=td.make<TH1D>("mu1ecalRelIso",title.c_str(),50,0.,5.);  
  title=std::string("caloRelIso(#mu_{1}) ")+post;
  mu1caloRelIso=td.make<TH1D>("mu1caloRelIso",title.c_str(),50,0.,5.);  
  title=std::string("Dxy(#mu_{1}) ")+post;

  title=std::string("trackRelIso(#mu_{2}) ")+post;
  mu2trackRelIso=td.make<TH1D>("mu2trackRelIso",title.c_str(),50,0.,5.);  
  title=std::string("hcalRelIso(#mu_{2}) ")+post;
  mu2hcalRelIso=td.make<TH1D>("mu2hcalRelIso",title.c_str(),50,0.,5.);  
  title=std::string("ecalRelIso(#mu_{2}) ")+post;
  mu2ecalRelIso=td.make<TH1D>("mu2ecalRelIso",title.c_str(),50,0.,5.);  
  title=std::string("caloRelIso(#mu_{2}) ")+post;
  mu2caloRelIso=td.make<TH1D>("mu2caloRelIso",title.c_str(),50,0.,5.);  
  title=std::string("Dxy(#mu_{2}) ")+post;

  // Jet histograms 
  title=std::string("p_{T}(j_{1}) ")+post;
  ptJet1=td.make<TH1D>("ptJet1",title.c_str(),50,0.,500.);  
  title=std::string("p_{T}(j_{2}) ")+post;
  ptJet2=td.make<TH1D>("ptJet2",title.c_str(),50,0.,500.);  
  title=std::string("#eta(j_{1}) ")+post;
  etaJet1=td.make<TH1D>("etaJet1",title.c_str(),40,-5,5);  
  title=std::string("#eta(j_{2}) ")+post;
  etaJet2=td.make<TH1D>("etaJet2",title.c_str(),40,-5,5);  
  title=std::string("#phi(j_{1}) ")+post;
  phiJet1=td.make<TH1D>("phiJet1",title.c_str(),30,-3.14159,3.14159);  
  title=std::string("#phi(j_{2}) ")+post;
  phiJet2=td.make<TH1D>("phiJet2",title.c_str(),30,-3.14159,3.14159);  
  title=std::string("#Delta#eta(j_{1},j_{2}) ")+post;
  dEtaJet=td.make<TH1D>("dEtaJet",title.c_str(),40,0,5);  
  title=std::string("#Delta#phi(j_{1},j_{2}) ")+post;
  dPhiJet=td.make<TH1D>("dPhiJet",title.c_str(),30,0,3.14159);  
  title=std::string("Jet #Delta#eta vs. #Delta#phi ")+post;
  dEtaPhiJet=td.make<TH2D>("dEtaPhiJet",title.c_str(),50,0,5,30,0,3.14159);  

  // Mu/Jet histograms
  title=std::string("Minimum #Delta R(#mu_{1},jet) ")+post;
  dRminMu1jet=td.make<TH1D>("dRminMu1jet",title.c_str(),50,0,5.);  
  title=std::string("Minimum #Delta R(#mu_{2},jet) ")+post;
  dRminMu2jet=td.make<TH1D>("dRminMu2jet",title.c_str(),50,0,5.);  

  title=std::string("p_{T,rel}(#mu_{1},jet)")+post;
  hptrelMu1=td.make<TH1D>("ptrelMu1",title.c_str(),50,0,1000.);
  title=std::string("p_{T,rel}(#mu_{2},jet)")+post;
  hptrelMu2=td.make<TH1D>("ptrelMu2",title.c_str(),50,0,1000.);

  title=std::string("p_{T,rel}(#mu_{1},jet) vs #Delta R(#mu_{1},jet)")+post;
  title+="; #Delta R(#mu_{1},jet); p_{T,rel}(#mu_{1},jet)";
  ptrelVsdRminMu1jet = td.make<TH2D>("ptrelVsdRminMu1jet",title.c_str(),
				     50,0,5.,50,0,1000);
  title=std::string("p_{T,rel}(#mu_{2},jet) vs #Delta R(#mu_{2},jet)")+post;
  title+="; #Delta R(#mu_{2},jet); p_{T,rel}(#mu_{2},jet)";
  ptrelVsdRminMu2jet = td.make<TH2D>("ptrelVsdRminMu2jet",title.c_str(),
				     50,0,5.,50,0,1000);
  
  // Composite histograms 
  title=std::string("M(W_{R}) ")+post;
  mWR=td.make<TH1D>("mWR",title.c_str(),50,0,2000);  
  title=std::string("M(N_{R}) with #mu_{1} ")+post;
  mNuR1=td.make<TH1D>("mNuR1",title.c_str(),50,0,2000);  
  title=std::string("M(N_{R}) with #mu_{2} ")+post;
  mNuR2=td.make<TH1D>("mNuR2",title.c_str(),50,0,1000);  
  title=std::string("M(N_{R}) #mu_{1} vs. #mu_{2} ")+post;
  mNuR2D=td.make<TH2D>("mNuR2D",title.c_str(),50,0,2000,50,0,1000);  

  title=std::string("M(#mu #mu)")+post;
  mMuMu=td.make<TH1D>("mMuMu",title.c_str(),50,0,2000);  
  title=std::string("M(#mu #mu)")+post;
  mMuMuZoom=td.make<TH1D>("mMuMuZoom",title.c_str(),50,0,200);  
  title=std::string("M(jj)")+post;
  mJJ=td.make<TH1D>("mJJ",title.c_str(),50,0,2000);  

  // crazy angles
  title=std::string("cT(mumu)")+post;
  ctheta_mumu=td.make<TH1D>("ctMM",title.c_str(),50,0,1);  
  title=std::string("cT(jj)")+post;
  ctheta_jj=td.make<TH1D>("ctJJ",title.c_str(),50,0,1);  
  title=std::string("cT(mu1-jj)")+post;
  ctheta_mu1_jj=td.make<TH1D>("ctM1JJ",title.c_str(),50,0,1);  
  title=std::string("cT(mu2-jj)")+post;
  ctheta_mu2_jj=td.make<TH1D>("ctM2JJ",title.c_str(),50,0,1);  
  title=std::string("cTz(mumu)")+post;
  cthetaz_mumu=td.make<TH1D>("ctzMM",title.c_str(),50,0,1);  
  title=std::string("cTz(jj)")+post;
  cthetaz_jj=td.make<TH1D>("ctzJJ",title.c_str(),50,0,1);  
  title=std::string("cTz(mu1-jj)")+post;
  cthetaz_mu1_jj=td.make<TH1D>("ctzM1JJ",title.c_str(),50,0,1);  
  title=std::string("cTz(mu2-jj)")+post;
  cthetaz_mu2_jj=td.make<TH1D>("ctzM2JJ",title.c_str(),50,0,1);  
  
}

void HeavyNu::HistPerDef::fill(pat::MuonCollection muons,
			       pat::JetCollection  jets,
			       bool isMC)
{  
  std::sort(muons.begin(),muons.end(),compare()) ; 
  std::sort(jets.begin(),jets.end(),compare()) ; 

  reco::Particle::LorentzVector WR ; 

  // Muons 
  ptMu1->Fill(muons.at(0).pt()) ; 
  ptMu2->Fill(muons.at(1).pt()) ; 

  etaMu1->Fill(muons.at(0).eta()) ; 
  etaMu2->Fill(muons.at(1).eta()) ; 

  phiMu1->Fill(muons.at(0).phi()) ; 
  phiMu2->Fill(muons.at(1).phi()) ; 

  dPhiMu->Fill( fabs(deltaPhi(muons.at(0).phi(),muons.at(1).phi())) ) ; 
  dEtaMu->Fill( fabs(muons.at(0).eta() - muons.at(1).eta()) ) ; 
  dEtaPhiMu->Fill(fabs(muons.at(0).eta()-muons.at(1).eta()),
		  fabs(deltaPhi(muons.at(0).phi(),muons.at(1).phi()))) ; 
  
  mu1trackIso->Fill(muons.at(0).trackIso());
  mu1hcalIso ->Fill(muons.at(0).hcalIso());
  mu1ecalIso ->Fill(muons.at(0).ecalIso());
  mu1caloIso ->Fill(muons.at(0).caloIso());
  mu1dB      ->Fill(muons.at(0).dB());
  mu2trackIso->Fill(muons.at(1).trackIso());
  mu2hcalIso ->Fill(muons.at(1).hcalIso());
  mu2ecalIso ->Fill(muons.at(1).ecalIso());
  mu2caloIso ->Fill(muons.at(1).caloIso());
  mu2dB      ->Fill(muons.at(1).dB());
  
  mu1trackRelIso->Fill(muons.at(0).trackIso()/muons.at(0).pt());
  mu1hcalRelIso ->Fill(muons.at(0).hcalIso()/muons.at(0).pt());
  mu1ecalRelIso ->Fill(muons.at(0).ecalIso()/muons.at(0).pt());
  mu1caloRelIso ->Fill(muons.at(0).caloIso()/muons.at(0).pt());
  mu2trackRelIso->Fill(muons.at(1).trackIso()/muons.at(1).pt());
  mu2hcalRelIso ->Fill(muons.at(1).hcalIso()/muons.at(1).pt());
  mu2ecalRelIso ->Fill(muons.at(1).ecalIso()/muons.at(1).pt());
  mu2caloRelIso ->Fill(muons.at(1).caloIso()/muons.at(1).pt());

  if (isMC) {
    for (unsigned int i=0; i<2; i++) { 
      if ( muons.at(i).genLepton() != 0 ) {
	float dpt = muons.at(i).pt()-muons.at(i).genLepton()->pt() ; 
	float dR = deltaR(muons.at(i).eta(),muons.at(i).phi(),
			  muons.at(i).genLepton()->eta(),muons.at(i).genLepton()->phi()) ; 
	if ( i == 0 ) { 
	  dptMu1gen->Fill(dpt/muons.at(i).genLepton()->pt()) ; 
	  dRMu1gen->Fill(dR) ; 
	} else { 
	  dptMu2gen->Fill(dpt/muons.at(i).genLepton()->pt()) ; 
	  dRMu2gen->Fill(dR) ; 
	}
      }
    }
  }
  for (int i=0; i<muonQualityFlags; i++) { 
    if (muons.at(0).muonID(muonQuality[i])) qualMu1->Fill( i ) ; 
    if (muons.at(1).muonID(muonQuality[i])) qualMu2->Fill( i ) ; 
  }

  // Jets 
  ptJet1->Fill(jets.at(0).pt()) ; 
  ptJet2->Fill(jets.at(1).pt()) ; 

  etaJet1->Fill(jets.at(0).eta()) ; 
  etaJet2->Fill(jets.at(1).eta()) ; 

  phiJet1->Fill(jets.at(0).phi()) ; 
  phiJet2->Fill(jets.at(1).phi()) ; 

  dPhiJet->Fill( fabs(deltaPhi(jets.at(0).phi(),jets.at(1).phi())) ) ; 
  dEtaJet->Fill( fabs(jets.at(0).eta() - jets.at(1).eta()) ) ; 
  dEtaPhiJet->Fill(fabs(jets.at(0).eta()-jets.at(1).eta()),
		   fabs(deltaPhi(jets.at(0).phi(),jets.at(1).phi()))) ;  

  // Muon-Jet plots
  float dRmu1jet1 = deltaR(muons.at(0).eta(),muons.at(0).phi(),jets.at(0).eta(),jets.at(0).phi()) ; 
  float dRmu1jet2 = deltaR(muons.at(0).eta(),muons.at(0).phi(),jets.at(1).eta(),jets.at(1).phi()) ; 
  float dRmu2jet1 = deltaR(muons.at(1).eta(),muons.at(1).phi(),jets.at(0).eta(),jets.at(0).phi()) ; 
  float dRmu2jet2 = deltaR(muons.at(1).eta(),muons.at(1).phi(),jets.at(1).eta(),jets.at(1).phi()) ; 

  const pat::Jet&  j4mu1 = (dRmu1jet1 < dRmu1jet2) ? jets.at(0) : jets.at(1);
  const pat::Jet&  j4mu2 = (dRmu2jet1 < dRmu2jet2) ? jets.at(0) : jets.at(1);

  TVector3 mu1vec(muons.at(0).momentum().X(), muons.at(0).momentum().Y(), muons.at(0).momentum().Z());
  TVector3 mu2vec(muons.at(1).momentum().X(), muons.at(1).momentum().Y(), muons.at(1).momentum().Z());

  TVector3 jt1vec(j4mu1.p4().Vect().X(), j4mu1.p4().Vect().Y(), j4mu1.p4().Vect().Z() );
  TVector3 jt2vec(j4mu2.p4().Vect().X(), j4mu2.p4().Vect().Y(), j4mu2.p4().Vect().Z() );

  double ptrelMu1 = mu1vec.Perp(jt1vec);
  double ptrelMu2 = mu2vec.Perp(jt2vec);

  dRminMu1jet->Fill( min(dRmu1jet1,dRmu1jet2) ) ; 
  dRminMu2jet->Fill( min(dRmu2jet1,dRmu2jet2) ) ; 

  hptrelMu1->Fill( ptrelMu1 );
  hptrelMu2->Fill( ptrelMu2 );

  ptrelVsdRminMu1jet->Fill(min(dRmu1jet1,dRmu1jet2),ptrelMu1);
  ptrelVsdRminMu2jet->Fill(min(dRmu2jet1,dRmu1jet2),ptrelMu2);

  // Composite objects
  reco::Particle::LorentzVector vWR = jets.at(0).p4() + jets.at(1).p4() ; 
  WR = vWR + muons.at(0).p4() + muons.at(1).p4() ; 

  mWR->Fill(WR.M()) ; 
  mNuR1->Fill( (vWR + muons.at(0).p4()).M() ) ; 
  mNuR2->Fill( (vWR + muons.at(1).p4()).M() ) ; 
  mNuR2D->Fill( (vWR + muons.at(0).p4()).M(),(vWR + muons.at(1).p4()).M() ) ; 

  reco::Particle::LorentzVector mumu=muons.at(0).p4()+muons.at(1).p4();
  reco::Particle::LorentzVector jj=jets.at(0).p4()+jets.at(1).p4();

  mMuMu->Fill(mumu.M() );
  mMuMuZoom->Fill(mumu.M() );
  mJJ->Fill(jj.M() );

}// end of fill()

void HeavyNu::HistPerDef::fill(const HeavyNuEvent& hne) {  

  // Muons 
  ptMu1->Fill(hne.mu1->pt()) ; 
  ptMu2->Fill(hne.mu2->pt()) ; 

  etaMu1->Fill(hne.mu1->eta()) ; 
  etaMu2->Fill(hne.mu2->eta()) ; 

  phiMu1->Fill(hne.mu1->phi()) ; 
  phiMu2->Fill(hne.mu2->phi()) ; 

  dPhiMu->Fill( fabs(deltaPhi(hne.mu1->phi(),hne.mu2->phi())) ) ; 
  dEtaMu->Fill( fabs(hne.mu1->eta() - hne.mu2->eta()) ) ; 
  dEtaPhiMu->Fill(fabs(hne.mu1->eta()-hne.mu2->eta()),
		  fabs(deltaPhi(hne.mu1->phi(),hne.mu2->phi()))) ; 

  mu1trackIso->Fill(hne.mu1->trackIso());
  mu1hcalIso ->Fill(hne.mu1->hcalIso());
  mu1ecalIso ->Fill(hne.mu1->ecalIso());
  mu1caloIso ->Fill(hne.mu1->caloIso());
  mu1dB      ->Fill(hne.mu1->dB());
  mu2trackIso->Fill(hne.mu2->trackIso());
  mu2hcalIso ->Fill(hne.mu2->hcalIso());
  mu2ecalIso ->Fill(hne.mu2->ecalIso());
  mu2caloIso ->Fill(hne.mu2->caloIso());
  mu2dB      ->Fill(hne.mu2->dB());
  
  mu1trackRelIso->Fill(hne.mu1->trackIso()/hne.mu1->pt());
  mu1hcalRelIso ->Fill(hne.mu1->hcalIso()/hne.mu1->pt());
  mu1ecalRelIso ->Fill(hne.mu1->ecalIso()/hne.mu1->pt());
  mu1caloRelIso ->Fill(hne.mu1->caloIso()/hne.mu1->pt());
  mu2trackRelIso->Fill(hne.mu2->trackIso()/hne.mu2->pt());
  mu2hcalRelIso ->Fill(hne.mu2->hcalIso()/hne.mu2->pt());
  mu2ecalRelIso ->Fill(hne.mu2->ecalIso()/hne.mu2->pt());
  mu2caloRelIso ->Fill(hne.mu2->caloIso()/hne.mu2->pt());
  
  if (hne.isMC) {
    for (unsigned int i=0; i<2; i++) { 
      if ( hne.mu[i]->genLepton() != 0 ) {
	float dpt = hne.mu[i]->pt()-hne.mu[i]->genLepton()->pt() ; 
	float dR  = deltaR(hne.mu[i]->eta(),hne.mu[i]->phi(),
			   hne.mu[i]->genLepton()->eta(),hne.mu[i]->genLepton()->phi()) ; 
	if ( i == 0 ) { 
	  dptMu1gen->Fill(dpt/hne.mu[i]->genLepton()->pt()) ; 
	  dRMu1gen->Fill(dR) ; 
	} else { 
	  dptMu2gen->Fill(dpt/hne.mu[i]->genLepton()->pt()) ; 
	  dRMu2gen->Fill(dR) ; 
	}
      }
    }
  }
  for (int i=0; i<muonQualityFlags; i++) { 
    if (hne.mu1->muonID(muonQuality[i])) qualMu1->Fill( i ) ; 
    if (hne.mu2->muonID(muonQuality[i])) qualMu2->Fill( i ) ; 
  }

  // Jets 
  if (hne.j1) {
    ptJet1->Fill(hne.j1->pt()) ; 
    etaJet1->Fill(hne.j1->eta()) ; 
    phiJet1->Fill(hne.j1->phi()) ; 

    if (hne.j2) {
      ptJet2->Fill(hne.j2->pt()) ; 
      etaJet2->Fill(hne.j2->eta()) ; 
      phiJet2->Fill(hne.j2->phi()) ; 

      dPhiJet->Fill( fabs(deltaPhi(hne.j1->phi(),hne.j2->phi())) ) ; 
      dEtaJet->Fill( fabs(hne.j1->eta() - hne.j2->eta()) ) ; 
      dEtaPhiJet->Fill(fabs(hne.j1->eta()-hne.j2->eta()),
		       fabs(deltaPhi(hne.j1->phi(),hne.j2->phi()))) ;  

      mWR->Fill   ( hne.mWR   ) ; 
      mNuR1->Fill ( hne.mNuR1 ) ; 
      mNuR2->Fill ( hne.mNuR2 ) ; 
      mNuR2D->Fill( hne.mNuR1, hne.mNuR2 );
      mJJ->Fill   ( hne.mJJ   );

      ctheta_jj->Fill(hne.ctheta_jj);
      ctheta_mu1_jj->Fill(hne.ctheta_mu1_jj);
      ctheta_mu2_jj->Fill(hne.ctheta_mu2_jj);
      cthetaz_jj->Fill(hne.cthetaz_jj);
      cthetaz_mu1_jj->Fill(hne.cthetaz_mu1_jj);
      cthetaz_mu2_jj->Fill(hne.cthetaz_mu2_jj);
    }

    dRminMu1jet->Fill(hne.dRminMu1jet);
    dRminMu2jet->Fill(hne.dRminMu2jet);

    hptrelMu1->Fill(hne.ptrelMu1);
    hptrelMu2->Fill(hne.ptrelMu2);

    ptrelVsdRminMu1jet->Fill(hne.dRminMu1jet,hne.ptrelMu1);
    ptrelVsdRminMu2jet->Fill(hne.dRminMu2jet,hne.ptrelMu2);
  }
  mMuMu->Fill( hne.mMuMu );
  mMuMuZoom->Fill( hne.mMuMu );

  ctheta_mumu->Fill(hne.ctheta_mumu);
  cthetaz_mumu->Fill(hne.cthetaz_mumu);

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
HeavyNu::HeavyNu(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  dolog_=iConfig.getParameter<bool>("DoLog");

  nnif_ = new HeavyNu_NNIF(iConfig);

  edm::Service<TFileService> fs;
  hists.nelec    = fs->make<TH1D>("nelec","N(e^{#pm})",10,-0.5,9.5);
  hists.nmu      = fs->make<TH1D>("nmu","N(#mu^{#pm})",10,-0.5,9.5);
  hists.njet     = fs->make<TH1D>("njet","N(Jet)",50,-0.5,49.5);
  hists.muPt     = fs->make<TH1D>("muPt","#mu p_{T} distribution",100,0,1000) ; 
  hists.muEta    = fs->make<TH1D>("muEta","#mu #eta distribution",50,-2.5,2.5) ; 
  hists.muPhi    = fs->make<TH1D>("muPhi","#mu #phi distribution",60,-3.14159,3.14159) ; 
  hists.jetPt    = fs->make<TH1D>("jetPt","jet p_{T} distribution",100,0,1000) ; 
  hists.jetEta   = fs->make<TH1D>("jetEta","jet #eta distribution",50,-5,5) ; 
  hists.jetPhi   = fs->make<TH1D>("jetPhi","jet #phi distribution",60,-3.14159,3.14159) ; 
  hists.jetPtvsNum=fs->make<TH2D>("jetPtvsNum","Jet P_{T} vs. Jet # ",11,-0.5,10.5,200,0.,1000.);
  hists.noCuts.book(fs->mkdir("noCuts"),"(no cuts)");
  hists.LLptCuts.book(fs->mkdir("LLptcuts"),"(dileptons with ptcuts:1)");
  hists.LLJJptCuts.book(fs->mkdir("LLJJptcuts"),"(4objects with ptcuts:1)");
  hists.massCut.book(fs->mkdir("masscut"),"(mumu mass cut:2)");
  init_=false;

  cuts.minimum_mu1_pt      = iConfig.getParameter<double>("minMu1pt");
  cuts.minimum_mu2_pt      = iConfig.getParameter<double>("minMu2pt");
  cuts.minimum_jet_pt      = iConfig.getParameter<double>("minJetPt");
  cuts.maximum_jet_abseta  = iConfig.getParameter<double>("maxJetAbsEta");
  cuts.minimum_mumu_mass   = iConfig.getParameter<double>("minMuMuMass");
  cuts.minimum_muon_jet_dR = iConfig.getParameter<double>("minMuonJetdR");

  // For the record...
  std::cout << "Configurable cut values applied:" << std::endl;
  std::cout << "minMu1pt      = " << cuts.minimum_mu1_pt      << std::endl;
  std::cout << "minMu2pt      = " << cuts.minimum_mu2_pt      << std::endl;
  std::cout << "minJetPt      = " << cuts.minimum_jet_pt      << std::endl;
  std::cout << "maxJetAbsEta  = " << cuts.maximum_jet_abseta  << std::endl;
  std::cout << "minMuonJetdR  = " << cuts.minimum_muon_jet_dR << std::endl;
  std::cout << "minMuMuMass   = " << cuts.minimum_mumu_mass   << std::endl;
}
  
HeavyNu::~HeavyNu()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
bool
HeavyNu::isVBTFloose(const pat::Muon& m)
{
  return (m.muonID("AllGlobalMuons") &&
	  (m.numberOfValidHits() > 10));
}

bool
HeavyNu::isVBTFtight(const pat::Muon& m)
{
  assert(isVBTFloose(m));

  reco::TrackRef gt = m.globalTrack();
  if (gt.isNull()) {
    std::cerr << "Mu1 global track reference is NULL" << std::endl;
    return false;
  }
  return (m.muonID("AllTrackerMuons") &&
	  (m.dB() < 0.2) &&
	  (m.normChi2() < 10) &&
	  (m.numberOfMatches() > 1) &&
	  (gt->hitPattern().numberOfValidMuonHits()>0) &&
	  (gt->hitPattern().numberOfValidPixelHits()>0) );
}

// ------------ method called to for each event  ------------
bool
HeavyNu::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  HeavyNuEvent hnuEvent;

  hnuEvent.isMC = !iEvent.isRealData();

  edm::Handle<pat::MuonCollection> patMuonCollection ; 
  iEvent.getByLabel("patMuons",patMuonCollection) ; 

  edm::Handle<pat::ElectronCollection> patElectronCollection ;
  iEvent.getByLabel("patElectrons", patElectronCollection) ;

  edm::Handle<pat::JetCollection> patJetCollection ;
  iEvent.getByLabel("patJets", patJetCollection) ;

  if ( !patElectronCollection.isValid() || 
       !patMuonCollection.isValid() || 
       !patJetCollection.isValid() ) {
    std::cout << "Exiting as valid PAT objects not found" << std::endl ;
    return false; 
  }

  const pat::ElectronCollection& pElecs = *(patElectronCollection.product()) ;
  const pat::MuonCollection& pMuons = *(patMuonCollection.product()) ;
  const pat::JetCollection& pJets = *(patJetCollection.product()) ;

  pat::ElectronCollection::const_iterator iE ;
  pat::MuonCollection::const_iterator iM ;
  pat::JetCollection::const_iterator iJ ;
  
  hists.nelec->Fill(pElecs.size()) ;
  hists.nmu->Fill(pMuons.size()) ;
  hists.njet->Fill(pJets.size()) ;

  int njet=1;
  for (iJ=pJets.begin(); iJ!=pJets.end(); iJ++) { 
    hists.jetPt->Fill( (*iJ).pt() ) ; 
    hists.jetEta->Fill( (*iJ).eta() ) ; 
    hists.jetPhi->Fill( (*iJ).phi() ) ; 
    hists.jetPtvsNum->Fill(njet++, (*iJ).pt() ) ; 
  }
  for (iM=pMuons.begin(); iM!=pMuons.end(); iM++) { 
    hists.muPt->Fill( (*iM).pt() ) ; 
    hists.muEta->Fill( (*iM).eta() ) ; 
    hists.muPhi->Fill( (*iM).phi() ) ; 
  }

  // Basic selection requirements: Require at least two muons, two jets
  if ( pMuons.size() >= 2 && pJets.size() >= 2 ) {
    hists.noCuts.fill( pMuons,pJets, hnuEvent.isMC ) ; 
  } else return false;

  
  // next, we look for valid muons and jets and put them into the Event
  for (iJ=pJets.begin(); iJ!=pJets.end(); iJ++) { 
    if (((*iJ).pt()        >cuts.minimum_jet_pt)   && // more later!
	(fabs((*iJ).eta())<=cuts.maximum_jet_abseta) ) {
      if (hnuEvent.j1==0 || hnuEvent.j1->pt()<(*iJ).pt()) {
	hnuEvent.j2=hnuEvent.j1;
	hnuEvent.j1=&(*iJ);
      } else 	if (hnuEvent.j2==0 || hnuEvent.j2->pt()<(*iJ).pt()) {
	hnuEvent.j2=&(*iJ);
      }
    }
  }

  for (iM=pMuons.begin(); iM!=pMuons.end(); iM++) { 
    double dr1=(hnuEvent.j1==0)?(10.0):(deltaR((*iM).eta(),(*iM).phi(),hnuEvent.j1->eta(),hnuEvent.j1->phi()));
    double dr2=(hnuEvent.j2==0)?(10.0):(deltaR((*iM).eta(),(*iM).phi(),hnuEvent.j2->eta(),hnuEvent.j2->phi()));

    if ((*iM).pt()>cuts.minimum_mu2_pt
	&& isVBTFloose(*iM)
	&& std::min(dr1,dr2)>cuts.minimum_muon_jet_dR) {
      if ( (hnuEvent.mu1==0) ||
	   (hnuEvent.mu1->pt()<(*iM).pt()) ) {
	hnuEvent.mu2=hnuEvent.mu1;
	hnuEvent.mu1=&(*iM);
      } else 	if (hnuEvent.mu2==0 ||
		    hnuEvent.mu2->pt()<(*iM).pt()) {
	hnuEvent.mu2=&(*iM);
      }
    }
  }

  // require two "loose" muons first
  if (hnuEvent.mu2==0)            return false;

  hnuEvent.regularize(); // assign internal standards
  hnuEvent.calculateMuMu();
  hists.LLptCuts.fill(hnuEvent);
  
  // require four objects
  if (hnuEvent.mu2==0 || hnuEvent.j2==0) return false;

  // apply the mu1 cut
  if ( !isVBTFtight(*(hnuEvent.mu1)) ||
       (hnuEvent.mu1->pt()<=cuts.minimum_mu1_pt))
      return false;

  hnuEvent.calculate(); // calculate various details
  hists.LLJJptCuts.fill(hnuEvent);

  if (hnuEvent.mMuMu<cuts.minimum_mumu_mass) return false;  // dimuon mass cut
  hists.massCut.fill(hnuEvent);

  nnif_->fillvector(hnuEvent);
  nnif_->print();

  if (iEvent.isRealData())
    std::cout<<iEvent.id() << std::endl;
  return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNu::beginJob() {
  nnif_->beginJob();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNu::endJob() {
  nnif_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNu);
