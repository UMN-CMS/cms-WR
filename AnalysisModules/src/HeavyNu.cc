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
// $Id: HeavyNu.cc,v 1.3 2010/09/17 16:47:02 mansj Exp $
//
//


// system include files
#include <memory>

#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"


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

class HeavyNu : public edm::EDAnalyzer {
public:
  explicit HeavyNu(const edm::ParameterSet&);
  ~HeavyNu();


private:
  virtual void respondToOpenInputFile(edm::FileBlock const& fb) {
    currentFile_=fb.fileName();
  }
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  std::string currentFile_;
  bool dolog_;

  // ----------member data ---------------------------
  static const std::string muonQuality[] ; 
  static const int muonQualityFlags ; 

  struct HistPerDef {
    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory td, const std::string&) ;
    // fill all histos of the set with the two electron candidates
    void fill(pat::MuonCollection muons, pat::JetCollection jets) ;
    // fill all histos of the set with the two electron candidates
    void fill(const HeavyNuEvent& hne) ;

    TH1 *ptMu1, *ptMu2, *ptJet1, *ptJet2 ;  
    TH1 *etaMu1, *etaMu2, *etaJet1, *etaJet2 ;  
    TH1 *phiMu1, *phiMu2, *phiJet1, *phiJet2 ;  
    TH1 *dEtaMu, *dPhiMu, *dEtaJet, *dPhiJet ;
    TH1 *dEtaPhiMu, *dEtaPhiJet ; 
    TH1 *dRminMu1jet, *dRminMu2jet ; 
    TH1 *dptMu1gen, *dptMu2gen ; 
    TH1 *dRMu1gen, *dRMu2gen ; 
    TH1 *qualMu1, *qualMu2 ; 

    TH1 *mWR, *mNuR1, *mNuR2, *mMuMu, *mJJ ; 
    TH2 *mNuR2D ; 
  };

  bool init_;

  // gf set of histo for all Z definitios in a stack
  struct HistStruct {
    TH1 *nelec, *nmu, *njet ;
    TH1 *muPt, *muEta, *muPhi ; 
    TH1 *jetPt, *jetEta, *jetPhi ; 
    HistPerDef noCuts ; 
    HistPerDef ptCuts;
  } hists;

  struct CutsStruct {
    double minimum_muon_pt;
    double minimum_jet_pt;
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
  title=std::string("M(jj)")+post;
  mJJ=td.make<TH1D>("mJJ",title.c_str(),50,0,2000);  


}

void HeavyNu::HistPerDef::fill(pat::MuonCollection muons, pat::JetCollection jets) {  

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

  dRminMu1jet->Fill( min(dRmu1jet1,dRmu1jet2) ) ; 
  dRminMu2jet->Fill( min(dRmu2jet1,dRmu2jet2) ) ; 

  // Composite objects
  reco::Particle::LorentzVector vWR = jets.at(0).p4() + jets.at(1).p4() ; 
  WR = vWR + muons.at(0).p4() + muons.at(1).p4() ; 

  mWR->Fill(WR.M()) ; 
  mNuR1->Fill( (vWR + muons.at(0).p4()).M() ) ; 
  mNuR2->Fill( (vWR + muons.at(1).p4()).M() ) ; 
  mNuR2D->Fill( (vWR + muons.at(0).p4()).M(),(vWR + muons.at(1).p4()).M() ) ; 

  mMuMu->Fill((muons.at(0).p4()+muons.at(1).p4()).M() );
  mJJ->Fill((jets.at(0).p4()+jets.at(1).p4()).M() );


}// end of fill()

void HeavyNu::HistPerDef::fill(const HeavyNuEvent& hne) {  

  reco::Particle::LorentzVector WR ; 

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
  

  for (unsigned int i=0; i<2; i++) { 
    if ( hne.mu[i]->genLepton() != 0 ) {
      float dpt = hne.mu[i]->pt()-hne.mu[i]->genLepton()->pt() ; 
      float dR = deltaR(hne.mu[i]->eta(),hne.mu[i]->phi(),
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
  for (int i=0; i<muonQualityFlags; i++) { 
    if (hne.mu1->muonID(muonQuality[i])) qualMu1->Fill( i ) ; 
    if (hne.mu2->muonID(muonQuality[i])) qualMu2->Fill( i ) ; 
  }

  // Jets 
  ptJet1->Fill(hne.j1->pt()) ; 
  ptJet2->Fill(hne.j2->pt()) ; 

  etaJet1->Fill(hne.j1->eta()) ; 
  etaJet2->Fill(hne.j2->eta()) ; 

  phiJet1->Fill(hne.j1->phi()) ; 
  phiJet2->Fill(hne.j2->phi()) ; 

  dPhiJet->Fill( fabs(deltaPhi(hne.j1->phi(),hne.j2->phi())) ) ; 
  dEtaJet->Fill( fabs(hne.j1->eta() - hne.j2->eta()) ) ; 
  dEtaPhiJet->Fill(fabs(hne.j1->eta()-hne.j2->eta()),
		   fabs(deltaPhi(hne.j1->phi(),hne.j2->phi()))) ;  

  // Muon-Jet plots
  float dRmu1jet1 = deltaR(hne.mu1->eta(),hne.mu1->phi(),hne.j1->eta(),hne.j1->phi()) ; 
  float dRmu1jet2 = deltaR(hne.mu1->eta(),hne.mu1->phi(),hne.j2->eta(),hne.j2->phi()) ; 
  float dRmu2jet1 = deltaR(hne.mu2->eta(),hne.mu2->phi(),hne.j1->eta(),hne.j1->phi()) ; 
  float dRmu2jet2 = deltaR(hne.mu2->eta(),hne.mu2->phi(),hne.j2->eta(),hne.j2->phi()) ; 

  dRminMu1jet->Fill( min(dRmu1jet1,dRmu1jet2) ) ; 
  dRminMu2jet->Fill( min(dRmu2jet1,dRmu2jet2) ) ; 

  // Composite objects
  reco::Particle::LorentzVector vWR = hne.j1->p4() + hne.j2->p4() ; 
  WR = vWR + hne.mu1->p4() + hne.mu2->p4() ; 

  mWR->Fill(WR.M()) ; 
  mNuR1->Fill( (vWR + hne.mu1->p4()).M() ) ; 
  mNuR2->Fill( (vWR + hne.mu2->p4()).M() ) ; 
  mNuR2D->Fill( (vWR + hne.mu1->p4()).M(),(vWR + hne.mu2->p4()).M() ) ; 

  mMuMu->Fill((hne.mu1->p4()+hne.mu2->p4()).M() );
  mJJ->Fill((hne.j1->p4()+hne.j2->p4()).M() );

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
  hists.noCuts.book(fs->mkdir("noCuts"),"(no cuts)");
  hists.ptCuts.book(fs->mkdir("ptcuts"),"(diobjects with ptcuts)");
  init_=false;


  cuts.minimum_muon_pt=20;
  cuts.minimum_jet_pt=40;

}
  
HeavyNu::~HeavyNu()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HeavyNu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  HeavyNuEvent hnuEvent;


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
    return ; 
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

  for (iJ=pJets.begin(); iJ!=pJets.end(); iJ++) { 
    hists.jetPt->Fill( (*iJ).pt() ) ; 
    hists.jetEta->Fill( (*iJ).eta() ) ; 
    hists.jetPhi->Fill( (*iJ).phi() ) ; 
  }
  for (iM=pMuons.begin(); iM!=pMuons.end(); iM++) { 
    hists.muPt->Fill( (*iM).pt() ) ; 
    hists.muEta->Fill( (*iM).eta() ) ; 
    hists.muPhi->Fill( (*iM).phi() ) ; 
  }

  // Basic selection requirements: Require at least two muons, two jets
  if ( pMuons.size() >= 2 && pJets.size() >= 2 ) {
    hists.noCuts.fill( pMuons,pJets ) ; 
  } else return;

  
  // next, we look for valid muons and jets and put them into the Event
  for (iM=pMuons.begin(); iM!=pMuons.end(); iM++) { 
    if ((*iM).pt()>cuts.minimum_muon_pt && (*iM).muonID("AllGlobalMuons")) {
      if (hnuEvent.mu1==0 || hnuEvent.mu1->pt()<(*iM).pt()) {
	hnuEvent.mu2=hnuEvent.mu1;
	hnuEvent.mu1=&(*iM);
      } else 	if (hnuEvent.mu2==0 || hnuEvent.mu2->pt()<(*iM).pt()) {
	hnuEvent.mu2=&(*iM);
      }
    }
  }
  
  for (iJ=pJets.begin(); iJ!=pJets.end(); iJ++) { 
    if ((*iJ).pt()>cuts.minimum_jet_pt) { // more later!
      if (hnuEvent.j1==0 || hnuEvent.j1->pt()<(*iJ).pt()) {
	hnuEvent.j2=hnuEvent.j1;
	hnuEvent.j1=&(*iJ);
      } else 	if (hnuEvent.j2==0 || hnuEvent.j2->pt()<(*iJ).pt()) {
	hnuEvent.j2=&(*iJ);
      }
    }
  } 
  // require four objects
  if (hnuEvent.mu2==0 || hnuEvent.j2==0) return;
  hnuEvent.regularize(); // calculate various details
  hists.ptCuts.fill(hnuEvent);

}

// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNu::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNu::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNu);
