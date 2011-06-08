// -*- C++ -*-
//
// Package:    HeavyNuGenLevel
// Class:      HeavyNuGenLevel
// 
/**\class HeavyNuGenLevel HeavyNuGenLevel.cc JetCheck/HeavyNuGenLevel/src/HeavyNuGenLevel.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Alexander Gude
//         Created:  Thu May 12 11:15:22 CDT 2011
// $Id: HeavyNuGenLevel.cc,v 1.1 2011/06/03 03:35:31 mansj Exp $
//
//


// system include files
#include <memory>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"


#include "Math/VectorUtil.h"

#include <iostream>

//
// class declaration
//

class HeavyNuGenLevel : public edm::EDAnalyzer {
    public:
        explicit HeavyNuGenLevel(const edm::ParameterSet&);
        ~HeavyNuGenLevel();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

  struct CutStruct {
    double minJetPT;
    double minL1PT;
    double minL2PT;
    double minLJdR;
    double minLLMass;
    double min4objMass;
  } cuts;

  struct HistStruct {
    TH1 *cutProgress;
    TH1 *ptl1, *ptl2;
    TH1 *ptj1, *ptj2;
    TH1 *ljdR, *min_ljdR;
    TH1 *mll, *m4obj;
  } hists;

        // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HeavyNuGenLevel::HeavyNuGenLevel(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  cuts.minJetPT=iConfig.getParameter<double>("minJetPt");
  cuts.minL1PT=iConfig.getParameter<double>("minMu1pt");
  cuts.minL2PT=iConfig.getParameter<double>("minMu2pt");
  cuts.minLJdR=iConfig.getParameter<double>("minMuonJetdR");
  cuts.minLLMass=iConfig.getParameter<double>("minMuMuMass");
  cuts.min4objMass=iConfig.getParameter<double>("min4objMass");
  
  edm::Service<TFileService> fileService;

  hists.cutProgress = fileService->make<TH1F>("cutProgress","Progress through cut series",20,-0.5,19.5);
  hists.ptl1 = fileService->make<TH1F>("ptl1","Lepton 1 PT",50,0,500);
  hists.ptl2 = fileService->make<TH1F>("ptl2","Lepton 2 PT",50,0,400);
  hists.ptj1 = fileService->make<TH1F>("ptj1","Jet 1 PT",50,0,500);
  hists.ptj2 = fileService->make<TH1F>("ptj2","Jet 2 PT",50,0,500);
  
  hists.ljdR = fileService->make<TH1F>("ljdR","Lepton-Jet dR",50,0,5);
  hists.min_ljdR = fileService->make<TH1F>("min_ljdR","Minimum Lepton-Jet dR",50,0,5);

  hists.mll = fileService->make<TH1F>("mll","MLL",50,0,700);
  hists.m4obj = fileService->make<TH1F>("m4obj","MWR",50,0,2500);
}


HeavyNuGenLevel::~HeavyNuGenLevel()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // Need recoCaloJets_ak5CaloJets_RECO
    // recoJPTJets_JetPlusTrackZSPCorJetAntiKt5_RECO
}


//
// member functions
//

// ------------ method called to for each event  ------------
void HeavyNuGenLevel::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //using namespace edm;
    //using namespace ROOT::Math::VectorUtil; // DeltaR

  edm::Handle<edm::HepMCProduct> hepMCEvt;
  edm::Handle<reco::GenJetCollection> genJets;
  hists.cutProgress->Fill(0);
 
  iEvent.getByLabel("generator",hepMCEvt);
  iEvent.getByLabel("ak5GenJetsNoMuNoNu",genJets);

  const HepMC::GenEvent* genEvt=hepMCEvt->GetEvent(); 
  const HepMC::GenEvent* genE=genEvt;

  HepMC::GenEvent::vertex_const_iterator vtex;
  HepMC::GenVertex::particles_out_const_iterator Pout;
  HepMC::GenParticle* theLep=0;
  HepMC::GenParticle* theLep2=0;
  reco::GenJetCollection::const_iterator ji, jet1=genJets->end(), jet2=genJets->end();
  
  for (vtex=genE->vertices_begin();vtex!=genE->vertices_end();vtex++){
    for(Pout=(*vtex)->particles_out_const_begin();Pout!=(*vtex)->particles_out_const_end();Pout++){
      if (abs((*Pout)->pdg_id())==13 && (*Pout)->status()==1) {
	if (theLep==0 || theLep->momentum().perp()<(*Pout)->momentum().perp()) {
	  theLep2=theLep;
	  theLep=*Pout;
	} else if (theLep2==0 || theLep2->momentum().perp()<(*Pout)->momentum().perp()) {
	  theLep2=*Pout;
	}
      }
    }
  }
  if (theLep2==0) {
    std::cout << "Got less than two!\n";
    return;
  }


  hists.cutProgress->Fill(1);

  hists.ptl1->Fill(theLep->momentum().perp());
  hists.ptl2->Fill(theLep2->momentum().perp());

  if (theLep->momentum().perp()>cuts.minL1PT &&
      theLep2->momentum().perp()>cuts.minL2PT &&
      fabs(theLep->momentum().eta())<2.4 &&
      fabs(theLep2->momentum().eta())<2.4 &&
      (fabs(theLep->momentum().eta())<2.1 ||fabs(theLep2->momentum().eta())<2.1)
      ) {
    
    hists.cutProgress->Fill(2);

    for (ji=genJets->begin(); ji!=genJets->end(); ji++) {
      if (fabs(ji->eta())>2.5) continue;
      if (jet1==genJets->end() || jet1->pt()<ji->pt()) {
	jet2=jet1;
	jet1=ji;
      } else if (jet2==genJets->end() || jet2->pt()<ji->pt()) {
	jet2=ji;
      }
    }

    if (jet1!=genJets->end() && jet2!=genJets->end()) {
    
      hists.ptj1->Fill(jet1->pt());
      hists.ptj2->Fill(jet2->pt());

      if (jet1->pt()>cuts.minJetPT && 
	  jet2->pt()>cuts.minJetPT) {
	hists.cutProgress->Fill(3);
	
	double dR11=deltaR(theLep->momentum().eta(),theLep->momentum().phi(),
			   jet1->eta(),jet1->phi());
	double dR12=deltaR(theLep->momentum().eta(),theLep->momentum().phi(),
			   jet2->eta(),jet2->phi());
	double dR21=deltaR(theLep2->momentum().eta(),theLep2->momentum().phi(),
			   jet1->eta(),jet1->phi());
	double dR22=deltaR(theLep2->momentum().eta(),theLep2->momentum().phi(),
			   jet2->eta(),jet2->phi());

	hists.ljdR->Fill(dR11);
	hists.ljdR->Fill(dR12);
	hists.ljdR->Fill(dR21);
	hists.ljdR->Fill(dR22);
	hists.ljdR->Fill(dR11);
	hists.min_ljdR->Fill(std::min(std::min(dR11,dR12),std::min(dR21,dR22)));

	if (dR11 > cuts.minLJdR && dR12 > cuts.minLJdR && 
	    dR21 > cuts.minLJdR && dR22 > cuts.minLJdR) {
	  hists.cutProgress->Fill(4);

	  reco::Particle::LorentzVector lvl1(theLep->momentum().px(),theLep->momentum().py(),theLep->momentum().pz(),theLep->momentum().e());
	  reco::Particle::LorentzVector lvl2(theLep2->momentum().px(),theLep2->momentum().py(),theLep2->momentum().pz(),theLep2->momentum().e());

	  reco::Particle::LorentzVector ll=lvl1+lvl2;
	  
	  hists.mll->Fill(ll.M());
	  
	  if (ll.M()>cuts.minLLMass) {
	    hists.cutProgress->Fill(5);

	    reco::Particle::LorentzVector jj=jet1->p4()+jet2->p4();
	    reco::Particle::LorentzVector wr=jj+ll;

	    hists.m4obj->Fill(wr.M());

	    if (wr.M()>cuts.min4objMass) {
	      hists.cutProgress->Fill(6);
	    }

	  }


	}
      }

    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
    void 
HeavyNuGenLevel::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNuGenLevel::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuGenLevel);
