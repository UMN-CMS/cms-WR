// -*- C++ -*-
//
// Package:    JetChecker
// Class:      JetChecker
// 
/**\class JetChecker JetChecker.cc JetCheck/JetChecker/src/JetChecker.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Alexander Gude
//         Created:  Thu May 12 11:15:22 CDT 2011
// $Id$
//
//


// system include files
#include <memory>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"

#include "Math/VectorUtil.h"

#include <iostream>

//
// class declaration
//

class JetChecker : public edm::EDAnalyzer {
    public:
        explicit JetChecker(const edm::ParameterSet&);
        ~JetChecker();

    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        struct HistStruct {
            TH1F *c1pt, *c2pt;
            TH1F *c1t1dR, *c2t2dR;
            TH1F *ctdz, *ttdz;
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
JetChecker::JetChecker(const edm::ParameterSet& iConfig)
{
    //now do what ever initialization is needed
    std::vector<double> maxDeltaVRs = iConfig.getParameter< std::vector<double> >("maxDeltaVRs");

    edm::Service<TFileService> fileService;

    for ( std::vector<double>::const_iterator i = maxDeltaVRs.begin(); i != maxDeltaVRs.end(); i++){
        std::cout << *i << std::endl;
    }
    hists.c1pt = fileService->make<TH1F>("c1pt","Calo Jet 1 PT",50,0,500);
    hists.c2pt = fileService->make<TH1F>("c2pt","Calo Jet 2 PT",50,0,300);
    hists.c1t1dR = fileService->make<TH1F>("c1t1dR","Jet 1 Delta R",50,0,2);
    hists.c2t2dR = fileService->make<TH1F>("c2t2dR","Jet 2 Delta R",50,0,2);
    hists.ctdz = fileService->make<TH1F>("ctdz","Calo - Track Vertex Z",50,-20,20);
    hists.ttdz = fileService->make<TH1F>("ttdz","Track 1 - Track 2 Vertex Z",200,-10,10);

}


JetChecker::~JetChecker()
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
void JetChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //using namespace edm;
    //using namespace ROOT::Math::VectorUtil; // DeltaR

    // Open our jets
    edm::Handle<reco::JPTJetCollection> jptjets;
    iEvent.getByLabel("JetPlusTrackZSPCorJetAntiKt5", jptjets);

    edm::Handle<reco::CaloJetCollection> calojets;
    iEvent.getByLabel("ak5CaloJets", calojets);

    // Loop through Calojets to find 2 leading
    reco::CaloJetCollection::const_iterator cjet1=calojets->end();
    reco::CaloJetCollection::const_iterator cjet2=calojets->end();
    for (reco::CaloJetCollection::const_iterator i=calojets->begin(); i!=calojets->end(); i++) {
        if (cjet1==calojets->end() || i->pt()>cjet1->pt()) {
            cjet2=cjet1;
            cjet1=i;
        } else if (cjet2==calojets->end() || i->pt()>cjet2->pt()) 
            cjet2=i;
    }
    if (cjet1==calojets->end() || cjet2==calojets->end()) return; // We don't _have_ two jets...

    hists.c1pt->Fill(cjet1->pt());
    hists.c2pt->Fill(cjet2->pt());

    if (cjet1->pt()<30 || cjet2->pt()<30) return; // need more umph 

    // Find best match jptjets
    reco::JPTJetCollection::const_iterator tjet1=jptjets->end();
    reco::JPTJetCollection::const_iterator tjet2=jptjets->end();
    for (reco::JPTJetCollection::const_iterator i=jptjets->begin(); i!=jptjets->end(); i++) {
        if (tjet1==jptjets->end() || ROOT::Math::VectorUtil::DeltaR(cjet1->p4(),i->p4()) < ROOT::Math::VectorUtil::DeltaR(cjet1->p4(),tjet1->p4())){
            tjet1=i;
        }
        if (tjet2==jptjets->end() || ROOT::Math::VectorUtil::DeltaR(cjet2->p4(),i->p4()) < ROOT::Math::VectorUtil::DeltaR(cjet2->p4(),tjet2->p4())){
            tjet2=i;
        }
    }
    hists.c1t1dR->Fill(ROOT::Math::VectorUtil::DeltaR(cjet1->p4(),tjet1->p4()));
    hists.c2t2dR->Fill(ROOT::Math::VectorUtil::DeltaR(cjet2->p4(),tjet2->p4()));

    // Calculating average z of track jets
    double jt1z=hnu::avgVertex(*tjet1, 1.);
    double jt2z=hnu::avgVertex(*tjet2, 1.);

    // Getting vertex from calojets
    double c1z=cjet1->vz();
    
    hists.ctdz->Fill(c1z-jt1z);
    hists.ctdz->Fill(c1z-jt2z);
    hists.ttdz->Fill(jt1z-jt2z);
}


// ------------ method called once each job just before starting event loop  ------------
    void 
JetChecker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetChecker::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(JetChecker);
