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
// $Id: HeavyNuFilter.cc,v 1.6 2012/12/01 00:35:07 pastika Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>
//#include <c++/4.1.2/bits/stl_vector.h>

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
#include "TH3.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuMuHist.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

class HeavyNuFilter : public edm::EDFilter
{
public:
    explicit HeavyNuFilter(const edm::ParameterSet&);
    ~HeavyNuFilter();


private:

    virtual void beginJob();
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob();

    bool Overlap(double eleEta, double elePhi, double xEta, double xPhi)
    {
        return ( deltaR(eleEta, elePhi, xEta, xPhi) <  0.1) ;
    }

    edm::InputTag muonTag_;
    edm::InputTag jetTag_;
    edm::InputTag elecTag_;
    edm::InputTag elecRhoTag_ ; 

    int evtCounter;
    bool isPFJets_; // true if PFJets are used (turns off jet id requirement)

    int mode_;

    int heepVersion_;

    // ----------member data ---------------------------
    bool init_;
    
    TRandom3 *tr3;

    // gf set of histo for all Z definitions in a stack

    struct CutsStruct
    {
        double minimum_mu1_pt;
        double minimum_mu2_pt;
        double minimum_jet_pt;
        double maximum_mu_abseta;
        double maximum_elec_abseta;
        double maximum_jet_abseta;
        double minimum_mumu_mass;
        double minimum_mWR_mass;
        double minimum_muon_jet_dR;
        double muon_trackiso_limit;
        double maxVertexZsep;
        double maxJetVZsepCM;
    } cuts;

} ;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HeavyNuFilter::HeavyNuFilter(const edm::ParameterSet& iConfig)
{
    // ==================== Get all parameters ====================
    //

    muonTag_ = iConfig.getParameter< edm::InputTag > ("muonTag");
    jetTag_ = iConfig.getParameter< edm::InputTag > ("jetTag");
    elecTag_ = iConfig.getParameter< edm::InputTag > ("electronTag");

    mode_ = iConfig.getParameter<int>("mode");

    cuts.minimum_mu1_pt = iConfig.getParameter<double>("minMu1pt");
    cuts.minimum_mu2_pt = iConfig.getParameter<double>("minMu2pt");
    cuts.minimum_jet_pt = iConfig.getParameter<double>("minJetPt");
    cuts.maximum_mu_abseta = iConfig.getParameter<double>("maxMuAbsEta");
    cuts.maximum_elec_abseta = iConfig.getParameter<double>("maxElecAbsEta");
    cuts.maximum_jet_abseta = iConfig.getParameter<double>("maxJetAbsEta");
    cuts.minimum_muon_jet_dR = iConfig.getParameter<double>("minMuonJetdR");
    cuts.muon_trackiso_limit = iConfig.getParameter<double>("muonTrackRelIsoLimit");
    cuts.maxVertexZsep = iConfig.getParameter<double>("maxVertexZsepCM");
    cuts.maxJetVZsepCM = iConfig.getParameter<double>("maxJetVZsepCM");

    isPFJets_ = iConfig.getParameter<bool>("isPFJets");

    heepVersion_ = iConfig.getUntrackedParameter<int>("heepVersion");

    elecRhoTag_ = iConfig.getParameter< edm::InputTag > ("electronRho");


    // For the record...
    std::cout << "Configurable cut values applied:" << std::endl;
    std::cout << "muonTag           = " << muonTag_ << std::endl;
    std::cout << "jetTag            = " << jetTag_ << std::endl;
    std::cout << "electronTag       = " << elecTag_ << std::endl;
    std::cout << "minMu1pt          = " << cuts.minimum_mu1_pt << " GeV" << std::endl;
    std::cout << "minMu2pt          = " << cuts.minimum_mu2_pt << " GeV" << std::endl;
    std::cout << "minJetPt          = " << cuts.minimum_jet_pt << " GeV" << std::endl;
    std::cout << "maxMuAbsEta       = " << cuts.maximum_mu_abseta << std::endl;
    std::cout << "maxElecAbsEta     = " << cuts.maximum_mu_abseta << std::endl;
    std::cout << "maxJetAbsEta      = " << cuts.maximum_jet_abseta << std::endl;
    std::cout << "minMuonJetdR      = " << cuts.minimum_muon_jet_dR << std::endl;
    std::cout << "muonTrackRelIso   = " << cuts.muon_trackiso_limit << std::endl;
    std::cout << "isPFJets          = " << isPFJets_ << std::endl;
    
    tr3 = new TRandom3(16);

}

HeavyNuFilter::~HeavyNuFilter(){

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to for each event  ------------

bool HeavyNuFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    HeavyNuEvent hnuEvent(HeavyNuEvent::HNUMU);

    evtCounter++;

    hnuEvent.isMC = !iEvent.isRealData();
    hnuEvent.pfJets = isPFJets_;
    hnuEvent.scaleMuE(1.0, 1.0);

    edm::Handle<pat::MuonCollection> pMuons;
    iEvent.getByLabel(muonTag_, pMuons);

    edm::Handle<pat::ElectronCollection> pElecs;
    iEvent.getByLabel(elecTag_, pElecs);

    edm::Handle<pat::JetCollection> pJets;
    iEvent.getByLabel(jetTag_, pJets);

    edm::Handle<reco::MuonCollection> tevMuons;
    iEvent.getByLabel("refitMuons", tevMuons);

    //Shirpa reweighting info
    edm::Handle<GenEventInfoProduct> geneventinfo;
    iEvent.getByLabel("generator", geneventinfo);

    edm::Handle<double> electronRhoHandle ;
    iEvent.getByLabel(elecRhoTag_, electronRhoHandle) ;

    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    
    if(mode_ == 1)
    { // This mode is to single out 0 Jet events from the inclusive DY sample
    	edm::Handle<LHEEventProduct> LHEinfo;
        iEvent.getByLabel("source", LHEinfo);
        if(LHEinfo->hepeup().NUP != 5) return false;
    }

    if(!pElecs.isValid() || !pMuons.isValid() || !pJets.isValid())
    {
        std::cout << "Exiting as valid PAT objects not found" << std::endl;
        std::cout << "Electrons: " << pElecs.isValid() << std::endl;
        std::cout << "Muons:     " << pMuons.isValid() << std::endl;
        std::cout << "Jets:      " << pJets.isValid() << std::endl;
        return false;
    }

    // Look for valid jets and put them in the event
    std::vector< std::pair<pat::Jet, float> > jetCands =
            hnu::getJetList(pJets, NULL, cuts.minimum_jet_pt, cuts.maximum_jet_abseta, 0, 0, false, 0); //the isMC = false is intentional eventhough this is a MC skim
    hnuEvent.nJets = jetCands.size();

    // Look for valid muons
    std::vector<pat::Muon> muCands =
            hnu::getMuonList(pMuons, pvHandle, 0, cuts.minimum_mu2_pt, cuts.maximum_mu_abseta, 1.0, false, false);

    // Look for valid electrons
    std::vector< std::pair<pat::Electron, float> > eCands =
            hnu::getElectronList(pElecs, cuts.maximum_elec_abseta, cuts.minimum_mu2_pt, cuts.minimum_mu2_pt, 0, pvHandle, ((electronRhoHandle.isValid()) ? (*(electronRhoHandle.product())) : 0.));


    if(hnuEvent.nJets < 2) return false;

    hnuEvent.j1 = jetCands.at(0).first;
    hnuEvent.j2 = jetCands.at(1).first;
    hnuEvent.j1scale = jetCands.at(0).second;
    hnuEvent.j2scale = jetCands.at(1).second;

    if(muCands.size() >= 1)
    {
        hnuEvent.mu1 = muCands.at(0);
        hnuEvent.nMuons++;
    }
    if(muCands.size() >= 2)
    {
        hnuEvent.mu2 = muCands.at(1);
        hnuEvent.nMuons++;
    }

    if(eCands.size() >= 1)
    {
        hnuEvent.e1 = eCands.at(0).first;
        hnuEvent.nElectrons++;
    }
    if(eCands.size() >= 2)
    {
        hnuEvent.e2 = eCands.at(1).first;
        hnuEvent.nElectrons++;
    }

    if(hnuEvent.j1.pt() < cuts.minimum_jet_pt || hnuEvent.j2.pt() < cuts.minimum_jet_pt) return false;

  
    if(tr3->Uniform() > 0.9)
    {
            for(std::vector<pat::Muon>::const_iterator iM = muCands.begin(); iM != muCands.end(); ++iM)
            {
                double dRlj1 = deltaR(iM->p4(), hnuEvent.j1.p4());
                double dRlj2 = deltaR(iM->p4(), hnuEvent.j2.p4());
                if(iM->pt() > cuts.minimum_mu2_pt && std::min(dRlj1, dRlj2) > cuts.minimum_muon_jet_dR) return true;
            }

            if(hnuEvent.nElectrons >= 1 &&  hnu::getElectronEt(hnuEvent.e1, false) > cuts.minimum_mu2_pt)
            {
                if(hnuEvent.nElectrons >= 2 && hnu::getElectronEt(hnuEvent.e2, false) > cuts.minimum_mu2_pt) return true;

                edm::Handle<reco::SuperClusterCollection> ebSCCollection ;
                iEvent.getByLabel("correctedHybridSuperClusters", ebSCCollection) ;
                const reco::SuperClusterCollection ebSCs = *(ebSCCollection.product()) ;

                edm::Handle<reco::SuperClusterCollection> eeSCCollection ;
                iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower", eeSCCollection) ;
                const reco::SuperClusterCollection eeSCs = *(eeSCCollection.product()) ;

                // EE sc's
                for(unsigned int eesc = 0; eesc < eeSCs.size(); eesc++ )
                {
                    if ( eeSCs.at(eesc).energy() / cosh(eeSCs.at(eesc).eta()) < cuts.minimum_mu2_pt) continue;

                    // make sure this SC is not already among the selected electrons
                    bool reject = false;
                    for(unsigned int sEl = 0; sEl < eCands.size(); sEl++)
                    {
                        if((reject |= Overlap(eeSCs.at(eesc).eta(), eeSCs.at(eesc).phi(), eCands.at(sEl).first.eta(), eCands.at(sEl).first.phi()))) break;
                    }// loop over selected electrons
                    if(reject) continue;

                    if ( fabs( eeSCs.at(eesc).eta() ) < cuts.maximum_mu_abseta) return true;
                }// loop over EE SC

                // EB sc's
                for(unsigned int ebsc = 0; ebsc < ebSCs.size(); ebsc++ )
                {
                    if ( ebSCs.at(ebsc).energy() / cosh(ebSCs.at(ebsc).eta()) < cuts.minimum_mu2_pt) continue;

                    // make sure this SC is not already among the selected electrons
                    bool reject = false;
                    for(unsigned int sEl = 0; sEl < eCands.size(); sEl++)
                    {
                        if((reject |= Overlap(ebSCs.at(ebsc).eta(), ebSCs.at(ebsc).phi(), eCands.at(sEl).first.eta(), eCands.at(sEl).first.phi()))) break;
                    }// loop over selected electrons
                    if(reject) continue;

                    if ( fabs( ebSCs.at(ebsc).eta() ) < cuts.maximum_mu_abseta) return true;
                }
            }
            if(hnuEvent.nMuons + hnuEvent.nElectrons < 2) return false;
    }

    if(hnuEvent.nMuons >= 2     && hnuEvent.mu1.pt() > cuts.minimum_mu2_pt && hnuEvent.mu2.pt() > cuts.minimum_mu2_pt) return true;
    if(hnuEvent.nElectrons >= 2 &&  hnu::getElectronEt(hnuEvent.e1, false) > cuts.minimum_mu2_pt &&  hnu::getElectronEt(hnuEvent.e2, false) > cuts.minimum_mu2_pt) return true;
    if(hnuEvent.nMuons >= 1 && hnuEvent.nElectrons >= 1)
    {
        if(hnu::getElectronEt(hnuEvent.e1, false)  > cuts.minimum_mu2_pt && hnuEvent.mu1.pt() > cuts.minimum_mu2_pt) return true;
    }
    
    return false;
}

// ------------ method called once each job just before starting event loop  ------------

void HeavyNuFilter::beginJob(){ }

// ------------ method called once each job just after ending the event loop  ------------

void HeavyNuFilter::endJob(){ }

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuFilter);
