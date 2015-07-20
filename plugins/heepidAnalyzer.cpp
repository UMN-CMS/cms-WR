
// -*- C++ -*-
//
// Package:    ExoAnalysis/cmsWR
// Class:      heepIDAnalyzer
// 
/**\class heepIDAnalyzer heepIDAnalyzer.cc ExoAnalysis/cmsWR/plugins/heepIDAnalyzer.cc
 Description: Test analyzer to demonstrate use of heepID
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Peter Hansen
//         Created:  Mon Jul 20 16:03:38 2015
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

class heepIDAnalyzer : public edm::EDAnalyzer {
   public:
      explicit heepIDAnalyzer(const edm::ParameterSet&);
      ~heepIDAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
				 const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

      void findFirstNonElectronMother(const reco::Candidate *particle,
					 int &ancestorPID, int &ancestorStatus);

      // ----------member data ---------------------------

      // Data members that are the same for AOD and miniAOD
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

      // AOD case data members
      edm::EDGetToken electronsToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleHEEPIdMapCFRToken_;

  TTree *electronTree_;

  // all variables for the output tree
  Int_t nElectrons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;

  std::vector<Float_t> dz_;
  std::vector<Int_t> passConversionVeto_;

  std::vector<Int_t> passHEEPId_;

  std::vector<Int_t> isTrue_;

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
heepIDAnalyzer::heepIDAnalyzer(const edm::ParameterSet& iConfig):
  eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"))),
  eleHEEPIdMapCFRToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
{

  //
  // Prepare tokens for all input collections and objects
  //

  beamSpotToken_    = consumes<reco::BeamSpot> 
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));

  // AOD tokens
  electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));

  // MiniAOD tokens
  // For electrons, use the fact that pat::Electron can be cast into 
  // GsfElectron
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

  vtxMiniAODToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("verticesMiniAOD"));

  conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversionsMiniAOD"));

  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");
  
  electronTree_->Branch("nEle",  &nElectrons_ , "nEle/I");
  electronTree_->Branch("pt"  ,  &pt_    );
  electronTree_->Branch("eta" ,  &eta_ );
  electronTree_->Branch("phi" ,  &phi_ );

  electronTree_->Branch("dz" ,  &dz_ );
  electronTree_->Branch("passHEEPId"  ,  &passHEEPId_ );

  electronTree_->Branch("isTrue"             , &isTrue_);

}


heepIDAnalyzer::~heepIDAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
heepIDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  
  
  // Retrieve the collection of electrons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name.
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Electron objects can be recast as reco::GsfElectron objects.
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  bool isAOD = true;
  iEvent.getByToken(electronsToken_, electrons);
  if( !electrons.isValid() ){
    isAOD = false;
    iEvent.getByToken(electronsMiniAODToken_,electrons);
  }
  
  // Get the MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  if( isAOD )
    iEvent.getByToken(vtxToken_, vertices);
  else
    iEvent.getByToken(vtxMiniAODToken_, vertices);
  
  if (vertices->empty()) return; // skip the event if no PV found
  
  // Find the first vertex in the collection that passes
  // good quality criteria
  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = vtx->isFake();
    if( !isAOD )
      isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if ( !isFake
	 &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() )
    return; // skip event if there are no good PVs

  // Get the conversions collection
  edm::Handle<reco::ConversionCollection> conversions;
  if(isAOD)
    iEvent.getByToken(conversionsToken_, conversions);
  else
    iEvent.getByToken(conversionsMiniAODToken_, conversions);

  // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  iEvent.getByToken(eleHEEPIdMapToken_ ,heep_id_decisions);

  edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep_id_CFR;
  iEvent.getByToken(eleHEEPIdMapCFRToken_ ,heep_id_CFR);

  // Clear vectors
  nElectrons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  //
  dz_.clear();
  //
  passHEEPId_ .clear();
  //
  isTrue_.clear();

  // Loop over electrons
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);

    // Kinematics
    if( el->pt() < 10 ) // keep only electrons above 10 GeV
      continue;
    
    nElectrons_++;

    //
    // Save electron kinematics
    //
    pt_  .push_back( el->pt() );
    eta_ .push_back( el->superCluster()->eta() );
    phi_ .push_back( el->superCluster()->phi() );

    // Impact parameter
    reco::GsfTrackRef theTrack = el->gsfTrack();
    dz_.push_back( theTrack->dz( firstGoodVertex->position() ) );
    
    // Conversion rejection
    bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, 
									 conversions,
									 theBeamSpot->position());
    passConversionVeto_.push_back( (int) passConvVeto );

    //
    // Look up and save the ID decisions
    // 
    bool isPassHEEP = (*heep_id_decisions)[el];
    passHEEPId_.push_back ( (int)isPassHEEP );

	 cout << "[DEBUG]" << (*heep_id_CFR)[el].cutFlowName() << std::endl;
	 for (size_t i = 0; i < (*heep_id_CFR)[el].cutFlowSize(); i++)
	 {
		 cout << "[DEBUG]" << i << ' ' << (*heep_id_CFR)[el].getNameAtIndex(i) << std::endl;
		 cout << "[DEBUG]" << i << ' ' << (*heep_id_CFR)[el].getCutResultByIndex(i) << std::endl;
		 cout << "[DEBUG]" << i << ' ' << (*heep_id_CFR)[el].getValueCutUpon(i) << std::endl;
	 }

    // Save MC truth match
    isTrue_.push_back( matchToTruth( el, genParticles) );

   }
   
  // Save the info
  electronTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
heepIDAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
heepIDAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
heepIDAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
heepIDAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
heepIDAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
heepIDAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
heepIDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int heepIDAnalyzer::matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
				  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void heepIDAnalyzer::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(heepIDAnalyzer);
