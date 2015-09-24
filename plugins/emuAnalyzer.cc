// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/emuAnalyzer
// Class:      emuAnalyzer
// 
/**\class emuAnalyzer emuAnalyzer.cc doubleElectronTracklessTrigger/emuAnalyzer/plugins/emuAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sean Kalafut
//         Created:  Wed, 15 April 2015 
//
//


// system include files
#include <memory>
#include <map>
#include <utility>
#include <cstring>
#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>



// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include "TTree.h"
//#include "TLorentzVector.h"
#include <TFile.h>
#include <TBranch.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include "TCollection.h"

//#define DEBUG

//
// class declaration
//

class emuAnalyzer : public edm::EDAnalyzer {
   public:
      explicit emuAnalyzer(const edm::ParameterSet&);
      ~emuAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	  ///calculate the dilepton mass using two const_iterators to reco::Candidate objects
	  double getDileptonMass(edm::OwnVector<reco::Candidate>::const_iterator& one, edm::OwnVector<reco::Candidate>::const_iterator& two){
		  double mass = TMath::Sqrt( 2*(one->pt())*(two->pt())*( TMath::CosH( (one->eta())-(two->eta()) ) - TMath::Cos( (one->phi())-(two->phi()) ) ) );
		  return mass;
	  }///end getDileptonMass()

	  bool isCloseToMuon(edm::OwnVector<reco::Candidate>::const_iterator objIt,edm::Handle<std::vector<pat::Muon> > particlesToAvoid){
		  for(std::vector<pat::Muon>::const_iterator avoidIt=particlesToAvoid->begin(); avoidIt!=particlesToAvoid->end(); avoidIt++){
			  if(avoidIt->pt() < 5) continue;
			  if(deltaR(objIt->eta(),objIt->phi(),avoidIt->eta(),avoidIt->phi()) < 0.1) return true;
		  }///end loop over objects in particlesToAvoid collection
		  return false;
	  }///end isCloseToMuon()

	  void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::Handle<edm::OwnVector<reco::Candidate> > collectionOne, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collectionTwo,bool doDileptonMassCut, edm::Handle<std::vector<pat::Muon> > musToAvoid){

#ifdef DEBUG
		  std::cout<<"checking pt of particles in two handled findLeadingAndSubleading fxn"<<std::endl;
#endif

		  ///find the highest pT object in collectionOne by looping over all contents in collectionOne
		  for(edm::OwnVector<reco::Candidate>::const_iterator genItOne = collectionOne->begin(); genItOne != collectionOne->end(); genItOne++){
#ifdef DEBUG
			  std::cout<<"a particle from collectionOne has pT = \t"<< genItOne->pt() << std::endl;
#endif

			  if(isCloseToMuon(genItOne, musToAvoid) ) continue;	///skip genItOne if it is close to an object in musToAvoid with pt>5
			  if(first==collectionOne->end()) first=genItOne;
			  else if(genItOne->pt() > first->pt() ) first = genItOne;
		  }//end loop over reco::Candidate objects in collectionOne

		  if(first==collectionOne->end()) return;	///don't look for a second lepton if a candidate for the first lepton is not found!

		  ///make sure the candidate chosen from collectionTwo has the opposite charge of the reco candidate ref named first
		  if(!doDileptonMassCut){
			  ///now find the highest pT object in collectionTwo which is separated from the highest pT object chosen from collectionOne
			  for(edm::OwnVector<reco::Candidate>::const_iterator genItTwo = collectionTwo->begin(); genItTwo != collectionTwo->end(); genItTwo++){
#ifdef DEBUG
				  std::cout<<"a particle from collectionTwo has pT = \t"<< genItTwo->pt() << std::endl;
#endif

				  if(second==collectionTwo->end() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > minDrSep && genItTwo->charge() != first->charge() ) second = genItTwo;

				  if(second!=collectionTwo->end()){

					  if(genItTwo->pt() > second->pt() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > minDrSep && genItTwo->charge() != first->charge()) second = genItTwo;
	
				  }///check that second points to a real reco::Candidate
			  
			  }///end loop over objects in collectionTwo

		  }///end if(!doDileptonMassCut)

		  if(doDileptonMassCut){
			  ///now find the highest pT object in collectionTwo which is separated from the highest pT object chosen from collectionOne
			  ///and whose dilepton mass, when combined with the object chosen from collectionOne, is above the dileptonMass threshold value
#ifdef DEBUG
			  std::cout<<"inside if(doDileptonMassCut) of two handled findLeadingAndSubleading() fxn"<<std::endl;
#endif

			  for(edm::OwnVector<reco::Candidate>::const_iterator genItTwo = collectionTwo->begin(); genItTwo != collectionTwo->end(); genItTwo++){
#ifdef DEBUG
				  std::cout<<"a particle from collectionTwo has pT = \t"<< genItTwo->pt() << std::endl;
#endif
			  
				  if(second==collectionTwo->end() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > minDrSep && getDileptonMass(first,genItTwo) > minDileptonMass && genItTwo->charge() != first->charge() ) second = genItTwo;

				  if(second!=collectionTwo->end() ){

					  if(genItTwo->pt() > second->pt() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > minDrSep && getDileptonMass(first,genItTwo) > minDileptonMass && genItTwo->charge() != first->charge() ) second = genItTwo;
	
				  }///check that second points to a real reco::Candidate
			  
			  }///end loop over objects in collectionTwo with dilepton mass cut applied

		  }///end if(doDileptonMassCut)

#ifdef DEBUG
		  std::cout<<"leaving two handled findLeadingAndSubleading fxn"<<std::endl;
#endif

	  }///end two handled findLeadingAndSubleading()


	  void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collection,bool doDileptonMassCut){

#ifdef DEBUG
		  std::cout<<"checking pt of particles in one handled findLeadingAndSubleading fxn"<<std::endl;
#endif

		  if(!doDileptonMassCut){

			  for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++){
#ifdef DEBUG
				  std::cout<<"pT = \t"<< genIt->pt() << std::endl;
#endif
				  if(first==collection->end()) first=genIt;
				  else{
					  if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > minDrSep ){
						  second = first;
						  first = genIt;
					  }
					  else if( (second==collection->end() || genIt->pt() > second->pt() ) && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > minDrSep  ) second = genIt;
				  }
			  }//end loop over reco::Candidate collection

		  }///end if(!doDileptonMassCut)

		  if(doDileptonMassCut){

			  for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++){
#ifdef DEBUG
				  std::cout<<"pT = \t"<< genIt->pt() << std::endl;
#endif
				  if(first==collection->end()) first=genIt;
				  else{
					  if(genIt->pt() > first->pt() && getDileptonMass(first,genIt) > minDileptonMass && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > minDrSep ){
						  second = first;
						  first = genIt;
					  }
					  else if( (second==collection->end() || genIt->pt() > second->pt()) && getDileptonMass(first,genIt) > minDileptonMass && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > minDrSep ) second = genIt;
				  }
			  }//end loop over reco::Candidate collection

		  }///end if(doDileptonMassCut)

#ifdef DEBUG
		  std::cout<<"leaving one handled findLeadingAndSubleading fxn"<<std::endl;
#endif

	  }///end one handled findLeadingAndSubleading()
	 

private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;



// ----------member data ---------------------------

std::string tName;
bool applyDileptonMassCut;
double minDileptonMass;
double minDrSep;

///Handles to RECO object collections
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > leptonsOne;
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > leptonsTwo;
edm::Handle<GenEventInfoProduct> genEvtInfo;
edm::Handle<std::vector<reco::Vertex> > vertices;
edm::Handle<std::vector<pat::Muon> > muonsToAvoid;

///tokens to input collections
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > leptonsOneToken;
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > leptonsTwoToken;
edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken;
edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;
edm::EDGetTokenT<std::vector<pat::Muon> > muonsToAvoidToken;


TTree * tree;

Int_t runNumber;
ULong64_t evtNumber;

Int_t nLeptons;
Int_t nLeptonsOne;
Int_t nLeptonsTwo;
Int_t nVertices;

///charge of the two selected leptons
Int_t chargeLeptonOne;
Int_t chargeLeptonTwo;

//first element is leading (highest pT) lepton 
//second element is subleading lepton 
Float_t etaEle[2];
Float_t ptEle[2];
Float_t phiEle[2];
Float_t dileptonMass;

Float_t dR_leadingLeptonSubleadingLepton;

///leadingIsHardest = 1 when leading lepton pT > subleading lepton pT
///leadingIsHardest = 0 when leading lepton pT < subleading lepton pT
Int_t leadingIsHardest;

Float_t evWeight;	///< weight of the event.  defaults to 1, only changes for bkgnd MC
Float_t evWeightSign;	///< if the sign of evWeight is negative, then the evt should not be plotted in a histo

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

emuAnalyzer::emuAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName")),
	applyDileptonMassCut(iConfig.getParameter<bool>("doDileptonMassCut")),
	minDileptonMass(iConfig.getParameter<double>("minDileptonMass")),
	minDrSep(iConfig.getParameter<double>("minDr"))

{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName.c_str(),"event kinematic info");

   tree->Branch("etaEle",etaEle,"etaEle[2]/F");
   tree->Branch("ptEle",ptEle,"ptEle[2]/F");
   tree->Branch("phiEle",phiEle,"phiEle[2]/F");
   tree->Branch("dileptonMass",&dileptonMass,"dileptonMass/F");
   tree->Branch("chargeLeptonOne",&chargeLeptonOne,"chargeLeptonOne/I");
   tree->Branch("chargeLeptonTwo",&chargeLeptonTwo,"chargeLeptonTwo/I");
   
   tree->Branch("leadingIsHardest",&leadingIsHardest,"leadingIsHardest/I");

   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
   tree->Branch("runNumber",&runNumber,"runNumber/I");

   tree->Branch("nLeptons",&nLeptons,"nLeptons/I");
   tree->Branch("nLeptonsOne",&nLeptonsOne,"nLeptonsOne/I");
   tree->Branch("nLeptonsTwo",&nLeptonsTwo,"nLeptonsTwo/I");
   tree->Branch("nVertices",&nVertices,"nVertices/I");

   tree->Branch("dR_leadingLeptonSubleadingLepton",&dR_leadingLeptonSubleadingLepton,"dR_leadingLeptonSubleadingLepton/F");

   tree->Branch("evWeight",&evWeight,"evWeight/F");
   tree->Branch("evWeightSign",&evWeightSign,"evWeightSign/F");
 
   leptonsOneToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("leptonsOneCollection"));///electrons
   leptonsTwoToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("leptonsTwoCollection"));///muons
   genEventInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   verticesToken = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));
   muonsToAvoidToken = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));

}


emuAnalyzer::~emuAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
emuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

#ifdef DEBUG
	std::cout<<"in analyze method of emuAnalyzer class"<<std::endl;
#endif

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();
	evWeight = 1.0;
	evWeightSign = 1.0;

	iEvent.getByToken(leptonsOneToken, leptonsOne);	///electrons
	iEvent.getByToken(leptonsTwoToken, leptonsTwo);	///muons
	iEvent.getByToken(verticesToken, vertices);
	iEvent.getByToken(muonsToAvoidToken, muonsToAvoid);

	iEvent.getByToken(genEventInfoToken, genEvtInfo);	///< get evt weights if analyzing MC

	if(genEvtInfo.isValid() ){
		///real data does not have gen lvl event weights
		evWeight = genEvtInfo->weight();
		if(evWeight < 0) evWeightSign = -1.0;
	}

	nLeptonsOne = leptonsOne->size();
	nLeptonsTwo = leptonsTwo->size();
	nLeptons = nLeptonsOne+nLeptonsTwo;
	nVertices = vertices->size();
	
	///assign iterators to both input leptons collections, and use these iterators to find the hardest and second hardest leptons in the evt
	///if one lepton is not found in both collections, then skip this evt 
	edm::OwnVector<reco::Candidate>::const_iterator leadingLeptonOne = leptonsOne->end(), leadingLeptonTwo = leptonsTwo->end();
	
	///now use the iterators to find the hardest and second hardest leptons
	///avoid leptons from leptonsOne which are within dR < 0.1 of an object in muonsToAvoid with pt>5
	findLeadingAndSubleading(leadingLeptonOne, leptonsOne, leadingLeptonTwo, leptonsTwo, applyDileptonMassCut, muonsToAvoid);
	
	if(leadingLeptonOne==leptonsOne->end() || leadingLeptonTwo==leptonsTwo->end()) return;	///< skip this evt if two leptons are not found

	///now declare the leadingLepton and subleadingLepton iterators which will be used below, and
	///assign them to the appropriate leadingLeptonOne and leadingLeptonTwo iterators
	///NOTE it is imperative that the electron and muon identities be preserved
	///for this reason the leadingLepton iterator will always be assigned to leadingLeptonOne
	///and the subleadingLepton iterator will always be assigned to leadingLeptonTwo

	edm::OwnVector<reco::Candidate>::const_iterator leadingLepton = leadingLeptonOne, subleadingLepton = leadingLeptonTwo;

	///now that the leading and subleading leptons have been found, fill all of the arrays and single Float_ values
	///which will be saved into the tree
#ifdef DEBUG
	std::cout<<"leading lepton pt: \t"<<leadingLepton->pt()<<std::endl;
	std::cout<<"subleading lepton pt: \t"<<subleadingLepton->pt()<<std::endl;
#endif

	chargeLeptonOne = leadingLepton->charge();
	chargeLeptonTwo = subleadingLepton->charge();
	etaEle[0] = leadingLepton->eta();
	ptEle[0] = leadingLepton->pt();
	phiEle[0] = leadingLepton->phi();
	etaEle[1] = subleadingLepton->eta();
	ptEle[1] = subleadingLepton->pt();
	phiEle[1] = subleadingLepton->phi();
	if(ptEle[0] > ptEle[1]){
		leadingIsHardest = 1;
	}
	else leadingIsHardest = 0;

	///use LorentzVector objects to calculate invariant mass
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > l1 = leadingLepton->p4(), l2 = subleadingLepton->p4();
	dileptonMass = (l1+l2).M();
	
#ifdef DEBUG
	std::cout<<"dilepton mass = \t"<< dileptonMass << std::endl;
	std::cout<<"\t"<<std::endl;
#endif

	///now use the individual object eta and phi values to calculate dR between different (lepton, jet) pairs
	dR_leadingLeptonSubleadingLepton = deltaR(etaEle[0],phiEle[0],etaEle[1],phiEle[1]);
	
	tree->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
emuAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
emuAnalyzer::endJob() 
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
emuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(emuAnalyzer);
