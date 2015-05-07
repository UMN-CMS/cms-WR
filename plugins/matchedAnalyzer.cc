// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/matchedAnalyzer
// Class:      matchedAnalyzer
// 
/**\class matchedAnalyzer matchedAnalyzer.cc doubleElectronTracklessTrigger/matchedAnalyzer/plugins/matchedAnalyzer.cc

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
#include "TLorentzVector.h"
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

class matchedAnalyzer : public edm::EDAnalyzer {
   public:
      explicit matchedAnalyzer(const edm::ParameterSet&);
      ~matchedAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	  ///use this fxn to find the two highest pT jets in the evt which are both at least dR >= 0.4 away from the two leading leptons in the evt
	  ///input params: collection of hadrons (jets, quarks, GEN or RECO lvl), const_iterators to the two leading leptons, and const_iterators
	  ///to two objects in the collection of hadrons
	  ///the two hadron iterators will be updated by this fxn
	  void findFarHadrons(edm::Handle<edm::OwnVector<reco::Candidate> > hadronColl, edm::OwnVector<reco::Candidate>::const_iterator& leadingLept, edm::OwnVector<reco::Candidate>::const_iterator& subleadingLept, edm::OwnVector<reco::Candidate>::const_iterator& hadronOne, edm::OwnVector<reco::Candidate>::const_iterator& hadronTwo){
#ifdef DEBUG
		  std::cout<<"entered findFarHadrons function"<<std::endl;
#endif
		  for(edm::OwnVector<reco::Candidate>::const_iterator hadIt = hadronColl->begin(); hadIt!=hadronColl->end(); hadIt++){
			  double drOne = deltaR(leadingLept->eta(),leadingLept->phi(),hadIt->eta(),hadIt->phi());
			  double drTwo = deltaR(subleadingLept->eta(),subleadingLept->phi(),hadIt->eta(),hadIt->phi());
			  if(drOne >= minDr && drTwo >= minDr){
				  if(hadronOne==hadronColl->end()) hadronOne=hadIt;
				  else{
					  if(hadIt->pt() > hadronOne->pt()){
						  hadronTwo = hadronOne;
						  hadronOne = hadIt;
					  }
					  else if(hadronTwo==hadronColl->end() || hadIt->pt() > hadronTwo->pt()) hadronTwo = hadIt;
				  }///end else
			  }///end if(current iterator hadIt points to a hadron which is well separated from the two leading leptons
		  }///end loop over objects in hadron collection
#ifdef DEBUG
		  std::cout<<"leaving findFarHadrons fxn"<<std::endl;
#endif

	  }

	  void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collection){

#ifdef DEBUG
		  std::cout<<"checking pt of all gen particles in findLeadingAndSubleading fxn"<<std::endl;
#endif

		  for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++){
#ifdef DEBUG
			  std::cout<<"pT = \t"<< genIt->pt() << std::endl;
#endif
			  if(first==collection->end()) first=genIt;
			  else{
				  if(genIt->pt() > first->pt()){
					  second = first;
					  first = genIt;
				  }
				  else if(second==collection->end() || genIt->pt() > second->pt()) second = genIt;
			  }
		  }//end loop over reco::Candidate collection

#ifdef DEBUG
		  std::cout<<"leaving findLeadingAndSubleading fxn"<<std::endl;
#endif

	  }///end findLeadingAndSubleading()

	  /** this fxn sifts through a collection of reco::Candidate objects, finds the highest pT object in the collection, and assigns
	   * a pointer to this object to the iterator named iter
	   * a reference to iter is input to this fxn
	   */
	  void findHighestPt(edm::OwnVector<reco::Candidate>::const_iterator& iter, edm::Handle<edm::OwnVector<reco::Candidate> > coll){
		  for(edm::OwnVector<reco::Candidate>::const_iterator it = coll->begin(); it!=coll->end(); it++){
			  if(iter==coll->end()) iter=it;
			  else{
				  if(it->pt() > iter->pt()) iter=it;
			  }
		  }///end loop over reco::Candidate objects in the collection named coll

	  }///end findHighestPt()
    

private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;


//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

// ----------member data ---------------------------

std::string tName;
bool applyDrCut;
bool applyFourObjMassCut;
double fourObjMassCutVal;
double minDr;

edm::InputTag genLeadingLeptonCollTag;
edm::InputTag genSubleadingLeptonCollTag;
edm::InputTag genQuarkCollTag;

//Handles to GEN object collections
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > leadingLeptons;
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > subleadingLeptons;
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > quarks;

TTree * tree;

Int_t runNumber;
ULong64_t evtNumber;

//first element is leading (highest pT) electron
//second element is subleading electron
Float_t etaGenEle[2];
Float_t ptGenEle[2];
Float_t phiGenEle[2];
Float_t dileptonMassGen;

Float_t etaGenJet[2];
Float_t ptGenJet[2];
Float_t phiGenJet[2];
Float_t dijetMassGen;

///three and four object masses
Float_t fourObjectMassGen;
Float_t leadLeptonThreeObjMassGen;
Float_t subleadingLeptonThreeObjMassGen;

///deltaR between each lepton and both jets
Float_t dR_leadingLeptonLeadingJetGen;
Float_t dR_leadingLeptonSubleadingJetGen;
Float_t dR_subleadingLeptonLeadingJetGen;
Float_t dR_subleadingLeptonSubleadingJetGen;
Float_t dR_leadingLeptonSubleadingLeptonGen;
Float_t dR_leadingJetSubleadingJetGen;

///leadingIsHardest = 1 when leading ele pT > subleading ele pT
///leadingIsHardest = 0 when leading ele pT < subleading ele pT
Int_t leadingIsHardest;

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

matchedAnalyzer::matchedAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName")),
	applyDrCut(iConfig.getParameter<bool>("doDeltaRcut")),
	applyFourObjMassCut(iConfig.getParameter<bool>("doFourObjMassCut")),
	fourObjMassCutVal(iConfig.getParameter<double>("minFourObjMass")),
	minDr(iConfig.getParameter<double>("minDeltaRforLeptonJetExclusion")),
	genLeadingLeptonCollTag(iConfig.getParameter<edm::InputTag>("genLeadingLeptonCollection")),
	genSubleadingLeptonCollTag(iConfig.getParameter<edm::InputTag>("genSubleadingLeptonCollection")),
	genQuarkCollTag(iConfig.getParameter<edm::InputTag>("genQuarkCollection"))

{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName.c_str(),"event kinematic info");

   tree->Branch("etaGenEle",etaGenEle,"etaGenEle[2]/F");
   tree->Branch("ptGenEle",ptGenEle,"ptGenEle[2]/F");
   tree->Branch("phiGenEle",phiGenEle,"phiGenEle[2]/F");
   tree->Branch("dileptonMassGen",&dileptonMassGen,"dileptonMassGen/F");
   
   tree->Branch("leadingIsHardest",&leadingIsHardest,"leadingIsHardest/I");

   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
   tree->Branch("runNumber",&runNumber,"runNumber/I");

   tree->Branch("etaGenJet",etaGenJet,"etaGenJet[2]/F");
   tree->Branch("ptGenJet",ptGenJet,"ptGenJet[2]/F");
   tree->Branch("phiGenJet",phiGenJet,"phiGenJet[2]/F");
   tree->Branch("dijetMassGen",&dijetMassGen,"dijetMassGen/F");

   tree->Branch("fourObjectMassGen",&fourObjectMassGen,"fourObjectMassGen/F");
   tree->Branch("leadLeptonThreeObjMassGen",&leadLeptonThreeObjMassGen,"leadLeptonThreeObjMassGen/F");
   tree->Branch("subleadingLeptonThreeObjMassGen",&subleadingLeptonThreeObjMassGen,"subleadingLeptonThreeObjMassGen/F");

   tree->Branch("dR_leadingLeptonLeadingJetGen",&dR_leadingLeptonLeadingJetGen,"dR_leadingLeptonLeadingJetGen/F");
   tree->Branch("dR_leadingLeptonSubleadingJetGen",&dR_leadingLeptonSubleadingJetGen,"dR_leadingLeptonSubleadingJetGen/F");
   tree->Branch("dR_subleadingLeptonLeadingJetGen",&dR_subleadingLeptonLeadingJetGen,"dR_subleadingLeptonLeadingJetGen/F");
   tree->Branch("dR_subleadingLeptonSubleadingJetGen",&dR_subleadingLeptonSubleadingJetGen,"dR_subleadingLeptonSubleadingJetGen/F");
   tree->Branch("dR_leadingLeptonSubleadingLeptonGen",&dR_leadingLeptonSubleadingLeptonGen,"dR_leadingLeptonSubleadingLeptonGen/F");
   tree->Branch("dR_leadingJetSubleadingJetGen",&dR_leadingJetSubleadingJetGen,"dR_leadingJetSubleadingJetGen/F");
 
}


matchedAnalyzer::~matchedAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
matchedAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

#ifdef DEBUG
	std::cout<<"in analyze method of matchedAnalyzer class"<<std::endl;
#endif

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();

	iEvent.getByLabel(genLeadingLeptonCollTag, leadingLeptons);
	iEvent.getByLabel(genSubleadingLeptonCollTag, subleadingLeptons);
	iEvent.getByLabel(genQuarkCollTag, quarks);

	edm::OwnVector<reco::Candidate>::const_iterator leadingLepton = leadingLeptons->end(), subleadingLepton = subleadingLeptons->end();
	edm::OwnVector<reco::Candidate>::const_iterator leadingJet = quarks->end(), subleadingJet = quarks->end();
	   

	findHighestPt(leadingLepton, leadingLeptons);
	findHighestPt(subleadingLepton, subleadingLeptons);
	if(!applyDrCut) findLeadingAndSubleading(leadingJet, subleadingJet, quarks);

	if(applyDrCut){
		findFarHadrons(quarks,leadingLepton,subleadingLepton,leadingJet,subleadingJet);
		///skip this evt if the dR cut is applied and no two jet candidates
		///are found outside dR=minDr away from the two leading leptons
		if(leadingJet == quarks->end() || subleadingJet == quarks->end()) return;
	}///end if(apply dR cut)

	///now that the leading and subleading leptons and quarks have been found, fill all of the arrays and single Float_ values
	///which will be saved into the tree
#ifdef DEBUG
	std::cout<<"leading lepton pt: \t"<<leadingLepton->pt()<<std::endl;
	std::cout<<"subleading lepton pt: \t"<<subleadingLepton->pt()<<std::endl;
#endif

	etaGenEle[0] = leadingLepton->eta();
	ptGenEle[0] = leadingLepton->pt();
	phiGenEle[0] = leadingLepton->phi();
	etaGenEle[1] = subleadingLepton->eta();
	ptGenEle[1] = subleadingLepton->pt();
	phiGenEle[1] = subleadingLepton->phi();
	if(ptGenEle[0] > ptGenEle[1]){
		leadingIsHardest = 1;
	}
	else leadingIsHardest = 0;


#ifdef DEBUG
	std::cout<<"leading jet pt: \t"<<leadingJet->pt()<<std::endl;
	std::cout<<"subleading jet pt: \t"<<subleadingJet->pt()<<std::endl;
#endif

	etaGenJet[0] = leadingJet->eta();
	ptGenJet[0] = leadingJet->pt();
	phiGenJet[0] = leadingJet->phi();
	etaGenJet[1] = subleadingJet->eta();
	ptGenJet[1] = subleadingJet->pt();
	phiGenJet[1] = subleadingJet->phi();

	///make TLorentzVector objects for the four GEN objects.  Then use these LorentzVectors to calculate three and four object
	///invariant mass values.
	TLorentzVector l1(ptGenEle[0]*TMath::Cos(phiGenEle[0]),ptGenEle[0]*TMath::Sin(phiGenEle[0]),ptGenEle[0]*TMath::SinH(etaGenEle[0]),ptGenEle[0]*TMath::CosH(etaGenEle[0]));	///leading lepton lorentz vector (px, py, pz, E)
	TLorentzVector l2(ptGenEle[1]*TMath::Cos(phiGenEle[1]),ptGenEle[1]*TMath::Sin(phiGenEle[1]),ptGenEle[1]*TMath::SinH(etaGenEle[1]),ptGenEle[1]*TMath::CosH(etaGenEle[1]));
	TLorentzVector j1(ptGenJet[0]*TMath::Cos(phiGenJet[0]),ptGenJet[0]*TMath::Sin(phiGenJet[0]),ptGenJet[0]*TMath::SinH(etaGenJet[0]),ptGenJet[0]*TMath::CosH(etaGenJet[0]));	///leading jet lorentz vector (px, py, pz, E)
	TLorentzVector j2(ptGenJet[1]*TMath::Cos(phiGenJet[1]),ptGenJet[1]*TMath::Sin(phiGenJet[1]),ptGenJet[1]*TMath::SinH(etaGenJet[1]),ptGenJet[1]*TMath::CosH(etaGenJet[1]));

	fourObjectMassGen = (l1+l2+j1+j2).M();
	if(applyFourObjMassCut && (fourObjectMassGen < fourObjMassCutVal) ) return;	///don't add an entry to the tree if this evt fails the cut
	
	dileptonMassGen = (l1+l2).M();
	dijetMassGen = (j1+j2).M();
	leadLeptonThreeObjMassGen = (l1+j1+j2).M();
	subleadingLeptonThreeObjMassGen = (l2+j1+j2).M();

#ifdef DEBUG
	std::cout<<"dilepton mass = \t"<< dileptonMassGen << std::endl;
	std::cout<<"dijet mass = \t"<< dijetMassGen << std::endl;
	std::cout<<"\t"<<std::endl;
#endif

	///now use the individual GEN object eta and phi values to calculate dR between different (lepton, jet) pairs
	dR_leadingLeptonLeadingJetGen = deltaR(etaGenEle[0],phiGenEle[0],etaGenJet[0],phiGenJet[0]);
	dR_leadingLeptonSubleadingJetGen = deltaR(etaGenEle[0],phiGenEle[0],etaGenJet[1],phiGenJet[1]);
	dR_subleadingLeptonLeadingJetGen = deltaR(etaGenEle[1],phiGenEle[1],etaGenJet[0],phiGenJet[0]);
	dR_subleadingLeptonSubleadingJetGen = deltaR(etaGenEle[1],phiGenEle[1],etaGenJet[1],phiGenJet[1]);
	dR_leadingLeptonSubleadingLeptonGen = deltaR(etaGenEle[0],phiGenEle[0],etaGenEle[1],phiGenEle[1]);
	dR_leadingJetSubleadingJetGen = deltaR(etaGenJet[0],phiGenJet[0],etaGenJet[1],phiGenJet[1]);
	
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
matchedAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
matchedAnalyzer::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
/*
void 
matchedAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
matchedAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
matchedAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
matchedAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
matchedAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(matchedAnalyzer);
