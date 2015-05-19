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
bool saveGenMatchedInfo;


///Matching = GEN objects
//Handles to GEN object collections
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > matchingLeadingLeptons;
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > matchingSubleadingLeptons;
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > matchingQuarks;

///Handles to RECO object collections
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > leadingLeptons;
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > subleadingLeptons;
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > quarks;

///tokens to input collections
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > matchingLeadingLeptonsToken;
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > matchingSubleadingLeptonsToken;
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > matchingQuarksToken;
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > leadingLeptonsToken;
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > subleadingLeptonsToken;
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > quarksToken;


TTree * tree;

Int_t runNumber;
ULong64_t evtNumber;

//first element is leading (highest pT) electron
//second element is subleading electron
Float_t etaMatchingEle[2];
Float_t ptMatchingEle[2];
Float_t phiMatchingEle[2];
Float_t dileptonMassMatching;

Float_t etaEle[2];
Float_t ptEle[2];
Float_t phiEle[2];
Float_t dileptonMass;

Float_t etaMatchingJet[2];
Float_t ptMatchingJet[2];
Float_t phiMatchingJet[2];
Float_t dijetMassMatching;

Float_t etaJet[2];
Float_t ptJet[2];
Float_t phiJet[2];
Float_t dijetMass;

///three and four object masses
Float_t fourObjectMassMatching;
Float_t leadLeptonThreeObjMassMatching;
Float_t subleadingLeptonThreeObjMassMatching;

Float_t fourObjectMass;
Float_t leadLeptonThreeObjMass;
Float_t subleadingLeptonThreeObjMass;

///deltaR between each lepton and both jets
Float_t dR_leadingLeptonLeadingJetMatching;
Float_t dR_leadingLeptonSubleadingJetMatching;
Float_t dR_subleadingLeptonLeadingJetMatching;
Float_t dR_subleadingLeptonSubleadingJetMatching;
Float_t dR_leadingLeptonSubleadingLeptonMatching;
Float_t dR_leadingJetSubleadingJetMatching;

Float_t dR_leadingLeptonLeadingJet;
Float_t dR_leadingLeptonSubleadingJet;
Float_t dR_subleadingLeptonLeadingJet;
Float_t dR_subleadingLeptonSubleadingJet;
Float_t dR_leadingLeptonSubleadingLepton;
Float_t dR_leadingJetSubleadingJet;


///leadingIsHardest = 1 when leading ele pT > subleading ele pT
///leadingIsHardest = 0 when leading ele pT < subleading ele pT
///leadingIsHardestMatching is for GEN electrons
Int_t leadingIsHardestMatching;
Int_t leadingIsHardest;

///pT, eta, and phi of the heavy Nu and WR
Float_t etaHvyNu;
Float_t ptHvyNu;
Float_t phiHvyNu;
Float_t etaWr;
Float_t ptWr;
Float_t phiWr;

Float_t etaHvyNuMatching;
Float_t ptHvyNuMatching;
Float_t phiHvyNuMatching;
Float_t etaWrMatching;
Float_t ptWrMatching;
Float_t phiWrMatching;


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
	saveGenMatchedInfo(iConfig.getParameter<bool>("saveGenMatched"))

{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName.c_str(),"event kinematic info");

   tree->Branch("etaMatchingEle",etaMatchingEle,"etaMatchingEle[2]/F");
   tree->Branch("ptMatchingEle",ptMatchingEle,"ptMatchingEle[2]/F");
   tree->Branch("phiMatchingEle",phiMatchingEle,"phiMatchingEle[2]/F");
   tree->Branch("dileptonMassMatching",&dileptonMassMatching,"dileptonMassMatching/F");
   
   tree->Branch("leadingIsHardestMatching",&leadingIsHardestMatching,"leadingIsHardestMatching/I");

   tree->Branch("etaEle",etaEle,"etaEle[2]/F");
   tree->Branch("ptEle",ptEle,"ptEle[2]/F");
   tree->Branch("phiEle",phiEle,"phiEle[2]/F");
   tree->Branch("dileptonMass",&dileptonMass,"dileptonMass/F");
   
   tree->Branch("leadingIsHardest",&leadingIsHardest,"leadingIsHardest/I");

   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
   tree->Branch("runNumber",&runNumber,"runNumber/I");

   tree->Branch("etaMatchingJet",etaMatchingJet,"etaMatchingJet[2]/F");
   tree->Branch("ptMatchingJet",ptMatchingJet,"ptMatchingJet[2]/F");
   tree->Branch("phiMatchingJet",phiMatchingJet,"phiMatchingJet[2]/F");
   tree->Branch("dijetMassMatching",&dijetMassMatching,"dijetMassMatching/F");

   tree->Branch("fourObjectMassMatching",&fourObjectMassMatching,"fourObjectMassMatching/F");
   tree->Branch("leadLeptonThreeObjMassMatching",&leadLeptonThreeObjMassMatching,"leadLeptonThreeObjMassMatching/F");
   tree->Branch("subleadingLeptonThreeObjMassMatching",&subleadingLeptonThreeObjMassMatching,"subleadingLeptonThreeObjMassMatching/F");
 
   tree->Branch("etaJet",etaJet,"etaJet[2]/F");
   tree->Branch("ptJet",ptJet,"ptJet[2]/F");
   tree->Branch("phiJet",phiJet,"phiJet[2]/F");
   tree->Branch("dijetMass",&dijetMass,"dijetMass/F");

   tree->Branch("fourObjectMass",&fourObjectMass,"fourObjectMass/F");
   tree->Branch("leadLeptonThreeObjMass",&leadLeptonThreeObjMass,"leadLeptonThreeObjMass/F");
   tree->Branch("subleadingLeptonThreeObjMass",&subleadingLeptonThreeObjMass,"subleadingLeptonThreeObjMass/F");
   
   tree->Branch("etaHvyNu",&etaHvyNu,"etaHvyNu/F");
   tree->Branch("ptHvyNu",&ptHvyNu,"ptHvyNu/F");
   tree->Branch("phiHvyNu",&phiHvyNu,"phiHvyNu/F");

   tree->Branch("etaWr",&etaWr,"etaWr/F");
   tree->Branch("ptWr",&ptWr,"ptWr/F");
   tree->Branch("phiWr",&phiWr,"phiWr/F");
   
   tree->Branch("etaHvyNuMatching",&etaHvyNuMatching,"etaHvyNuMatching/F");
   tree->Branch("ptHvyNuMatching",&ptHvyNuMatching,"ptHvyNuMatching/F");
   tree->Branch("phiHvyNuMatching",&phiHvyNuMatching,"phiHvyNuMatching/F");

   tree->Branch("etaWrMatching",&etaWrMatching,"etaWrMatching/F");
   tree->Branch("ptWrMatching",&ptWrMatching,"ptWrMatching/F");
   tree->Branch("phiWrMatching",&phiWrMatching,"phiWrMatching/F");

   tree->Branch("dR_leadingLeptonLeadingJetMatching",&dR_leadingLeptonLeadingJetMatching,"dR_leadingLeptonLeadingJetMatching/F");
   tree->Branch("dR_leadingLeptonSubleadingJetMatching",&dR_leadingLeptonSubleadingJetMatching,"dR_leadingLeptonSubleadingJetMatching/F");
   tree->Branch("dR_subleadingLeptonLeadingJetMatching",&dR_subleadingLeptonLeadingJetMatching,"dR_subleadingLeptonLeadingJetMatching/F");
   tree->Branch("dR_subleadingLeptonSubleadingJetMatching",&dR_subleadingLeptonSubleadingJetMatching,"dR_subleadingLeptonSubleadingJetMatching/F");
   tree->Branch("dR_leadingLeptonSubleadingLeptonMatching",&dR_leadingLeptonSubleadingLeptonMatching,"dR_leadingLeptonSubleadingLeptonMatching/F");
   tree->Branch("dR_leadingJetSubleadingJetMatching",&dR_leadingJetSubleadingJetMatching,"dR_leadingJetSubleadingJetMatching/F");
 
   tree->Branch("dR_leadingLeptonLeadingJet",&dR_leadingLeptonLeadingJet,"dR_leadingLeptonLeadingJet/F");
   tree->Branch("dR_leadingLeptonSubleadingJet",&dR_leadingLeptonSubleadingJet,"dR_leadingLeptonSubleadingJet/F");
   tree->Branch("dR_subleadingLeptonLeadingJet",&dR_subleadingLeptonLeadingJet,"dR_subleadingLeptonLeadingJet/F");
   tree->Branch("dR_subleadingLeptonSubleadingJet",&dR_subleadingLeptonSubleadingJet,"dR_subleadingLeptonSubleadingJet/F");
   tree->Branch("dR_leadingLeptonSubleadingLepton",&dR_leadingLeptonSubleadingLepton,"dR_leadingLeptonSubleadingLepton/F");
   tree->Branch("dR_leadingJetSubleadingJet",&dR_leadingJetSubleadingJet,"dR_leadingJetSubleadingJet/F");

   if(saveGenMatchedInfo){
	   ///setup the tokens to the matching (GEN) collections IF matching object info (pT, eta, etc) should be saved
	   matchingLeadingLeptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("matchingLeadingLeptonCollection")); 
	   matchingSubleadingLeptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("matchingSubleadingLeptonCollection"));
	   matchingQuarksToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("matchingQuarkCollection")); 
   }///end if(saveGenMatchedInfo)
   
   leadingLeptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("leadingLeptonCollection")); 
   subleadingLeptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("subleadingLeptonCollection"));
   quarksToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("quarkCollection")); 
 
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

	iEvent.getByToken(leadingLeptonsToken, leadingLeptons);
	iEvent.getByToken(subleadingLeptonsToken, subleadingLeptons);
	iEvent.getByToken(quarksToken, quarks);

	if(saveGenMatchedInfo){
		iEvent.getByToken(matchingLeadingLeptonsToken, matchingLeadingLeptons);
		iEvent.getByToken(matchingSubleadingLeptonsToken, matchingSubleadingLeptons);
		iEvent.getByToken(matchingQuarksToken, matchingQuarks);

		edm::OwnVector<reco::Candidate>::const_iterator matchingLeadingLepton = matchingLeadingLeptons->end(), matchingSubleadingLepton = matchingSubleadingLeptons->end();
		edm::OwnVector<reco::Candidate>::const_iterator matchingLeadingJet = matchingQuarks->end(), matchingSubleadingJet = matchingQuarks->end();

		findHighestPt(matchingLeadingLepton, matchingLeadingLeptons);
		findHighestPt(matchingSubleadingLepton, matchingSubleadingLeptons);
		findLeadingAndSubleading(matchingLeadingJet, matchingSubleadingJet, matchingQuarks);

		etaMatchingEle[0] = matchingLeadingLepton->eta();
		ptMatchingEle[0] = matchingLeadingLepton->pt();
		phiMatchingEle[0] = matchingLeadingLepton->phi();
		etaMatchingEle[1] = matchingSubleadingLepton->eta();
		ptMatchingEle[1] = matchingSubleadingLepton->pt();
		phiMatchingEle[1] = matchingSubleadingLepton->phi();
		if(ptMatchingEle[0] > ptMatchingEle[1]){
			leadingIsHardestMatching = 1;
		}
		else leadingIsHardestMatching = 0;


#ifdef DEBUG
		std::cout<<"matching leading jet pt: \t"<< matchingLeadingJet->pt()<<std::endl;
		std::cout<<"matchingSubleading jet pt: \t"<<matchingSubleadingJet->pt()<<std::endl;
#endif

		etaMatchingJet[0] = matchingLeadingJet->eta();
		ptMatchingJet[0] = matchingLeadingJet->pt();
		phiMatchingJet[0] = matchingLeadingJet->phi();
		etaMatchingJet[1] = matchingSubleadingJet->eta();
		ptMatchingJet[1] = matchingSubleadingJet->pt();
		phiMatchingJet[1] = matchingSubleadingJet->phi();

		///make TLorentzVector objects for the four GEN objects.  Then use these LorentzVectors to calculate three and four object
		///invariant mass values.
		TLorentzVector l1(ptMatchingEle[0]*TMath::Cos(phiMatchingEle[0]),ptMatchingEle[0]*TMath::Sin(phiMatchingEle[0]),ptMatchingEle[0]*TMath::SinH(etaMatchingEle[0]),ptMatchingEle[0]*TMath::CosH(etaMatchingEle[0]));	///leading lepton lorentz vector (px, py, pz, E)
		TLorentzVector l2(ptMatchingEle[1]*TMath::Cos(phiMatchingEle[1]),ptMatchingEle[1]*TMath::Sin(phiMatchingEle[1]),ptMatchingEle[1]*TMath::SinH(etaMatchingEle[1]),ptMatchingEle[1]*TMath::CosH(etaMatchingEle[1]));
		TLorentzVector j1(ptMatchingJet[0]*TMath::Cos(phiMatchingJet[0]),ptMatchingJet[0]*TMath::Sin(phiMatchingJet[0]),ptMatchingJet[0]*TMath::SinH(etaMatchingJet[0]),ptMatchingJet[0]*TMath::CosH(etaMatchingJet[0]));	///leading jet lorentz vector (px, py, pz, E)
		TLorentzVector j2(ptMatchingJet[1]*TMath::Cos(phiMatchingJet[1]),ptMatchingJet[1]*TMath::Sin(phiMatchingJet[1]),ptMatchingJet[1]*TMath::SinH(etaMatchingJet[1]),ptMatchingJet[1]*TMath::CosH(etaMatchingJet[1]));

		fourObjectMassMatching = (l1+l2+j1+j2).M();
		if(applyFourObjMassCut && (fourObjectMassMatching < fourObjMassCutVal) ) return;	///don't add an entry to the tree if this evt fails the cut

		dileptonMassMatching = (l1+l2).M();
		dijetMassMatching = (j1+j2).M();
		leadLeptonThreeObjMassMatching = (l1+j1+j2).M();
		subleadingLeptonThreeObjMassMatching = (l2+j1+j2).M();

		etaHvyNuMatching = (l2+j1+j2).Eta();
		ptHvyNuMatching = (l2+j1+j2).Pt();
		phiHvyNuMatching = (l2+j1+j2).Phi();
		etaWrMatching = (l1+l2+j1+j2).Eta();
		ptWrMatching = (l1+l2+j1+j2).Pt();
		phiWrMatching = (l1+l2+j1+j2).Phi();

#ifdef DEBUG
		std::cout<<"matching dilepton mass = \t"<< dileptonMassMatching << std::endl;
		std::cout<<"matching dijet mass = \t"<< dijetMassMatching << std::endl;
		std::cout<<"\t"<<std::endl;
#endif

		///now use the individual GEN object eta and phi values to calculate dR between different (lepton, jet) pairs
		dR_leadingLeptonLeadingJetMatching = deltaR(etaMatchingEle[0],phiMatchingEle[0],etaMatchingJet[0],phiMatchingJet[0]);
		dR_leadingLeptonSubleadingJetMatching = deltaR(etaMatchingEle[0],phiMatchingEle[0],etaMatchingJet[1],phiMatchingJet[1]);
		dR_subleadingLeptonLeadingJetMatching = deltaR(etaMatchingEle[1],phiMatchingEle[1],etaMatchingJet[0],phiMatchingJet[0]);
		dR_subleadingLeptonSubleadingJetMatching = deltaR(etaMatchingEle[1],phiMatchingEle[1],etaMatchingJet[1],phiMatchingJet[1]);
		dR_leadingLeptonSubleadingLeptonMatching = deltaR(etaMatchingEle[0],phiMatchingEle[0],etaMatchingEle[1],phiMatchingEle[1]);
		dR_leadingJetSubleadingJetMatching = deltaR(etaMatchingJet[0],phiMatchingJet[0],etaMatchingJet[1],phiMatchingJet[1]);

	}///end if(saveGenMatchedInfo)

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


#ifdef DEBUG
	std::cout<<"leading jet pt: \t"<<leadingJet->pt()<<std::endl;
	std::cout<<"subleading jet pt: \t"<<subleadingJet->pt()<<std::endl;
#endif

	etaJet[0] = leadingJet->eta();
	ptJet[0] = leadingJet->pt();
	phiJet[0] = leadingJet->phi();
	etaJet[1] = subleadingJet->eta();
	ptJet[1] = subleadingJet->pt();
	phiJet[1] = subleadingJet->phi();

	///make TLorentzVector objects for the four GEN objects.  Then use these LorentzVectors to calculate three and four object
	///invariant mass values.
	TLorentzVector l1(ptEle[0]*TMath::Cos(phiEle[0]),ptEle[0]*TMath::Sin(phiEle[0]),ptEle[0]*TMath::SinH(etaEle[0]),ptEle[0]*TMath::CosH(etaEle[0]));	///leading lepton lorentz vector (px, py, pz, E)
	TLorentzVector l2(ptEle[1]*TMath::Cos(phiEle[1]),ptEle[1]*TMath::Sin(phiEle[1]),ptEle[1]*TMath::SinH(etaEle[1]),ptEle[1]*TMath::CosH(etaEle[1]));
	TLorentzVector j1(ptJet[0]*TMath::Cos(phiJet[0]),ptJet[0]*TMath::Sin(phiJet[0]),ptJet[0]*TMath::SinH(etaJet[0]),ptJet[0]*TMath::CosH(etaJet[0]));	///leading jet lorentz vector (px, py, pz, E)
	TLorentzVector j2(ptJet[1]*TMath::Cos(phiJet[1]),ptJet[1]*TMath::Sin(phiJet[1]),ptJet[1]*TMath::SinH(etaJet[1]),ptJet[1]*TMath::CosH(etaJet[1]));

	fourObjectMass = (l1+l2+j1+j2).M();
	if(applyFourObjMassCut && (fourObjectMass < fourObjMassCutVal) ) return;	///don't add an entry to the tree if this evt fails the cut
	
	dileptonMass = (l1+l2).M();
	dijetMass = (j1+j2).M();
	leadLeptonThreeObjMass = (l1+j1+j2).M();
	subleadingLeptonThreeObjMass = (l2+j1+j2).M();

	etaHvyNu = (l2+j1+j2).Eta();
	ptHvyNu = (l2+j1+j2).Pt();
	phiHvyNu = (l2+j1+j2).Phi();
	etaWr = (l1+l2+j1+j2).Eta();
	ptWr = (l1+l2+j1+j2).Pt();
	phiWr = (l1+l2+j1+j2).Phi();

#ifdef DEBUG
	std::cout<<"dilepton mass = \t"<< dileptonMass << std::endl;
	std::cout<<"dijet mass = \t"<< dijetMass << std::endl;
	std::cout<<"\t"<<std::endl;
#endif

	///now use the individual GEN object eta and phi values to calculate dR between different (lepton, jet) pairs
	dR_leadingLeptonLeadingJet = deltaR(etaEle[0],phiEle[0],etaJet[0],phiJet[0]);
	dR_leadingLeptonSubleadingJet = deltaR(etaEle[0],phiEle[0],etaJet[1],phiJet[1]);
	dR_subleadingLeptonLeadingJet = deltaR(etaEle[1],phiEle[1],etaJet[0],phiJet[0]);
	dR_subleadingLeptonSubleadingJet = deltaR(etaEle[1],phiEle[1],etaJet[1],phiJet[1]);
	dR_leadingLeptonSubleadingLepton = deltaR(etaEle[0],phiEle[0],etaEle[1],phiEle[1]);
	dR_leadingJetSubleadingJet = deltaR(etaJet[0],phiJet[0],etaJet[1],phiJet[1]);
	
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
