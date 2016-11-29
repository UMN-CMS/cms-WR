// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/generalPurposeGenAndRecoAnalyzer
// Class:      generalPurposeGenAndRecoAnalyzer
//
/**\class generalPurposeGenAndRecoAnalyzer generalPurposeGenAndRecoAnalyzer.cc doubleElectronTracklessTrigger/generalPurposeGenAndRecoAnalyzer/plugins/generalPurposeGenAndRecoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
	 
	 this analyzer is used for both GEN leptons, GEN quarks, GEN jets, and reco leptons and jets
	 any branch or variable denoted as "Quark" corresponds to a hadron, which could be at GEN or RECO lvl
	 no dR, PDGID, or status based matching between particles is done here
	 this analyzer simply takes one lepton collection, one hadron collection (qrks or jets), finds the
	 leading two leptons and hadrons (qrks or jets), and saves the kinematic info of these individual
	 objects, and the kinematic properties of multi particle systems (dilepton, dilepton dihadron, etc)
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
//#include "DataFormats/EgammaCandidates/interface/Electron.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
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

class generalPurposeGenAndRecoAnalyzer : public edm::EDAnalyzer
{
public:
	explicit generalPurposeGenAndRecoAnalyzer(const edm::ParameterSet&);
	~generalPurposeGenAndRecoAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	/**calculate diobject mass, return zero if result is 0 or imaginary
	 */
	Float_t diobjectMass(edm::OwnVector<reco::Candidate>::const_iterator objOne, edm::OwnVector<reco::Candidate>::const_iterator objTwo){
		Float_t result = 0.;
		Float_t massSqd = 2*(objOne->pt())*(objTwo->pt())*(TMath::CosH( (objOne->eta()) - (objTwo->eta()) ) - TMath::Cos( (objOne->phi()) - (objTwo->phi()) ) );
		if(massSqd > 0.) result = TMath::Sqrt(massSqd);
		return result;
	}

	/** this fxn sifts through a collection of reco::Candidate objects, finds the highest pT object in the collection, and assigns
	 * a pointer to this object to the iterator named iter
	 * a reference to iter is input to this fxn
	 */
	void findHighestPt(edm::OwnVector<reco::Candidate>::const_iterator& iter, edm::Handle<edm::OwnVector<reco::Candidate> > coll)
	{
		for(edm::OwnVector<reco::Candidate>::const_iterator it = coll->begin(); it != coll->end(); it++) {
			if(iter == coll->end()) iter = it;
			else {
				if(it->pt() > iter->pt()) iter = it;
			}
		}///end loop over reco::Candidate objects in the collection named coll

	}///end findHighestPt()

	void findLeadingAndSubleadingVectorInput(edm::OwnVector<reco::Candidate>::iterator& first, edm::OwnVector<reco::Candidate>::iterator& second, edm::OwnVector<reco::Candidate> collection)
	{
#ifdef DEBUG
		std::cout << "in findLeadingAndSubleadingVectorInput()" << std::endl;
		std::cout << "num elements in input collection =\t" << collection.size() << std::endl;
#endif

	}///end findLeadingAndSubleadingVectorInput()

	void findLeadingAndSubleadingJet(std::vector<pat::Jet>::const_iterator& first, std::vector<pat::Jet>::const_iterator& second, edm::Handle<std::vector<pat::Jet> > collection)
	{

		for(std::vector<pat::Jet>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++) {
			if(first == collection->end()) first = genIt;
			else {
				if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1 ) {
					second = first;
					first = genIt;
				} else if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) <= 0.1 ) first = genIt;
				else if( (second == collection->end() || genIt->pt() > second->pt() ) && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1  ) second = genIt;
			}
		}//end loop over reco::Candidate collection

	}///end findLeadingAndSubleadingJet()


	void findLeadingAndSubleadingLepton(std::vector<pat::Electron>::const_iterator& first, std::vector<pat::Electron>::const_iterator& second, std::vector<pat::Electron> collection)
	{

		for(std::vector<pat::Electron>::const_iterator genIt = collection.begin(); genIt != collection.end(); genIt++) {
			if(first == collection.end()) first = genIt;
			else {
				if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1 ) {
					second = first;
					first = genIt;
				} else if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) <= 0.1 ) first = genIt;
				else if( (second == collection.end() || genIt->pt() > second->pt() ) && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1  ) second = genIt;
			}
		}//end loop over reco::Candidate collection

	}///end findLeadingAndSubleading()



	void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collection)
	{

		for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++) {
			if(first == collection->end()) first = genIt;
			else {
				if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1 ) {
					second = first;
					first = genIt;
				} else if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) <= 0.1 ) first = genIt;
				else if( (second == collection->end() || genIt->pt() > second->pt() ) && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1  ) second = genIt;
			}
		}//end loop over reco::Candidate collection

	}///end findLeadingAndSubleading()

	void resetTreeVars()
	{
		runNumber = -1, evtNumber = -1, numLeptons = -1, numQuarks = -1;
		evWeight = 1, evWeightSign = 1;

		motherPdgIdLeptOne = -1, motherStatusLeptOne = -1, motherPdgIdLeptTwo = -1, motherStatusLeptTwo = -1, motherPdgIdQuarkOne = -1, motherStatusQuarkOne = -1, motherPdgIdQuarkTwo = -1, motherStatusQuarkTwo = -1;

		etaLeptOne = -9, ptLeptOne = -9, phiLeptOne = -9, etaLeptTwo = -9, ptLeptTwo = -9, phiLeptTwo = -9, etaQuarkOne = -9, ptQuarkOne = -9, phiQuarkOne = -9, etaQuarkTwo = -9, ptQuarkTwo = -9, phiQuarkTwo = -9;

		dileptonMass = -9;
		dileptonPt = -9;
		dileptonEta = -9;
		dileptonPhi = -9;
		
		subleadLeptonBothHadronsMass = -9, subleadLeptonBothHadronsPt = -9, subleadLeptonBothHadronsEta = -9, subleadLeptonBothHadronsPhi = -9;
		leadLeptonBothHadronsMass = -9, leadLeptonBothHadronsPt = -9, leadLeptonBothHadronsEta = -9, leadLeptonBothHadronsPhi = -9;
		dileptonDihadronMass = -9, dileptonDihadronPt = -9, dileptonDihadronEta = -9, dileptonDihadronPhi = -9;

		dRleptonOneQuarkOne = -9;
		dRleptonOneQuarkTwo = -9;
		dRleptonTwoQuarkOne = -9;
		dRleptonTwoQuarkTwo = -9;

	}///end resetTreeVars()

private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

// ----------member data ---------------------------

	std::string _tName;
	bool processingGen, processingQuarks;

///Handles to MINIAOD object collections
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > leptonCollection;
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > hadronCollection;
	edm::Handle<GenEventInfoProduct> genEvtInfo;	///<gen event weight

///tokens to input collections
	edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > leptonCollectionToken;
	edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > hadronCollectionToken;
	edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken;

	TTree * tree;

	Int_t runNumber;
	ULong64_t evtNumber;

	Float_t evWeight;	///< gen weight of the event.  defaults to 1
	Float_t evWeightSign;	///< evtWeightSign used to determine weight of the event when filling a histo

///pt, eta, phi of leptons and jets or quarks
	Float_t etaLeptOne;
	Float_t ptLeptOne;
	Float_t phiLeptOne;
	Float_t etaLeptTwo;
	Float_t ptLeptTwo;
	Float_t phiLeptTwo;
	Float_t etaQuarkOne;
	Float_t ptQuarkOne;
	Float_t phiQuarkOne;
	Float_t etaQuarkTwo;
	Float_t ptQuarkTwo;
	Float_t phiQuarkTwo;

///only for GEN leptons and quarks
	Int_t motherPdgIdLeptOne;
	Int_t motherStatusLeptOne;
	Int_t motherPdgIdLeptTwo;
	Int_t motherStatusLeptTwo;
	Int_t motherPdgIdQuarkOne;
	Int_t motherStatusQuarkOne;
	Int_t motherPdgIdQuarkTwo;
	Int_t motherStatusQuarkTwo;

///pt, eta, phi, and mass of objects made from two or more GEN leptons or jets or quarks	
///dilepton system
	Float_t dileptonMass;
	Float_t dileptonPt;
	Float_t dileptonEta;
	Float_t dileptonPhi;

///subleading lepton and two leading hadron system
	Float_t subleadLeptonBothHadronsMass;
	Float_t subleadLeptonBothHadronsPt;
	Float_t subleadLeptonBothHadronsEta;
	Float_t subleadLeptonBothHadronsPhi;

///leading lepton and two leading hadron system
	Float_t leadLeptonBothHadronsMass;
	Float_t leadLeptonBothHadronsPt;
	Float_t leadLeptonBothHadronsEta;
	Float_t leadLeptonBothHadronsPhi;

///dilepton dihadron system
	Float_t dileptonDihadronMass;
	Float_t dileptonDihadronPt;
	Float_t dileptonDihadronEta;
	Float_t dileptonDihadronPhi;

///dR between leptons and hadrons
	Float_t dRleptonOneQuarkOne;
	Float_t dRleptonOneQuarkTwo;
	Float_t dRleptonTwoQuarkOne;
	Float_t dRleptonTwoQuarkTwo;

///num leptons and hadrons
	Int_t numLeptons, numQuarks;

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

generalPurposeGenAndRecoAnalyzer::generalPurposeGenAndRecoAnalyzer(const edm::ParameterSet& iConfig):
	_tName(iConfig.getParameter<std::string>("treeName")),
	processingGen(iConfig.getParameter<bool>("inputParticleCollsAreGEN")),
	processingQuarks(iConfig.getParameter<bool>("inputHadronsAreQuarks"))

{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>(_tName.c_str(), "event kinematic info");

	tree->Branch("runNumber", &runNumber, "runNumber/I");
	tree->Branch("evtNumber", &evtNumber, "evtNumber/l");

	tree->Branch("evWeight", &evWeight, "evWeight/F");
	tree->Branch("evWeightSign", &evWeightSign, "evWeightSign/F");

	///branches for GEN quantities
	tree->Branch("etaLeptOne", &etaLeptOne, "etaLeptOne/F");
	tree->Branch("ptLeptOne", &ptLeptOne, "ptLeptOne/F");
	tree->Branch("phiLeptOne", &phiLeptOne, "phiLeptOne/F");
	tree->Branch("etaLeptTwo", &etaLeptTwo, "etaLeptTwo/F");
	tree->Branch("ptLeptTwo", &ptLeptTwo, "ptLeptTwo/F");
	tree->Branch("phiLeptTwo", &phiLeptTwo, "phiLeptTwo/F");
	tree->Branch("motherPdgIdLeptOne", &motherPdgIdLeptOne, "motherPdgIdLeptOne/I");
	tree->Branch("motherStatusLeptOne", &motherStatusLeptOne, "motherStatusLeptOne/I");
	tree->Branch("motherPdgIdLeptTwo", &motherPdgIdLeptTwo, "motherPdgIdLeptTwo/I");
	tree->Branch("motherStatusLeptTwo", &motherStatusLeptTwo, "motherStatusLeptTwo/I");
	tree->Branch("motherPdgIdQuarkOne", &motherPdgIdQuarkOne, "motherPdgIdQuarkOne/I");
	tree->Branch("motherStatusQuarkOne", &motherStatusQuarkOne, "motherStatusQuarkOne/I");
	tree->Branch("motherPdgIdQuarkTwo", &motherPdgIdQuarkTwo, "motherPdgIdQuarkTwo/I");
	tree->Branch("motherStatusQuarkTwo", &motherStatusQuarkTwo, "motherStatusQuarkTwo/I");
	tree->Branch("etaQuarkOne", &etaQuarkOne, "etaQuarkOne/F");
	tree->Branch("ptQuarkOne", &ptQuarkOne, "ptQuarkOne/F");
	tree->Branch("phiQuarkOne", &phiQuarkOne, "phiQuarkOne/F");
	tree->Branch("etaQuarkTwo", &etaQuarkTwo, "etaQuarkTwo/F");
	tree->Branch("ptQuarkTwo", &ptQuarkTwo, "ptQuarkTwo/F");
	tree->Branch("phiQuarkTwo", &phiQuarkTwo, "phiQuarkTwo/F");

	//dilepton system
	tree->Branch("dileptonMass", &dileptonMass, "dileptonMass/F");
	tree->Branch("dileptonPt", &dileptonPt, "dileptonPt/F");
	tree->Branch("dileptonEta", &dileptonEta, "dileptonEta/F");
	tree->Branch("dileptonPhi", &dileptonPhi, "dileptonPhi/F");

	//sublead lepton + dihadron system
	tree->Branch("subleadLeptonBothHadronsPt", &subleadLeptonBothHadronsPt, "subleadLeptonBothHadronsPt/F");
	tree->Branch("subleadLeptonBothHadronsMass", &subleadLeptonBothHadronsMass, "subleadLeptonBothHadronsMass/F");
	tree->Branch("subleadLeptonBothHadronsEta", &subleadLeptonBothHadronsEta, "subleadLeptonBothHadronsEta/F");
	tree->Branch("subleadLeptonBothHadronsPhi", &subleadLeptonBothHadronsPhi, "subleadLeptonBothHadronsPhi/F");
	
	//lead lepton + dihadron system
	tree->Branch("leadLeptonBothHadronsPt", &leadLeptonBothHadronsPt, "leadLeptonBothHadronsPt/F");
	tree->Branch("leadLeptonBothHadronsMass", &leadLeptonBothHadronsMass, "leadLeptonBothHadronsMass/F");
	tree->Branch("leadLeptonBothHadronsEta", &leadLeptonBothHadronsEta, "leadLeptonBothHadronsEta/F");
	tree->Branch("leadLeptonBothHadronsPhi", &leadLeptonBothHadronsPhi, "leadLeptonBothHadronsPhi/F");

	//dilepton + dihadron system
	tree->Branch("dileptonDihadronPt", &dileptonDihadronPt, "dileptonDihadronPt/F");
	tree->Branch("dileptonDihadronMass", &dileptonDihadronMass, "dileptonDihadronMass/F");
	tree->Branch("dileptonDihadronEta", &dileptonDihadronEta, "dileptonDihadronEta/F");
	tree->Branch("dileptonDihadronPhi", &dileptonDihadronPhi, "dileptonDihadronPhi/F");

	tree->Branch("dRleptonOneQuarkOne", &dRleptonOneQuarkOne, "dRleptonOneQuarkOne/F");
	tree->Branch("dRleptonOneQuarkTwo", &dRleptonOneQuarkTwo, "dRleptonOneQuarkTwo/F");
	tree->Branch("dRleptonTwoQuarkOne", &dRleptonTwoQuarkOne, "dRleptonTwoQuarkOne/F");
	tree->Branch("dRleptonTwoQuarkTwo", &dRleptonTwoQuarkTwo, "dRleptonTwoQuarkTwo/F");

	tree->Branch("numLeptons", &numLeptons, "numLeptons/I");
	tree->Branch("numQuarks", &numQuarks, "numQuarks/I");

	///tokens to input collections
	leptonCollectionToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("leptonCollection"));
	hadronCollectionToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("hadronCollection"));
	genEventInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

}


generalPurposeGenAndRecoAnalyzer::~generalPurposeGenAndRecoAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
generalPurposeGenAndRecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	resetTreeVars();
#ifdef DEBUG
	std::cout << "in analyze method of generalPurposeGenAndRecoAnalyzer class" << std::endl;
#endif

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();
	evWeight = 1.0;
	evWeightSign = 1.0;

	iEvent.getByToken(leptonCollectionToken, leptonCollection);
	iEvent.getByToken(hadronCollectionToken, hadronCollection);
	iEvent.getByToken(genEventInfoToken, genEvtInfo);	///< get evt weights if analyzing MC

#ifdef DEBUG
	std::cout << "called getByToken to fill leptonCollection and genEvtInfo handles" << std::endl;
#endif

	if(genEvtInfo.isValid() ) {
		///real data does not have gen lvl event weights
		evWeight = genEvtInfo->weight();
		if(evWeight < 0) evWeightSign = -1.0;
	}

	if(leptonCollection->size() == 0 || hadronCollection->size() == 0) {
		tree->Fill();
		return;
	}

	edm::OwnVector<reco::Candidate>::const_iterator leadHadronPtr = hadronCollection->end(), subleadHadronPtr = hadronCollection->end(), leadLept = leptonCollection->end(), subleadLept = leptonCollection->end();

	for(edm::OwnVector<reco::Candidate>::const_iterator it = leptonCollection->begin(); it != leptonCollection->end(); it++) {
		//find the leading and subleading leptons
		if(std::fabs(it->eta()) <= 10 && it->pt() > 0.01) {
			numLeptons++;
#ifdef DEBUG
			std::cout << " " << std::endl;
			std::cout << "lepton pt=\t" << it->pt() << std::endl;
			std::cout << "lepton eta=\t" << it->eta() << std::endl;
			std::cout << "lepton phi=\t" << it->phi() << std::endl;
			std::cout << " " << std::endl;
#endif

			if(leadLept == leptonCollection->end()) leadLept = it;
			if(leadLept->pt() == it->pt() || diobjectMass(leadLept, it) == 0.) continue;	///<avoid selecting lepton pairs with dilepton mass zero
			if(subleadLept != leptonCollection->end()){
				if(subleadLept->pt() == it->pt()) continue;
			}
			else {
				if(it->pt() > leadLept->pt() && deltaR(it->eta(), it->phi(), leadLept->eta(), leadLept->phi() ) > 0.1) {
					subleadLept = leadLept;
					leadLept = it;
				} else if(it->pt() > leadLept->pt() && deltaR(it->eta(), it->phi(), leadLept->eta(), leadLept->phi() ) <= 0.1) leadLept = it;
				else if((subleadLept == leptonCollection->end() || it->pt() > subleadLept->pt()) && deltaR(it->eta(), it->phi(), leadLept->eta(), leadLept->phi() ) > 0.1) subleadLept = it;

			}//end reassigning lepton iterators
		}///end if(iterator points to a lepton with sensible pt and eta)
	}///end loop over leptons collection

	for(edm::OwnVector<reco::Candidate>::const_iterator it = hadronCollection->begin(); it != hadronCollection->end(); it++) {
		//find the leading and subleading hadrons
		if(std::fabs(it->eta()) < 10 && it->pt() > 0.01) {
			//if a leading and or subleading lepton has already been found, make sure a selected hadron does not overlap with a selected lepton
			if(subleadLept != leptonCollection->end()){
				if(deltaR(it->eta(), it->phi(), subleadLept->eta(), subleadLept->phi())<=0.01 ) continue;
			}
			if(leadLept != leptonCollection->end()){
				if(deltaR(it->eta(), it->phi(), leadLept->eta(), leadLept->phi())<=0.01 ) continue;
			}

#ifdef DEBUG
			std::cout << " " << std::endl;
			std::cout << "hadron pt=\t" << it->pt() << std::endl;
			std::cout << "hadron eta=\t" << it->eta() << std::endl;
			std::cout << "hadron phi=\t" << it->phi() << std::endl;
			std::cout << " " << std::endl;
#endif
			
			numQuarks++;
			if(leadHadronPtr == hadronCollection->end()) leadHadronPtr = it;
			if(diobjectMass(leadHadronPtr, it) == 0.) continue;
			else {
				if(it->pt() > leadHadronPtr->pt() && deltaR(it->eta(), it->phi(), leadHadronPtr->eta(), leadHadronPtr->phi() ) > 0.1) {
					subleadHadronPtr = leadHadronPtr;
					leadHadronPtr = it;
				} else if(it->pt() > leadHadronPtr->pt() && deltaR(it->eta(), it->phi(), leadHadronPtr->eta(), leadHadronPtr->phi() ) <= 0.1) leadHadronPtr = it;
				else if((subleadHadronPtr == hadronCollection->end() || it->pt() > subleadHadronPtr->pt()) && deltaR(it->eta(), it->phi(), leadHadronPtr->eta(), leadHadronPtr->phi() ) > 0.1) subleadHadronPtr = it;

			}//end reassigning hadron iterators
		}///end if(iterator points to a hadron with sensible pt and eta)
	}///end loop over hadrons collection

	if(leadLept == leptonCollection->end() || subleadLept == leptonCollection->end() || leadHadronPtr == hadronCollection->end() || subleadHadronPtr == hadronCollection->end() ) {
		tree->Fill();
		return;
	}
#ifdef DEBUG
	std::cout << "found leading leptons and hadrons" << std::endl;
#endif
	
	///now calculate the desired kinematic quantities of all particles, and save them
	etaLeptOne = leadLept->eta(), ptLeptOne = leadLept->pt(), phiLeptOne = leadLept->phi();
	etaLeptTwo = subleadLept->eta(), ptLeptTwo = subleadLept->pt(), phiLeptTwo = subleadLept->phi();
	etaQuarkTwo = subleadHadronPtr->eta(), ptQuarkTwo = subleadHadronPtr->pt(), phiQuarkTwo = subleadHadronPtr->phi();
	etaQuarkOne = leadHadronPtr->eta(), ptQuarkOne = leadHadronPtr->pt(), phiQuarkOne = leadHadronPtr->phi();

#ifdef DEBUG
	std::cout <<"saved pt, eta and phi of two leading leptons and hadrons into tree"<<std::endl;
	std::cout<<"leptOne pt eta phi mass:\t"<< ptLeptOne <<"  "<< etaLeptOne <<"  "<< phiLeptOne <<"  " << leadLept->mass() <<std::endl;
	std::cout<<"leptTwo pt eta phi mass:\t"<< ptLeptTwo <<"  "<< etaLeptTwo <<"  "<< phiLeptTwo <<"  " << subleadLept->mass() <<std::endl;
	std::cout<<"hadronOne pt eta phi mass:\t"<< ptQuarkOne <<"  "<< etaQuarkOne <<"  "<< phiQuarkOne << "  " << leadHadronPtr->mass() <<std::endl;
	std::cout<<"hadronTwo pt eta phi mass:\t"<< ptQuarkTwo <<"  "<< etaQuarkTwo <<"  "<< phiQuarkTwo <<"  " << subleadHadronPtr->mass() <<std::endl;
#endif

	TLorentzVector fourMomLeadLept, fourMomSubleadLept, fourMomLeadHadron, fourMomSubleadHadron;
	fourMomLeadLept.SetPtEtaPhiM(ptLeptOne, etaLeptOne, phiLeptOne, leadLept->mass());
	fourMomSubleadLept.SetPtEtaPhiM(ptLeptTwo, etaLeptTwo, phiLeptTwo, subleadLept->mass());
	fourMomLeadHadron.SetPtEtaPhiM(ptQuarkOne, etaQuarkOne, phiQuarkOne, leadHadronPtr->mass());
	fourMomSubleadHadron.SetPtEtaPhiM(ptQuarkTwo, etaQuarkTwo, phiQuarkTwo, subleadHadronPtr->mass());

#ifdef DEBUG
	std::cout <<"built TLorentzVector objects to represent leading leptons and hadrons"<<std::endl;
#endif

	TLorentzVector dileptonVect = fourMomLeadLept + fourMomSubleadLept;
	if(dileptonVect.M2() > 0. && dileptonVect.Perp2() > 0. && dileptonVect.Et2() > 0. && std::fabs(dileptonVect.CosTheta()) < 1.0){
		dileptonMass = dileptonVect.M();
		dileptonPt = dileptonVect.Pt();
		dileptonEta = dileptonVect.Eta();
		dileptonPhi = dileptonVect.Phi();
	}

#ifdef DEBUG
	std::cout <<"calculated and saved dilepton system kinematics"<<std::endl;
#endif

	TLorentzVector subleadLJJVect = fourMomSubleadLept + fourMomLeadHadron + fourMomSubleadHadron;
	if(subleadLJJVect.M2() > 0. && subleadLJJVect.Perp2() > 0. && subleadLJJVect.Et2() > 0. && std::fabs(subleadLJJVect.CosTheta()) < 1.0){
		subleadLeptonBothHadronsMass = subleadLJJVect.M();
		subleadLeptonBothHadronsPt = subleadLJJVect.Pt();
		subleadLeptonBothHadronsEta = subleadLJJVect.Eta();
		subleadLeptonBothHadronsPhi = subleadLJJVect.Phi();
	}

#ifdef DEBUG
	std::cout <<"calculated and saved sublead lepton plus dihadron system kinematics"<<std::endl;
#endif

	TLorentzVector leadLJJVect = fourMomLeadLept + fourMomLeadHadron + fourMomSubleadHadron;
	if(leadLJJVect.M2() > 0. && leadLJJVect.Perp2() > 0. && leadLJJVect.Et2() > 0. && std::fabs(leadLJJVect.CosTheta()) < 1.0){
		leadLeptonBothHadronsMass = leadLJJVect.M();
		leadLeptonBothHadronsPt = leadLJJVect.Pt();
		leadLeptonBothHadronsEta = leadLJJVect.Eta();
		leadLeptonBothHadronsPhi = leadLJJVect.Phi();
	}

#ifdef DEBUG
	std::cout <<"calculated and saved lead lepton plus dihadron system kinematics"<<std::endl;
#endif

	TLorentzVector lljjFourVect = fourMomLeadLept + fourMomSubleadLept + fourMomLeadHadron + fourMomSubleadHadron;
#ifdef DEBUG
	std::cout <<"made four vector representing dilepton plus dihadron system"<<std::endl;
#endif
	if(lljjFourVect.M2() > 0. && lljjFourVect.Perp2() > 0. && lljjFourVect.Et2() > 0. && std::fabs(lljjFourVect.CosTheta()) < 1.0){
#ifdef DEBUG
	std::cout <<"dilepton plus dihadron system mass, squared pt, and squared et are all greater than zero"<<std::endl;
	std::cout<<"CosTheta =\t"<< lljjFourVect.CosTheta() <<"  "<<"fZ=\t"<< lljjFourVect(2) << std::endl;
#endif
		dileptonDihadronMass = (lljjFourVect).M();
		dileptonDihadronPt = (lljjFourVect).Pt();
		dileptonDihadronEta = (lljjFourVect).Eta();
		dileptonDihadronPhi = (lljjFourVect).Phi();
	}

#ifdef DEBUG
	std::cout <<"calculated and saved dilepton plus dihadron system kinematics"<<std::endl;
#endif

	dRleptonOneQuarkOne = deltaR(fourMomLeadLept.Eta(), fourMomLeadLept.Phi(), fourMomLeadHadron.Eta(), fourMomLeadHadron.Phi());
	dRleptonOneQuarkTwo = deltaR(fourMomLeadLept.Eta(), fourMomLeadLept.Phi(), fourMomSubleadHadron.Eta(), fourMomSubleadHadron.Phi());
	dRleptonTwoQuarkOne = deltaR(fourMomSubleadLept.Eta(), fourMomSubleadLept.Phi(), fourMomLeadHadron.Eta(), fourMomLeadHadron.Phi());
	dRleptonTwoQuarkTwo = deltaR(fourMomSubleadLept.Eta(), fourMomSubleadLept.Phi(), fourMomSubleadHadron.Eta(), fourMomSubleadHadron.Phi());

#ifdef DEBUG
	std::cout << "saved kinematic info from leptons and hadrons, and multi particle systems built from them" << std::endl;
#endif

	if(processingGen) {
		motherPdgIdLeptOne = (leadLept->mother(0))->pdgId(), motherStatusLeptOne = (leadLept->mother(0))->status();
		motherPdgIdLeptTwo = (subleadLept->mother(0))->pdgId(), motherStatusLeptTwo = (subleadLept->mother(0))->status();
		if(processingQuarks) motherPdgIdQuarkTwo = (subleadHadronPtr->mother(0))->pdgId(), motherStatusQuarkTwo = (subleadHadronPtr->mother(0))->status();
		if(processingQuarks) motherPdgIdQuarkOne = (leadHadronPtr->mother(0))->pdgId(), motherStatusQuarkOne = (leadHadronPtr->mother(0))->status();
	}///end filling GEN mother pdgId and status branches
	
	tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
generalPurposeGenAndRecoAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
generalPurposeGenAndRecoAnalyzer::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void generalPurposeGenAndRecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(generalPurposeGenAndRecoAnalyzer);
