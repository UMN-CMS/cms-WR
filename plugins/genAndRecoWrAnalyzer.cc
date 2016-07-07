// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/genAndRecoWrAnalyzer
// Class:      genAndRecoWrAnalyzer
//
/**\class genAndRecoWrAnalyzer genAndRecoWrAnalyzer.cc doubleElectronTracklessTrigger/genAndRecoWrAnalyzer/plugins/genAndRecoWrAnalyzer.cc

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

class genAndRecoWrAnalyzer : public edm::EDAnalyzer
{
public:
	explicit genAndRecoWrAnalyzer(const edm::ParameterSet&);
	~genAndRecoWrAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

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

		/*
		for(edm::OwnVector<reco::Candidate>::iterator genIt = collection.begin(); genIt != collection.end(); ++genIt){
		  if(first==collection.end()) first=genIt;
		  else{
		#ifdef DEBUG
			  std::cout<<"first iterator in findLeadingAndSubleadingVectorInput is no longer null"<<std::endl;
		#endif
			  if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1 ){
		#ifdef DEBUG
				  std::cout<<"reassigning second iterator and first iterator"<<std::endl;
		#endif
				  second = first;
				  first = genIt;
		#ifdef DEBUG
				  std::cout<<"second iterator in findLeadingAndSubleadingVectorInput is no longer null"<<std::endl;
		#endif
			  }
			  else if( (second==collection.end() || genIt->pt() > second->pt() ) && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.1  ){
		#ifdef DEBUG
				  std::cout<<"second iterator is about to be reset, but first iterator will not change"<<std::endl;
		#endif
				  second = genIt;
			  }
		  }
		}//end loop over reco::Candidate collection
		*/

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
		runNumber = -1, evtNumber = -1, numGenFstHvyPtcl = -1, numGenScdHvyPtcl = -1, numGenLeptons = -1, numGenQuarks = -1, leadGenLeptonNotFromFstHvyPtcl = -1, subleadGenLeptonNotFromScdHvyPtcl = -1, leadGenQuarkNotFromScdHvyPtcl = -1, subleadGenQuarkNotFromScdHvyPtcl = -1, wrLeptonDidRadiatePhoton = -1;
	
		motherPdgIdGenLeptFromFstHvyPtcl = -1, motherStatusGenLeptFromFstHvyPtcl = -1, motherPdgIdGenLeptFromScdHvyPtcl = -1, motherStatusGenLeptFromScdHvyPtcl = -1, motherPdgIdGenQuarkOneFromScdHvyPtcl = -1, motherStatusGenQuarkOneFromScdHvyPtcl = -1, motherPdgIdGenQuarkTwoFromScdHvyPtcl = -1, motherStatusGenQuarkTwoFromScdHvyPtcl = -1;
		
		etaGenFstHvyPtcl = -9, ptGenFstHvyPtcl = -9, massGenFstHvyPtcl = -9, etaGenScdHvyPtcl = -9, ptGenScdHvyPtcl = -9, massGenScdHvyPtcl = -9, etaGenLeptFromFstHvyPtcl = -9, ptGenLeptFromFstHvyPtcl = -9, phiGenLeptFromFstHvyPtcl = -9, etaGenLeptFromScdHvyPtcl = -9, ptGenLeptFromScdHvyPtcl = -9, phiGenLeptFromScdHvyPtcl = -9, etaGenQuarkOneFromScdHvyPtcl = -9, ptGenQuarkOneFromScdHvyPtcl = -9, phiGenQuarkOneFromScdHvyPtcl = -9, etaGenQuarkTwoFromScdHvyPtcl = -9, ptGenQuarkTwoFromScdHvyPtcl = -9, phiGenQuarkTwoFromScdHvyPtcl = -9;

		ptLeadGenLepton = -9, etaLeadGenLepton = -9, phiLeadGenLepton = -9, ptSubleadGenLepton = -9, etaSubleadGenLepton = -9, phiSubleadGenLepton = -9, ptLeadGenQuark = -9, etaLeadGenQuark = -9, phiLeadGenQuark = -9, ptSubleadGenQuark = -9, etaSubleadGenQuark = -9, phiSubleadGenQuark = -9, evWeight = 1, evWeightSign = 1;

		ptRecoLeptMatchedToWrDau = -9, etaRecoLeptMatchedToWrDau = -9, phiRecoLeptMatchedToWrDau = -9, ptRecoLeptMatchedToNuDau = -9, etaRecoLeptMatchedToNuDau = -9, phiRecoLeptMatchedToNuDau = -9, ptRecoJetOneMatchedToNuDau = -9, etaRecoJetOneMatchedToNuDau = -9, phiRecoJetOneMatchedToNuDau = -9, ptRecoJetTwoMatchedToNuDau = -9, etaRecoJetTwoMatchedToNuDau = -9, phiRecoJetTwoMatchedToNuDau = -9, ptGenJetFromMatchedRecoJetOne = -9, etaGenJetFromMatchedRecoJetOne = -9, phiGenJetFromMatchedRecoJetOne = -9, ptGenJetFromMatchedRecoJetTwo = -9, etaGenJetFromMatchedRecoJetTwo = -9, phiGenJetFromMatchedRecoJetTwo = -9, ptLeadRecoLept = -9, etaLeadRecoLept = -9, phiLeadRecoLept = -9, ptSubleadRecoLept = -9, etaSubleadRecoLept = -9, phiSubleadRecoLept = -9, ptLeadRecoJet = -9, etaLeadRecoJet = -9, phiLeadRecoJet = -9, ptSubleadRecoJet = -9, etaSubleadRecoJet = -9, phiSubleadRecoJet = -9;
		
		motherPdgIdLeadGenLepton = -1, motherStatusLeadGenLepton = -1;
		motherPdgIdSubleadGenLepton = -1, motherStatusSubleadGenLepton = -1;
		motherPdgIdLeadGenQuark = -1, motherStatusLeadGenQuark = -1;
		motherPdgIdSubleadGenQuark = -1, motherStatusSubleadGenQuark = -1;

		fourObjMassFromGenObjsFromFstAndScdHvyPtcl = -9;
	
	}///end resetTreeVars()

private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

// ----------member data ---------------------------

	std::string _tName;
	double _leptonPdgId;
	double _minQuarkPdgId, _maxQuarkPdgId;
	double _firstHeavyParticlePdgId, _secondHeavyParticlePdgId;
	double _firstHeavyParticleStatus;
///set _cutOnStatusCode to true to require that the GEN firstHeavyParticle has a specific status code
	bool _saveMatchedRecoInfo, _cutOnStatusCode;

///Handles to MINIAOD object collections
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > genParticleCollection;
	edm::Handle<std::vector<pat::Jet> > recoJetCollection;
	edm::Handle<reco::CandidateView> recoLeptonCollection;

	edm::Handle<GenEventInfoProduct> genEvtInfo;	///<gen event weight

///tokens to input collections
	edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > genParticleCollectionToken;
	edm::EDGetTokenT<std::vector<pat::Jet> > recoJetCollectionToken;
	edm::EDGetTokenT<reco::CandidateView> recoLeptonCollectionToken;
	edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken;

	TTree * tree;

	Int_t runNumber;
	ULong64_t evtNumber;


///branches for GEN quantities
///mass, pt, eta of GEN lvl heavy, unstable particles
	Int_t numGenFstHvyPtcl;
	Float_t etaGenFstHvyPtcl;
	Float_t ptGenFstHvyPtcl;
	Float_t massGenFstHvyPtcl;

	Int_t numGenScdHvyPtcl;
	Float_t etaGenScdHvyPtcl;
	Float_t ptGenScdHvyPtcl;
	Float_t massGenScdHvyPtcl;

///pt and eta of GEN leptons and quarks whose mothers are the heavy, unstable particles
///it is assumed that the quarks only come from the second heavy particle
	Float_t etaGenLeptFromFstHvyPtcl;
	Float_t ptGenLeptFromFstHvyPtcl;
	Float_t phiGenLeptFromFstHvyPtcl;
	Float_t etaGenLeptFromScdHvyPtcl;
	Float_t ptGenLeptFromScdHvyPtcl;
	Float_t phiGenLeptFromScdHvyPtcl;
	Float_t etaGenQuarkOneFromScdHvyPtcl;
	Float_t ptGenQuarkOneFromScdHvyPtcl;
	Float_t phiGenQuarkOneFromScdHvyPtcl;
	Float_t etaGenQuarkTwoFromScdHvyPtcl;
	Float_t ptGenQuarkTwoFromScdHvyPtcl;
	Float_t phiGenQuarkTwoFromScdHvyPtcl;
	Int_t wrLeptonDidRadiatePhoton;	///<1 if true, 0 if false
	Int_t motherPdgIdGenLeptFromFstHvyPtcl;
	Int_t motherStatusGenLeptFromFstHvyPtcl;
	Int_t motherPdgIdGenLeptFromScdHvyPtcl;
	Int_t motherStatusGenLeptFromScdHvyPtcl;
	Int_t motherPdgIdGenQuarkOneFromScdHvyPtcl;
	Int_t motherStatusGenQuarkOneFromScdHvyPtcl;
	Int_t motherPdgIdGenQuarkTwoFromScdHvyPtcl;
	Int_t motherStatusGenQuarkTwoFromScdHvyPtcl;
	Float_t fourObjMassFromGenObjsFromFstAndScdHvyPtcl;
	


///pt and eta of the two leading GEN leptons and quarks in each event, no requirements on mother pdgId
///count how often the lead lepton does not come from the first heavy particle (use mother info), how often
///the sublead lepton does not come from the second heavy particle, how often the lead
///quark does not come from the second heavy particle, and how often the sublead quark
///does not come from the second heavy particle
	Int_t numGenLeptons, numGenQuarks;
	Float_t ptLeadGenLepton, etaLeadGenLepton, phiLeadGenLepton;
	Int_t leadGenLeptonNotFromFstHvyPtcl;	///< 1 if true, 0 if false
	Int_t motherPdgIdLeadGenLepton, motherStatusLeadGenLepton;
	Float_t ptSubleadGenLepton, etaSubleadGenLepton, phiSubleadGenLepton;
	Int_t subleadGenLeptonNotFromScdHvyPtcl;	///< 1 if true, 0 if false
	Int_t motherPdgIdSubleadGenLepton, motherStatusSubleadGenLepton;
	Float_t ptLeadGenQuark, etaLeadGenQuark, phiLeadGenQuark;
	Int_t leadGenQuarkNotFromScdHvyPtcl;	///< 1 if true, 0 if false
	Int_t motherPdgIdLeadGenQuark, motherStatusLeadGenQuark;
	Float_t ptSubleadGenQuark, etaSubleadGenQuark, phiSubleadGenQuark;
	Int_t subleadGenQuarkNotFromScdHvyPtcl;	///< 1 if true, 0 if false
	Int_t motherPdgIdSubleadGenQuark, motherStatusSubleadGenQuark;
	
	Float_t evWeight;	///< gen weight of the event.  defaults to 1, only changes for bkgnd MC
	Float_t evWeightSign;	///< evtWeightSign used to determine weight of the event when filling a histo

///branches for RECO quantities
///filled if _saveMatchedRecoInfo is true

///pt and eta of reco leptons and jets matched to GEN Wr and Nu daughter leptons and quarks
	Float_t ptRecoLeptMatchedToWrDau, etaRecoLeptMatchedToWrDau, phiRecoLeptMatchedToWrDau;
	Float_t ptRecoLeptMatchedToNuDau, etaRecoLeptMatchedToNuDau, phiRecoLeptMatchedToNuDau;
	Float_t ptRecoJetOneMatchedToNuDau, etaRecoJetOneMatchedToNuDau, phiRecoJetOneMatchedToNuDau;
	Float_t ptRecoJetTwoMatchedToNuDau, etaRecoJetTwoMatchedToNuDau, phiRecoJetTwoMatchedToNuDau;
	Float_t ptGenJetFromMatchedRecoJetOne, etaGenJetFromMatchedRecoJetOne, phiGenJetFromMatchedRecoJetOne;
	Float_t ptGenJetFromMatchedRecoJetTwo, etaGenJetFromMatchedRecoJetTwo, phiGenJetFromMatchedRecoJetTwo;

///pt and eta of reco leading leptons and jets, no matching to GEN Wr or Nu daughter leptons and quarks
	Float_t ptLeadRecoLept, etaLeadRecoLept, phiLeadRecoLept;
	Float_t ptSubleadRecoLept, etaSubleadRecoLept, phiSubleadRecoLept;
	Float_t ptLeadRecoJet, etaLeadRecoJet, phiLeadRecoJet;
	Float_t ptSubleadRecoJet, etaSubleadRecoJet, phiSubleadRecoJet;


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

genAndRecoWrAnalyzer::genAndRecoWrAnalyzer(const edm::ParameterSet& iConfig):
	_tName(iConfig.getParameter<std::string>("treeName")),
	_leptonPdgId(iConfig.getParameter<double>("leptonPdgId")),
	_minQuarkPdgId(iConfig.getParameter<double>("minQuarkPdgId")),
	_maxQuarkPdgId(iConfig.getParameter<double>("maxQuarkPdgId")),
	_firstHeavyParticlePdgId(iConfig.getParameter<double>("firstHeavyParticlePdgId")),
	_secondHeavyParticlePdgId(iConfig.getParameter<double>("secondHeavyParticlePdgId")),
	_firstHeavyParticleStatus(iConfig.getParameter<double>("firstHeavyParticleStatus")),
	_saveMatchedRecoInfo(iConfig.getParameter<bool>("saveMatchedRecoInfo")),
	_cutOnStatusCode(iConfig.getParameter<bool>("cutOnStatusCode"))

{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>(_tName.c_str(), "event kinematic info");

	tree->Branch("runNumber", &runNumber, "runNumber/I");
	tree->Branch("evtNumber", &evtNumber, "evtNumber/l");

	///branches for GEN quantities
	tree->Branch("etaGenFstHvyPtcl", &etaGenFstHvyPtcl, "etaGenFstHvyPtcl/F");
	tree->Branch("ptGenFstHvyPtcl", &ptGenFstHvyPtcl, "ptGenFstHvyPtcl/F");
	tree->Branch("massGenFstHvyPtcl", &massGenFstHvyPtcl, "massGenFstHvyPtcl/F");

	tree->Branch("etaGenScdHvyPtcl", &etaGenScdHvyPtcl, "etaGenScdHvyPtcl/F");
	tree->Branch("ptGenScdHvyPtcl", &ptGenScdHvyPtcl, "ptGenScdHvyPtcl/F");
	tree->Branch("massGenScdHvyPtcl", &massGenScdHvyPtcl, "massGenScdHvyPtcl/F");

	tree->Branch("etaGenLeptFromFstHvyPtcl", &etaGenLeptFromFstHvyPtcl, "etaGenLeptFromFstHvyPtcl/F");
	tree->Branch("ptGenLeptFromFstHvyPtcl", &ptGenLeptFromFstHvyPtcl, "ptGenLeptFromFstHvyPtcl/F");
	tree->Branch("phiGenLeptFromFstHvyPtcl", &phiGenLeptFromFstHvyPtcl, "phiGenLeptFromFstHvyPtcl/F");
	tree->Branch("etaGenLeptFromScdHvyPtcl", &etaGenLeptFromScdHvyPtcl, "etaGenLeptFromScdHvyPtcl/F");
	tree->Branch("ptGenLeptFromScdHvyPtcl", &ptGenLeptFromScdHvyPtcl, "ptGenLeptFromScdHvyPtcl/F");
	tree->Branch("phiGenLeptFromScdHvyPtcl", &phiGenLeptFromScdHvyPtcl, "phiGenLeptFromScdHvyPtcl/F");
	tree->Branch("wrLeptonDidRadiatePhoton", &wrLeptonDidRadiatePhoton, "wrLeptonDidRadiatePhoton/I");
	tree->Branch("motherPdgIdGenLeptFromFstHvyPtcl", &motherPdgIdGenLeptFromFstHvyPtcl, "motherPdgIdGenLeptFromFstHvyPtcl/I");
	tree->Branch("motherStatusGenLeptFromFstHvyPtcl", &motherStatusGenLeptFromFstHvyPtcl, "motherStatusGenLeptFromFstHvyPtcl/I");
	tree->Branch("motherPdgIdGenLeptFromScdHvyPtcl", &motherPdgIdGenLeptFromScdHvyPtcl, "motherPdgIdGenLeptFromScdHvyPtcl/I");
	tree->Branch("motherStatusGenLeptFromScdHvyPtcl", &motherStatusGenLeptFromScdHvyPtcl, "motherStatusGenLeptFromScdHvyPtcl/I");
	tree->Branch("motherPdgIdGenQuarkOneFromScdHvyPtcl", &motherPdgIdGenQuarkOneFromScdHvyPtcl, "motherPdgIdGenQuarkOneFromScdHvyPtcl/I");
	tree->Branch("motherStatusGenQuarkOneFromScdHvyPtcl", &motherStatusGenQuarkOneFromScdHvyPtcl, "motherStatusGenQuarkOneFromScdHvyPtcl/I");
	tree->Branch("motherPdgIdGenQuarkTwoFromScdHvyPtcl", &motherPdgIdGenQuarkTwoFromScdHvyPtcl, "motherPdgIdGenQuarkTwoFromScdHvyPtcl/I");
	tree->Branch("motherStatusGenQuarkTwoFromScdHvyPtcl", &motherStatusGenQuarkTwoFromScdHvyPtcl, "motherStatusGenQuarkTwoFromScdHvyPtcl/I");
	tree->Branch("etaGenQuarkOneFromScdHvyPtcl", &etaGenQuarkOneFromScdHvyPtcl, "etaGenQuarkOneFromScdHvyPtcl/F");
	tree->Branch("ptGenQuarkOneFromScdHvyPtcl", &ptGenQuarkOneFromScdHvyPtcl, "ptGenQuarkOneFromScdHvyPtcl/F");
	tree->Branch("phiGenQuarkOneFromScdHvyPtcl", &phiGenQuarkOneFromScdHvyPtcl, "phiGenQuarkOneFromScdHvyPtcl/F");
	tree->Branch("etaGenQuarkTwoFromScdHvyPtcl", &etaGenQuarkTwoFromScdHvyPtcl, "etaGenQuarkTwoFromScdHvyPtcl/F");
	tree->Branch("ptGenQuarkTwoFromScdHvyPtcl", &ptGenQuarkTwoFromScdHvyPtcl, "ptGenQuarkTwoFromScdHvyPtcl/F");
	tree->Branch("phiGenQuarkTwoFromScdHvyPtcl", &phiGenQuarkTwoFromScdHvyPtcl, "phiGenQuarkTwoFromScdHvyPtcl/F");
	tree->Branch("fourObjMassFromGenObjsFromFstAndScdHvyPtcl", &fourObjMassFromGenObjsFromFstAndScdHvyPtcl, "fourObjMassFromGenObjsFromFstAndScdHvyPtcl/F");
	

	tree->Branch("etaLeadGenLepton", &etaLeadGenLepton, "etaLeadGenLepton/F");
	tree->Branch("ptLeadGenLepton", &ptLeadGenLepton, "ptLeadGenLepton/F");
	tree->Branch("phiLeadGenLepton", &phiLeadGenLepton, "phiLeadGenLepton/F");
	tree->Branch("leadGenLeptonNotFromFstHvyPtcl", &leadGenLeptonNotFromFstHvyPtcl, "leadGenLeptonNotFromFstHvyPtcl/I");
	tree->Branch("motherPdgIdLeadGenLepton", &motherPdgIdLeadGenLepton, "motherPdgIdLeadGenLepton/I");
	tree->Branch("motherStatusLeadGenLepton", &motherStatusLeadGenLepton, "motherStatusLeadGenLepton/I");
	
	tree->Branch("etaSubleadGenLepton", &etaSubleadGenLepton, "etaSubleadGenLepton/F");
	tree->Branch("ptSubleadGenLepton", &ptSubleadGenLepton, "ptSubleadGenLepton/F");
	tree->Branch("phiSubleadGenLepton", &phiSubleadGenLepton, "phiSubleadGenLepton/F");
	tree->Branch("subleadGenLeptonNotFromScdHvyPtcl", &subleadGenLeptonNotFromScdHvyPtcl, "subleadGenLeptonNotFromScdHvyPtcl/I");
	tree->Branch("motherPdgIdSubleadGenLepton", &motherPdgIdSubleadGenLepton, "motherPdgIdSubleadGenLepton/I");
	tree->Branch("motherStatusSubleadGenLepton", &motherStatusSubleadGenLepton, "motherStatusSubleadGenLepton/I");
	
	tree->Branch("etaLeadGenQuark", &etaLeadGenQuark, "etaLeadGenQuark/F");
	tree->Branch("ptLeadGenQuark", &ptLeadGenQuark, "ptLeadGenQuark/F");
	tree->Branch("phiLeadGenQuark", &phiLeadGenQuark, "phiLeadGenQuark/F");
	tree->Branch("leadGenQuarkNotFromScdHvyPtcl", &leadGenQuarkNotFromScdHvyPtcl, "leadGenQuarkNotFromScdHvyPtcl/I");
	tree->Branch("motherPdgIdLeadGenQuark", &motherPdgIdLeadGenQuark, "motherPdgIdLeadGenQuark/I");
	tree->Branch("motherStatusLeadGenQuark", &motherStatusLeadGenQuark, "motherStatusLeadGenQuark/I");
	
	tree->Branch("etaSubleadGenQuark", &etaSubleadGenQuark, "etaSubleadGenQuark/F");
	tree->Branch("ptSubleadGenQuark", &ptSubleadGenQuark, "ptSubleadGenQuark/F");
	tree->Branch("phiSubleadGenQuark", &phiSubleadGenQuark, "phiSubleadGenQuark/F");
	tree->Branch("subleadGenQuarkNotFromScdHvyPtcl", &subleadGenQuarkNotFromScdHvyPtcl, "subleadGenQuarkNotFromScdHvyPtcl/I");
	tree->Branch("motherPdgIdSubleadGenQuark", &motherPdgIdSubleadGenQuark, "motherPdgIdSubleadGenQuark/I");
	tree->Branch("motherStatusSubleadGenQuark", &motherStatusSubleadGenQuark, "motherStatusSubleadGenQuark/I");
	

	tree->Branch("numGenFstHvyPtcl", &numGenFstHvyPtcl, "numGenFstHvyPtcl/I");
	tree->Branch("numGenScdHvyPtcl", &numGenScdHvyPtcl, "numGenScdHvyPtcl/I");
	tree->Branch("numGenLeptons", &numGenLeptons, "numGenLeptons/I");
	tree->Branch("numGenQuarks", &numGenQuarks, "numGenQuarks/I");

	tree->Branch("evWeight", &evWeight, "evWeight/F");
	tree->Branch("evWeightSign", &evWeightSign, "evWeightSign/F");

	///branches for RECO quantities
	if(_saveMatchedRecoInfo) {
		tree->Branch("ptRecoLeptMatchedToWrDau", &ptRecoLeptMatchedToWrDau, "ptRecoLeptMatchedToWrDau/F");
		tree->Branch("etaRecoLeptMatchedToWrDau", &etaRecoLeptMatchedToWrDau, "etaRecoLeptMatchedToWrDau/F");
		tree->Branch("phiRecoLeptMatchedToWrDau", &phiRecoLeptMatchedToWrDau, "phiRecoLeptMatchedToWrDau/F");
		tree->Branch("ptRecoLeptMatchedToNuDau", &ptRecoLeptMatchedToNuDau, "ptRecoLeptMatchedToNuDau/F");
		tree->Branch("etaRecoLeptMatchedToNuDau", &etaRecoLeptMatchedToNuDau, "etaRecoLeptMatchedToNuDau/F");
		tree->Branch("phiRecoLeptMatchedToNuDau", &phiRecoLeptMatchedToNuDau, "phiRecoLeptMatchedToNuDau/F");
		tree->Branch("ptRecoJetOneMatchedToNuDau", &ptRecoJetOneMatchedToNuDau, "ptRecoJetOneMatchedToNuDau/F");
		tree->Branch("etaRecoJetOneMatchedToNuDau", &etaRecoJetOneMatchedToNuDau, "etaRecoJetOneMatchedToNuDau/F");
		tree->Branch("phiRecoJetOneMatchedToNuDau", &phiRecoJetOneMatchedToNuDau, "phiRecoJetOneMatchedToNuDau/F");
		tree->Branch("ptRecoJetTwoMatchedToNuDau", &ptRecoJetTwoMatchedToNuDau, "ptRecoJetTwoMatchedToNuDau/F");
		tree->Branch("etaRecoJetTwoMatchedToNuDau", &etaRecoJetTwoMatchedToNuDau, "etaRecoJetTwoMatchedToNuDau/F");
		tree->Branch("phiRecoJetTwoMatchedToNuDau", &phiRecoJetTwoMatchedToNuDau, "phiRecoJetTwoMatchedToNuDau/F");

		///kinematics of gen jets linked with reco jets which were matched with gen quarks from WR and Nu decays
		tree->Branch("ptGenJetFromMatchedRecoJetOne", &ptGenJetFromMatchedRecoJetOne, "ptGenJetFromMatchedRecoJetOne/F");
		tree->Branch("etaGenJetFromMatchedRecoJetOne", &etaGenJetFromMatchedRecoJetOne, "etaGenJetFromMatchedRecoJetOne/F");
		tree->Branch("phiGenJetFromMatchedRecoJetOne", &phiGenJetFromMatchedRecoJetOne, "phiGenJetFromMatchedRecoJetOne/F");
		tree->Branch("ptGenJetFromMatchedRecoJetTwo", &ptGenJetFromMatchedRecoJetTwo, "ptGenJetFromMatchedRecoJetTwo/F");
		tree->Branch("etaGenJetFromMatchedRecoJetTwo", &etaGenJetFromMatchedRecoJetTwo, "etaGenJetFromMatchedRecoJetTwo/F");
		tree->Branch("phiGenJetFromMatchedRecoJetTwo", &phiGenJetFromMatchedRecoJetTwo, "phiGenJetFromMatchedRecoJetTwo/F");

		///kinematics of leading reco jets and leptons, no matching required
		tree->Branch("ptLeadRecoLept", &ptLeadRecoLept, "ptLeadRecoLept/F");
		tree->Branch("etaLeadRecoLept", &etaLeadRecoLept, "etaLeadRecoLept/F");
		tree->Branch("phiLeadRecoLept", &phiLeadRecoLept, "phiLeadRecoLept/F");
		tree->Branch("ptSubleadRecoLept", &ptSubleadRecoLept, "ptSubleadRecoLept/F");
		tree->Branch("etaSubleadRecoLept", &etaSubleadRecoLept, "etaSubleadRecoLept/F");
		tree->Branch("phiSubleadRecoLept", &phiSubleadRecoLept, "phiSubleadRecoLept/F");
		tree->Branch("ptLeadRecoJet", &ptLeadRecoJet, "ptLeadRecoJet/F");
		tree->Branch("etaLeadRecoJet", &etaLeadRecoJet, "etaLeadRecoJet/F");
		tree->Branch("phiLeadRecoJet", &phiLeadRecoJet, "phiLeadRecoJet/F");
		tree->Branch("ptSubleadRecoJet", &ptSubleadRecoJet, "ptSubleadRecoJet/F");
		tree->Branch("etaSubleadRecoJet", &etaSubleadRecoJet, "etaSubleadRecoJet/F");
		tree->Branch("phiSubleadRecoJet", &phiSubleadRecoJet, "phiSubleadRecoJet/F");

		recoJetCollectionToken = consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetCollection"));
		recoLeptonCollectionToken = consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("recoLeptonCollection"));
	}


	///tokens to input collections
	genParticleCollectionToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("genParticlesCollection"));
	genEventInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

}


genAndRecoWrAnalyzer::~genAndRecoWrAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
genAndRecoWrAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	resetTreeVars();
#ifdef DEBUG
	std::cout << "in analyze method of genAndRecoWrAnalyzer class" << std::endl;
#endif

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();
	evWeight = 1.0;
	evWeightSign = 1.0;

	iEvent.getByToken(genParticleCollectionToken, genParticleCollection);
	iEvent.getByToken(genEventInfoToken, genEvtInfo);	///< get evt weights if analyzing MC

#ifdef DEBUG
	std::cout << "called getByToken to fill genParticleCollection and genEvtInfo handles" << std::endl;
#endif

	if(genEvtInfo.isValid() ) {
		///real data does not have gen lvl event weights
		evWeight = genEvtInfo->weight();
		if(evWeight < 0) evWeightSign = -1.0;
	}

	if(genParticleCollection->size() == 0 ) {
		tree->Fill();
		return;
	}

	edm::OwnVector<reco::Candidate>::const_iterator genFstHvyPtcl = genParticleCollection->end(), genScdHvyPtcl = genParticleCollection->end(), genLeptFromFstHvyPtcl = genParticleCollection->end(), genLeptFromScdHvyPtcl = genParticleCollection->end(), genQuarkOneFromScdHvyPtcl = genParticleCollection->end(), genQuarkTwoFromScdHvyPtcl = genParticleCollection->end(), leadGenQrk = genParticleCollection->end(), subleadGenQrk = genParticleCollection->end(), leadGenLept = genParticleCollection->end(), subleadGenLept = genParticleCollection->end();

	bool photonRadiatedFromWrLepton = false;
	for(edm::OwnVector<reco::Candidate>::const_iterator it = genParticleCollection->begin(); it != genParticleCollection->end(); it++) {
		//find the two heavy particles in the gen decay chain
		if(std::fabs(it->pdgId()) == (Int_t) _firstHeavyParticlePdgId) {
			if(!_cutOnStatusCode) genFstHvyPtcl = it, numGenFstHvyPtcl++;

			else {
				if(std::fabs(it->status()) == (Int_t) _firstHeavyParticleStatus) {
					genFstHvyPtcl = it;
					numGenFstHvyPtcl++;
				}//end status match requirement
			}//end else
		}
		if(std::fabs(it->pdgId()) == (Int_t) _secondHeavyParticlePdgId && std::fabs( (it->mother(0))->pdgId() ) == (Int_t) _firstHeavyParticlePdgId) {
			if(!photonRadiatedFromWrLepton) genScdHvyPtcl = it, wrLeptonDidRadiatePhoton = 0;
			numGenScdHvyPtcl++;
		}

		if(std::fabs(it->pdgId()) == (Int_t) _secondHeavyParticlePdgId && std::fabs(it->status()) == 52 && std::fabs( (it->mother(0))->pdgId() ) == (Int_t) _secondHeavyParticlePdgId) {
			genScdHvyPtcl = it;
			numGenScdHvyPtcl++;
			photonRadiatedFromWrLepton = true;	//if this is reached, then dont change the genScdHvyPtcl iterator again
			wrLeptonDidRadiatePhoton = 1;	///<1 if true, 0 if false
		}
	


		//find the gen quarks which are decay products of the second heavy particle, and the leading and subleading gen quarks
		if(std::fabs(it->pdgId()) >= (Int_t) _minQuarkPdgId && std::fabs(it->pdgId()) <= (Int_t) _maxQuarkPdgId && std::fabs( (it->mother(0))->pdgId() ) == (Int_t) _secondHeavyParticlePdgId) {
#ifdef DEBUG
			std::cout << "found a quark whose mother is the secondHeavyParticle" << std::endl;
#endif
			//if genQuarkOne is not empty but genQuarkTwo is empty, then assign the ref to genQuarkTwo
			if(genQuarkOneFromScdHvyPtcl != genParticleCollection->end() && genQuarkTwoFromScdHvyPtcl == genParticleCollection->end()) genQuarkTwoFromScdHvyPtcl = it;

			//if both genQuark refs are empty, then assign the reference to genQuarkOne
			if(genQuarkOneFromScdHvyPtcl == genParticleCollection->end() && genQuarkTwoFromScdHvyPtcl == genParticleCollection->end()) genQuarkOneFromScdHvyPtcl = it;

		}///end gen quark selection

		//find the leading and subleading quarks
		if(std::fabs(it->pdgId()) >= (Int_t) _minQuarkPdgId && std::fabs(it->pdgId()) <= (Int_t) _maxQuarkPdgId && std::fabs(it->eta()) < 10 && it->pt() > 0.01) {
#ifdef DEBUG
			std::cout << " " << std::endl;
			std::cout << "gen quark pt=\t" << it->pt() << std::endl;
			std::cout << "gen quark eta=\t" << it->eta() << std::endl;
			std::cout << "gen quark phi=\t" << it->phi() << std::endl;
			std::cout << " " << std::endl;
#endif
			numGenQuarks++;
			if(leadGenQrk == genParticleCollection->end()) leadGenQrk = it;
			else {
				if(it->pt() > leadGenQrk->pt() && deltaR(it->eta(), it->phi(), leadGenQrk->eta(), leadGenQrk->phi() ) > 0.1) {
					subleadGenQrk = leadGenQrk;
					leadGenQrk = it;
				} else if(it->pt() > leadGenQrk->pt() && deltaR(it->eta(), it->phi(), leadGenQrk->eta(), leadGenQrk->phi() ) <= 0.1) leadGenQrk = it;
				else if((subleadGenQrk == genParticleCollection->end() || it->pt() > subleadGenQrk->pt()) && deltaR(it->eta(), it->phi(), leadGenQrk->eta(), leadGenQrk->phi() ) > 0.1) subleadGenQrk = it;

			}//end reassigning gen quark iterators
		}///end if(iterator points to a gen quark)

		//find the gen leptons which are decay products of the first and second heavy particles
		//and add all gen lepts to allGenLeptons vector
		if(std::fabs(it->pdgId()) == (Int_t) _leptonPdgId && std::fabs( (it->mother(0))->pdgId() ) == (Int_t) _secondHeavyParticlePdgId) genLeptFromScdHvyPtcl = it;
		if(std::fabs(it->pdgId()) == (Int_t) _leptonPdgId && std::fabs( (it->mother(0))->pdgId() ) == (Int_t) _firstHeavyParticlePdgId) genLeptFromFstHvyPtcl = it;
		if(std::fabs(it->pdgId()) == (Int_t) _leptonPdgId) {
			numGenLeptons++;
			if(leadGenLept == genParticleCollection->end()) leadGenLept = it;
			else {
				if(it->pt() > leadGenLept->pt() && deltaR(it->eta(), it->phi(), leadGenLept->eta(), leadGenLept->phi() ) > 0.1) {
					subleadGenLept = leadGenLept;
					leadGenLept = it;
				} else if(it->pt() > leadGenLept->pt() && deltaR(it->eta(), it->phi(), leadGenLept->eta(), leadGenLept->phi() ) <= 0.1) leadGenLept = it;
				else if((subleadGenLept == genParticleCollection->end() || it->pt() > subleadGenLept->pt()) && deltaR(it->eta(), it->phi(), leadGenLept->eta(), leadGenLept->phi() ) > 0.1) subleadGenLept = it;

			}//end reassigning gen lepton iterators
		}///end if(iterator points to a gen lepton)

	}///end loop over reco::GenParticle collection

	if(genFstHvyPtcl == genParticleCollection->end() || genScdHvyPtcl == genParticleCollection->end() ) {
		tree->Fill();
		return;
	}
#ifdef DEBUG
	std::cout << "found hvy unstable particles, and their daughter GEN leptons and quarks" << std::endl;
#endif
	///now calculate the desired kinematic quantities of all GEN objects, and save them to variables tied to tree branches
	etaGenFstHvyPtcl = genFstHvyPtcl->eta(), ptGenFstHvyPtcl = genFstHvyPtcl->pt(), massGenFstHvyPtcl = genFstHvyPtcl->mass();
	etaGenScdHvyPtcl = genScdHvyPtcl->eta(), ptGenScdHvyPtcl = genScdHvyPtcl->pt(), massGenScdHvyPtcl = genScdHvyPtcl->mass();

#ifdef DEBUG
	std::cout << "found wr and nu mass" << std::endl;
	std::cout << "wr mass=\t" << massGenFstHvyPtcl << std::endl;
	std::cout << "nu mass=\t" << massGenScdHvyPtcl << std::endl;
#endif

	if(genLeptFromFstHvyPtcl == genParticleCollection->end() || genLeptFromScdHvyPtcl == genParticleCollection->end() || genQuarkTwoFromScdHvyPtcl == genParticleCollection->end() || genQuarkOneFromScdHvyPtcl == genParticleCollection->end()) {
		std::cout << "a ref to one of the gen leptons or quarks coming from the WR decay points to a NULL object" << std::endl;
		tree->Fill();
		return;
	}
	//gen leptons and quarks with hvy particle mothers
	etaGenLeptFromFstHvyPtcl = genLeptFromFstHvyPtcl->eta(), ptGenLeptFromFstHvyPtcl = genLeptFromFstHvyPtcl->pt(), phiGenLeptFromFstHvyPtcl = genLeptFromFstHvyPtcl->phi();
	motherPdgIdGenLeptFromFstHvyPtcl = (genLeptFromFstHvyPtcl->mother(0))->pdgId(), motherStatusGenLeptFromFstHvyPtcl = (genLeptFromFstHvyPtcl->mother(0))->status();
	etaGenLeptFromScdHvyPtcl = genLeptFromScdHvyPtcl->eta(), ptGenLeptFromScdHvyPtcl = genLeptFromScdHvyPtcl->pt(), phiGenLeptFromScdHvyPtcl = genLeptFromScdHvyPtcl->phi();
	motherPdgIdGenLeptFromScdHvyPtcl = (genLeptFromScdHvyPtcl->mother(0))->pdgId(), motherStatusGenLeptFromScdHvyPtcl = (genLeptFromScdHvyPtcl->mother(0))->status();
	etaGenQuarkTwoFromScdHvyPtcl = genQuarkTwoFromScdHvyPtcl->eta(), ptGenQuarkTwoFromScdHvyPtcl = genQuarkTwoFromScdHvyPtcl->pt(), phiGenQuarkTwoFromScdHvyPtcl = genQuarkTwoFromScdHvyPtcl->phi();
	motherPdgIdGenQuarkTwoFromScdHvyPtcl = (genQuarkTwoFromScdHvyPtcl->mother(0))->pdgId(), motherStatusGenQuarkTwoFromScdHvyPtcl = (genQuarkTwoFromScdHvyPtcl->mother(0))->status();
	etaGenQuarkOneFromScdHvyPtcl = genQuarkOneFromScdHvyPtcl->eta(), ptGenQuarkOneFromScdHvyPtcl = genQuarkOneFromScdHvyPtcl->pt(), phiGenQuarkOneFromScdHvyPtcl = genQuarkOneFromScdHvyPtcl->phi();
	motherPdgIdGenQuarkOneFromScdHvyPtcl = (genQuarkOneFromScdHvyPtcl->mother(0))->pdgId(), motherStatusGenQuarkOneFromScdHvyPtcl = (genQuarkOneFromScdHvyPtcl->mother(0))->status();

	TLorentzVector genLeptFromFstHvyPtclFourMom, genLeptFromScdHvyPtclFourMom, genQuarkOneFromScdHvyPtclFourMom, genQuarkTwoFromScdHvyPtclFourMom;
	genLeptFromFstHvyPtclFourMom.SetPtEtaPhiM(ptGenLeptFromFstHvyPtcl, etaGenLeptFromFstHvyPtcl, phiGenLeptFromFstHvyPtcl, genLeptFromFstHvyPtcl->mass());
	genLeptFromScdHvyPtclFourMom.SetPtEtaPhiM(ptGenLeptFromScdHvyPtcl, etaGenLeptFromScdHvyPtcl, phiGenLeptFromScdHvyPtcl, genLeptFromScdHvyPtcl->mass());
	genQuarkOneFromScdHvyPtclFourMom.SetPtEtaPhiM(ptGenQuarkOneFromScdHvyPtcl, etaGenQuarkOneFromScdHvyPtcl, phiGenQuarkOneFromScdHvyPtcl, genQuarkOneFromScdHvyPtcl->mass());
	genQuarkTwoFromScdHvyPtclFourMom.SetPtEtaPhiM(ptGenQuarkTwoFromScdHvyPtcl, etaGenQuarkTwoFromScdHvyPtcl, phiGenQuarkTwoFromScdHvyPtcl, genQuarkTwoFromScdHvyPtcl->mass());
	
	fourObjMassFromGenObjsFromFstAndScdHvyPtcl = (genLeptFromFstHvyPtclFourMom + genLeptFromScdHvyPtclFourMom + genQuarkOneFromScdHvyPtclFourMom + genQuarkTwoFromScdHvyPtclFourMom).M();

#ifdef DEBUG
	std::cout << "saved kinematic info from GEN lepton and quark daughters of the two hvy particles" << std::endl;
#endif

	//////////////////////
	//gen leptons and quarks with no mother requirements
	leadGenLeptonNotFromFstHvyPtcl = 0, subleadGenLeptonNotFromScdHvyPtcl = 0;
	leadGenQuarkNotFromScdHvyPtcl = 0, subleadGenQuarkNotFromScdHvyPtcl = 0;

	///one or more of the four const_iterators could point to a null reference if the dR > 0.1 cut is failed in findLeadingAndSubleading()
	if(leadGenQrk == genParticleCollection->end() || subleadGenQrk == genParticleCollection->end() || leadGenLept == genParticleCollection->end() || subleadGenLept == genParticleCollection->end() ) {
		std::cout << "a ref to one of the leading gen leptons or quarks points to a NULL object" << std::endl;
		tree->Fill();
		return;
	}

	//now save the kinematic info from the leading quarks and leptons into the tree
	ptLeadGenLepton = leadGenLept->pt(), etaLeadGenLepton = leadGenLept->eta(), phiLeadGenLepton = leadGenLept->phi();
	motherPdgIdLeadGenLepton = (leadGenLept->mother(0))->pdgId(), motherStatusLeadGenLepton = (leadGenLept->mother(0))->status();
	if(leadGenLept != genLeptFromFstHvyPtcl) leadGenLeptonNotFromFstHvyPtcl = 1;
	
	ptSubleadGenLepton = subleadGenLept->pt(), etaSubleadGenLepton = subleadGenLept->eta(), phiSubleadGenLepton = subleadGenLept->phi();
	motherPdgIdSubleadGenLepton = (subleadGenLept->mother(0))->pdgId(), motherStatusSubleadGenLepton = (subleadGenLept->mother(0))->status();
	if(subleadGenLept != genLeptFromScdHvyPtcl) subleadGenLeptonNotFromScdHvyPtcl = 1;

	ptLeadGenQuark = leadGenQrk->pt(), etaLeadGenQuark = leadGenQrk->eta(), phiLeadGenQuark = leadGenQrk->phi();
	motherPdgIdLeadGenQuark = (leadGenQrk->mother(0))->pdgId(), motherStatusLeadGenQuark = (leadGenQrk->mother(0))->status();
	if(leadGenQrk != genQuarkOneFromScdHvyPtcl && leadGenQrk != genQuarkTwoFromScdHvyPtcl) leadGenQuarkNotFromScdHvyPtcl = 1;

	ptSubleadGenQuark = subleadGenQrk->pt(), etaSubleadGenQuark = subleadGenQrk->eta(), phiSubleadGenQuark = subleadGenQrk->phi();
	motherPdgIdSubleadGenQuark = (subleadGenQrk->mother(0))->pdgId(), motherStatusSubleadGenQuark = (subleadGenQrk->mother(0))->status();
	if(subleadGenQrk != genQuarkOneFromScdHvyPtcl && subleadGenQrk != genQuarkTwoFromScdHvyPtcl) subleadGenQuarkNotFromScdHvyPtcl = 1;

#ifdef DEBUG
	std::cout << "saved kinematic info from leading gen quarks and leptons" << std::endl;
#endif
	if(_saveMatchedRecoInfo) {
#ifdef DEBUG
		std::cout << "in RECO matching if statement" << std::endl;
#endif

		iEvent.getByToken(recoJetCollectionToken, recoJetCollection);
		iEvent.getByToken(recoLeptonCollectionToken, recoLeptonCollection);

		std::vector<pat::Electron> patLepts;
		unsigned int maxLepts = recoLeptonCollection->size();
		for(unsigned int it = 0; it < maxLepts; it++) {
			reco::CandidateBaseRef tempRef = recoLeptonCollection->refAt(it);
			const pat::Electron& eleTemp = dynamic_cast<const pat::Electron&>(*tempRef);
			patLepts.push_back(eleTemp);
		}

#ifdef DEBUG
		std::cout << "there are\t" << patLepts.size() << "\treco leptons" << std::endl;
#endif
		if(patLepts.size() == 0) {
			tree->Fill();
			return;
		}

		///refs to leptons and jets which are dR matched to WR and Nu daughter leptons and quarks
		//std::vector<pat::Jet>::const_iterator recoJetOneFromScdHvyPtcl = recoJetCollection->end(), recoJetTwoFromScdHvyPtcl = recoJetCollection->end();
		Bool_t filledFirstMatchedRecoJet = false;

		for(std::vector<pat::Jet>::const_iterator it = recoJetCollection->begin(); it != recoJetCollection->end(); it++) {
#ifdef DEBUG
			std::cout << "in loop over recoJetCollection objs" << std::endl;
#endif

			if( it->genParton() == NULL ) continue;
			if(std::fabs( (it->genParton()->mother(0))->pdgId() ) == (Int_t) _secondHeavyParticlePdgId && std::fabs( it->genParton()->pdgId() ) >= (Int_t) _minQuarkPdgId && std::fabs( it->genParton()->pdgId() ) <= (Int_t) _maxQuarkPdgId) {
#ifdef DEBUG
				std::cout << "found pat::Jet matched to a gen quark whose mother is secondHeavyParticle" << std::endl;
				std::cout << "reco jet pt\t" << it->pt() << "\tjet eta\t" << it->eta() << std::endl;
				std::cout << "gen quark pt\t" << it->genParton()->pt() << "\tjet eta\t" << it->genParton()->eta() << std::endl;
#endif
				if(filledFirstMatchedRecoJet) {
					ptRecoJetTwoMatchedToNuDau = it->pt(), etaRecoJetTwoMatchedToNuDau = it->eta(), phiRecoJetTwoMatchedToNuDau = it->phi();
#ifdef DEBUG
					std::cout << "about to access gen jet matched to second reco jet and gen quark" << std::endl;
#endif
					if(it->genJet() != NULL) {
						ptGenJetFromMatchedRecoJetTwo = it->genJet()->pt(), etaGenJetFromMatchedRecoJetTwo = it->genJet()->eta(), phiGenJetFromMatchedRecoJetTwo = it->genJet()->phi();
					}
				}//end saving kinematic info from second matched reco jet

				if(!filledFirstMatchedRecoJet) {
					filledFirstMatchedRecoJet = true;
					ptRecoJetOneMatchedToNuDau = it->pt(), etaRecoJetOneMatchedToNuDau = it->eta(), phiRecoJetOneMatchedToNuDau = it->phi();
#ifdef DEBUG
					std::cout << "about to access gen jet matched to first reco jet and gen quark" << std::endl;
#endif
					if(it->genJet() != NULL) ptGenJetFromMatchedRecoJetOne = it->genJet()->pt(), etaGenJetFromMatchedRecoJetOne = it->genJet()->eta(), phiGenJetFromMatchedRecoJetOne = it->genJet()->phi();
				}//end saving kinematic info from first matched reco jet

			}//end saving info from reco jets and gen jets matched to gen quarks

		}///end loop over reco jet collection
#ifdef DEBUG
		std::cout << "finished saving matched reco jet kinematic info" << std::endl;
#endif

		for(std::vector<pat::Electron>::const_iterator it = patLepts.cbegin(); it != patLepts.cend(); it++) {
			if(it->genLepton() == NULL) continue;
			if(std::fabs( (it->genLepton()->mother(0))->pdgId() ) == (Int_t) _secondHeavyParticlePdgId && std::fabs( it->genLepton()->pdgId() ) == (Int_t) _leptonPdgId) {
				ptRecoLeptMatchedToNuDau = it->pt(), etaRecoLeptMatchedToNuDau = it->eta(), phiRecoLeptMatchedToNuDau = it->phi();
			}//end reco lepton matched to Nu daughter lepton

			if(std::fabs( (it->genLepton()->mother(0))->pdgId() ) == (Int_t) _firstHeavyParticlePdgId && std::fabs( it->genLepton()->pdgId() ) == (Int_t) _leptonPdgId) {
				ptRecoLeptMatchedToWrDau = it->pt(), etaRecoLeptMatchedToWrDau = it->eta(), phiRecoLeptMatchedToWrDau = it->phi();
			}//end reco lepton matched to WR daughter lepton
		}///end loop over reco lepton collection

#ifdef DEBUG
		std::cout << "saved kinematic info of reco leptons matched to GEN Nu and WR" << std::endl;
#endif
		std::vector<pat::Electron>::const_iterator recoLeadLept = patLepts.cend(), recoSubleadLept = patLepts.cend();
		std::vector<pat::Jet>::const_iterator recoLeadJet = recoJetCollection->end(), recoSubleadJet = recoJetCollection->end();

#ifdef DEBUG
		std::cout << "about to call findLeadingAndSubleadingLepton()" << std::endl;
#endif
		findLeadingAndSubleadingLepton(recoLeadLept, recoSubleadLept, patLepts);

#ifdef DEBUG
		std::cout << "about to call findLeadingAndSubleadingJet()" << std::endl;
#endif
		findLeadingAndSubleadingJet(recoLeadJet, recoSubleadJet, recoJetCollection);
		if(recoLeadJet == recoJetCollection->end() || recoSubleadJet == recoJetCollection->end() || recoLeadLept == patLepts.end() || recoSubleadLept == patLepts.end() ) {
			std::cout << "a ref to one of the leading reco leptons or quarks points to a NULL object" << std::endl;
			tree->Fill();
			return;
		}
		ptLeadRecoLept = recoLeadLept->pt(), etaLeadRecoLept = recoLeadLept->eta(), phiLeadRecoLept = recoLeadLept->phi();
		ptSubleadRecoLept = recoSubleadLept->pt(), etaSubleadRecoLept = recoSubleadLept->eta(), phiSubleadRecoLept = recoSubleadLept->phi();
		ptLeadRecoJet = recoLeadJet->pt(), etaLeadRecoJet = recoLeadJet->eta(), phiLeadRecoJet = recoLeadJet->phi();
		ptSubleadRecoJet = recoSubleadJet->pt(), etaSubleadRecoJet = recoSubleadJet->eta(), phiSubleadRecoJet = recoSubleadJet->phi();

	}///end filling matched RECO branches

	tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
genAndRecoWrAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
genAndRecoWrAnalyzer::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void genAndRecoWrAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(genAndRecoWrAnalyzer);
