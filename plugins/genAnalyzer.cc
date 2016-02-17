// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/genAnalyzer
// Class:      genAnalyzer
//
/**\class genAnalyzer genAnalyzer.cc doubleElectronTracklessTrigger/genAnalyzer/plugins/genAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sean Kalafut
//         Created:  Thu, 06 Nov 2014 23:16:33 GMT
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
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include "TTree.h"
#include "TLorentzVector.h"

//#define DEBUG

//
// class declaration
//

class genAnalyzer : public edm::EDAnalyzer
{
public:
	explicit genAnalyzer(const edm::ParameterSet&);
	~genAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collection)
	{

#ifdef DEBUG
		std::cout << "checking pt of all gen particles in findLeadingAndSubleading fxn" << std::endl;
#endif

		for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++) {
#ifdef DEBUG
			std::cout << "pT = \t" << genIt->pt() << std::endl;
#endif
			if(first == collection->end()) first = genIt;
			else {
				if(genIt->pt() > first->pt()) {
					second = first;
					first = genIt;
				} else if(second == collection->end() || genIt->pt() > second->pt()) second = genIt;
			}
		}//end loop over reco::Candidate collection

#ifdef DEBUG
		std::cout << "leaving findLeadingAndSubleading fxn" << std::endl;
#endif

	}///end findLeadingAndSubleading()


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

	edm::InputTag genLeptonCollTag;
	edm::InputTag genJetCollTag;

//Handles to GEN object collections
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > leptons;
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > jets;

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

genAnalyzer::genAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName")),
	genLeptonCollTag(iConfig.getParameter<edm::InputTag>("genLeptonCollection")),
	genJetCollTag(iConfig.getParameter<edm::InputTag>("genJetCollection"))

{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>(tName.c_str(), "GEN electron kinematic info");

	tree->Branch("etaGenEle", etaGenEle, "etaGenEle[2]/F");
	tree->Branch("ptGenEle", ptGenEle, "ptGenEle[2]/F");
	tree->Branch("phiGenEle", phiGenEle, "phiGenEle[2]/F");
	tree->Branch("dileptonMassGen", &dileptonMassGen, "dileptonMassGen/F");

	tree->Branch("evtNumber", &evtNumber, "evtNumber/l");
	tree->Branch("runNumber", &runNumber, "runNumber/I");

	tree->Branch("etaGenJet", etaGenJet, "etaGenJet[2]/F");
	tree->Branch("ptGenJet", ptGenJet, "ptGenJet[2]/F");
	tree->Branch("phiGenJet", phiGenJet, "phiGenJet[2]/F");
	tree->Branch("dijetMassGen", &dijetMassGen, "dijetMassGen/F");

	tree->Branch("fourObjectMassGen", &fourObjectMassGen, "fourObjectMassGen/F");
	tree->Branch("leadLeptonThreeObjMassGen", &leadLeptonThreeObjMassGen, "leadLeptonThreeObjMassGen/F");
	tree->Branch("subleadingLeptonThreeObjMassGen", &subleadingLeptonThreeObjMassGen, "subleadingLeptonThreeObjMassGen/F");

	tree->Branch("dR_leadingLeptonLeadingJetGen", &dR_leadingLeptonLeadingJetGen, "dR_leadingLeptonLeadingJetGen/F");
	tree->Branch("dR_leadingLeptonSubleadingJetGen", &dR_leadingLeptonSubleadingJetGen, "dR_leadingLeptonSubleadingJetGen/F");
	tree->Branch("dR_subleadingLeptonLeadingJetGen", &dR_subleadingLeptonLeadingJetGen, "dR_subleadingLeptonLeadingJetGen/F");
	tree->Branch("dR_subleadingLeptonSubleadingJetGen", &dR_subleadingLeptonSubleadingJetGen, "dR_subleadingLeptonSubleadingJetGen/F");

}


genAnalyzer::~genAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
genAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();

	iEvent.getByLabel(genLeptonCollTag, leptons);
	iEvent.getByLabel(genJetCollTag, jets);

	edm::OwnVector<reco::Candidate>::const_iterator leadingLepton = leptons->end(), subleadingLepton = leptons->end();
	edm::OwnVector<reco::Candidate>::const_iterator leadingJet = jets->end(), subleadingJet = jets->end();

	findLeadingAndSubleading(leadingLepton, subleadingLepton, leptons);
	findLeadingAndSubleading(leadingJet, subleadingJet, jets);

	///now that the leading and subleading leptons and jets have been found, fill all of the arrays and single Float_ values
	///which will be saved into the tree
#ifdef DEBUG
	std::cout << "leading lepton pt: \t" << leadingLepton->pt() << std::endl;
	std::cout << "subleading lepton pt: \t" << subleadingLepton->pt() << std::endl;
#endif

	etaGenEle[0] = leadingLepton->eta();
	ptGenEle[0] = leadingLepton->pt();
	phiGenEle[0] = leadingLepton->phi();
	etaGenEle[1] = subleadingLepton->eta();
	ptGenEle[1] = subleadingLepton->pt();
	phiGenEle[1] = subleadingLepton->phi();

#ifdef DEBUG
	std::cout << "leading jet pt: \t" << leadingJet->pt() << std::endl;
	std::cout << "subleading jet pt: \t" << subleadingJet->pt() << std::endl;
#endif

	etaGenJet[0] = leadingJet->eta();
	ptGenJet[0] = leadingJet->pt();
	phiGenJet[0] = leadingJet->phi();
	etaGenJet[1] = subleadingJet->eta();
	ptGenJet[1] = subleadingJet->pt();
	phiGenJet[1] = subleadingJet->phi();

	///make TLorentzVector objects for the four GEN objects.  Then use these LorentzVectors to calculate three and four object
	///invariant mass values.
	TLorentzVector l1(ptGenEle[0]*TMath::Cos(phiGenEle[0]), ptGenEle[0]*TMath::Sin(phiGenEle[0]), ptGenEle[0]*TMath::SinH(etaGenEle[0]), ptGenEle[0]*TMath::CosH(etaGenEle[0]));	///leading lepton lorentz vector (px, py, pz, E)
	TLorentzVector l2(ptGenEle[1]*TMath::Cos(phiGenEle[1]), ptGenEle[1]*TMath::Sin(phiGenEle[1]), ptGenEle[1]*TMath::SinH(etaGenEle[1]), ptGenEle[1]*TMath::CosH(etaGenEle[1]));
	TLorentzVector j1(ptGenJet[0]*TMath::Cos(phiGenJet[0]), ptGenJet[0]*TMath::Sin(phiGenJet[0]), ptGenJet[0]*TMath::SinH(etaGenJet[0]), ptGenJet[0]*TMath::CosH(etaGenJet[0]));	///leading jet lorentz vector (px, py, pz, E)
	TLorentzVector j2(ptGenJet[1]*TMath::Cos(phiGenJet[1]), ptGenJet[1]*TMath::Sin(phiGenJet[1]), ptGenJet[1]*TMath::SinH(etaGenJet[1]), ptGenJet[1]*TMath::CosH(etaGenJet[1]));

	dileptonMassGen = (l1 + l2).M();
	dijetMassGen = (j1 + j2).M();
	fourObjectMassGen = (l1 + l2 + j1 + j2).M();
	leadLeptonThreeObjMassGen = (l1 + j1 + j2).M();
	subleadingLeptonThreeObjMassGen = (l2 + j1 + j2).M();

	///now use the individual GEN object eta and phi values to calculate dR between different (lepton, jet) pairs
	dR_leadingLeptonLeadingJetGen = deltaR(etaGenEle[0], phiGenEle[0], etaGenJet[0], phiGenJet[0]);
	dR_leadingLeptonSubleadingJetGen = deltaR(etaGenEle[0], phiGenEle[0], etaGenJet[1], phiGenJet[1]);
	dR_subleadingLeptonLeadingJetGen = deltaR(etaGenEle[1], phiGenEle[1], etaGenJet[0], phiGenJet[0]);
	dR_subleadingLeptonSubleadingJetGen = deltaR(etaGenEle[1], phiGenEle[1], etaGenJet[1], phiGenJet[1]);

	tree->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example", pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
genAnalyzer::beginJob()
{
	/*
	  tree_file = new TFile(foutName.c_str(), "recreate");
	  if(tree_file->IsZombie()){
	    throw cms::Exception("OutputError") <<  "Output tree not created (Zombie): " << foutName;
	    return;
	  }
	  tree_file->cd();

	  //now do what ever initialization is needed
	  tree = new TTree("selected","selected");
	  tree->SetDirectory(tree_file);

	  //InitNewTree();
	*/

}

// ------------ method called once each job just after ending the event loop  ------------
void
genAnalyzer::endJob()
{
	//loop over bins of "EventFraction", divide each bin content by totalNumEvents, then reset the bin content to the old content divided by totalNumEvents

	/*
	for(int i=1; i<=getXBins("EventFraction"); i++){
		if( getXBins("EventFraction") < 3) break;	//shouldn't need this, but just in case

		std::cout<<"bin # "<< i <<" content equals "<< get1DBinContents("EventFraction",i) <<std::endl;
		set1DBinContents("EventFraction",i, (get1DBinContents("EventFraction",i)/getTotalNumEvents() ) );
		std::cout<<"bin # "<< i <<" content equals "<< get1DBinContents("EventFraction",i) <<std::endl;

	}

	std::cout<< "the trackless leg of trigger fired on "<< getNumTriggeredEvents() << " events out of "<< getEfficiencyDenominator() << " total events which should have fired trackless leg of trigger" <<std::endl;
	set1DBinContents("HLTRecoEff",1, getNumTriggeredEvents()/getEfficiencyDenominator());

	*/


	/*
	//save the tree into the file
	tree_file->cd();
	tree->Write();
	tree_file->Close();
	*/

}

// ------------ method called when starting to processes a run  ------------
/*
void
genAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
genAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
genAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
genAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
genAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

/*
void genAnalyzer::InitNewTree(void){

  //make one branch for each unique variable I want to track - ecal iso, lepton pT, invariant mass of dilepton system, etc

  std::cout << "[STATUS] InitNewTree" << std::endl;
  if(tree==NULL) return;
  tree->Branch("runNumber",     &runNumber,     "runNumber/I");
  tree->Branch("eventNumber",   &eventNumber, "eventNumber/l");
  tree->Branch("lumiBlock",     &lumiBlock,     "lumiBlock/I");
  tree->Branch("runTime",       &runTime,         "runTime/i");

  tree->Branch("mcGenWeight",   &mcGenWeight, "mcGenWeight/F");

  tree->Branch("nPU", nPU, "nPU[1]/I");
  tree->Branch("rho", &rho, "rho/F");
  tree->Branch("nPV", &nPV, "nPV/I");


  tree->Branch("chargeEle",   chargeEle,    "chargeEle[2]/I");	//[nEle]
  tree->Branch("etaSCEle",      etaSCEle,       "etaSCEle[2]/F");	//[nSCEle]
  tree->Branch("phiSCEle",      phiSCEle,       "phiSCEle[2]/F");	//[nSCEle]

  tree->Branch("PtEle",       PtEle,        "PtEle[2]/F");

  tree->Branch("seedXSCEle",           seedXSCEle,      "seedXSCEle[2]/F");
  tree->Branch("seedYSCEle",           seedYSCEle,      "seedYSCEle[2]/F");
  tree->Branch("seedEnergySCEle", seedEnergySCEle, "seedEnergySCEle[2]/F");

  tree->Branch("gainEle", gainEle, "gainEle[2]/b");

  tree->Branch("energyMCEle", energyMCEle, "energyMCEle[2]/F");
  tree->Branch("energySCEle", energySCEle, "energySCEle[2]/F");
  tree->Branch("rawEnergySCEle", rawEnergySCEle, "rawEnergySCEle[2]/F");
  tree->Branch("esEnergySCEle", esEnergySCEle, "esEnergySCEle[2]/F");


  tree->Branch("R9Ele", R9Ele, "R9Ele[2]/F");

  tree->Branch("e5x5SCEle", e5x5SCEle, "e5x5SCEle[2]/F");

  tree->Branch("invMass",    &invMass,      "invMass/F");   // invariant mass ele+SC
  tree->Branch("invMass_SC", &invMass_SC,   "invMass_SC/F"); // invariant mass SC+SC


  tree->Branch("invMass_MC", &invMass_MC, "invMass_MC/F");

  tree->Branch("etaMCEle",      etaMCEle,       "etaMCEle[2]/F");	//[nEle]
  tree->Branch("phiMCEle",      phiMCEle,       "phiMCEle[2]/F");	//[nEle]

  tree->Branch("nHitsSCEle", nHitsSCEle, "nHitsSCEle[2]/I");

  tree->Branch("sigmaIEtaIEtaSCEle", sigmaIEtaIEtaSCEle, "sigmaIEtaIEtaSCEle[2]/F");
  tree->Branch("sigmaIEtaIEtaSCEle", sigmaIEtaIEtaSCEle, "sigmaIEtaIEtaSCEle[2]/F");

  return;
}

//negative index means the corresponding electron does not exist
void genAnalyzer::TreeSetSingleElectronVar(const pat::Electron& electron1, int index){

  if(index<0){
    PtEle[-index] 	  = 0;
    chargeEle[-index] = 0;
    etaEle[-index]    = 0;
    phiEle[-index]    = 0;
    return;
  }

  PtEle[index]     = electron1.et();
  chargeEle[index] = electron1.charge();
  etaEle[index]    = electron1.eta();
  phiEle[index]    = electron1.phi();
}

void genAnalyzer::TreeSetSingleElectronVar(const reco::SuperCluster& electron1, int index){

  if(index<0){
    PtEle[-index] 	  = 0;
    chargeEle[-index] = 0;
    etaEle[-index]    = 0;
    phiEle[-index]    = 0;
    return;
  }

//checks

  PtEle[index]     = electron1.energy()/cosh(electron1.eta());
  chargeEle[index] = -100; // dont know the charge for SC
  etaEle[index]    = electron1.eta(); // eta = etaSC
  phiEle[index]    = electron1.phi();
}

void genAnalyzer::TreeSetDiElectronVar(const pat::Electron& electron1, const reco::SuperCluster& electron2){

  TreeSetSingleElectronVar(electron1, 0);
  TreeSetSingleElectronVar(electron2, 1);

  double t1=TMath::Exp(-etaEle[0]);
  double t2=TMath::Exp(-etaEle[1]);
  double t1q = t1*t1;
  double t2q = t2*t2;

  double angle=1- ( (1-t1q)*(1-t2q)+4*t1*t2*cos(phiEle[0]-phiEle[1]))/( (1+t1q)*(1+t2q) );


  invMass = sqrt(2*electron1.energy()*electron2.energy() *angle);
  invMass_e5x5   = sqrt(2*electron1.e5x5()*(clustertools->e5x5(*electron2.seed())) * angle);

}
*/


//define this as a plug-in
DEFINE_FWK_MODULE(genAnalyzer);
