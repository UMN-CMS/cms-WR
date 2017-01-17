// -*- C++ -*-
//
// Package:    testFilter/htFilter
// Class:      htFilter
//
/**\class htFilter htFilter.cc testFilter/htFilter/plugins/htFilter.cc

 Description: [one line class summary]

 filter events based on an HT threshold
 where HT = sum of pt of all GEN quarks and gluons leaving the hard interaction

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  sean kalafut
//         Created:  Tue, 04 Aug 2015 22:09:29 GMT
//
//



// system include files
#include <memory>
#include <cstring>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

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
#include "DataFormats/Math/interface/LorentzVector.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "TTree.h"
#include <TFile.h>
#include <TBranch.h>


//#define DEBUG
//
// class declaration
//

class htFilter : public edm::EDFilter
{
public:
	explicit htFilter(const edm::ParameterSet&);
	~htFilter();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	virtual bool filter(edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

	// ----------member data ---------------------------
	///Handles to input object collections
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > inputParticles;

	///tokens to input collections
	edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputParticlesToken;
	
	double _threshold;	//threshold for cut
	bool _thresholdIsLowerBound;	//threshold indicates lower bound

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
htFilter::htFilter(const edm::ParameterSet& iConfig):
	_threshold(iConfig.getParameter<double>("cutThreshold")),
	_thresholdIsLowerBound(iConfig.getParameter<double>("isLowerBound"))

{
	//now do what ever initialization is needed
	inputParticlesToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputCollection"));

}


htFilter::~htFilter()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
htFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

#ifdef DEBUG
	std::cout << "in filter method of htFilter class" << std::endl;
#endif

	double sumPt = 0;
	iEvent.getByToken(inputParticlesToken, inputParticles);
	
	for(edm::OwnVector<reco::Candidate>::const_iterator it = inputParticles->begin(); it!=inputParticles->end(); it++){
		sumPt += it->pt();
	}//end loop over particles in input collection

	if(_thresholdIsLowerBound){
		if(sumPt < _threshold) return false;
	} else{
		if(sumPt > _threshold) return false;
	}

	return true;
}

// ------------ method called once each job just before starting event loop  ------------
void
htFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
htFilter::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
htFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
htFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
htFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
htFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
htFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(htFilter);
