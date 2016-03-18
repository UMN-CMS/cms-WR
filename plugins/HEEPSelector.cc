
// -*- C++ -*-
//
// Package:    ExoAnalysis/cmsWR
// Class:      HEEPIDSelector
//
/**\class HEEPIDSelector HEEPIDSelector.cc ExoAnalysis/cmsWR/plugins/HEEPIDSelector.cc
 Description: Producer that creats a collection that passes the HEEPID
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Peter Hansen
//
//         Created:  Thu Aug  6 14:17:04 CDT 2015
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

class HEEPIDSelector : public edm::EDProducer
{
public:
	explicit HEEPIDSelector(const edm::ParameterSet&);
	~HEEPIDSelector();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void produce(edm::Event&, const edm::EventSetup&) override;

	// ----------member data ---------------------------

	edm::EDGetToken electronsToken_;
	edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;


};

//
// constructors and destructor
//
HEEPIDSelector::HEEPIDSelector(const edm::ParameterSet& iConfig):
	eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap")))
{

	electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
	                     (iConfig.getParameter<edm::InputTag>
	                      ("electrons"));

	produces<pat::ElectronCollection>();
}


HEEPIDSelector::~HEEPIDSelector()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HEEPIDSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// Retrieve the collection of electrons from the event.
	edm::Handle<edm::View<reco::GsfElectron> > electrons;
	iEvent.getByToken(electronsToken_, electrons);

	// Get the electron ID data from the event stream.
	edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
	iEvent.getByToken(eleHEEPIdMapToken_ , heep_id_decisions);

	std::unique_ptr<pat::ElectronCollection> elec_out(new pat::ElectronCollection());

	// Loop over electrons
	for (size_t i = 0; i < electrons->size(); ++i) {
		const auto el = electrons->ptrAt(i);

		if ((*heep_id_decisions)[el])
			elec_out->push_back(el);

	}
	iEvent.put(std::move(elec_out));

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HEEPIDSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.add<edm::InputTag>("electrons", edm::InputTag("slimmedElectrons"));
	desc.add<edm::InputTag>("eleHEEPIdMap", edm::InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"));
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HEEPIDSelector);
