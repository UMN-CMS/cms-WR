
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

class heepIDAnalyzer : public edm::EDAnalyzer
{
public:
	explicit heepIDAnalyzer(const edm::ParameterSet&);
	~heepIDAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

	// ----------member data ---------------------------

	// MiniAOD case data members
	edm::EDGetToken electronsMiniAODToken_;

	// ID decisions objects
	edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
	edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleHEEPIdMapCFRToken_;


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
	// MiniAOD tokens
	// For electrons, use the fact that pat::Electron can be cast into
	// GsfElectron
	electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
	                            (iConfig.getParameter<edm::InputTag>
	                             ("electronsMiniAOD"));
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


	// Retrieve the collection of electrons from the event.
	edm::Handle<edm::View<reco::GsfElectron> > electrons;
	iEvent.getByToken(electronsMiniAODToken_, electrons);

	// Get the electron ID data from the event stream.
	// Note: this implies that the VID ID modules have been run upstream.
	// If you need more info, check with the EGM group.
	edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
	iEvent.getByToken(eleHEEPIdMapToken_ , heep_id_decisions);

	edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep_id_CFR;
	iEvent.getByToken(eleHEEPIdMapCFRToken_ , heep_id_CFR);

	// Loop over electrons

	std::cout << "elec: " << electrons->size() << std::endl;
	for (size_t i = 0; i < electrons->size(); ++i) {
		const auto el = electrons->ptrAt(i);

		// Kinematics
		if( el->pt() < 10 ) // keep only electrons above 10 GeV
			continue;

		cout << "[DEBUG]" << (*heep_id_CFR)[el].cutFlowName() << std::endl;
		for (size_t i = 0; i < (*heep_id_CFR)[el].cutFlowSize(); i++) {
			cout << "[DEBUG]" << i << ' ' << (*heep_id_CFR)[el].getNameAtIndex(i) << std::endl;
			cout << "[DEBUG]" << i << ' ' << (*heep_id_CFR)[el].getCutResultByIndex(i) << std::endl;
			cout << "[DEBUG]" << i << ' ' << (*heep_id_CFR)[el].getValueCutUpon(i) << std::endl;
		}

	}

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
heepIDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(heepIDAnalyzer);
