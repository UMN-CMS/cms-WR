// -*- C++ -*-
//
// Package:    ExoAnalysis/produceJECUnc
// Class:      produceJECUnc
// 
/**\class produceJECUnc produceJECUnc.cc ExoAnalysis/produceJECUnc/plugins/produceJECUnc.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Peter Hansen
//         Created:  Fri, 29 Jan 2016 04:32:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/Common/interface/ValueMap.h"

typedef edm::ValueMap<float> JECUnc_Map;

//
// class declaration
//

class produceJECUnc : public edm::stream::EDProducer<> {
   public:
      explicit produceJECUnc(const edm::ParameterSet&);
      ~produceJECUnc();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::InputTag src;        ///<input particle objects
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
produceJECUnc::produceJECUnc(const edm::ParameterSet& iConfig)
        : src(iConfig.getParameter<edm::InputTag>("src"))
{
        produces<JECUnc_Map>("JECUnc");
}


produceJECUnc::~produceJECUnc()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
produceJECUnc::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	ESHandle<JetCorrectorParametersCollection> JetCorParColl;
	iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
	JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
	JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

	float eta = 0.0;
	float ptCor = 80.0;
	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(ptCor); // here you must use the CORRECTED jet pt
	double unc = jecUnc->getUncertainty(true);
	double ptCor_shifted = ptCor*(1+unc);

	std::cout << "eta: " << eta << " pt: " << ptCor << " unc: " << unc << " cor: " << ptCor_shifted << std::endl;

 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
produceJECUnc::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
produceJECUnc::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
produceJECUnc::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
produceJECUnc::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
produceJECUnc::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
produceJECUnc::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
produceJECUnc::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(produceJECUnc);
