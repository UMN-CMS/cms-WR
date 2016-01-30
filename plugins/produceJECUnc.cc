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

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/Common/interface/ValueMap.h"

typedef double JECUnc_t;
typedef edm::ValueMap<JECUnc_t> JECUnc_Map;

//
// class declaration
//

class produceJECUnc : public edm::stream::EDProducer<> {
   public:
      explicit produceJECUnc(const edm::ParameterSet&);
      ~produceJECUnc();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&) override;

      // ----------member data ---------------------------
      const edm::InputTag src;        ///<input particle objects
		const std::string jetUncOutput;
		const std::string jetType;
};

//
// constructors and destructor
//
produceJECUnc::produceJECUnc(const edm::ParameterSet& iConfig)
   : src(iConfig.getParameter<edm::InputTag>("src")),
   jetUncOutput(iConfig.getParameter<std::string>("jetUncOutput")),
   jetType(iConfig.getParameter<std::string>("jetType"))
{
	consumes<pat::JetCollection>(src);
   produces<JECUnc_Map>(jetUncOutput);
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
	iSetup.get<JetCorrectionsRecord>().get(jetType,JetCorParColl); 
	JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
	JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

   edm::Handle<pat::JetCollection> jets;
   iEvent.getByLabel(src, jets);
	std::vector<JECUnc_t> JECUncertainties;
   for(auto jet : *jets) {
		jecUnc->setJetEta(jet.eta());
		jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
		double unc = jecUnc->getUncertainty(true);
		JECUncertainties.push_back(unc);
		
   }
   std::auto_ptr<JECUnc_Map> JEC_uncertainty_Map(new JECUnc_Map());
   JECUnc_Map::Filler JEC_uncertainty_filler(*JEC_uncertainty_Map);
   JEC_uncertainty_filler.insert(jets,JECUncertainties.begin(),JECUncertainties.end());
   JEC_uncertainty_filler.fill();

   iEvent.put(JEC_uncertainty_Map, jetUncOutput);
 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
produceJECUnc::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<std::string>("jetUncOutput", "JECUncertainty");
  desc.add<std::string>("jetType", "AK4PF");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(produceJECUnc);
