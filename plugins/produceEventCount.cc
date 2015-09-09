#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <string>

class produceEventCount : public edm::EDProducer
{
public:
	produceEventCount(const edm::ParameterSet&);

	virtual void produce(edm::Event&, const edm::EventSetup&);
private:
	std::string outputCollName;	///<user defined

};

produceEventCount::produceEventCount(const edm::ParameterSet& cfg):
	outputCollName(cfg.getParameter<std::string>("outputCollectionName"))
{
	produces <Float_t>(outputCollName);
}

void produceEventCount::produce(edm::Event& event, const edm::EventSetup&)
{
	std::auto_ptr<Float_t> evtCount(new Float_t(1.0));
	event.put(evtCount, outputCollName);
}

DEFINE_FWK_MODULE(produceEventCount);
