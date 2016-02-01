#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <string>

class produceStringTag : public edm::EDProducer
{
public:
	produceStringTag(const edm::ParameterSet&);

	virtual void produce(edm::Event&, const edm::EventSetup&);
private:
	std::string _outputCollName;	///<user defined
	std::string _stringTagName;		///<user defined

};

produceStringTag::produceStringTag(const edm::ParameterSet& cfg):
	_outputCollName(cfg.getParameter<std::string>("outputCollectionName")),
	_stringTagName(cfg.getParameter<std::string>("stringStoredInOutputCollection"))
{
	produces <std::string>(_outputCollName);
}

void produceStringTag::produce(edm::Event& event, const edm::EventSetup&)
{
	std::auto_ptr<std::string> strgTg(new std::string(_stringTagName));
	event.put(strgTg, _outputCollName);
}

DEFINE_FWK_MODULE(produceStringTag);
