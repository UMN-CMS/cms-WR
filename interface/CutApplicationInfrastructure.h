///including C++ std lib header files, and header files I wrote
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>

///cmssw include files
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/OwnVector.h"



/**
 * the cutAndProduceOutputCollection producer needs a flexible interface to accommodate
 * a variable number of input object collections (# of output collections = # of input
 * collections), to apply cuts to these input object collections, and to put the objects
 * passing the cuts into the event.  This class is the flexible interface. 
 *
 */

class CutApplicationInfrastructure {
	public:
		CutApplicationInfrastructure(){};

		///place objects passing cuts into event record
		void putObjectsIntoEvent(edm::Event & anEvent);

	private:
		///use maps with string keys to link tokens, input tags, and object collections
		typedef std::map<std::string, edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > > tokenMap;
		tokenMap inputTokens;
		typedef std::map<std::string, edm::InputTag> inputTagMap;
		inputTagMap inputTags;
		typedef std::map<std::string, edm::OwnVector<reco::Candidate> > inputTagMap;
	
		std::vector<edm::OwnVector<reco::Candidate> > l;

};	///end CutApplicationInfrastructure class declaration

