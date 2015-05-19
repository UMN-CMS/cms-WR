///including C++ std lib header files, and header files I wrote
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include "SlimCutVar.h" 

///cmssw include files
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


/**
 * the cutAndProduceOutputCollection producer needs a flexible interface to accommodate
 * a variable number of input object collections (# of output collections = # of input
 * collections), to apply cuts to these input object collections, and to put the objects
 * passing the cuts into the event.  This class is the flexible interface. 
 * 
 */

class CutApplicationInfrastructure {
	typedef std::map<std::string, edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > > tokenMap;
	typedef std::map<std::string, edm::InputTag> inputTagMap;
	typedef std::map<std::string, edm::OwnVector<reco::Candidate> > collectionMap;

	public:
		CutApplicationInfrastructure(){};

		///fill member var maps with string keys, tokens, inputtag objects, and collections
		void initializeMaps(tokenMap tokens, inputTagMap tags, collectionMap colls);

		///use this method to apply a cut to objects in input collections, and add the objects
		///passing the cut into output collections
		void applyCutFillOutputColls(SlimCutVar & aCutObject);

		///place objects in outputCollections into event record
		void putObjectsIntoEvent(edm::Event & anEvent);

	private:
		///use maps with string keys to link tokens, input tags, and object collections
		///the same keys will be used in all maps
		tokenMap inputTokens;
		inputTagMap inputTags;
		collectionMap inputCollections;
		collectionMap outputCollections;	///< the two collectionMap objects have the same # of map entries

};	///end CutApplicationInfrastructure class declaration

