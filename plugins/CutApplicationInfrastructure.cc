#include "../interface/CutApplicationInfrastructure.h"

using namespace std;

void CutApplicationInfrastructure::initializeMaps(tokenMap tokens, inputTagMap tags, collectionMap colls){

}///end initializeMaps()

void CutApplicationInfrastructure::applyCutFillOutputColls(SlimCutVar & aCutObject){

}///end applyCutFillOutputColls()

void CutApplicationInfrastructure::putObjectsIntoEvent(edm::Event & anEvent){
	/*
	for(collectionMap::iterator collMapIt = outputCollections->begin(); collMapIt!=outputCollections->end(); collMapIt++){
		anEvent.put(collMapIt->second,collMapIt->first);
	}///end loop over entries in outputCollections map
	*/
}///end putObjectsIntoEvent()

