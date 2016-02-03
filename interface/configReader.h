#ifndef configreader_h
#define configreader_h

#include <map>
#include <fstream>
#include <iostream>
#include <cassert>

class configLine{
public:
	typedef std::string datasetName_t;
	typedef long long int numEvents_t;
	typedef double crossSection_t;
	typedef double norm_t;

	datasetName_t primaryDatasetPath;
	datasetName_t skimmedDatasetPath;
	crossSection_t crossSection, crossSection_error;
	numEvents_t primaryDatasetEvents;
	numEvents_t skimmedDatasetEvents;

	friend std::istream& operator >>(std::istream& os, configLine& line);
};


class configReader{

	configReader(std::string filename);
	
	std::map<std::string, configLine> configMap;

	inline bool checkCategory(std::string datasetName) const{
		if(configMap.count(datasetName)==0){
			std::cerr << "[WARNING] datasetName not found" << std::endl;
			return false;
		}
		return true;
	}
		
	inline const configLine& getConfigLine(std::string datasetName) const{
		assert(checkCategory(datasetName)==true);
		//return configMap[datasetName];
		return configMap.find(datasetName)->second; // configMap[datasetName] does not return const objects, so it won't compile
	}
		
	inline configLine::datasetName_t getDatasetPath(std::string datasetName) const{ 
		if(checkCategory(datasetName)==true) return getConfigLine(datasetName).primaryDatasetPath;
		else return configLine::datasetName_t();
	};

	inline configLine::norm_t getNorm1fb(std::string datasetName) const {
		if(checkCategory(datasetName)==true){
			const configLine& c = getConfigLine(datasetName);
			return c.crossSection/c.primaryDatasetEvents;
		} 
		return -1;
	}
			
};



#endif
