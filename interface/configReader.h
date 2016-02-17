#ifndef configreader_h
#define configreader_h

#include <map>
#include <fstream>
#include <iostream>
#include <cassert>

#include <TChain.h>

class configLine
{
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
	friend std::ostream& operator <<(std::ostream&os, configLine& line);
};



class configReader
{
public:
	configReader(std::string filename);

	friend std::ostream& operator <<(std::ostream& os, configReader& r);

	inline configLine::norm_t getNorm1fb(std::string datasetName) const
	{
		if(checkCategory(datasetName) == true) {
			const configLine& c = getConfigLine(datasetName);
			return c.crossSection / c.primaryDatasetEvents;
		}
		return -1;
	}

	TChain *getMiniTreeChain(std::string datasetName, std::string tag);
	TChain *getMiniTreeChain(std::vector<std::string> datasetNames, std::string tag);


private:
	std::map<std::string, configLine> configMap;
	std::map<std::string, std::string> configFile;

	void datasetsFileReader(std::string filename);

	inline bool checkCategory(std::string datasetName) const
	{
		if(configMap.count(datasetName) == 0) {
			std::cerr << "[WARNING] datasetName " << datasetName << " not found" << std::endl;
			exit(1);
			return false;
		}
		return true;
	}

	inline const configLine& getConfigLine(std::string datasetName) const
	{
		assert(checkCategory(datasetName) == true);
		//return configMap[datasetName];
		return configMap.find(datasetName)->second; // configMap[datasetName] does not return const objects, so it won't compile
	}

	inline configLine::datasetName_t getDatasetPath(std::string datasetName) const
	{
		if(checkCategory(datasetName) == true) return getConfigLine(datasetName).primaryDatasetPath;
		else return configLine::datasetName_t();
	};

	inline configLine::datasetName_t getSkimmedDatasetPath(std::string datasetName) const
	{
		if(checkCategory(datasetName) == true) return getConfigLine(datasetName).skimmedDatasetPath;
		else return configLine::datasetName_t();
	}

};



#endif
