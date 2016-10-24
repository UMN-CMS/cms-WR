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
	std::vector<std::string> getDatasetNames();

	void setupDyMllScaleFactor(std::string inputFile)
	{
		std::string ch = "", dataset = "";
		double scalefactor = 1.0, scalefactorUncertainty = 0.0;
		std::ifstream readScaleFactors;
		readScaleFactors.open(inputFile);
		if(readScaleFactors.is_open() == false) {
			std::cerr << "[ERROR] File " << inputFile << " not opened. Check if the file exists" << std::endl;
		}

		while(readScaleFactors.peek() != EOF && readScaleFactors.good()) {
			readScaleFactors >> ch >> dataset >> scalefactor >> scalefactorUncertainty;
			DyMllScaleFactorMap[ch + dataset] = scalefactor;
		}
	}
	/// return the dataOverMC dilepton mass correction factor for one DY sample
	double getDyMllScaleFactor(std::string chnl, std::string mode_str, std::string inputFile = "")
	{
		if(DyMllScaleFactorMap.size() == 0) {
			setupDyMllScaleFactor(inputFile);
		}

		if (DyMllScaleFactorMap.find(chnl + mode_str) != DyMllScaleFactorMap.end()) {
			return DyMllScaleFactorMap[chnl + mode_str];
		}
		return 1.0;	///< for DYPOWINCL
	}///end getDyMllScaleFactor()

	inline configLine::numEvents_t getPrimaryEvents(std::string datasetName) const
	{
		if(checkCategory(datasetName) == true) return getConfigLine(datasetName).primaryDatasetEvents;
		else return -1;
	}

	inline std::string rochesterFile(void) const
	{
		return configFile["rochester"];
	};

private:
	std::map<std::string, configLine> configMap;
	std::map<std::string, std::string> configFile;
	std::map<std::string, double> DyMllScaleFactorMap;

	void datasetsFileReader(std::string filename);

	std::string unblindTag(void);

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
