#include "ExoAnalysis/cmsWR/interface/configReader.h"
#include <iostream>

std::string configReader::unblindTag(void)
{
	if(configFile["unblind"] == "true") return "_unblind";
	return "";
}


configReader::configReader(std::string filename)
{
	std::ifstream f_in(filename);
	char  var[100], varValue[100];

	if(!f_in.good()) {
		std::cerr << "[ERROR] file " << filename << " not readable" << std::endl;
		exit(1);
		return;
	}

	while(f_in.peek() != EOF && f_in.good()) {

		if(f_in.peek() == 10) { // 10 = \n
			f_in.get();
			continue;
		}

		if(f_in.peek() == 35) { // 35 = #
			f_in.ignore(1000, 10); // ignore the rest of the line until \n
			continue;
		}

		f_in.get(var, 100, '=');
		f_in.get(); //reading the delim char
		f_in.get(varValue, 100);
		//sscanf(line, "%s=%s", var, varValue);
		configFile[var] = varValue;
		if(std::string(var) == "datasetFile") {
			std::cout << var << "\t" << varValue << std::endl;
			datasetsFileReader(varValue);
		}
	}

}

// datasetFile=configs/datasets.dat
// jsonFile=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver.txt
// jsonName=246908-260627-Prompt_25ns-golden_silver-v1
// productionTAG=_SHv3

// KEY: TDirectoryFile   miniTree_flavoursideband;1      miniTree_flavoursideband
// KEY: TDirectoryFile   miniTree_lowdileptonsideband;1  miniTree_lowdileptonsideband
// KEY: TDirectoryFile   miniTree_signal_ee;1    miniTree_signal_ee
// KEY: TDirectoryFile   miniTree_signal_mumu;1  miniTree_signal_mumu
// KEY: TDirectoryFile   miniTree_dytagandprobe;1        miniTree_dytagandprobe
// KEY: TDirectoryFile   zToEEAnalyzer;1 zToEEAnalyzer
// KEY: TDirectoryFile   zToMuMuAnalyzer;1       zToMuMuAnalyzer

TChain *configReader::getMiniTreeChain(std::string datasetName, std::string tag)
{
	TChain *chain = new TChain((tag + "/t").c_str(), "");
	//////////standard minitrees stored in shervin's user area
	//chain->Add(("root://eoscms.cern.ch//store/user/shervin/ntuples/" + datasetName + configFile["productionTAG"] + unblindTag() + "/unmerged-allRange.root").c_str());
	
	///////////other bkgnd minitrees (WW, single top, t+W, QCD), DYMadIncl gen HT < 100 with pdf (uncert) weights, and DYMadHT and DYAMC with pdf (uncert) weights all stored in WR group space on EOS
	///////////also WR minitrees with pdf (uncertainty) weights
	chain->Add(("root://eoscms.cern.ch//store/group/phys_exotica/leptonsPlusJets/WR/tuples/" + datasetName + configFile["productionTAG"] + unblindTag() + "/unmerged-allRange.root").c_str());
	chain->GetEntries();
	return chain;
}

TChain *configReader::getMiniTreeChain(std::vector<std::string> datasetNames, std::string tag)
{

	TChain *chain = new TChain((tag + "/t").c_str(), "");
	for(auto s : datasetNames) {
		std::cout << "[INFO] Reading chain for: " << s << "\t" << tag << std::endl;
		TChain *c = getMiniTreeChain(s, tag);
		chain->Add(c);
		delete c;
	}

	return chain;
}

std::vector<std::string> configReader::getDatasetNames()
{
	std::vector<std::string> names;
	for(auto d : configMap) {
		names.push_back(d.first);
	}
	return names;
}

void configReader::datasetsFileReader(std::string filename)
{
	std::ifstream f_in(filename);
	std::string datasetName;
	configLine line;

	if(!f_in.good()) {
		std::cerr << "[ERROR] file " << filename << " not readable" << std::endl;
		exit(1);
		return;
	}

	while(f_in.peek() != EOF && f_in.good()) {

		if(f_in.peek() == 10) { // 10 = \n
			f_in.get();
			continue;
		}

		if(f_in.peek() == 35) { // 35 = #
			f_in.ignore(1000, 10); // ignore the rest of the line until \n
			continue;
		}

		f_in >> datasetName >> line;
		configMap[datasetName] = line;
#ifdef DEBUG
		std::cout << "##" << datasetName << "##" << "\t" << configMap.size() << std::endl;
		std::cout << configMap.begin()->first << "\t" << configMap.count(datasetName) << std::endl;
		//assert(checkCategory(datasetName)==false);
#endif
	}


}




std::ostream& operator<<(std::ostream& os, configReader& r)
{
	os << r.configMap.size() << std::endl;
	for(auto line : r.configMap) {
		os << line.first << "\t" << line.second << std::endl;
	}
	return os;
}

std::istream& operator >>(std::istream& os, configLine& line)
{
	char tmp[100];
	os >> line.primaryDatasetPath;
	os >> tmp;
	if(tmp[0] != '-') sscanf(tmp, "%lf", &line.crossSection);
	else line.crossSection = 0.;

	os >> tmp;
	if(tmp[0] != '-') sscanf(tmp, "%lf", &line.crossSection_error);
	else line.crossSection_error = 0.;

	os >> line.primaryDatasetEvents >> line.skimmedDatasetPath;

	os >> tmp;
	if(tmp[0] != '-') sscanf(tmp, "%lld", &line.skimmedDatasetEvents);
	else line.skimmedDatasetEvents = 0.;

	return os;
}

std::ostream& operator <<(std::ostream& os, configLine& line)
{
	os << line.primaryDatasetPath << "\t" <<  line.crossSection << "\t" <<  line.crossSection_error << "\t" <<  line.primaryDatasetEvents << "\t" <<  line.skimmedDatasetPath << "\t" <<  line.skimmedDatasetEvents;
	return os;
}
