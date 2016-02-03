#include "ExoAnalysis/cmsWR/interface/configReader.h"
#include <iostream>

configReader::configReader(std::string filename){
	std::ifstream f_in(filename);
	std::string datasetName; 
	configLine line;

	if(!f_in.good()){
		std::cerr << "[ERROR] file " << filename << " not readable" << std::endl;
		exit(1);
		return;
	}

	while(f_in.peek()!=EOF && f_in.good()){
		if(f_in.peek()==10){ // 10 = \n
			f_in.get(); 
			continue;
		}

		if(f_in.peek() == 35){ // 35 = #
			f_in.ignore(1000,10); // ignore the rest of the line until \n
			continue;
		}

		f_in >> datasetName >> line;
	}
	

}



std::istream& operator >>(std::istream& os, configLine& line){
	os >> line.primaryDatasetPath >> line.crossSection >> line.crossSection_error >> line.primaryDatasetEvents >> line.skimmedDatasetPath >> line.skimmedDatasetEvents;
	return os;
};
