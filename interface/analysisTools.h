#include <map>
#include <utility>
#include <fstream>
#include <string>
#include "Selector.h"

#include "RooAbsPdf.h"
#include "RooRealVar.h"
using namespace RooFit;

typedef std::map< std::pair<Selector::tag_t, int>, std::pair<int, int> > mass_cut_map_t;

mass_cut_map_t getMassCutMap()
{
	mass_cut_map_t mass_cut;
	std::ifstream ifs;
	ifs.open("configs/mass_cuts.txt", std::ifstream::in);
	std::string line;
	std::string mass, ch, low, high;
	while (std::getline(ifs, line)) {
		if (line[0] == '#' || line.size() == 0) continue;
		std::stringstream ss(line);
		std::getline(ss, ch, ' ');
		std::getline(ss, mass, ' ');
		std::getline(ss, low, ' ');
		std::getline(ss, high, ' ');
		mass_cut[std::make_pair(Selector::getTag(ch), std::stoi(mass))] = std::make_pair(std::stoi(low), std::stoi(high));
	}

	return mass_cut;
}

std::vector<int> getMassVec()
{
	std::vector<int> mass;
	std::ifstream ifs;
	ifs.open("configs/mass_cuts.txt", std::ifstream::in);
	std::string line;
	std::string word;
	while (std::getline(ifs, line)) {
		if (line[0] == '#' || line.size() == 0) continue;
		std::stringstream ss(line);
		// check channel, only use ee.
		std::getline(ss, word, ' ');
		if (word == "EE") {
			std::getline(ss, word, ' ');
			mass.push_back(std::stoi(word));
		}
	}

	return mass;
}

double NormalizedIntegral(RooAbsPdf * function, RooRealVar& integrationVar, double lowerLimit, double upperLimit)
{
	integrationVar.setRange("integralRange", lowerLimit, upperLimit) ;
	RooAbsReal* integral = (*function).createIntegral(integrationVar, NormSet(integrationVar), Range("integralRange")) ;
	double normalizedIntegralValue = integral->getVal();
	return normalizedIntegralValue;
}
