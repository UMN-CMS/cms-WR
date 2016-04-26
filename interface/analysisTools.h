#include <map>
#include <utility>
#include <fstream>
#include <string>

#include "RooAbsPdf.h"
#include "RooRealVar.h"
using namespace RooFit;

std::map<int, std::pair<int,int> >  getMassCutMap()
{
	std::map<int, std::pair<int, int> > mass_cut;
	std::ifstream ifs;
	ifs.open("configs/mass_cuts.txt", std::ifstream::in);
	std::string line;
	std::string mass, low, high;
	while (std::getline(ifs, line)) {
		if (line[0] == '#' || line.size() == 0) continue;
		std::stringstream ss(line);
		std::getline(ss, mass, ' ');
		std::getline(ss, low, ' ');
		std::getline(ss, high, ' ');
		mass_cut[std::stoi(mass)] = std::make_pair(std::stoi(low), std::stoi(high));
	}

	return mass_cut;
}

std::vector<int> getMassVec()
{
	std::vector<int> mass;
	std::ifstream ifs;
	ifs.open("configs/mass_cuts.txt", std::ifstream::in);
	std::string line;
	std::string m;
	while (std::getline(ifs, line)) {
		if (line[0] == '#' || line.size() == 0) continue;
		std::stringstream ss(line);
		std::getline(ss, m, ' ');
		mass.push_back(std::stoi(m));
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
