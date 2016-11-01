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
  std::string m;
  while (std::getline(ifs, line)) {
    if (line[0] == '#' || line.size() == 0) continue;
    std::stringstream ss(line);
    std::getline(ss, m, ' ');
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

std::map<Double_t, float> PU_map = {
  {0,0.0},
  {1,0.00037479523852},
  {2,0.00723350728499},
  {3,0.0149819481392},
  {4,0.026772485431},
  {5,0.0469226006499},
  {6,0.036530283623},
  {7,0.0378144988499},
  {8,0.119371164953},
  {9,0.285872999831},
  {10,0.491813793768},
  {11,0.689632013241},
  {12,0.959594440995},
  {13,1.170233541},
  {14,1.27517863083},
  {15,1.4401376652},
  {16,1.44213795678},
  {17,1.35998258979},
  {18,1.41063135333},
  {19,1.34139222043},
  {20,1.4238085196},
  {21,1.22503213772},
  {22,1.10774533938},
  {23,1.04844682795},
  {24,1.01400723977},
  {25,0.905968755792},
  {26,0.897884712575},
  {27,0.693673170731},
  {28,0.619167741211},
  {29,0.434955151335},
  {30,0.353393947441},
  {31,0.246646716677},
  {32,0.16923169631},
  {33,0.14852029452},
  {34,0.11676760677},
  {35,0.12839888127},
  {36,0.166653218986},
  {37,0.188896371236},
  {38,0.238628610164},
  {39,0.0},
  {40,0.0},
  {41,0.0},
  {42,0.0},
  {43,0.0},
  {44,0.0},
  {45,0.0},
  {46,0.0},
  {47,0.0},
  {48,0.0},
  {49,0.0},
  {50,0.0}
};
float PUreweight(Double_t nPU)
{
  return PU_map[nPU];
  //return nPU;
}
