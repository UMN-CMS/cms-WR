#ifndef toythoer_h
#define toythoer_h

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"

#include "TRandom3.h"

#include <string>
#include <vector>

void ToyThrower(miniTreeEvent myEvent, TRandom3& rand, int random_seed, std::vector<std::string> list, int isData);

#endif
