#ifndef toythoer_h
#define toythoer_h

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"

#include "TRandom3.h"

#include <string>
#include <vector>

void ToyThrower(miniTreeEvent* myEvent,  float *rand_smear, float *rand_up_down, int random_seed, std::vector<std::string> list, bool isData);

#endif
