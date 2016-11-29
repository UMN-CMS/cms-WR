#ifndef jetres_h
#define jetres_h

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"

#include "TRandom3.h"

#include <string>
#include <vector>

void JetResolution(miniTreeEvent *myEvent, TRandom3 rand, bool isData);

#endif
