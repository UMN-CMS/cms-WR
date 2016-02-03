#include "ExoAnalysis/cmsWR/interface/Selector.h"


Selector::Selector(const miniTreeEvent& myEvent) : 
	miniTreeEvent(myEvent),
	WR_mass(-1),
	dilepton_mass(-1),
	pass_selection(false)
{
    clear();
}

void Selector::Clear() {
  WR_mass = dilepton_mass = 0.0;
  pass_selection = false;
  clear();
};
