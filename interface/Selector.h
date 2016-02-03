#ifndef selector_h
#define selector_h

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"

class Selector: public miniTreeEvent{
public:
    Float_t WR_mass; // this is of Float_t because want to save it into a tree
    Float_t dilepton_mass;
    Bool_t pass_selection;
	
	Selector(const miniTreeEvent& myEvent);

private:
	
    void Clear();

};



#endif
