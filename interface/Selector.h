#include "../interface/miniTreeEvent.h"

class Selector: public miniTreeEvent{
 public:
    Float_t WR_mass;
    Float_t dilepton_mass;
    Bool_t pass_selection;
    miniTreeEvent event;

    void clear();

 Selector(miniTreeEvent myEvent)
   : miniTreeEvent()
  {
    clear();
    event = myEvent;
  }
};

void Selector::clear() {
  WR_mass = dilepton_mass = 0.0;
  pass_selection = false;
  event.clear();
};
