/** -- C++ -- **/
#ifndef HeavyNuEvent_h_included
#define HeavyNuEvent_h_included

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

/** The purpose of this class is contain the key items for
    a HeavyNuEvent and provide a simple way to pass this information
    between sections of the code.
*/
class HeavyNuEvent {
public:
  HeavyNuEvent() { reset(); }
  void reset();
  void regularize();

  const pat::Muon *mu1, *mu2, *mu[2];
  const pat::Jet  *j1, *j2, *j[2];

};

#endif
