#ifndef HeavyNuNetworks_h_included
#define HeavyNuNetworks_h_included 1

#include <vector>
#include <utility>

class HeavyNuNetworks {
 public:
  static float evaluate(int mw, int mnu, const std::vector<float>& invalstd);
  static std::vector<std::pair<int,int> > getMassPoints();
};

#endif // HeavyNuNetworks_h_included
