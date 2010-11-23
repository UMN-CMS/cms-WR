#ifndef HeavyNuNetworks_h_included
#define HeavyNuNetworks_h_included 1

#include <vector>

class HeavyNuNetworks {
 public:
  static float evaluate(int mw, int mnu, const std::vector<float>& invalstd);
};

#endif // HeavyNuNetworks_h_included
