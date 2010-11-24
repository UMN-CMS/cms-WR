/** -- C++ -- **/
#ifndef HeavyNu_NNIF_h_included
#define HeavyNu_NNIF_h_included

#include <vector>
#include <fstream>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HeavyNuEvent.h"

/** The purpose of this class is provide a self-contained
    interface to neural net training/operation for Heavy Nu analysis
*/
class HeavyNu_NNIF {
 public:
  explicit HeavyNu_NNIF(const edm::ParameterSet&);

  void fillvector(const HeavyNuEvent& e);
  void print();
  void beginJob();
  void endJob();

 private:
  bool isSignal_;
  std::string trainingFileName_;
  std::ofstream trainingFile_;
  bool trainingMode_;
  double mNuRnorm_;

  std::vector<float> netinputs_;
  std::vector<float> netoutputs_;
};

#endif
