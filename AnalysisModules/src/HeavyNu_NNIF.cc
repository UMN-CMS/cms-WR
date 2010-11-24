#include <string>
#include <iostream>
#include "HeavyNu_NNIF.h"

HeavyNu_NNIF::HeavyNu_NNIF(const edm::ParameterSet& iConfig) :
  isSignal_        (iConfig.getParameter<bool>("isSignal")),
  trainingFileName_(iConfig.getUntrackedParameter<std::string>("trainingFileName","")),
  mNuRnorm_(iConfig.getParameter<double>("mNuRnormalization"))
{
  if (mNuRnorm_ == 0.0)
    throw cms::Exception("mNuRnormalization parameter cannot be 0!!!");

  if (trainingFileName_.length()) {
    trainingMode_ = true;
    netoutputs_.push_back((float)isSignal_);
  }
}

//======================================================================

void
HeavyNu_NNIF::fillvector(const HeavyNuEvent& hne)
{
  netinputs_.clear();

  // This method maps the variables of interest into the neural net input vector
  netinputs_.push_back(hne.mNuR1/mNuRnorm_);
  netinputs_.push_back(hne.mNuR2/mNuRnorm_);
}

//======================================================================

void
HeavyNu_NNIF::print()
{
  using namespace std;

  if (!trainingMode_) return;

  for (size_t i=0; i<netinputs_.size(); i++) {
    if (i) trainingFile_ << " ";
    trainingFile_ << netinputs_[i];
  }
  trainingFile_ << endl;

  for (size_t i=0; i<netoutputs_.size(); i++) {
    if (i) trainingFile_ << " ";
    trainingFile_ << netoutputs_[i];
  }
  trainingFile_ << endl;
}

//======================================================================

void
HeavyNu_NNIF::beginJob()
{
  if (trainingMode_)
    trainingFile_.open(trainingFileName_.c_str());
}

//======================================================================

void
HeavyNu_NNIF::endJob()
{
  if (trainingMode_)
    trainingFile_.close();
}

//======================================================================
