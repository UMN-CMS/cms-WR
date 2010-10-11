#include <iostream>
#include "HeavyNu_NNIF.h"

HeavyNu_NNIF::HeavyNu_NNIF(const edm::ParameterSet& iConfig) :
  trainingMode_(iConfig.getParameter<bool>("trainingMode")),
  isSignal_    (iConfig.getParameter<bool>("isSignal"))
{
  netoutputs.push_back((float)isSignal_);
}

//======================================================================

void
HeavyNu_NNIF::fillvector(const HeavyNuEvent& hne)
{
  netinputs.clear();

  // This method maps the variables of interest into the neural net input vector
  netinputs.push_back(hne.mu1->pt());
  netinputs.push_back(hne.j1->pt());
  netinputs.push_back(deltaPhi(hne.j1->phi(),hne.j2->phi()));
  netinputs.push_back(hne.dRminMu1jet);
  netinputs.push_back(hne.dRminMu2jet);
  netinputs.push_back(hne.mWR);
  netinputs.push_back(hne.mNuR1);
  netinputs.push_back(hne.mNuR2);
  netinputs.push_back(hne.mMuMu);
  netinputs.push_back(hne.ctheta_mu1_jj);
}

//======================================================================

void
HeavyNu_NNIF::print()
{
  using namespace std;

  if (!trainingMode_) return;

  for (size_t i=0; i<netinputs.size(); i++) {
    if (i) trainingFile_ << " ";
    trainingFile_ << netinputs[i];
  }
  trainingFile_ << endl;

  for (size_t i=0; i<netoutputs.size(); i++) {
    if (i) trainingFile_ << " ";
    trainingFile_ << netoutputs[i];
  }
  trainingFile_ << endl;
}

//======================================================================

void
HeavyNu_NNIF::beginJob()
{
  if (trainingMode_)
    trainingFile_.open("heavynu_nntraining.txt");
}

//======================================================================

void
HeavyNu_NNIF::endJob()
{
  if (trainingMode_)
    trainingFile_.close();
}

//======================================================================
