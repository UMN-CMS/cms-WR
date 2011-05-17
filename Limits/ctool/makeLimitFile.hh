#ifndef MAKE_LIMIT_FILE
#define MAKE_LIMIT_FILE 1

static const int i_SIGNAL=0;
static const int i_TT=1;
static const int i_ZJ=2;
static const int i_OTHER=3;

#include "TFile.h"
#include <vector>

struct PerBinInfo {
  double lumi;
  double signal;
  int data;
  double lowEdge, highEdge;
  int year;
};


void formatLimitFile(const std::vector<PerBinInfo>& pbi, const char* limitFileName);

std::vector<PerBinInfo> makeLimitContent2010(TFile* dataf, TFile* signalf);

void makeLimitFile2010(TFile* dataf, TFile* signalf, const char* limitFileName);

std::vector<PerBinInfo> makeLimitContent2011(double lumi, TFile* dataf, TFile* signalf);

void makeLimitFile2011(double lumi, TFile* dataf, TFile* signalf, const char* limitFileName);

void makeLimitFileTwoYear(double lumi11, TFile* dataf11, TFile* signalf11, TFile* dataf10, TFile* signalf10, const char* limitFileName);


#endif // MAKE_LIMIT_FILE 1
