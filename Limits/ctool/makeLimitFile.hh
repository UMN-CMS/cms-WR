#ifndef MAKE_LIMIT_FILE
#define MAKE_LIMIT_FILE 1

static const int i_SIGNAL=0;
static const int i_TT=1;
static const int i_ZJ=2;
static const int i_OTHER=3;

class TFile;
class SystematicsDB;

#include <vector>
#include <string>

struct PerBinInfo {
  double lumi;
  double signal;
  double bkgd;
  int data;
  int sourceBin;
  double lowEdge, highEdge;
  int year;
  std::string binName;
};

struct MassPoint {
  int mwr;
  int mnr;
};

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const MassPoint& pt, const char* limitFileName, const SystematicsDB& syst);

std::vector<double> extractBins(TFile* f, const std::string& histName);

std::vector<PerBinInfo> makeLimitContent(double lumi, double xsec, const MassPoint& pt, TFile* dataf, TFile* signalf);

void makeLimitFile(double lumi, double xsec, const MassPoint& pt, TFile* dataf, TFile* signalf, const char* limitFileName, const SystematicsDB& syst);


#endif // MAKE_LIMIT_FILE 1
