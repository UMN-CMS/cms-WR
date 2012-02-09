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

struct LimitPoint {
  double lumi;
  double xsec;
  int mwr, mwr_syst;
  int mnr, mnr_syst;
};

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const LimitPoint& pt, const char* limitFileName, const SystematicsDB& syst);

std::vector<double> extractBins(TFile* f, const std::string& histName);

std::vector<PerBinInfo> makeLimitContent(const LimitPoint& pt, TFile* dataf, TFile* signalf, bool fullRange=false);

void makeLimitFile(const LimitPoint& pt, TFile* dataf, TFile* signalf, const char* limitFileName, const SystematicsDB& syst);

void makeLimitFileInterpolate(const LimitPoint& pt, TFile* dataf, 
			      TFile* signalf1, const LimitPoint& signalp1, 
			      TFile* signalf2, const LimitPoint& signalp2, 
			      const char* limitFileName, const SystematicsDB& syst);


#endif // MAKE_LIMIT_FILE 1
