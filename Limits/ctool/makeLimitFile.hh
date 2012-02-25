#ifndef MAKE_LIMIT_FILE
#define MAKE_LIMIT_FILE 1

static const int i_SIGNAL=0;
static const int i_TT=1;
static const int i_ZJ=2;
static const int i_OTHER=3;

class TFile;
class SystematicsDB;
class RateDB;

#include <vector>
#include <string>

struct PerBinInfo {
  double lumi;
  double signal;
  double bkgd[3];
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

std::vector<PerBinInfo> makeLimitContent(const LimitPoint& pt, TFile* dataf, const RateDB& dbr, bool fullRange=false);

void makeLimitFile(const LimitPoint& pt, TFile* dataf, const RateDB& dbr, const char* limitFileName, const SystematicsDB& syst);

void makeLimitFileInterpolate(const LimitPoint& pt, TFile* dataf, const RateDB& dbf,
			      const LimitPoint& signalp1, 
			      const LimitPoint& signalp2, 
			      const char* limitFileName, const SystematicsDB& syst);

std::string whichSyst();


#endif // MAKE_LIMIT_FILE 1
