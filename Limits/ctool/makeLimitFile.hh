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
#include <map>
#include <string>

struct PerBinSystematic {
  double signal;
  double bkgd[3];
};

struct PerBinInfo {
  double lumi;
  double signal;
  double bkgd[3];

  std::map<std::string,PerBinSystematic> perBinSyst;
  
  int data;
  int sourceBin;
  double lowEdge, highEdge;
  int year;
  std::string binName;
};

struct LimitPoint {
  double lumi;
  double xsec;
  int year;
  enum { lp_Muon1ECM, lp_Elec1ECM, lp_Muon2ECM, lp_MuonElec } mode;
  double rebin_above_mlljj;
  int mwr, mwr_syst;
  int mnr, mnr_syst;
};

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const LimitPoint& pt, const char* limitFileName);

std::vector<double> extractBins(TFile* f, const std::string& histName);

std::vector<PerBinInfo> makeLimitContent(const LimitPoint& pt, TFile* dataf, const RateDB& dbr, 
					 const SystematicsDB& syst, char bin_prefix='b',bool fullRange=false);

void makeLimitFile(const LimitPoint& pt, TFile* dataf, const RateDB& dbr, const char* limitFileName, const SystematicsDB& syst);

void makeLimitFileInterpolate(const LimitPoint& pt, TFile* dataf, const RateDB& dbf,
			      const LimitPoint& signalp1, 
			      const LimitPoint& signalp2, 
			      const char* limitFileName, const SystematicsDB& syst);

std::string whichSyst();


#endif // MAKE_LIMIT_FILE 1
