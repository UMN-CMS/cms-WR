#ifndef doLimitSetting_hh_included
#define doLimitSetting_hh_included 1

class TFile;

struct CLInfo {
  float obs, exp;
  float exp_p1s, exp_m1s;
  float exp_p2s, exp_m2s;
};
struct BaseInfo {
  float lumi;
  float mwr;
  float mnu;    
  float xsec;
  float signal, background, data;
  float signal_eff;
};
struct LimitPointStruct {
  struct BaseInfo base;
  struct CLInfo cl_sb, cl_b, cls;
};

void doLimitSetting(TFile* dataf, TFile* signal, int ntoys, LimitPointStruct& info, bool doSyst);


#endif // doLimitSetting_hh_included
