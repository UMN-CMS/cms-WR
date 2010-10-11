#ifndef doLimitSetting_hh_included
#define doLimitSetting_hh_included 1

class TFile;

struct LimitPointStruct {
  float lumi;
  float mwr;
  float mnu;    
  float xsec;
  float cl_sb_obs, cl_sb_exp;
  float cl_b_obs, cl_b_exp;
  float cl_sb_exp_p1s, cl_sb_exp_m1s;
  float cl_sb_exp_p2s, cl_sb_exp_m2s;
  float cl_b_exp_p1s, cl_b_exp_m1s;
  float cl_b_exp_p2s, cl_b_exp_m2s;

};

void doLimitSetting(TFile* dataf, TFile* signal, int ntoys, LimitPointStruct& info);


#endif // doLimitSetting_hh_included
