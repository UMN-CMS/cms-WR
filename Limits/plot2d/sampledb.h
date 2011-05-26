#ifndef sampledb_h
#define sampledb_h

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

// originally, signal crossections in the table is LO cross sections
// we use the same k-factor for 7, 10 and 14 TeV samples
//const float k = 1.13; // k-factor for the signal (old version)
// const float k = 1.3; // k-factor for the signal for 1 < MWR < 1.3
//const float kZ = 1.37; // old
const float kZ = 1.25; // k-factor for Alpgen Z+jets Fall 2010 production (tot = 3048 pb)
const float kW = 1.236;//// k-factor for Alpgen W+jets Fall 2010 production (tot = 31314 pb)

// record of sample properties
struct Sample {
  int type;            // 0 - background (MC), 1 - signal (MC), 2 - data, 3 - unknown
  int ECM;             // beam energy [TeV]
  int MW;              // W_R mass [GeV]
  int MNu;             // Nu_R mass [GeV]
  int channel;         // 0 - all, 1 - electron, 2 - muon
  float lumi;          // Integrated luminosity [pb^-1], the field valid only for data samples
  int stat;            // Number of events
  float CS;            // Cross-section [pb]
  std::string name;    // name
  int pid;             // Process ID
  std::string fname;   // file name
  std::string AlCa;    // frontier conditions (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
  std::string HLTmenu; // HLT Menu            (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideGlobalHLT)
  std::string CMSSW;   // version
  std::string DBS;     // sample path in DBS  (https://cmsweb.cern.ch/dbs_discovery/)
  
  // function return weight of event in sample
  //   L - integrated luminosity [pb^-1]
  float weight(const float L) const
  {
    // return unit weight for data or unknown samples
    if (type == 2 || type == 3) return 1;
    
    // return -1 if event weight stored in tree
    if (isIntegratedWeightPid()) return -1;
    
    return L*CS/stat;
  }
  
  // return true if real values of Process ID and Event Weight
  // should be read from the tree for this sample
  // this is necessary for several CSA07 MC background samples
  // (because events from several processes was mixed in one dataset)
  bool isIntegratedWeightPid() const
  {
    return ((type == 0) && (CS < 0) && (pid < 0));
  }
  
  // the function check, if "par" == "val" for this sample
  bool check(const std::string par, const int val) const
  {
         if (par == "type")    return (val == type);
    else if (par == "ECM")     return (val == ECM);
    else if (par == "MW")      return (val == MW);
    else if (par == "MNu")     return (val == MNu);
    else if (par == "channel") return (val == channel);
    else if (par == "stat")    return (val == stat);
    else if (par == "pid")     return (val == pid);
    else return false;
  }
  
  bool check(const std::string par, const std::string val) const
  {
         if (par == "name")    return (val == name);
    else if (par == "fname") {
      // get file name from full path
      std::string fn = val;
      const std::string::size_type pos = fn.find_last_of('/');
      if (pos != std::string::npos)
        fn.erase(0, pos + 1);
      
      return (fn == fname);
    }
    else if (par == "AlCa")    return (val == AlCa);
    else if (par == "HLTmenu") return (val == HLTmenu);
    else if (par == "CMSSW")   return (val == CMSSW);
    else if (par == "DBS")     return (val == DBS);
    else return false;
  }
  
  // equality operator, check only by file name (have to be uniq)
  bool operator == (const Sample& sa) const
  {
    return (fname == sa.fname);
  }
  
  // overload of << operator to allow easy printing of Sample
  friend std::ostream& operator<<(std::ostream& os, const Sample& sa)
  {
    const char* type[4] = {"background", "signal", "data", "unknown"};
    const char* channel[3] = {"all", "electron", "muon"};
    
    os << "[Sample]: " << sa.name << " (" << type[sa.type] << ")\n";
    
    if (sa.type == 0) // background
      os << "  ECM = "   << sa.ECM << " TeV\n"
         << "  stat = "  << sa.stat << " events\n"
         << "  CS = "    << sa.CS << " pb\n"
         << "  fname = " << sa.fname << "\n"
         << "  AlCa = "  << sa.AlCa << "\n"
         << "  HLTmenu = "  << sa.HLTmenu << "\n"
         << "  CMSSW = " << sa.CMSSW << "\n"
         << "  DBS = "   << sa.DBS << "\n";
    
    else if (sa.type == 1) // signal
      os << "  ECM = "   << sa.ECM << " TeV\n"
         << "  channel = "  << channel[sa.channel] << "\n"
         << "  MW / MNu = "  << sa.MW << " / " << sa.MNu << " GeV\n"
         << "  stat = "  << sa.stat << " events\n"
         << "  CS = "    << sa.CS << " pb\n"
         << "  fname = " << sa.fname << "\n"
         << "  AlCa = "  << sa.AlCa << "\n"
         << "  HLTmenu = "  << sa.HLTmenu << "\n"
         << "  CMSSW = " << sa.CMSSW << "\n";
    
    else if (sa.type == 2) // data
      os << "  ECM = "   << sa.ECM << " TeV\n"
         << "  integrated luminosity = "  << sa.lumi << " pb^-1\n"
         << "  stat = "  << sa.stat << " events\n"
         << "  fname = " << sa.fname << "\n"
         << "  AlCa = "  << sa.AlCa << "\n"
         << "  HLTmenu = "  << sa.HLTmenu << "\n"
         << "  CMSSW = " << sa.CMSSW << "\n"
         << "  DBS = "   << sa.DBS << "\n";
    else // unknown
      os << "  fname = " << sa.fname << "\n";
    
    return os;
  }
};

// array of all samples
const Sample AllSamples[] = {
  // Background
  // --- 7 TeV samples ---
// type  ECM MW  MNu  channel  lumi   stat            CS      name        pid        fname                                    AlCa       HLTmenu  CMSSW       DBS
 
  // ttbar is simulated with tauola and pythia, LO cs = 94.3
  // ttbar cs is NLO taken from: http://alcaraz.web.cern.ch/alcaraz/CROSS_SECTIONS.txt
  // Zee and Zmumu cs are NLO Mll>20 taken from: http://alcaraz.web.cern.ch/alcaraz/CROSS_SECTIONS.txt, LO cs is 1300
  // Wenu and Wmunu cs are NLO taken from: http://alcaraz.web.cern.ch/alcaraz/CROSS_SECTIONS.txt, LO cs is ~6000
  {0,    7,  0,  0,   0,       0,    626610,       1.62E+02, "ttbar",  1040000,  "7TeV-TTbar.root",                             "MC_31X_V3", "8E29", "3_1_2", "/TTbar/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,    626610, 100 * 1.62E+02, "ttbar",  1040000,  "7TeV-TTbar-SD_Ele15.root",                    "MC_31X_V3", "8E29", "", "/TTbar/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,    626610, 100 * 1.62E+02, "ttbar",  1040000,  "7TeV-TTbar-SD_Mu9.root",                      "MC_31X_V3", "8E29", "", "/TTbar/Summer09-MC_31X_V3_7TeV_SD_Mu9-v1/GEN-SIM-RECO"},
  
  {0,    7,  0,  0,   0,       0,   2078361,       1.03E+04, "Wenu",   1054000,  "7TeV-Wenu.root",                              "MC_31X_V3", "8E29", "3_1_2", "/Wenu/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,   2078361,  20 * 1.03E+04, "Wenu",   1054000,  "7TeV-Wenu-SD_Ele15.root",                     "MC_31X_V3", "8E29", "", "/Wenu/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,   2073721,       1.03E+04, "Wmunu",  1053000,  "7TeV-Wmunu.root",                             "MC_31X_V3", "8E29", "", "/Wmunu/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,   2073721,  20 * 1.03E+04, "Wmunu",  1053000,  "7TeV-Wmunu-SD_Mu9.root",                      "MC_31X_V3", "8E29", "", "/Wmunu/Summer09-MC_31X_V3_7TeV_SD_Mu9-v1/GEN-SIM-RECO"},
  
  {0,    7,  0,  0,   0,       0,   2538855,       1.67E+03, "Zee",    1071010,  "7TeV-Zee.root",                               "MC_31X_V3", "8E29", "", "/Zee/Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,   2538855, 100 * 1.67E+03, "Zee",    1071010,  "7TeV-Zee-SD_Ele15.root",                      "MC_31X_V3", "8E29", "", "/Zee/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO"},
  // part (30000 events) of 7TeV-Zmumu dataset are not read due to corrupted .root file
  // for details, see 7TeV-Zmumu/crab_0_091213_003103/res/CMSSW_58.stderr
  {0,    7,  0,  0,   0,       0,   2310156,       1.67E+03, "Zmumu",  1071020,  "7TeV-Zmumu.root",                             "MC_31X_V3", "8E29", "", "/Zmumu/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,   2310156, 100 * 1.67E+03, "Zmumu",  1071020,  "7TeV-Zmumu-SD_Mu9.root",                      "MC_31X_V3", "8E29", "", "/Zmumu/Summer09-MC_31X_V3_7TeV_SD_Mu9-v1/GEN-SIM-RECO"},
  
  // WW, WZ, ZZ cs are NLO taken from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CrossSections_3XSeries, LO cs are 28, 10.5, 4.3
  {0,    7,  0,  0,   0,       0,    120280,            43.,    "WW",  1072000,  "7TeV-WW.root",                                "MC_31X_V3", "8E29", "", "/WW/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,    114070,            18.,    "WZ",  1076000,  "7TeV-WZ.root",                                "MC_31X_V3", "8E29", "", "/WZ/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,    145368,            5.9,    "ZZ",  1073000,  "7TeV-ZZ.root",                                "MC_31X_V3", "8E29", "", "/ZZ/Summer09-MC_31X_V3_7TeV-v1/GEN-SIM-RECO"},

// === CSA10 ======================================================================================================================================================================

  // In fact trigger table name in the following CSA10 samples is "GRun".
  // This trigger table consists of the 8E29 physics triggers plus the online
  // commissioning triggers migrated from CMSSW_22X_HLT.
// type  ECM MW  MNu  channel  lumi   stat           CS      name        pid        fname                                         AlCa        HLTmenu  CMSSW       DBS
  // NLO ttbar inclusive cross section calculated with MCFM, initial (LO) was 95
  // NNLO Z+jets cross section for Mll > 50, l=e,mu,tau calculated with FEWZ, initial (LO) was 2224
  // NNLO W+jets cross section for l=e,mu,tau  calculated with FEWZ, initial (LO) was 25090
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MadGraphStandardModel
  // https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  {0,    7,  0,  0,   0,       0,   1047008,     3048,        "Z+jets", 4052001,   "ZJets_7TeV-madgraph-tauola.root",            "START36_V9",  "8E29", "3_6_2", "/ZJets_7TeV-madgraph-tauola/Summer10-START36_V9_S09-v2/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,  10218854,    31314,        "W+jets", 4051001,   "WJets_7TeV-madgraph-tauola-36.root",            "START36_V9",  "8E29", "3_6_2", "/WJets_7TeV-madgraph-tauola/Summer10-START36_V9_S09-v1/GEN-SIM-RECO"},
  {0,    7,  0,  0,   0,       0,   1463572,    157.5,        "ttbar",  4040001,   "TTbarJets_7TeV-madgraph-tauola.root",        "START36_V9",  "8E29", "3_6_2", "/TTbarJets_Tauola-madgraph/Summer10-START36_V9_S09-v1/GEN-SIM-RECO"},
  
  
  // https://twiki.cern.ch/twiki/bin/view/CMS/ProductionFall2010
  // NLO ttbar inclusive cross section calculated with MCFM, initial (LO) was 94
  // Total Z+jets CS (LO) = 2436. k factor kZ (calculated from the above LO and NNLO cs) 1.37
  // Total WZ CS: NLO, see above
  
  {0,     7,  0,  0,   0,       0,1435909 ,  kZ*1.929e+03,  "Z+jets",  8000022,   "Z0Jets_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z0Jets_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},   
  
  {0,     7,  0,  0,   0,       0, 317567 ,  kZ*3.808e+02,   "Z+jets",  8000023,   "Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 265406 ,  kZ*8.721e+00,   "Z+jets",  8000028,   "Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 110308 ,  kZ*7.386e-02,   "Z+jets",  8000033,   "Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  32549 ,  kZ*1.374e-04,   "Z+jets",  8000038,   "Z1Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z1Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  
  {0,     7,  0,  0,   0,       0, 118361 ,  kZ*1.035e+02,   "Z+jets",  8000024,   "Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 129106 ,  kZ*8.534e+00,   "Z+jets",  8000029,   "Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 109393 ,  kZ*1.151e-01,   "Z+jets",  8000034,   "Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  11040 ,  kZ*3.023e-04,   "Z+jets",  8000039,   "Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
 
  {0,     7,  0,  0,   0,       0,  55037 ,  kZ*2.289e+01,   "Z+jets",  8000025,   "Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  55058 ,  kZ*3.951e+00,   "Z+jets",  8000030,   "Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  54262 ,  kZ*8.344e-02,   "Z+jets",  8000035,   "Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  11328 ,  kZ*2.480e-04,   "Z+jets",  8000040,   "Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
 
  {0,     7,  0,  0,   0,       0,  44432 ,  kZ*4.619e+00,   "Z+jets",  8000026,   "Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  44276 ,  kZ*1.298e+00,   "Z+jets",  8000031,   "Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  10874 ,  kZ*3.935e-02,   "Z+jets",  8000036,   "Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  10622 ,  kZ*1.394e-04,   "Z+jets",  8000041,   "Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
 
  {0,     7,  0,  0,   0,       0,  10934 ,  kZ*1.135e+00,   "Z+jets",  8000027,   "Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  10563 ,  kZ*4.758e-01,   "Z+jets",  8000032,   "Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  11146 ,  kZ*1.946e-02,   "Z+jets",  8000037,   "Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  11172 ,  kZ*7.195e-05,   "Z+jets",  8000042,   "Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
 
  // {0,     7,  0,  0,   0,       0,1404230 ,  kZ*1.917e+03,   "Z+jets",  8000064,   "Z0Jets_TuneD6T_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z0Jets_TuneD6T_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  // {0,     7,  0,  0,   0,       0,  11094 ,  kZ*2.403e-04,   "Z+jets",  8000081,   "Z2Jets_ptZ-800to1600_TuneD6T_7TeV-alpgen-tauola.root",  "START38_V12::All",  "1E31",   "3_8_5", "/Z2Jets_ptZ-800to1600_TuneD6T_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  
  
  {0,    7,  0,  0,   0,        0,10218854,         31314,   "W+jets",  4051001,   "WJets_7TeV-madgraph-tauola.root",                          "START38_V14::All",  "1E31",   "3_8_7", "/WJets_7TeV-madgraph-tauola/Summer10-START36_V9_S09-v1/GEN-SIM-RECO"},
 
  {0,     7,  0,  0,   0,       0, 1099550,          165,   "ttbar",  1000454,   "TT_TuneZ2_7TeV-pythia6-tauola.root",                       "START38_V12::All",  "1E31",   "3_8_5", "/TT_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 1164732,          165,   "ttbar",  4000040,   "TTJets_TuneZ2_7TeV-madgraph-tauola.root",                  "START38_V12::All",  "1E31",   "3_8_5", "/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
 
  {0,     7,  0,  0,   0,       0, 2194752,           18.2,      "WZ",  1000031,   "WZtoAnything_TuneZ2_7TeV-pythia6-tauola.root",             "START38_V12::All",  "1E31",   "3_8_5", "/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 2113368,            5.9,      "ZZ",  1000032,   "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola.root",             "START38_V12::All",  "1E31",   "3_8_7", "/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 2061760,             43,      "WW",  1000030,   "WWtoAnything_TuneZ2_7TeV-pythia6-tauola.root",             "START38_V12::All",  "1E31",   "3_8_5", "/WWtoAnything_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
 
  {0,     7,  0,  0,   0,       0,  494961,           10.6,      "tW",  4000018,   "TToBLNu_TuneZ2_tW-channel_7TeV-madgraph.root",             "START38_V14::All",  "1E31",   "3_8_7", "/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  871720,           0.1835,    "tW",  9000374,   "TW_dr_7TeV-mcatnlo.root",                                  "START38_V14::All",  "1E31",   "3_8_7", "/TW_dr_7TeV-mcatnlo/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  861319,           0.1782,    "tW",  9000375,   "TW_ds_7TeV-mcatnlo.root",                                  "START38_V14::All",  "1E31",   "3_8_7", "/TW_dr_7TeV-mcatnlo/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  
  // background samples for systematics study:
  {0,     7,  0,  0,   0,       0,  1662884,      2205,     "Z+jets",   4000036,    "DYJetsToLL_TuneD6T_matchingdown_7TeV-madgraph-tauola.root",    "START38_V12::All",  "",   "3_8_5", "/DYJetsToLL_TuneD6T_matchingdown_7TeV-madgraph-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1667367,      2366,     "Z+jets",   4000037,    "DYJetsToLL_TuneD6T_matchingup_7TeV-madgraph-tauola.root",      "START38_V12::All",  "",   "3_8_5", "/DYJetsToLL_TuneD6T_matchingup_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1436150,      2541,     "Z+jets",   4000034,    "DYJetsToLL_TuneD6T_scaledown_7TeV-madgraph-tauola.root",       "START38_V12::All",  "",   "3_8_5", "/DYJetsToLL_TuneD6T_scaledown_7TeV-madgraph-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1329028,      2583,     "Z+jets",   4000035,    "DYJetsToLL_TuneD6T_scaleup_7TeV-madgraph-tauola.root",         "START38_V12::All",  "",   "3_8_5", "/DYJetsToLL_TuneD6T_scaleup_7TeV-madgraph-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  
  {0,     7,  0,  0,   0,       0,  1209681,      94.6,     "ttbar",    4000024,    "TTJets_TuneD6T_smallerISRFSR_7TeV-madgraph-tauola.root",       "START38_V12::All",  "",   "3_8_5", "/TTJets_TuneD6T_smallerISRFSR_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1394010,     100.5,     "ttbar",    4000025,    "TTJets_TuneD6T_largerISRFSR_7TeV-madgraph-tauola.root",        "START38_V12::All",  "",   "3_8_5", "/TTJets_TuneD6T_largerISRFSR_7TeV-madgraph-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,   938005,     111.1,     "ttbar",    4000028,    "TTJets_TuneD6T_matchingdown_7TeV-madgraph-tauola.root",        "START38_V12::All",  "",   "3_8_5", "/TTJets_TuneD6T_matchingdown_7TeV-madgraph-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1036968,     105.9,     "ttbar",    4000029,    "TTJets_TuneD6T_matchingup_7TeV-madgraph-tauola.root",          "START38_V12::All",  "",   "3_8_5", "/TTJets_TuneD6T_matchingup_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1096071,     186.6,     "ttbar",    4000026,    "TTJets_TuneD6T_scaledown_7TeV-madgraph-tauola.root",           "START38_V12::All",  "",   "3_8_5", "/TTJets_TuneD6T_scaledown_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1146546,      75.0,     "ttbar",    4000027,    "TTJets_TuneD6T_scaleup_7TeV-madgraph-tauola.root",             "START38_V12::All",  "",   "3_8_5", "/TTJets_TuneD6T_scaleup_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1098152,     106.4,     "ttbar",    1000481,    "TT_smallerISRFSR_TuneZ2_7TeV-pythia6-tauola.root",             "START38_V12::All",  "",   "3_8_5", "/TT_smallerISRFSR_TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  1100000,     108.3,     "ttbar",    1000480,    "TT_largerISRFSR-TuneZ2_7TeV-pythia6-tauola.root",              "START38_V12::All",  "",   "3_8_5", "/TT_largerISRFSR-TuneZ2_7TeV-pythia6-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  
  {0,     7,  0,  0,   0,       0,  6303252,     27230,     "W+jets",   4000031,    "WJets_TuneD6T_scaleup_7TeV-madgraph-tauola.root",              "START38_V12::All",  "",   "3_8_5", "/WJets_TuneD6T_scaleup_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  5165729,     26530,     "W+jets",   4000030,    "WJets_TuneD6T_scaledown_7TeV-madgraph-tauola.root",            "START38_V12::All",  "",   "3_8_5", "/WJets_TuneD6T_scaledown_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  9882119,     25950,     "W+jets",   4000033,    "WJets_TuneD6T_matchingup_7TeV-madgraph-tauola.root",           "START38_V12::All",  "",   "3_8_5", "/WJets_TuneD6T_matchingup_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  2689051,     24990,     "W+jets",   4000032,    "WJets_TuneD6T_matchingdown_7TeV-madgraph-tauola.root",         "START38_V12::All",  "",   "3_8_5", "/WJets_TuneD6T_matchingdown_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  

  {0,     7,  0,  0,   0,       0,  4199091, kW*2.024e+04,     "W+0 Jets",     8000001,    "W0Jets_TuneZ2_7TeV-alpgen-tauola.root",             "START38_V14::All",  "",   "3_8_7_patch2", "/W0Jets_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},

  {0,     7,  0,  0,   0,       0,   883023,      kW*3693,     "W+1 Jets",   8000002,    "W1Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W1Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 21060182,     kW*71.97,     "W+1 Jets",   8000007,    "W1Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W1Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,  5694804,    kW*0.5658,     "W+1 Jets",   8000012,    "W1Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W1Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,   317150,   kW*0.00109,     "W+1 Jets",   8000017,    "W1Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W1Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},

  {0,     7,  0,  0,   0,       0,   2065969,     kW*943.4,     "W+2 Jets",   8000003,    "W2Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W2Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,   1658420,     kW*67.18,     "W+2 Jets",   8000008,    "W2Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W2Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    745242,    kW*0.8553,     "W+2 Jets",   8000013,    "W2Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W2Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,     39900,    kW*0.0022,     "W+2 Jets",   8000018,    "W2Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W2Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},

  {0,     7,  0,  0,   0,       0,   783313,    kW*208.7,     "W+3 Jets",   8000004,    "W3Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W3Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    33974,     kW*32.43,     "W+3 Jets",   8000009,    "W3Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W3Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,   174115,    kW*0.6229,     "W+3 Jets",   8000014,    "W3Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W3Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    62501,  kW*0.001974,     "W+3 Jets",   8000019,    "W3Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W3Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},

  {0,     7,  0,  0,   0,       0,   154846,     kW*44.46,     "W+4 Jets",   8000005,    "W4Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W4Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    75536,     kW*11.38,     "W+4 Jets",   8000010,    "W4Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W4Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    26679,     kW*0.295,     "W+4 Jets",   8000015,    "W4Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W4Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    21932,  kW*0.001025,     "W+4 Jets",   8000020,    "W4Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W4Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
 
  {0,     7,  0,  0,   0,       0,    41761,     kW*11.11,     "W+5 Jets",   8000006,    "W5Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W5Jets_ptW-0to100_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    54565,     kW*3.789,     "W+5 Jets",   8000011,    "W5Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W5Jets_ptW-100to300_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    23959,    kW*0.1565,     "W+5 Jets",   8000016,    "W5Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W5Jets_ptW-300to800_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0,    15163, kW*0.0005882,     "W+5 Jets",   8000021,    "W5Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/W5Jets_ptW-800to1600_TuneZ2_7TeV-alpgen-tauola/Fall10-START38_V12-v2/GEN-SIM-RECO"},
  {0,     7,  0,  0,   0,       0, 14923740,     31314,"W + Jets -> Leptons",4000041, "WJetsToLNu_TuneZ2_7TeV-madgraph-tauola.root",         "START38_V14::All",  "",   "3_8_7_patch2", "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall10-START38_V12-v1/GEN-SIM-RECO"},
// === private background samples ================================================================================================================ 	 
// type  ECM  MW  MNu  channel  lumi   stat      CS        name     pid    fname                              AlCa               HLTmenu   CMSSW   DBS 	 
  // based on Zjets_pT20_MZ180 configuration from heavynu_cfi.py 	 
  {0,     7,  0,  0,   0,       0,    83650,    0.9277,   "Z+jets",   0,    "Zjets_pT20_MZ180.root",          "START36_V10::All",  "8E29",   "3_6_3_patch2", ""},
  // Private production, see directory /Alpgen , electron and muon decay modes
  {0,     7,  0,  0,   0,       0,     4000, 2.*0.1644,   "Z+jets",   0,    "ZjetsAlpgen_Mll180.root",        "START38_V12::All",  "1E31",   "3_8_5",        ""},
  // FastSim production by Minnesota (see config in FastSim/)
  {0,     7,  0,  0,   0,       0,  9900000,     157.5,    "ttbar",   0,    "38x-ttbar-tauola-fastsim.root",  "MC_38Y_V10::All",   "1E31",   "3_8_2",        ""},


// === signal =====================================================================================================================================================================
// type  ECM   MW    MNu  channel  lumi   stat     CS       name                 pid    fname                              AlCa                HLTmenu  CMSSW          DBS
// --- 7 TeV ---  
{1,7,700,100,2,0,10000,1.32*3.745,"WR700_NuRmu100",0,"heavyNuReco_WR700_nuRmu100_001-heavyNuAnalysis.root","","","",""},
{1,7,700,200,2,0,10000,1.32*3.332,"WR700_NuRmu200",0,"heavyNuReco_WR700_nuRmu200_001-heavyNuAnalysis.root","","","",""},
{1,7,700,300,2,0,10000,1.32*2.798,"WR700_NuRmu300",0,"heavyNuReco_WR700_nuRmu300_001-heavyNuAnalysis.root","","","",""},
{1,7,700,400,2,0,10000,1.32*2.026,"WR700_NuRmu400",0,"heavyNuReco_WR700_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,700,500,2,0,10000,1.32*1.169,"WR700_NuRmu500",0,"heavyNuReco_WR700_nuRmu500_001-heavyNuAnalysis.root","","","",""},
{1,7,700,600,2,0,10000,1.32*3.946E-1,"WR700_NuRmu600",0,"heavyNuReco_WR700_nuRmu600_001-heavyNuAnalysis.root","","","",""},

{1,7,800,100,2,0,10000,1.32*2.100,"WR800_NuRmu100",0,"heavyNuReco_WR800_nuRmu100_001-heavyNuAnalysis.root","","","",""},
{1,7,800,200,2,0,10000,1.32*1.923,"WR800_NuRmu200",0,"heavyNuReco_WR800_nuRmu200_001-heavyNuAnalysis.root","","","",""},
{1,7,800,300,2,0,10000,1.32*1.662,"WR800_NuRmu300",0,"heavyNuReco_WR800_nuRmu300_001-heavyNuAnalysis.root","","","",""},
{1,7,800,400,2,0,10000,1.32*1.348,"WR800_NuRmu400",0,"heavyNuReco_WR800_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,800,500,2,0,10000,1.32*9.480E-1,"WR800_NuRmu500",0,"heavyNuReco_WR800_nuRmu500_001-heavyNuAnalysis.root","","","",""},
{1,7,800,600,2,0,10000,1.32*5.358E-1,"WR800_NuRmu600",0,"heavyNuReco_WR800_nuRmu600_001-heavyNuAnalysis.root","","","",""},
{1,7,800,700,2,0,10000,1.32*1.737E-1,"WR800_NuRmu700",0,"heavyNuReco_WR800_nuRmu700_001-heavyNuAnalysis.root","","","",""},

{1,7,900,100,2,0,10000,1.32*1.238,"WR900_NuRmu100",0,"heavyNuReco_WR900_nuRmu100_001-heavyNuAnalysis.root","","","",""},
{1,7,900,200,2,0,10000,1.32*1.154,"WR900_NuRmu200",0,"heavyNuReco_WR900_nuRmu200_001-heavyNuAnalysis.root","","","",""},
{1,7,900,300,2,0,10000,1.32*1.024,"WR900_NuRmu300",0,"heavyNuReco_WR900_nuRmu300_001-heavyNuAnalysis.root","","","",""},
{1,7,900,400,2,0,10000,1.32*8.700E-1,"WR900_NuRmu400",0,"heavyNuReco_WR900_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,900,500,2,0,10000,1.32*6.825E-1,"WR900_NuRmu500",0,"heavyNuReco_WR900_nuRmu500_001-heavyNuAnalysis.root","","","",""},
{1,7,900,600,2,0,10000,1.32*4.759E-1,"WR900_NuRmu600",0,"heavyNuReco_WR900_nuRmu600_001-heavyNuAnalysis.root","","","",""},
{1,7,900,700,2,0,10000,1.32*2.552E-1,"WR900_NuRmu700",0,"heavyNuReco_WR900_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,900,800,2,0,10000,1.32*8.271E-2,"WR900_NuRmu800",0,"heavyNuReco_WR900_nuRmu800_001-heavyNuAnalysis.root","","","",""},

{1,7,1000,100,2,0,10000,1.31*7.394E-1,"WR1000_NuRmu100",0,"heavyNuReco_WR1000_nuRmu100_1-heavyNuAnalysis.root","","","",""},
{1,7,1000,200,2,0,10000,1.31*7.007E-1,"WR1000_NuRmu200",0,"heavyNuReco_WR1000_nuRmu200_1-heavyNuAnalysis.root","","","",""},
{1,7,1000,300,2,0,10000,1.31*6.359E-1,"WR1000_NuRmu300",0,"heavyNuReco_WR1000_nuRmu300_1-heavyNuAnalysis.root","","","",""},
{1,7,1000,400,2,0,10000,1.31*5.620E-1,"WR1000_NuRmu400",0,"heavyNuReco_WR1000_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,1000,500,2,0,10000,1.31*4.724E-1,"WR1000_NuRmu500",0,"heavyNuReco_WR1000_nuRmu500_1-heavyNuAnalysis.root","","","",""},
{1,7,1000,600,2,0,10000,1.31*3.583E-1,"WR1000_NuRmu600",0,"heavyNuReco_WR1000_nuRmu600_001-heavyNuAnalysis.root","","","",""},
{1,7,1000,700,2,0,10000,1.31*2.405E-1,"WR1000_NuRmu700",0,"heavyNuReco_WR1000_nuRmu700_1-heavyNuAnalysis.root","","","",""},
{1,7,1000,800,2,0,10000,1.31*1.281E-1,"WR1000_NuRmu800",0,"heavyNuReco_WR1000_nuRmu800_1-heavyNuAnalysis.root","","","",""},
{1,7,1000,900,2,0,10000,1.31*4.142E-2,"WR1000_NuRmu900",0,"heavyNuReco_WR1000_nuRmu900_1-heavyNuAnalysis.root","","","",""},

{1,7,1100,100,2,0,10000,1.30*4.624E-1,"WR1100_NuRmu100",0,"heavyNuReco_WR1100_nuRmu100_1-heavyNuAnalysis.root","","","",""},
{1,7,1100,200,2,0,10000,1.30*4.405E-1,"WR1100_NuRmu200",0,"heavyNuReco_WR1100_nuRmu200_1-heavyNuAnalysis.root","","","",""},
{1,7,1100,300,2,0,10000,1.30*4.081E-1,"WR1100_NuRmu300",0,"heavyNuReco_WR1100_nuRmu300_001-heavyNuAnalysis.root","","","",""},
{1,7,1100,400,2,0,10000,1.30*3.704E-1,"WR1100_NuRmu400",0,"heavyNuReco_WR1100_nuRmu400_1-heavyNuAnalysis.root","","","",""},
{1,7,1100,500,2,0,10000,1.30*3.184E-1,"WR1100_NuRmu500",0,"heavyNuReco_WR1100_nuRmu500_001-heavyNuAnalysis.root","","","",""},
{1,7,1100,600,2,0,10000,1.30*2.629E-1,"WR1100_NuRmu600",0,"heavyNuReco_WR1100_nuRmu600_1-heavyNuAnalysis.root","","","",""},
{1,7,1100,700,2,0,10000,1.30*1.967E-1,"WR1100_NuRmu700",0,"heavyNuReco_WR1100_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,1100,800,2,0,10000,1.30*1.296E-1,"WR1100_NuRmu800",0,"heavyNuReco_WR1100_nuRmu800_1-heavyNuAnalysis.root","","","",""},
{1,7,1100,900,2,0,10000,1.30*6.833E-2,"WR1100_NuRmu900",0,"heavyNuReco_WR1100_nuRmu900_1-heavyNuAnalysis.root","","","",""},
{1,7,1100,1000,2,0,10000,1.30*2.150E-2,"WR1100_NuRmu1000",0,"heavyNuReco_WR1100_nuRmu1000_1-heavyNuAnalysis.root","","","",""},

{1,7,1200,100,2,0,10000,1.30*2.968E-1,"WR1200_NuRmu100",0,"heavyNuReco_WR1200_nuRmu100_1-heavyNuAnalysis.root","","","",""},
{1,7,1200,200,2,0,10000,1.30*2.810E-1,"WR1200_NuRmu200",0,"heavyNuReco_WR1200_nuRmu200_1-heavyNuAnalysis.root","","","",""},
{1,7,1200,300,2,0,10000,1.30*2.639E-1,"WR1200_NuRmu300",0,"heavyNuReco_WR1200_nuRmu300_1-heavyNuAnalysis.root","","","",""},
{1,7,1200,400,2,0,10000,1.30*2.426E-1,"WR1200_NuRmu400",0,"heavyNuReco_WR1200_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,1200,500,2,0,10000,1.30*2.141E-1,"WR1200_NuRmu500",0,"heavyNuReco_WR1200_nuRmu500_1-heavyNuAnalysis.root","","","",""},
{1,7,1200,600,2,0,10000,1.30*1.828E-1,"WR1200_NuRmu600",0,"heavyNuReco_WR1200_nuRmu600_001-heavyNuAnalysis.root","","","",""},
{1,7,1200,700,2,0,10000,1.30*1.477E-1,"WR1200_NuRmu700",0,"heavyNuReco_WR1200_nuRmu700_1-heavyNuAnalysis.root","","","",""},
{1,7,1200,800,2,0,10000,1.30*1.106E-1,"WR1200_NuRmu800",0,"heavyNuReco_WR1200_nuRmu800_001-heavyNuAnalysis.root","","","",""},
{1,7,1200,900,2,0,10000,1.30*7.154E-2,"WR1200_NuRmu900",0,"heavyNuReco_WR1200_nuRmu900_1-heavyNuAnalysis.root","","","",""},
{1,7,1200,1000,2,0,10000,1.30*3.732E-2,"WR1200_NuRmu1000",0,"heavyNuReco_WR1200_nuRmu1000_1-heavyNuAnalysis.root","","","",""},
{1,7,1200,1100,2,0,10000,1.30*1.170E-2,"WR1200_NuRmu1100",0,"heavyNuReco_WR1200_nuRmu1100_1-heavyNuAnalysis.root","","","",""},

{1,7,1300,100,2,0,10000,1.29*1.890E-1,"WR1300_NuRmu100",0,"heavyNuReco_WR1300_nuRmu100_1-heavyNuAnalysis.root","","","",""},
{1,7,1300,200,2,0,10000,1.29*1.820E-1,"WR1300_NuRmu200",0,"heavyNuReco_WR1300_nuRmu200_1-heavyNuAnalysis.root","","","",""},
{1,7,1300,300,2,0,10000,1.29*1.711E-1,"WR1300_NuRmu300",0,"heavyNuReco_WR1300_nuRmu300_001-heavyNuAnalysis.root","","","",""},
{1,7,1300,400,2,0,10000,1.29*1.586E-1,"WR1300_NuRmu400",0,"heavyNuReco_WR1300_nuRmu400_1-heavyNuAnalysis.root","","","",""},
{1,7,1300,500,2,0,10000,1.29*1.442E-1,"WR1300_NuRmu500",0,"heavyNuReco_WR1300_nuRmu500_001-heavyNuAnalysis.root","","","",""},
{1,7,1300,600,2,0,10000,1.29*1.277E-1,"WR1300_NuRmu600",0,"heavyNuReco_WR1300_nuRmu600_1-heavyNuAnalysis.root","","","",""},
{1,7,1300,700,2,0,10000,1.29*1.060E-1,"WR1300_NuRmu700",0,"heavyNuReco_WR1300_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,1300,800,2,0,10000,1.29*8.495E-2,"WR1300_NuRmu800",0,"heavyNuReco_WR1300_nuRmu800_1-heavyNuAnalysis.root","","","",""},
{1,7,1300,900,2,0,10000,1.29*6.152E-2,"WR1300_NuRmu900",0,"heavyNuReco_WR1300_nuRmu900_001-heavyNuAnalysis.root","","","",""},
{1,7,1300,1000,2,0,10000,1.29*4.052E-2,"WR1300_NuRmu1000",0,"heavyNuReco_WR1300_nuRmu1000_1-heavyNuAnalysis.root","","","",""},
{1,7,1300,1100,2,0,10000,1.29*2.080E-2,"WR1300_NuRmu1100",0,"heavyNuReco_WR1300_nuRmu1100_1-heavyNuAnalysis.root","","","",""},
{1,7,1300,1200,2,0,10000,1.29*6.482E-3,"WR1300_NuRmu1200",0,"heavyNuReco_WR1300_nuRmu1200_1-heavyNuAnalysis.root","","","",""},

{1,7,1400,100,2,0,10000,1.28*1.254E-1,"WR1400_NuRmu100",0,"heavyNuReco_WR1400_nuRmu100_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,200,2,0,10000,1.28*1.209E-1,"WR1400_NuRmu200",0,"heavyNuReco_WR1400_nuRmu200_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,300,2,0,10000,1.28*1.154E-1,"WR1400_NuRmu300",0,"heavyNuReco_WR1400_nuRmu300_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,400,2,0,10000,1.28*1.059E-1,"WR1400_NuRmu400",0,"heavyNuReco_WR1400_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,1400,500,2,0,10000,1.28*9.789E-2,"WR1400_NuRmu500",0,"heavyNuReco_WR1400_nuRmu500_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,600,2,0,10000,1.28*8.614E-2,"WR1400_NuRmu600",0,"heavyNuReco_WR1400_nuRmu600_001-heavyNuAnalysis.root","","","",""},
{1,7,1400,700,2,0,10000,1.28*7.595E-2,"WR1400_NuRmu700",0,"heavyNuReco_WR1400_nuRmu700_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,800,2,0,10000,1.28*6.310E-2,"WR1400_NuRmu800",0,"heavyNuReco_WR1400_nuRmu800_001-heavyNuAnalysis.root","","","",""},
{1,7,1400,900,2,0,10000,1.28*4.999E-2,"WR1400_NuRmu900",0,"heavyNuReco_WR1400_nuRmu900_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,1000,2,0,10000,1.28*3.616E-2,"WR1400_NuRmu1000",0,"heavyNuReco_WR1400_nuRmu1000_001-heavyNuAnalysis.root","","","",""},
{1,7,1400,1100,2,0,10000,1.28*2.317E-2,"WR1400_NuRmu1100",0,"heavyNuReco_WR1400_nuRmu1100_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,1200,2,0,10000,1.28*1.169E-2,"WR1400_NuRmu1200",0,"heavyNuReco_WR1400_nuRmu1200_1-heavyNuAnalysis.root","","","",""},
{1,7,1400,1300,2,0,10000,1.28*3.658E-3,"WR1400_NuRmu1300",0,"heavyNuReco_WR1400_nuRmu1300_1-heavyNuAnalysis.root","","","",""},

{1,7,1500,100,2,0,10000,1.28*8.357E-2,"WR1500_NuRmu100",0,"heavyNuReco_WR1500_nuRmu100_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,200,2,0,10000,1.28*7.963E-2,"WR1500_NuRmu200",0,"heavyNuReco_WR1500_nuRmu200_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,300,2,0,10000,1.28*7.591E-2,"WR1500_NuRmu300",0,"heavyNuReco_WR1500_nuRmu300_001-heavyNuAnalysis.root","","","",""},
{1,7,1500,400,2,0,10000,1.28*7.111E-2,"WR1500_NuRmu400",0,"heavyNuReco_WR1500_nuRmu400_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,500,2,0,10000,1.28*6.571E-2,"WR1500_NuRmu500",0,"heavyNuReco_WR1500_nuRmu500_001-heavyNuAnalysis.root","","","",""},
{1,7,1500,600,2,0,10000,1.28*6.007E-2,"WR1500_NuRmu600",0,"heavyNuReco_WR1500_nuRmu600_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,700,2,0,10000,1.28*5.250E-2,"WR1500_NuRmu700",0,"heavyNuReco_WR1500_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,1500,800,2,0,10000,1.28*4.543E-2,"WR1500_NuRmu800",0,"heavyNuReco_WR1500_nuRmu800_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,900,2,0,10000,1.28*3.777E-2,"WR1500_NuRmu900",0,"heavyNuReco_WR1500_nuRmu900_001-heavyNuAnalysis.root","","","",""},
{1,7,1500,1000,2,0,10000,1.28*2.965E-2,"WR1500_NuRmu1000",0,"heavyNuReco_WR1500_nuRmu1000_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,1100,2,0,10000,1.28*2.129E-2,"WR1500_NuRmu1100",0,"heavyNuReco_WR1500_nuRmu1100_001-heavyNuAnalysis.root","","","",""},
{1,7,1500,1200,2,0,10000,1.28*1.350E-2,"WR1500_NuRmu1200",0,"heavyNuReco_WR1500_nuRmu1200_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,1300,2,0,10000,1.28*6.783E-3,"WR1500_NuRmu1300",0,"heavyNuReco_WR1500_nuRmu1300_1-heavyNuAnalysis.root","","","",""},
{1,7,1500,1400,2,0,10000,1.28*2.138E-3,"WR1500_NuRmu1400",0,"heavyNuReco_WR1500_nuRmu1400_1-heavyNuAnalysis.root","","","",""},

{1,7,1600,100,2,0,10000,1.28*5.563E-2,"WR1600_NuRmu100",0,"heavyNuReco_WR1600_nuRmu100_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,200,2,0,10000,1.28*5.365E-2,"WR1600_NuRmu200",0,"heavyNuReco_WR1600_nuRmu200_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,300,2,0,10000,1.28*5.080E-2,"WR1600_NuRmu300",0,"heavyNuReco_WR1600_nuRmu300_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,400,2,0,10000,1.28*4.837E-2,"WR1600_NuRmu400",0,"heavyNuReco_WR1600_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,1600,500,2,0,10000,1.28*4.487E-2,"WR1600_NuRmu500",0,"heavyNuReco_WR1600_nuRmu500_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,600,2,0,10000,1.28*4.143E-2,"WR1600_NuRmu600",0,"heavyNuReco_WR1600_nuRmu600_001-heavyNuAnalysis.root","","","",""},
{1,7,1600,700,2,0,10000,1.28*3.698E-2,"WR1600_NuRmu700",0,"heavyNuReco_WR1600_nuRmu700_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,800,2,0,10000,1.28*3.304E-2,"WR1600_NuRmu800",0,"heavyNuReco_WR1600_nuRmu800_001-heavyNuAnalysis.root","","","",""},
{1,7,1600,900,2,0,10000,1.28*2.812E-2,"WR1600_NuRmu900",0,"heavyNuReco_WR1600_nuRmu900_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,1000,2,0,10000,1.28*2.267E-2,"WR1600_NuRmu1000",0,"heavyNuReco_WR1600_nuRmu1000_001-heavyNuAnalysis.root","","","",""},
{1,7,1600,1100,2,0,10000,1.28*1.781E-2,"WR1600_NuRmu1100",0,"heavyNuReco_WR1600_nuRmu1100_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,1200,2,0,10000,1.28*1.259E-2,"WR1600_NuRmu1200",0,"heavyNuReco_WR1600_nuRmu1200_001-heavyNuAnalysis.root","","","",""},
{1,7,1600,1300,2,0,10000,1.28*7.799E-3,"WR1600_NuRmu1300",0,"heavyNuReco_WR1600_nuRmu1300_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,1400,2,0,10000,1.28*3.936E-3,"WR1600_NuRmu1400",0,"heavyNuReco_WR1600_nuRmu1400_1-heavyNuAnalysis.root","","","",""},
{1,7,1600,1500,2,0,10000,1.28*1.248E-3,"WR1600_NuRmu1500",0,"heavyNuReco_WR1600_nuRmu1500_1-heavyNuAnalysis.root","","","",""},

{1,7,1800,100,2,0,10000,1.27*0.02590000,"WR1800_nuRmu100",0,"heavyNuReco_WR1800_nuRmu100_001-heavyNuAnalysis.root","","","",""},
{1,7,1800,400,2,0,10000,1.27*0.02228000,"WR1800_nuRmu400",0,"heavyNuReco_WR1800_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,1800,700,2,0,10000,1.27*0.01811000,"WR1800_nuRmu700",0,"heavyNuReco_WR1800_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,1800,1000,2,0,10000,1.27*0.01241000,"WR1800_nuRmu1000",0,"heavyNuReco_WR1800_nuRmu1000_001-heavyNuAnalysis.root","","","",""},
{1,7,1800,1300,2,0,10000,1.27*0.00649100,"WR1800_nuRmu1300",0,"heavyNuReco_WR1800_nuRmu1300_001-heavyNuAnalysis.root","","","",""},
{1,7,1800,1600,2,0,10000,1.27*0.00141700,"WR1800_nuRmu1600",0,"heavyNuReco_WR1800_nuRmu1600_001-heavyNuAnalysis.root","","","",""},

{1,7,2000,100,2,0,10000,1.26*0.01229000,"WR2000_nuRmu100",0,"heavyNuReco_WR2000_nuRmu100_001-heavyNuAnalysis.root","","","",""},
{1,7,2000,400,2,0,10000,1.26*0.01044000,"WR2000_nuRmu400",0,"heavyNuReco_WR2000_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,2000,700,2,0,10000,1.26*0.00862500,"WR2000_nuRmu700",0,"heavyNuReco_WR2000_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,2000,1000,2,0,10000,1.26*0.00645900,"WR2000_nuRmu1000",0,"heavyNuReco_WR2000_nuRmu1000_001-heavyNuAnalysis.root","","","",""},
{1,7,2000,1300,2,0,10000,1.26*0.00404100,"WR2000_nuRmu1300",0,"heavyNuReco_WR2000_nuRmu1300_001-heavyNuAnalysis.root","","","",""},
{1,7,2000,1600,2,0,10000,1.26*0.00168900,"WR2000_nuRmu1600",0,"heavyNuReco_WR2000_nuRmu1600_001-heavyNuAnalysis.root","","","",""},
{1,7,2000,1900,2,0,10000,1.26*0.00016150,"WR2000_nuRmu1900",0,"heavyNuReco_WR2000_nuRmu1900_001-heavyNuAnalysis.root","","","",""},

{1,7,2200,100,2,0,10000,1.25*0.00606200,"WR2200_nuRmu100",0,"heavyNuReco_WR2200_nuRmu100_001-heavyNuAnalysis.root","","","",""},
{1,7,2200,400,2,0,10000,1.25*0.00497600,"WR2200_nuRmu400",0,"heavyNuReco_WR2200_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,2200,700,2,0,10000,1.25*0.00413300,"WR2200_nuRmu700",0,"heavyNuReco_WR2200_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,2200,1000,2,0,10000,1.25*0.00323000,"WR2200_nuRmu1000",0,"heavyNuReco_WR2200_nuRmu1000_001-heavyNuAnalysis.root","","","",""},
{1,7,2200,1300,2,0,10000,1.25*0.00227000,"WR2200_nuRmu1300",0,"heavyNuReco_WR2200_nuRmu1300_001-heavyNuAnalysis.root","","","",""},
{1,7,2200,1600,2,0,10000,1.25*0.00123400,"WR2200_nuRmu1600",0,"heavyNuReco_WR2200_nuRmu1600_001-heavyNuAnalysis.root","","","",""},
{1,7,2200,1900,2,0,10000,1.25*0.00038980,"WR2200_nuRmu1900",0,"heavyNuReco_WR2200_nuRmu1900_001-heavyNuAnalysis.root","","","",""},

{1,7,2400,100,2,0,10000,1.24*0.00315500,"WR2400_nuRmu100",0,"heavyNuReco_WR2400_nuRmu100_001-heavyNuAnalysis.root","","","",""},
{1,7,2400,400,2,0,10000,1.24*0.00241800,"WR2400_nuRmu400",0,"heavyNuReco_WR2400_nuRmu400_001-heavyNuAnalysis.root","","","",""},
{1,7,2400,700,2,0,10000,1.24*0.00196200,"WR2400_nuRmu700",0,"heavyNuReco_WR2400_nuRmu700_001-heavyNuAnalysis.root","","","",""},
{1,7,2400,1000,2,0,10000,1.24*0.00158700,"WR2400_nuRmu1000",0,"heavyNuReco_WR2400_nuRmu1000_001-heavyNuAnalysis.root","","","",""},
{1,7,2400,1300,2,0,10000,1.24*0.00117100,"WR2400_nuRmu1300",0,"heavyNuReco_WR2400_nuRmu1300_001-heavyNuAnalysis.root","","","",""},
{1,7,2400,1600,2,0,10000,1.24*0.00078450,"WR2400_nuRmu1600",0,"heavyNuReco_WR2400_nuRmu1600_001-heavyNuAnalysis.root","","","",""},
{1,7,2400,1900,2,0,10000,1.24*0.00035130,"WR2400_nuRmu1900",0,"heavyNuReco_WR2400_nuRmu1900_001-heavyNuAnalysis.root","","","",""},
{1,7,2400,2200,2,0,10000,1.24*0.00007056,"WR2400_nuRmu2200",0,"heavyNuReco_WR2400_nuRmu2200_001-heavyNuAnalysis.root","","","",""},

 };

// macro for array length calculation
#define DIM(a) (sizeof(a)/sizeof(a[0]))

// class for samples array managing
class SampleDB {
public:
  SampleDB()
  : info("All Samples"),
    samples(AllSamples, AllSamples + DIM(AllSamples))
  {
    // TODO: add a procedure to check samples
    // for example, all file names have to be uniq
  }
  
  // return i-th sample
  Sample operator [] (int i) const {return samples[i];}
  
  // number of samples in database
  int size() const {return samples.size();}
  
  // select samples for which "par" == "val"
  template<typename T>
  SampleDB find(const std::string par, const T val) const
  {
    std::vector<Sample> s;
    
    for (size_t i = 0; i < samples.size(); ++i)
      if (samples[i].check(par, val))
        s.push_back(samples[i]);
    
    return SampleDB(info + " / " + par + " = " + toString(val), s);
  }
  
  // union operation
  SampleDB operator + (const SampleDB& db) const
  {
    std::vector<Sample> s(samples);
    
    for (int i = 0; i < db.size(); ++i)
      if (std::find(s.begin(), s.end(), db[i]) == s.end())
        s.push_back(db[i]);
    
    return SampleDB("(" + info + ") + (" + db.info + ")", s);
  }
  
  std::string info;

private:
  SampleDB(const std::string i, const std::vector<Sample> s)
  : info(i), samples(s)
  {}
  
  std::vector<Sample> samples;
  
  // convert to string
  template<typename T>
  std::string toString(const T x) const
  {
    std::ostringstream ss;
    ss << x;
    return ss.str();
  }
  
  // overload of << operator to allow easy printing of SampleDB
  friend std::ostream& operator<<(std::ostream& os, const SampleDB& db)
  {
    os << "[SampleDB]: " << db.info << ", total = " << db.size() << " samples\n";
    
    for (int i = 0; i < db.size(); ++i) {
      const Sample s = db[i];
      os << "  " << i << ": "
           << s.name << "  " << s.fname << "\n";
    }
    
    return os;
  }
};

#endif
