/*
makeLimitFile

Some important constants are set at the top of the file.
*/

#include <math.h>
#include <stdio.h>
#include "makeLimitFile.hh"
#include "systematics.h"

struct BShape { BShape(double v, double s) : value(v),slope(s) {}
  double value;
  double slope;
};

// name of the histogram containing the observations
const char* data_hist_name = "hNu/cut6_mWRmass/mWR";
// name of the histogram containing the signal
const char* signal_hist_name = "hNuMu40/cut6_mWRmass/mWR";
const char* signal_norm_hist = "hNuMu40/mc_type";

const double bkgd_tt2011A[] = {25.462737,9.741921,3.770282,1.470655,0.575594,
			       0.227078,0.091100,0.037637,0.016665,0.009696};
const double bkgd_tt2011B[] = {27.368073,10.262338,3.887219,1.482203,0.566266,
			       0.217628,0.085133,0.034789,0.015650,0.009899};

const double bkgd_zj2011A[] = {14.641240,7.760743,4.087397,2.131959,1.107135,
			       0.575533,0.294940,0.150599,0.076322,0.037205};

const double bkgd_zj2011B[] = {18.616888,9.249375,4.569226,2.257969,1.113659,
			       0.544038,0.265238,0.129253,0.063021,0.030772};

const double bkgd_other2011A[] = {1.828816,0.827960,0.384629,0.182891,0.088731,
				  0.043773,0.021884,0.011056,0.005630,0.002885};

const double bkgd_other2011B[] = {1.874635,0.612362,0.226846,0.090535,0.037523,
				  0.015828,0.006729,0.002871,0.001226,0.000524};

// names
const char* jnames[]= {"WR","TT","ZJ","OT"};
const char* snames[]= {"SIGNAL_%d_%d","TTJETS","ZJETS","OTHER"};
const int jmax=3;
const double* jbkgd2011A[]={bkgd_tt2011A,bkgd_zj2011A,bkgd_other2011A};
const double* jbkgd2011B[]={bkgd_tt2011B,bkgd_zj2011B,bkgd_other2011B};

// Histogram manipulation
const double minimum_signal_content=0.01;
const double bkgd_norm_low=600.0;
const double bkgd_norm_high=2000.0;

#include <stdio.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const MassPoint& mp, const char* limitFileName, const SystematicsDB& syst) {
  double* bkgdh[100];
  const int nbins=int(pbi.size());
  char temp[128];
  std::vector<std::string> procName;
  std::vector<std::string> systematicsList=syst.getSystematicsList();

  // set up the official names for the systematics DB
  sprintf(temp,snames[0],mp.mwr,mp.mnr);
  procName.push_back(temp);
  for (int j=1; j<=jmax; j++) procName.push_back(snames[j]);
  
  //  TF1* f1=new TF1("f1","exp([0]*x)");
 
  for (int j=0; j<jmax; j++) {
    bkgdh[j]=new double[nbins];
    
    for (int ib=0; ib<nbins; ib++) {
      bkgdh[j][ib]=0;

      int srcBin=pbi[ib].sourceBin;

      bkgdh[j][ib]+=jbkgd2011A[j][srcBin];
      bkgdh[j][ib]+=jbkgd2011B[j][srcBin];

    }
  }
  
  FILE* limitFile=fopen(limitFileName,"wt");

  if (limitFile==0) {
    fprintf(stderr,"Unable to open '%s' for writing\n",limitFileName);
    return;
  }
  
  fprintf(limitFile,"imax %d\n",nbins);
  fprintf(limitFile,"jmax %d  # tt and zjets\n",jmax);
  fprintf(limitFile,"kmax %d \n",int(systematicsList.size()));  

  // these are the data loops
  fprintf(limitFile,"bin            ");
  for (int ibin=0; ibin<nbins; ibin++) fprintf(limitFile,"%4s ",pbi[ibin].binName.c_str());
  fprintf(limitFile,"\nobservation    ");
  for (int ibin=0; ibin<nbins; ibin++) fprintf(limitFile," %3d ",pbi[ibin].data);
  fprintf(limitFile,"\n\n");

  // these are the signal and background loops
  fprintf(limitFile,"bin            ");
  for (int ibin=0; ibin<nbins; ibin++)  
    for (int j=0; j<=jmax; j++) fprintf(limitFile," %4s ",pbi[ibin].binName.c_str());
  fprintf(limitFile,"\n");
  fprintf(limitFile,"process        ");
  for (int ibin=0; ibin<nbins; ibin++)  
    for (int j=0; j<=jmax; j++) fprintf(limitFile," %3s  ",jnames[j]);
  fprintf(limitFile,"\n");
  fprintf(limitFile,"process        ");
  for (int ibin=0; ibin<nbins; ibin++)  
    for (int j=0; j<=jmax; j++) fprintf(limitFile," %2d   ",j);
  fprintf(limitFile,"\n");
  fprintf(limitFile,"rate           ");
  for (int ibin=0; ibin<nbins; ibin++) {
    fprintf(limitFile,"%5.2f ", pbi[ibin].signal);
    for (int j=1; j<=jmax; j++) 
      if (bkgdh[j-1][ibin]<0.001)  fprintf(limitFile,"%5.3f ",0.001); //,bkgdh[j-1]->GetBinContent(ibin));
      else fprintf(limitFile,"%5.3f ",bkgdh[j-1][ibin]);
  }
  fprintf(limitFile,"\n");


  // systematics
  for (std::vector<std::string>::const_iterator i=systematicsList.begin();
       i!=systematicsList.end(); i++) {
    fprintf(limitFile,"%-10s lnN ",i->c_str());
    for (int ibin=0; ibin<nbins; ibin++) {
      int srcBin=pbi[ibin].sourceBin;
      for (int j=0; j<=jmax; j++) fprintf(limitFile,"%3.3f ",syst.getSystematic(*i,procName[j],srcBin));
    }
    fprintf(limitFile,"\n");
  }
  

  fclose(limitFile);

  
}

std::vector<double> extractBins(TFile* f, const std::string& histname) {
  std::vector<double> retval(10,0);
  TH1* h=(TH1*)(f->Get(histname.c_str()));
  if (h!=0) {
    for (int jbin=0; jbin<10; jbin++) 
      retval[jbin]=h->Integral(16+5*jbin,16+4+5*jbin);    
  }
  return retval;
}

std::vector<PerBinInfo> makeLimitContent(double lumi, double xsec, const MassPoint& mp, TFile* dataf, TFile* signalf) {
  //  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->GetBinContent(3);  //  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<double> vsignal, vdata;

  vsignal=extractBins(signalf,signal_hist_name);
  vdata=extractBins(dataf,data_hist_name);

  std::vector<PerBinInfo> pbi;

  int ilow=0;
  int ihigh=9;

  switch (mp.mwr) {
  case (700) : ihigh=4; break;
  case (800) : ihigh=4; break;
  case (900) : ihigh=4; break;
  case (1000) : ihigh=5; break;
  case (1100) : ihigh=5; break;
  case (1200) : ihigh=5; break;
  case (1300) : ihigh=6; break;
  case (1400) : ilow=1; ihigh=6; break;
  case (1500) : ilow=1; ihigh=6; break;
  case (1600) : ilow=1; ihigh=7; break;
  case (1700) :
  case (1800) : ilow=1; ihigh=7; break;
  case (1900) :
  case (2000) : ilow=2; ihigh=8; break;
  case (2100) :
  case (2200) : ilow=2; ihigh=9; break;
  case (2300) :
  case (2400) : ilow=3; ihigh=9; break;
  case (2500) : ilow=3; ihigh=9; break;
  };
    
  
  for (int ibin=ilow; ibin<=ihigh; ibin++) {
    PerBinInfo abin;
    abin.lowEdge=ibin*200+600;
    abin.highEdge=(ibin+1)*200+600;
    abin.signal=vsignal[ibin]*normSignal*lumi*xsec;
    abin.sourceBin=ibin;
    abin.lumi=lumi;
    abin.year=2011;
    abin.data=int(vdata[ibin]);
    char name[10];
    sprintf(name,"b%02d",ibin);
    abin.binName=name;
    if (abin.signal>0.01) 
      pbi.push_back(abin);
  }
  return pbi;
}

void makeLimitFile(double lumi, double xsec, const MassPoint& mp, TFile* dataf, TFile* signalf, const char* limitFileName, const SystematicsDB& syst) {
  std::vector<PerBinInfo> pbi = makeLimitContent(lumi,xsec,mp,dataf,signalf);
  formatLimitFile(pbi,mp,limitFileName,syst);
}

