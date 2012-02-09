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

const double bkgd_tt2011A[] = {26.925420,10.301448,3.986795,1.555096,0.608636,
			       0.240110,0.096327,0.039796,0.017621,0.010254};

const double bkgd_tt2011B[] = {28.000084,10.499279,3.976956,1.516418,0.579336,
			       0.222650,0.087097,0.035592,0.016012,0.010128};

const double bkgd_zj2011A[] = {14.491357,7.681296,4.045555,2.110133,1.095814,
			       0.569651,0.291925,0.149057,0.075542,0.036823};

const double bkgd_zj2011B[] = {18.538802,9.210629,4.550137,2.248559,1.109002,
			       0.541760,0.264128,0.128712,0.062757,0.030643};

const double bkgd_other2011A[] = {1.564294,0.681379,0.305788,0.140924,0.066417,
				  0.031874,0.015513,0.007631,0.003783,0.001887};

const double bkgd_other2011B[] = {1.652976,0.548266,0.194975,0.073102,0.028415,
				  0.011300,0.004557,0.001852,0.000757,0.000310};

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

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const LimitPoint& mp, const char* limitFileName, const SystematicsDB& syst) {
  double* bkgdh[100];
  const int nbins=int(pbi.size());
  char temp[128];
  std::vector<std::string> procName;
  std::vector<std::string> systematicsList=syst.getSystematicsList();

  // set up the official names for the systematics DB
  sprintf(temp,snames[0],mp.mwr_syst,mp.mnr_syst);
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
  fprintf(limitFile,"rate            ");
  for (int ibin=0; ibin<nbins; ibin++) {
    fprintf(limitFile,"%5.2f ", pbi[ibin].signal);
    for (int j=1; j<=jmax; j++) 
      if (bkgdh[j-1][ibin]<0.01)  fprintf(limitFile,"%5.2f ",0.01); //,bkgdh[j-1]->GetBinContent(ibin));
      else fprintf(limitFile,"%5.2f ",bkgdh[j-1][ibin]);
  }
  fprintf(limitFile,"\n");


  // systematics
  for (std::vector<std::string>::const_iterator i=systematicsList.begin();
       i!=systematicsList.end(); i++) {
    fprintf(limitFile,"%-10s lnN  ",i->c_str());
    for (int ibin=0; ibin<nbins; ibin++) {
      int srcBin=pbi[ibin].sourceBin;
      for (int j=0; j<=jmax; j++) {
	double systLevel=syst.getSystematic(*i,procName[j],srcBin);
	if (systLevel<0.001 || fabs(systLevel-1.0)<0.0015) fprintf(limitFile,"  -   ");
	else fprintf(limitFile,"%5.3f ", systLevel);
      }
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

static void binRanger(int mw, int& ilow, int& ihigh) {
  int mweff=((mw+50)/100);
  
  switch (mweff) {
  case (7) : ihigh=4; break;
  case (8) : ihigh=4; break;
  case (9) : ihigh=4; break;
  case (10) : ihigh=5; break;
  case (11) : ihigh=5; break;
  case (12) : ihigh=5; break;
  case (13) : ihigh=6; break;
  case (14) : ilow=1; ihigh=6; break;
  case (15) : ilow=1; ihigh=6; break;
  case (16) : ilow=1; ihigh=7; break;
  case (17) :
  case (18) : ilow=1; ihigh=7; break;
  case (19) :
  case (20) : ilow=2; ihigh=8; break;
  case (21) :
  case (22) : ilow=2; ihigh=9; break;
  case (23) :
  case (24) : ilow=3; ihigh=9; break;
  case (25) : ilow=3; ihigh=9; break;
  };
  
  
}


std::vector<PerBinInfo> makeLimitContent(const LimitPoint& mp, TFile* dataf, TFile* signalf, bool fullRange) {
  //  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->GetBinContent(3);  //  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<double> vsignal, vdata;

  vsignal=extractBins(signalf,signal_hist_name);
  vdata=extractBins(dataf,data_hist_name);

  std::vector<PerBinInfo> pbi;

  int ilow=0;
  int ihigh=9;

  if (!fullRange) binRanger(mp.mwr,ilow,ihigh);
    
  
  for (int ibin=ilow; ibin<=ihigh; ibin++) {
    PerBinInfo abin;
    abin.lowEdge=ibin*200+600;
    abin.highEdge=(ibin+1)*200+600;
    abin.signal=vsignal[ibin]*normSignal*mp.lumi*mp.xsec;
    abin.sourceBin=ibin;
    abin.lumi=mp.lumi;
    abin.year=2011;
    abin.data=int(vdata[ibin]);
    char name[10];
    sprintf(name,"b%02d",ibin);
    abin.binName=name;
    if (abin.signal>0.01 && !fullRange) 
      pbi.push_back(abin);
  }
  return pbi;
}

void makeLimitFile(const LimitPoint& mp, TFile* dataf, TFile* signalf, const char* limitFileName, const SystematicsDB& syst) {
  std::vector<PerBinInfo> pbi = makeLimitContent(mp,dataf,signalf);
  formatLimitFile(pbi,mp,limitFileName,syst);
}

void makeLimitFileInterpolate(const LimitPoint& pt, TFile* dataf, 
			      TFile* signalf1, const LimitPoint& signalp1, 
			      TFile* signalf2, const LimitPoint& signalp2, 
			      const char* limitFileName, const SystematicsDB& syst) {

  std::vector<PerBinInfo> pbi1=makeLimitContent(signalp1,dataf,signalf1,true);
  std::vector<PerBinInfo> pbi2=makeLimitContent(signalp2,dataf,signalf2,true);

  std::vector<PerBinInfo> pbiFinal;

  int ilow=0;
  int ihigh=9;

  binRanger(pt.mwr,ilow,ihigh);

  for (size_t i=0; i<pbi1.size(); i++) {
    // skip irrelevant mass bins
    if (i<size_t(ilow) || i>size_t(ihigh)) continue;
    
    PerBinInfo bin=pbi1[i];
    bin.signal=pbi1[i].signal+(pt.mwr-signalp1.mwr)*(pbi2[i].signal-pbi1[i].signal)/(signalp2.mwr-signalp1.mwr);
    pbiFinal.push_back(bin);
  }
  formatLimitFile(pbiFinal,pt,limitFileName,syst);

}
