/*
makeLimitFile

Some important constants are set at the top of the file.
*/

#include <math.h>
#include <stdio.h>
#include "makeLimitFile.hh"
#include "systematicsELEC.h"

struct BShape { BShape(double v, double s) : value(v),slope(s) {}
  double value;
  double slope;
};

// name of the histogram containing the observations
const char* data_hist_name = "Data";
// name of the histogram containing the signal
const char* signal_hist_name = "Signal";
const char* signal_norm_hist = 0;
// functional parameters for ttbar (total per ipb, exponential slope)
const BShape bkgd_tt2010(0.68/36.1, -5.8e-3);
const BShape bkgd_tt2011(7.49/191, -6.6e  -3);
// functional parameters for z+jets (total per ipb, exponential slope)
const BShape bkgd_zj2010(0.43/36.1, -4.3e-3);
const BShape bkgd_zj2011(4.34/191, -4.4e-3);
// functional parameters for other backgrounds (w+jets, VV, QCD, tW)
const BShape bkgd_other2010((0.12+0.04+0.03+0.2)/36.1, -4.6e-3);
const BShape bkgd_other2011((0.88+0.35+0.23+0.24)/191.0, -4.6e-3);

// names
const char* jnames[]= {"WR","TT","ZJ","OT"};
const int jmax=3;
const BShape jbkgd2010[]={bkgd_tt2010,bkgd_zj2010,bkgd_other2010};
const BShape jbkgd2011[]={bkgd_tt2011,bkgd_zj2011,bkgd_other2011};

// Histogram manipulation
const double minimum_signal_content=0.01;
const double bkgd_norm_low=520.0;
const double bkgd_norm_high=2000.0;

#include <stdio.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const char* limitFileName) {
  loadSystematics();

  double* bkgdh[100];
  const int nbins=int(pbi.size());
  
  TF1* f1=new TF1("f1","exp([0]*x)");
  
  for (int j=0; j<jmax; j++) {
    bkgdh[j]=new double[nbins];
    
    for (int ib=0; ib<nbins; ib++) {
      double v=(pbi[ib].year==2010)?(jbkgd2010[j].value):(jbkgd2011[j].value);
      double s=(pbi[ib].year==2010)?(jbkgd2010[j].slope):(jbkgd2011[j].slope);

      f1->SetParameter(0,s); 
   

      bkgdh[j][ib]=pbi[ib].lumi*v*
	f1->Integral(pbi[ib].lowEdge,pbi[ib].highEdge)/f1->Integral(bkgd_norm_low,bkgd_norm_high);
    }
  }
  
  FILE* limitFile=fopen(limitFileName,"wt");
  
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
      for (int j=0; j<=jmax; j++) fprintf(limitFile,"%s ",getSyst(j,pbi[ibin].year,*i)); 
    }
    fprintf(limitFile,"\n");
  }
  

  fclose(limitFile);

  
}

std::vector<PerBinInfo> makeLimitContent2010(int mwr,TFile* dataf, TFile* signalf) {
  const double lumi=36.1;
  
  TH1* datah=(TH1*)(dataf->Get(data_hist_name)->Clone("datah"));
  TH1* sigh=(TH1*)(signalf->Get(signal_hist_name)->Clone("sigh"));
  double normSignal=1.0;//((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
  //  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<PerBinInfo> pbi;

  datah->Rebin(5);
  sigh->Rebin(5);

  int ilow=3;
  int ihigh=10;

  switch (mwr) {
  case (700) : ihigh=7; break;
  case (800) : ihigh=7; break;
  case (900) : ihigh=8; break;
  case (1000) : ihigh=8; break;
  case (1100) : ihigh=8; break;
  case (1200) : ihigh=9; break;
  case (1300) : ihigh=9; break;
  case (1400) : ihigh=9; break;
  case (1500) : ilow=4; break;
  case (1600) : ilow=5; break;
  case (1700) : ilow=5; break;
  case (1800) : ilow=5; break;
  case (1900) : ilow=5; break;
  case (2000) : ilow=5; break;
  case (2100) : ilow=5; break;
  case (2200) : ilow=5; break;
  case (2300) : ilow=5; break;
  case (2400) : ilow=5; break;
  case (2500) : ilow=5; break;
  case (2600) : ilow=5; break;
  };
  
  for (int ibin=1; ibin<=sigh->GetNbinsX(); ibin++) {
    /*
      if (sigh->GetBinContent(ibin)<min_level_abs) {
      printf("  Dropping 2010 bin %d with %f (%f)\n",ibin,sigh->GetBinContent(ibin),min_level_abs);
      continue; 
    }
    */
    if (ibin<ilow || ibin>ihigh) continue;
    PerBinInfo abin;
    abin.lowEdge=std::max(520.0,sigh->GetXaxis()->GetBinLowEdge(ibin));
    abin.highEdge=sigh->GetXaxis()->GetBinUpEdge(ibin);
    abin.signal=sigh->GetBinContent(ibin)*normSignal*lumi;
    if (abin.signal<0.01) continue; // skip very empty bins (less than 0.01 event expected)  These cause formatting problems
    abin.lumi=lumi;
    abin.data=int(datah->GetBinContent(ibin));
    abin.year=2010;
    char name[10];
    sprintf(name,"c%02d",ibin);
    abin.binName=name;
    pbi.push_back(abin);
  }
  return pbi;
}

void makeLimitFile2010(int mwr, TFile* dataf, TFile* signalf, const char* limitFileName) {
  std::vector<PerBinInfo> pbi=makeLimitContent2010(mwr,dataf,signalf);
  formatLimitFile(pbi,limitFileName);
}  

std::vector<PerBinInfo> makeLimitContent2011(double lumi, int mwr, TFile* dataf, TFile* signalf) {
  TH1* datah=(TH1*)(dataf->Get(data_hist_name)->Clone("datah"));
  TH1* sigh=(TH1*)(signalf->Get(signal_hist_name)->Clone("sigh"));
  double normSignal=1.0;//((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
    //  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<PerBinInfo> pbi;

  datah->Rebin(5);
  sigh->Rebin(5);

  int ilow=3;
  int ihigh=10;

  switch (mwr) {
  case (700) : ihigh=6; break;
  case (800) : ihigh=7; break;
  case (900) : ihigh=7; break;
  case (1000) : ihigh=8; break;
  case (1100) : ihigh=8; break;
  case (1200) : ilow=4; ihigh=9; break;
  case (1300) : ilow=4; ihigh=9; break;
  case (1400) : ilow=4; ihigh=9; break;
  case (1500) : ilow=4; ihigh=10; break;
  case (1600) : ilow=4; ihigh=10; break;
  case (1700) : ilow=4; ihigh=10; break;
  case (1800) : ilow=4; ihigh=11; break;
  case (1900) : ilow=4; ihigh=11; break;
  case (2000) : ilow=5; ihigh=12; break;
  case (2100) : ilow=5; ihigh=12; break;
  case (2200) : ilow=6; ihigh=13; break;
  case (2300) : ilow=6; ihigh=13; break;
  case (2400) : ilow=7; ihigh=14; break;
  case (2500) : ilow=7; ihigh=14; break;
  };
    
  
  for (int ibin=1; ibin<=sigh->GetNbinsX(); ibin++) {
    /*
    if (sigh->GetBinContent(ibin)<min_level_abs) {
      printf("  Dropping 2011 bin %d with %f (%f)\n",ibin,sigh->GetBinContent(ibin),min_level_abs);
      continue;
    }
    */
    if (ibin<ilow || ibin>ihigh) continue;

    PerBinInfo abin;
    abin.lowEdge=std::max(520.0,sigh->GetXaxis()->GetBinLowEdge(ibin));
    abin.highEdge=sigh->GetXaxis()->GetBinUpEdge(ibin);
    abin.signal=sigh->GetBinContent(ibin)*normSignal*lumi;
    if (abin.signal<0.01) continue; // skip very empty bins (less
    abin.lumi=lumi;
    abin.year=2011;
    abin.data=int(datah->GetBinContent(ibin));
    char name[10];
    sprintf(name,"d%02d",ibin);
    abin.binName=name;
    pbi.push_back(abin);
  }
  return pbi;
}

void makeLimitFile2011(double lumi, int mwr, TFile* dataf, TFile* signalf, const char* limitFileName) {
  std::vector<PerBinInfo> pbi = makeLimitContent2011(lumi,mwr,dataf,signalf);
  formatLimitFile(pbi,limitFileName);
}

void makeLimitFileTwoYear(double lumi11, int mwr, TFile* dataf11, TFile* signalf11, TFile* dataf10, TFile* signalf10, const char* limitFileName) {

  std::vector<PerBinInfo> pbi2010=makeLimitContent2010(mwr,dataf10,signalf10);
  std::vector<PerBinInfo> pbi2011=makeLimitContent2011(lumi11,mwr,dataf11,signalf11);

  std::vector<PerBinInfo> pbi;
  std::vector<PerBinInfo>::const_iterator i;

  pbi.insert(pbi.end(),pbi2010.begin(),pbi2010.end());
  pbi.insert(pbi.end(),pbi2011.begin(),pbi2011.end());

  formatLimitFile(pbi,limitFileName);

}


