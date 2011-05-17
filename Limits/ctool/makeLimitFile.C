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
const char* data_hist_name = "hNu/cut8_mWRmass/mWR";
// name of the histogram containing the signal
const char* signal_hist_name = "hNu/cut8_mWRmass/mWR";
const char* signal_norm_hist = "hNu/nelec";
// functional parameters for ttbar (total per ipb, exponential slope)
const BShape bkgd_tt2010(1.03/36.1, -5.53e-3);
const BShape bkgd_tt2011(5.17/191, -4.94e-3);
// functional parameters for z+jets (total per ipb, exponential slope)
const BShape bkgd_zj2010(0.85/36.1, -4.63e-3);
const BShape bkgd_zj2011(5.26/191, -3.97e-3);
// functional parameters for other backgrounds (w+jets, VV, QCD, tW)
const BShape bkgd_other2010((0.03+0.037+0.022)/36.1+0.0026, -5e-3);
const BShape bkgd_other2011((0.19+0.14)/191.0+0.0026, -4.33e-3);

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
  for (int ibin=0; ibin<nbins; ibin++) fprintf(limitFile,"a%02d ",ibin);
  fprintf(limitFile,"\nobservation    ");
  for (int ibin=0; ibin<nbins; ibin++) fprintf(limitFile,"%3d ",pbi[ibin].data);
  fprintf(limitFile,"\n\n");

  // these are the signal and background loops
  fprintf(limitFile,"bin            ");
  for (int ibin=0; ibin<nbins; ibin++)  
    for (int j=0; j<=jmax; j++) fprintf(limitFile," a%02d  ",ibin);
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
      for (int j=0; j<=jmax; j++) fprintf(limitFile,"%s ",getSyst(j,*i)); 
    }
    fprintf(limitFile,"\n");
  }
  

  fclose(limitFile);

  
}

std::vector<PerBinInfo> makeLimitContent2010(TFile* dataf, TFile* signalf) {
  const double lumi=36.1;
  
  TH1* datah=(TH1*)(dataf->Get(data_hist_name)->Clone("datah"));
  TH1* sigh=(TH1*)(signalf->Get(signal_hist_name)->Clone("sigh"));
  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<PerBinInfo> pbi;

  datah->Rebin(5);
  sigh->Rebin(5);
  
  for (int ibin=1; ibin<=sigh->GetNbinsX(); ibin++) {
    if (sigh->GetBinContent(ibin)<min_level_abs) continue;
    PerBinInfo abin;
    abin.lowEdge=sigh->GetXaxis()->GetBinLowEdge(ibin);
    abin.highEdge=sigh->GetXaxis()->GetBinUpEdge(ibin);
    abin.signal=sigh->GetBinContent(ibin)*normSignal*lumi;
    abin.lumi=lumi;
    abin.data=int(datah->GetBinContent(ibin));
    abin.year=2010;
    pbi.push_back(abin);
  }
  return pbi;
}

void makeLimitFile2010(TFile* dataf, TFile* signalf, const char* limitFileName) {
  std::vector<PerBinInfo> pbi=makeLimitContent2010(dataf,signalf);
  formatLimitFile(pbi,limitFileName);
}  

std::vector<PerBinInfo> makeLimitContent2011(double lumi, TFile* dataf, TFile* signalf) {
  TH1* datah=(TH1*)(dataf->Get(data_hist_name)->Clone("datah"));
  TH1* sigh=(TH1*)(signalf->Get(signal_hist_name)->Clone("sigh"));
  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<PerBinInfo> pbi;

  datah->Rebin(5);
  sigh->Rebin(5);
  
  for (int ibin=1; ibin<=sigh->GetNbinsX(); ibin++) {
    if (sigh->GetBinContent(ibin)<min_level_abs) continue;
    PerBinInfo abin;
    abin.lowEdge=sigh->GetXaxis()->GetBinLowEdge(ibin);
    abin.highEdge=sigh->GetXaxis()->GetBinUpEdge(ibin);
    abin.signal=sigh->GetBinContent(ibin)*normSignal*lumi;
    abin.lumi=lumi;
    abin.year=2011;
    abin.data=int(datah->GetBinContent(ibin));
    pbi.push_back(abin);
  }
  return pbi;
}

void makeLimitFile2011(double lumi, TFile* dataf, TFile* signalf, const char* limitFileName) {
  std::vector<PerBinInfo> pbi = makeLimitContent2011(lumi,dataf,signalf);
  formatLimitFile(pbi,limitFileName);
}

void makeLimitFileTwoYear(double lumi11, TFile* dataf11, TFile* signalf11, TFile* dataf10, TFile* signalf10, const char* limitFileName) {

  std::vector<PerBinInfo> pbi2010=makeLimitContent2010(dataf10,signalf10);
  std::vector<PerBinInfo> pbi2011=makeLimitContent2011(lumi11,dataf11,signalf11);

  std::vector<PerBinInfo> pbi;
  std::vector<PerBinInfo>::const_iterator i;

  pbi.insert(pbi.end(),pbi2010.begin(),pbi2010.end());
  pbi.insert(pbi.end(),pbi2011.begin(),pbi2011.end());

  formatLimitFile(pbi,limitFileName);

}


