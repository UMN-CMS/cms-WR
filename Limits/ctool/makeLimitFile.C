/*
makeLimitFile

Some important constants are set at the top of the file.
*/

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
const BShape bkgd_tt(1.03/36.1, -5.53e-3);
// functional parameters for z+jets (total per ipb, exponential slope)
const BShape bkgd_zj(0.85/36.1, -4.63e-3);

// names
const char* jnames[]= {"WR","TT","ZJ"};
const int jmax=2;
const BShape jbkgd[]={bkgd_tt,bkgd_zj};

// Histogram manipulation
const double minimum_signal_content=0.01;
int rebinLevel = 5;
int firstBin = 3;
int lastBin = 10;

#include <stdio.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

TH1* bkgdh[100];

void makeLimitFile(double lumi, TFile* dataf, TFile* signalf, const char* limitFileName) {
  
  TH1* datah=(TH1*)(dataf->Get(data_hist_name)->Clone("datah"));
  TH1* sigh=(TH1*)(signalf->Get(signal_hist_name)->Clone("sigh"));
  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->Integral();

  datah->Rebin(rebinLevel);
  sigh->Rebin(rebinLevel);

  sigh->Scale(normSignal);

  TF1* f1=new TF1("f1","exp([0]*x)");
  for (int j=0; j<jmax; j++) {
    f1->SetParameter(0,jbkgd[j].slope); 
    bkgdh[j]=(TH1*)(sigh->Clone(jnames[j+1]));
    for (int ib=1; ib<=bkgdh[j]->GetNbinsX(); ib++) bkgdh[j]->SetBinContent(ib,0);
    double sum=0;
    for (int ib=firstBin; ib<=lastBin; ib++) {
      double val=(*f1)(bkgdh[j]->GetBinCenter(ib));
      sum+=val;
      bkgdh[j]->SetBinContent(ib,val);
    }
    bkgdh[j]->Scale(jbkgd[j].value*lumi/sum);    
  }
  
  bool useBin[1000];
  int nbins=0;

  for (int ibin=firstBin; ibin<=lastBin; ibin++) {
    useBin[ibin]=(sigh->GetBinContent(ibin)>minimum_signal_content*sigh->Integral());
    if (useBin[ibin]) nbins++;
  }

  sigh->Scale(lumi);

  FILE* limitFile=fopen(limitFileName,"wt");
  
  fprintf(limitFile,"imax %d\n",nbins);
  fprintf(limitFile,"jmax %d  # tt and zjets\n",jmax);
  fprintf(limitFile,"kmax 0  # no systematics yet\n");

  // these are the data loops
  fprintf(limitFile,"bin         ");
  for (int ibin=firstBin; ibin<=lastBin; ibin++)  if (useBin[ibin]) fprintf(limitFile,"a%02d ",ibin-firstBin);
  fprintf(limitFile,"\nobservation ");
  for (int ibin=firstBin; ibin<=lastBin; ibin++)  if (useBin[ibin]) fprintf(limitFile,"%3d ",int(datah->GetBinContent(ibin)));
  fprintf(limitFile,"\n\n");

  // these are the signal and background loops
  fprintf(limitFile,"bin         ");
  for (int ibin=firstBin; ibin<=lastBin; ibin++)  
    if (useBin[ibin])
      for (int j=0; j<=jmax; j++) fprintf(limitFile,"a%02d    ",ibin-firstBin);
  fprintf(limitFile,"\n");
  fprintf(limitFile,"process     ");
  for (int ibin=firstBin; ibin<=lastBin; ibin++)  
    if (useBin[ibin])
      for (int j=0; j<=jmax; j++) fprintf(limitFile," %3s   ",jnames[j]);
  fprintf(limitFile,"\n");
  fprintf(limitFile,"process     ");
  for (int ibin=firstBin; ibin<=lastBin; ibin++)  
    if (useBin[ibin])
      for (int j=0; j<=jmax; j++) fprintf(limitFile," %2d    ",j);
  fprintf(limitFile,"\n");
  fprintf(limitFile,"rate        ");
  for (int ibin=firstBin; ibin<=lastBin; ibin++) {
    if (!useBin[ibin]) continue;
    fprintf(limitFile,"%4.2f ",sigh->GetBinContent(ibin));
    for (int j=1; j<=jmax; j++) fprintf(limitFile,"%4.2e ",bkgdh[j-1]->GetBinContent(ibin));
  }
  fprintf(limitFile,"\n");
  

  fclose(limitFile);


}



