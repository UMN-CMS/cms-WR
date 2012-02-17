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
const char* data_hist_name[] = {"Data_A","Data_B",0};
// name of the histogram containing the signal
const char* signal_hist_name = "Signal";
const char* signal_norm_hist = 0;
// functional parameters for ttbar (total per ipb, exponential slope)
const BShape bkgd_tt2011(70.0/(2173+2511), -5.0e-3);
// functional parameters for z+jets (total per ipb, exponential slope)
const BShape bkgd_zj2011(48.0/(2173+2511), -3.5e-3);
// functional parameters for other backgrounds (w+jets, VV, QCD, tW)
const BShape bkgd_other2011(12.0/(2173+2511), -4.8e-3);

// names
const char* snames[]= {"SIGNAL_%d_%d","TTJETS","ZJETS","OTHER"};
const char* jnames[]= {"WR","TT","ZJ","OT"};
const int jmax=3;
const BShape jbkgd2011[]={bkgd_tt2011,bkgd_zj2011,bkgd_other2011};

// Histogram manipulation
const double minimum_signal_content=0.01;
const double bkgd_norm_low=600.0;
const double bkgd_norm_high=3000.0;

#include <stdio.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

std::string whichSyst() { return "ELEC"; }

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const LimitPoint& mp, const char* limitFileName, const SystematicsDB& syst) {

  double* bkgdh[100];
  char temp[128];
  const int nbins=int(pbi.size());
  std::vector<std::string> procName;
  std::vector<std::string> systematicsList=syst.getSystematicsList();

  // set up the official names for the systematics DB
  sprintf(temp,snames[0],mp.mwr_syst,mp.mnr_syst);
  procName.push_back(temp);
  for (int j=1; j<=jmax; j++) procName.push_back(snames[j]);
  
  TF1* f1=new TF1("f1","exp([0]*x)");
  
  for (int j=0; j<jmax; j++) {
    bkgdh[j]=new double[nbins];
    
    for (int ib=0; ib<nbins; ib++) {
      double v=jbkgd2011[j].value;
      double s=jbkgd2011[j].slope;

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
  TH1* datah=(TH1*)(dataf->Get(data_hist_name[0])->Clone("datah"));

  for (int i=1; data_hist_name[i]!=0; i++) {
    datah->Add((TH1*)(dataf->Get(data_hist_name[i])));
  }

  TH1* sigh=(TH1*)(signalf->Get(signal_hist_name)->Clone("sigh"));
  double normSignal=1.0;//((TH1*)(signalf->Get(signal_norm_hist)))->Integral();

  std::vector<PerBinInfo> pbi;

  int ilow=1;
  int ihigh=10;
  
  binRanger(mp.mwr,ilow,ihigh);

  
  for (int ibin=1; ibin<=sigh->GetNbinsX(); ibin++) {
    
    if (sigh->GetBinContent(ibin)<minimum_signal_content && !fullRange) {
      printf("  Dropping 2011 bin %d with %f (%f)\n",ibin,sigh->GetBinContent(ibin),minimum_signal_content);
      continue;
    }
    
    if ((ibin<ilow || ibin>ihigh) && !fullRange) continue;

    PerBinInfo abin;
    abin.lowEdge=std::max(600.0,sigh->GetXaxis()->GetBinLowEdge(ibin));
    abin.highEdge=sigh->GetXaxis()->GetBinUpEdge(ibin);
    abin.signal=sigh->GetBinContent(ibin)*normSignal*mp.lumi;
    if (abin.signal<0.01 && !fullRange) continue; // skip very empty bins (less
    abin.lumi=mp.lumi;
    abin.year=2011;
    abin.data=int(datah->GetBinContent(ibin));
    char name[10];
    sprintf(name,"d%02d",ibin);
    abin.binName=name;
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
    //    printf("%f\n",bin.signal);
    pbiFinal.push_back(bin);
  }
  formatLimitFile(pbiFinal,pt,limitFileName,syst);

}
