/*
makeLimitFile

Some important constants are set at the top of the file.
*/

#include <math.h>
#include <stdio.h>
#include "ratedb.hh"
#include "makeLimitFile.hh"
#include "systematics.h"

struct BShape { BShape(double v, double s) : value(v),slope(s) {}
  double value;
  double slope;
};

// name of the histogram containing the observations
const char* data_hist_name = "hNu/cut6_mWRmass/mWR";

// names
const char* jnames[]= {"WR","TT","ZJ","OT"};
const char* snames[]= {"SIGNAL_%d_%d","TTJETS","ZJETS","OTHER"};
const int jmax=3;

// Histogram manipulation
const double minimum_signal_content=0.01;

#include <stdio.h>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

std::string whichSyst() { return "MUON"; }

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const LimitPoint& mp, const char* limitFileName, const SystematicsDB& syst) {
  const int nbins=int(pbi.size());
  char temp[128];
  std::vector<std::string> procName;
  std::vector<std::string> systematicsList=syst.getSystematicsList();

  // set up the official names for the systematics DB
  sprintf(temp,snames[0],mp.mwr_syst,mp.mnr_syst);
  procName.push_back(temp);
  for (int j=1; j<=jmax; j++) procName.push_back(snames[j]);
   
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
      if (pbi[ibin].bkgd[j-1]<0.01)  fprintf(limitFile,"%5.2f ",0.01); //,bkgdh[j-1]->GetBinContent(ibin));
      else fprintf(limitFile,"%5.2f ",pbi[ibin].bkgd[j-1]);
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
  std::vector<double> retval(11,0);
  TH1* h=(TH1*)(f->Get(histname.c_str()));
  if (h!=0) {
    for (int jbin=0; jbin<11; jbin++) 
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
  case (11) : ilow=1; ihigh=5; break;
  case (12) : ilow=1; ihigh=5; break;
  case (13) : ilow=1; ihigh=6; break;
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
  case (24) : ilow=3; ihigh=10; break;
  case (25) : ilow=3; ihigh=10; break;
  case (26) : ilow=3; ihigh=10; break;
  case (27) : ilow=4; ihigh=10; break;
  case (28) : ilow=4; ihigh=10; break;
  case (29) : ilow=4; ihigh=10; break;
  case (30) : ilow=4; ihigh=10; break;
  };
  
  
}


std::vector<PerBinInfo> makeLimitContent(const LimitPoint& mp, TFile* dataf, const RateDB& db, bool fullRange) {
  //  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
  //  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->GetBinContent(3);  //  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<double> vdata;

  //  vsignal=extractBins(signalf,signal_hist_name);
  vdata=extractBins(dataf,data_hist_name);

  std::vector<PerBinInfo> pbi;

  int ilow=0;
  int ihigh=9;

  if (!fullRange) binRanger(mp.mwr,ilow,ihigh);
  char process[200];
  sprintf(process,snames[0],mp.mwr,mp.mnr);
  
  char signame[20];
  sprintf(signame,"eff_%d",mp.year);

  for (int ibin=ilow; ibin<=ihigh; ibin++) {
    PerBinInfo abin;
    abin.lowEdge=ibin*200+600;
    abin.highEdge=(ibin+1)*200+600;
    abin.signal=db.get(process,"sigeff",ibin)*mp.lumi*mp.xsec;
    for (int j=1; j<=3; j++) 
      if (mp.year==2011) {
	abin.bkgd[j-1]=db.get(snames[j],"2011A",ibin)+
	  db.get(snames[j],"2011B",ibin);    
      } else {
	abin.bkgd[j-1]=db.get(snames[j],"2012",ibin);
      }
    abin.sourceBin=ibin;
    abin.lumi=mp.lumi;
    abin.year=mp.year;
    abin.data=int(vdata[ibin]);
    char name[10];
    sprintf(name,"b%02d",ibin);
    abin.binName=name;
    if (abin.signal>0.01 || fullRange) 
      pbi.push_back(abin);
  }
  return pbi;
}

void makeLimitFile(const LimitPoint& mp, TFile* dataf, const RateDB& rates, const char* limitFileName, const SystematicsDB& syst) {
  std::vector<PerBinInfo> pbi = makeLimitContent(mp,dataf,rates);
  formatLimitFile(pbi,mp,limitFileName,syst);
}

void makeLimitFileInterpolate(const LimitPoint& pt, TFile* dataf,
			      const RateDB& ratedb,
			      const LimitPoint& signalp1, 
			      const LimitPoint& signalp2, 
			      const char* limitFileName, const SystematicsDB& syst) {

  std::vector<PerBinInfo> pbi1=makeLimitContent(signalp1,dataf,ratedb,true);
  std::vector<PerBinInfo> pbi2=makeLimitContent(signalp2,dataf,ratedb,true);

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
