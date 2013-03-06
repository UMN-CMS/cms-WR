/*
makeLimitFile

Some important constants are set at the top of the file.
*/

#include <math.h>
#include <stdio.h>
#include "ratedb.hh"
#include "makeLimitFile.hh"
#include "systematics.h"
#include <set>

struct BShape { BShape(double v, double s) : value(v),slope(s) {}
  double value;
  double slope;
};

// name of the histogram containing the observations
const char* data_hist_name_muon = "hNu/cut6_mWRmass/mWR";
const char* data_hist_name_elec = "hNuE/cut6_mWRmass/mWR";

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

void formatLimitFile(const std::vector<PerBinInfo>& pbi, const LimitPoint& mp, const char* limitFileName) {
  const int nbins=int(pbi.size());
  std::set<std::string> systematicList;
  int nextraSyst=0;

  // assemble the list of _all_ systematics
  for (int i=0; i<nbins; i++) {
    std::map<std::string,PerBinSystematic>::const_iterator jj;
    for (jj=pbi[i].perBinSyst.begin();jj!=pbi[i].perBinSyst.end();jj++) {
      if (jj->first==SystematicsDB::GAMMASTATS) {
	if (jj->second.signalN>0) nextraSyst++;
	for (int j=0; j<3; j++) if (jj->second.bkgdN[j]>0) nextraSyst++; 
      } 
      systematicList.insert(jj->first);
    }
  }
  if (nextraSyst>0) nextraSyst--; // because the name of the gammastats counts as one

  FILE* limitFile=fopen(limitFileName,"wt");

  if (limitFile==0) {
    fprintf(stderr,"Unable to open '%s' for writing\n",limitFileName);
    return;
  }
  
  fprintf(limitFile,"imax %d\n",nbins);
  fprintf(limitFile,"jmax %d  # tt and zjets\n",jmax);
  fprintf(limitFile,"kmax %d \n",int(systematicList.size())+nextraSyst);  

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
  for (std::set<std::string>::const_iterator i=systematicList.begin();
       i!=systematicList.end(); i++) {
    if (*i==SystematicsDB::GAMMASTATS) {
      for (int ibin=0; ibin<nbins; ibin++) {
	std::map<std::string,PerBinSystematic>::const_iterator pbsi=pbi[ibin].perBinSyst.find(*i);
	if (pbsi==pbi[ibin].perBinSyst.end()) continue;

	if (pbsi->second.signalN>0) {
	  fprintf(limitFile,"gss%s gmN %4d ",pbi[ibin].binName.c_str(),pbsi->second.signalN);
	  for (int iib=0; iib<ibin; iib++) fprintf(limitFile,"   -    -    -    -  "); // blanks
	  fprintf(limitFile,"%5.3f   -   -   - ",pbsi->second.signal);
	  for (int iib=ibin+1; iib<nbins; iib++) fprintf(limitFile,"  -    -    -    -  "); // blanks	  
	  fprintf(limitFile,"\n");
	}
	for (int j=0; j<3; j++) {
	  if (pbsi->second.bkgdN[j]<=0) continue; // no such systematic here

	  fprintf(limitFile,"gs%d%s gmN %4d ",j,pbi[ibin].binName.c_str(),pbsi->second.bkgdN[j]);
	  for (int iib=0; iib<ibin; iib++) fprintf(limitFile,"   -    -    -    -  "); // blanks
	  if (j==0) fprintf(limitFile,"-  %5.3f   -  - ",pbsi->second.bkgd[0]);
	  if (j==1) fprintf(limitFile,"-  -  %5.3f - ",pbsi->second.bkgd[1]);
	  if (j==2) fprintf(limitFile,"-  -   -  %5.3f  ",pbsi->second.bkgd[2]);
	  for (int iib=ibin+1; iib<nbins; iib++) fprintf(limitFile,"  -   -   -   -  "); // blanks	  
	  fprintf(limitFile,"\n");
	  
	}
	
      } 
    } else {
	
      fprintf(limitFile,"%-10s lnN  ",i->c_str());
      for (int ibin=0; ibin<nbins; ibin++) {
	std::map<std::string,PerBinSystematic>::const_iterator pbsi=pbi[ibin].perBinSyst.find(*i);
	if (pbsi==pbi[ibin].perBinSyst.end()) continue;
	
	if (pbsi->second.signal<0) fprintf(limitFile,"  -   ");
	else fprintf(limitFile,"%5.3f ", pbsi->second.signal);
	
	for (int j=1; j<=jmax; j++) {
	  double systLevel=pbsi->second.bkgd[j-1];
	  if (systLevel<0.001 || fabs(systLevel-1.0)<0.0015) fprintf(limitFile,"  -   ");
	  else fprintf(limitFile,"%5.3f ", systLevel);
	}
      }
      
      fprintf(limitFile,"\n");
    }
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

  ilow=0; ihigh=16;

  return;
  
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


std::vector<PerBinInfo> makeLimitContent(const LimitPoint& mp, TFile* dataf, const RateDB& db, const SystematicsDB& syst, char binprefix, bool fullRange) {
  //  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->Integral();
  //  double normSignal=1.0/((TH1*)(signalf->Get(signal_norm_hist)))->GetBinContent(3);  //  double min_level_abs=minimum_signal_content*sigh->Integral();

  std::vector<double> vdata;
  std::vector<std::string> systematicsList=syst.getSystematicsList();

  //  vsignal=extractBins(signalf,signal_hist_name);
  if (mp.mode==LimitPoint::lp_Muon1ECM || mp.mode==LimitPoint::lp_Muon2ECM) 
    vdata=extractBins(dataf,data_hist_name_muon);
  else
    vdata=extractBins(dataf,data_hist_name_elec);

  std::vector<PerBinInfo> pbi, pbi_alt;

  int ilow=0;
  int ihigh=16;


  if (!fullRange) binRanger(mp.mwr,ilow,ihigh);

  if (mp.rebin_above_mlljj>0) ihigh=16;

  int rebinabove=(int(mp.rebin_above_mlljj+0.5)-600)/200;

  char process[200];
  sprintf(process,snames[0],mp.mwr,mp.mnr);
  
  char signame[20];
  sprintf(signame,"eff_%d",mp.year);

  for (int ibin=ilow; ibin<=ihigh; ibin++) {
    PerBinInfo abin;
    if (ibin==0) abin.lowEdge=600;
    else abin.lowEdge=mp.bin_upper_edge[ibin-1];
    abin.highEdge=mp.bin_upper_edge[ibin];

    for (int jk=0; jk<16; jk++) {
      double bcenter=jk*200+700;
      if (bcenter<abin.lowEdge || bcenter>abin.highEdge) continue;
      
      double sigbineff=db.get(process,signame,ibin);
      abin.signal+=sigbineff*mp.lumi*mp.xsec;

      for (int j=1; j<=3; j++) 
	if (mp.year==2011) {
	  abin.bkgd[j-1]+=db.get(snames[j],"2011A",ibin)+
	    db.get(snames[j],"2011B",ibin);    
	} else {
	  abin.bkgd[j-1]+=db.get(snames[j],"2012",ibin);
	}
    }

    // Systematics
    for (std::vector<std::string>::const_iterator isyst=systematicsList.begin();
	 isyst!=systematicsList.end(); isyst++) {
      
      PerBinSystematic pbs;

      for (int jbin=0; jbin<16; jbin++) {
	double bcenter=jbin*200+700;
	if (bcenter<abin.lowEdge || bcenter>abin.highEdge) continue;

	if (*isyst==SystematicsDB::GAMMASTATS) {
	  double systLevel=syst.getSystematic(*isyst,process,jbin);
	  //	if (systLevel<0) {
	  pbs.signal=-1;
	  pbs.signalN=-1;
	  //} // really not ready for weighted signal.  Ignore this case.

	  for (int j=1; j<=3; j++) {
	    systLevel=syst.getSystematic(*isyst,snames[j],jbin);
	    if (systLevel<0.001 || fabs(systLevel-1)<0.0015) {
	      pbs.bkgd[j-1]=-1;
	      pbs.bkgdN[j-1]=0;
	  } else {
	    pbs.bkgdN[j-1]=int(abin.bkgd[j-1]/systLevel+0.51);
	    pbs.bkgd[j-1]=abin.bkgd[j-1]/pbs.bkgdN[j-1];
	  }
	}       
      } else {
	double systLevel=syst.getSystematic(*isyst,process,ibin);
	if (systLevel<0.001 || fabs(systLevel-1)<0.0015) systLevel=-1;
	pbs.signal=systLevel;
      
	for (int j=1; j<=3; j++) {
	  systLevel=syst.getSystematic(*isyst,snames[j],ibin);
	  if (systLevel<0.001 || fabs(systLevel-1)<0.0015) systLevel=-1;
	  pbs.bkgd[j-1]=systLevel;
	}
      }
      abin.perBinSyst.insert(std::pair<std::string,PerBinSystematic>(*isyst,pbs));
    }

    abin.sourceBin=ibin;
    abin.lumi=mp.lumi;
    abin.year=mp.year;
    abin.data=int(vdata[ibin]);
    char name[10];
    sprintf(name,"%c%02d",binprefix,ibin);
    abin.binName=name;
    if (sigbineff>0.01 || fullRange || (mp.rebin_above_mlljj>0 && ibin>=rebinabove)) 
      pbi.push_back(abin);
  }

  
  pbi_alt.swap(pbi);

    int srcabove=(int(mp.rebin_above_mlljj+0.5)-600)/200;

    PerBinInfo abin;
    abin.sourceBin=srcabove;
    abin.year=mp.year;
    abin.lumi=mp.lumi;
    abin.lowEdge=srcabove*200+600;
    abin.highEdge=4000;
    abin.data=0;
    char name[10];
    sprintf(name,"%c%02d",binprefix,srcabove);
    abin.binName=name;

    for (int j=0; j<3; j++) {
	abin.bkgd[j]=0;
    }
    abin.signal=0; 
    
    std::map<std::string,PerBinSystematic>::iterator ss;

    for (std::vector<PerBinInfo>::const_iterator jj=pbi_alt.begin();jj!=pbi_alt.end(); jj++) {
      if (jj->sourceBin<srcabove) pbi.push_back(*jj);
      else {
	abin.data+=jj->data;
	abin.signal+=jj->signal;
	for (int j=0; j<3; j++) abin.bkgd[j]+=jj->bkgd[j];

	std::map<std::string,PerBinSystematic>::const_iterator is;
	for (is=jj->perBinSyst.begin(); is!=jj->perBinSyst.end(); is++) {
	  ss=abin.perBinSyst.find(is->first);
	  if (ss==abin.perBinSyst.end()) {
	    PerBinSystematic blank;
	    blank.signal=0;  blank.signalN=0;
	    blank.bkgd[0]=0;	    blank.bkgd[1]=0;	    blank.bkgd[2]=0;
	    blank.bkgdN[0]=0;	    blank.bkgdN[1]=0;	    blank.bkgdN[2]=0;
	    abin.perBinSyst.insert(std::pair<std::string,PerBinSystematic>(is->first,blank));
	    ss=abin.perBinSyst.find(is->first);
	  }
	  // weighted sum...
	  ss->second.signal+=jj->signal*is->second.signal;
	  ss->second.signalN+=is->second.signalN;
	  for (int q=0; q<3; q++) {
	    ss->second.bkgd[q]+=jj->bkgd[q]*is->second.bkgd[q];
	    ss->second.bkgdN[q]+=is->second.bkgdN[q];
	  }
	}
      }
    }
    // normalize systematics...
    for (ss=abin.perBinSyst.begin(); ss!=abin.perBinSyst.end(); ss++) {
      ss->second.signal/=std::max(1e-9,abin.signal);
      for (int q=0; q<3; q++) 
	if (ss->first==SystematicsDB::GAMMASTATS) 
	  ss->second.bkgd[q]=abin.bkgd[q]/ss->second.bkgdN[q];
	else
	  ss->second.bkgd[q]/=std::max(1e-9,abin.bkgd[q]);      
      
      
    }
    double sigbineff=abin.signal/(mp.lumi*mp.xsec);
    if (sigbineff>0.01)
      pbi.push_back(abin); // combined bin
  }

  return pbi;
}

void makeLimitFile(const LimitPoint& mp, TFile* dataf, const RateDB& rates, const char* limitFileName, const SystematicsDB& syst) {
  std::vector<PerBinInfo> pbi = makeLimitContent(mp,dataf,rates,syst);
  formatLimitFile(pbi,mp,limitFileName);
}

void makeLimitFileInterpolate(const LimitPoint& pt, TFile* dataf,
			      const RateDB& ratedb,
			      const LimitPoint& signalp1, 
			      const LimitPoint& signalp2, 
			      const char* limitFileName, const SystematicsDB& syst) {

  std::vector<PerBinInfo> pbi1=makeLimitContent(signalp1,dataf,ratedb,syst,'b',true);
  std::vector<PerBinInfo> pbi2=makeLimitContent(signalp2,dataf,ratedb,syst,'b',true);

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
  formatLimitFile(pbiFinal,pt,limitFileName);

}
