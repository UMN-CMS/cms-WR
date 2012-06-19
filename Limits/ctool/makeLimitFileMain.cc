#include "makeLimitFile.hh"
#include "ratedb.hh"
#include "systematics.h"
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include <unistd.h>


int main(int argc, char* argv[]) {

  int opt;
  LimitPoint pt,pt1,pt2;
  std::string df_name[2], systf_name[2], of_name, ratedb_name[2];
  char temp1[100],temp2[100];
  std::string istring;
  bool interpolate=false;

  pt.rebin_above_mlljj=-1;
  pt1.rebin_above_mlljj=-1;
  pt2.rebin_above_mlljj=-1;
  pt.mode=LimitPoint::lp_Muon1ECM; // historical default

  while ((opt = getopt(argc, argv, "l:w:n:x:d:o:s:i:I:r:y:m:b:")) != -1) {
               switch (opt) {
	       case 'b':
		 pt.rebin_above_mlljj=atof(optarg);
		 pt1.rebin_above_mlljj=atof(optarg);
		 pt2.rebin_above_mlljj=atof(optarg);
	       case 'l':
		 if (pt.mode==LimitPoint::lp_MuonElec || pt.mode==LimitPoint::lp_Muon2ECM) {
		   sscanf(optarg,"%lf,%lf",&pt1.lumi,&pt2.lumi);
		 } else {
		   pt.lumi=atof(optarg);
		   pt1.lumi=atof(optarg);
		   pt2.lumi=atof(optarg);
		 }
		 break;
	       case 'x':
		 if (pt.mode==LimitPoint::lp_MuonElec || pt.mode==LimitPoint::lp_Muon2ECM) {
		   sscanf(optarg,"%lf,%lf",&pt1.xsec,&pt2.xsec);
		 } else {
		   pt.xsec=atof(optarg);
		   pt1.xsec=atof(optarg);
		   pt2.xsec=atof(optarg);
		 }
		 break;
	       case 'y':
		 if (pt.mode==LimitPoint::lp_Muon2ECM) {
		   sscanf(optarg,"%d,%d",&pt1.year,&pt2.year);
		 } else {
		   pt.year=atoi(optarg);
		   pt1.year=atoi(optarg);
		   pt2.year=atoi(optarg);
		 }
		 break;
	       case 'w':
		 pt.mwr=atoi(optarg);
		 pt.mwr_syst=pt.mwr;
		 pt1.mwr=pt.mwr;
		 pt2.mwr=pt.mwr;
		 break;
	       case 'r':
		 if (pt.mode==LimitPoint::lp_MuonElec || pt.mode==LimitPoint::lp_Muon2ECM) {
		   sscanf(optarg,"%[^,],%s",temp1,temp2);
		   ratedb_name[0]=temp1;
		   ratedb_name[1]=temp2;
		 } else {
		   ratedb_name[0]=optarg;
		 }
		 break;
	       case 'n':
		 if (pt.mode==LimitPoint::lp_Muon2ECM) {
		   sscanf(optarg,"%d,%d",&pt1.mnr,&pt2.mnr);
		   pt.mnr=pt1.mnr;
		 } else {
		   pt.mnr=atoi(optarg);
		   pt.mnr_syst=pt.mnr;
		   pt1.mnr=pt.mnr;
		   pt2.mnr=pt.mnr;
		 }
		 break;
	       case 'd':
		 if (pt.mode==LimitPoint::lp_MuonElec || pt.mode==LimitPoint::lp_Muon2ECM) {
		   sscanf(optarg,"%[^,],%s",temp1,temp2);
		   df_name[0]=temp1;
		   df_name[1]=temp2;
		 } else {
		   df_name[0]=optarg;
		 }
		 break;
	       case 'I':
		 istring=optarg;
		 interpolate=true;
		 break;
	       case 's':
		 if (pt.mode==LimitPoint::lp_MuonElec || pt.mode==LimitPoint::lp_Muon2ECM)  {
		   sscanf(optarg,"%[^,],%s",temp1,temp2);
		   systf_name[0]=temp1;
		   systf_name[1]=temp2;
		 } else {
		   systf_name[0]=optarg;
		 }
		 break;
	       case 'o':
		 of_name=optarg;
		 break;
	       case 'm':
		 if (!strcasecmp(optarg,"mu2ecm")) pt.mode=LimitPoint::lp_Muon2ECM;
		 if (!strcasecmp(optarg,"elec")) pt.mode=LimitPoint::lp_Elec1ECM;
		 if (!strcasecmp(optarg,"e-mu")) pt.mode=LimitPoint::lp_MuonElec;
		 break;
	       default:
		 break;
	       }
  };


		   /*
  if (argc<9) {
    printf("Usage: makeLimitFile [lumi] [mwr] [mn] [xsec] [data file] [signal file] [outfilename] [systematicsdb] [optional special syst mode]\n");
    return 1;
    }*/

  if (pt.mode==LimitPoint::lp_MuonElec) {
    TFile df1(df_name[0].c_str());
    TFile df2(df_name[1].c_str());

    SystematicsDB syst1,syst2;
    syst1.load(systf_name[0].c_str());
    syst1.standardSystematics("MUON"); // define the standard systematics
    syst2.load(systf_name[1].c_str());
    syst2.standardSystematics("ELEC"); // define the standard systematics

    RateDB rates1,rates2;
    rates1.load(ratedb_name[0].c_str());
    rates2.load(ratedb_name[1].c_str());

    pt1.mode=LimitPoint::lp_Muon1ECM;
    pt2.mode=LimitPoint::lp_Elec1ECM;
   
    std::vector<PerBinInfo> bins1=makeLimitContent(pt1, &df1, rates1, syst1, 'm');
    std::vector<PerBinInfo> bins2=makeLimitContent(pt2, &df2, rates2, syst2, 'e');

    std::vector<PerBinInfo> bins;
    for (std::vector<PerBinInfo>::const_iterator i=bins1.begin(); i!=bins1.end(); i++)
      bins.push_back(*i);
    for (std::vector<PerBinInfo>::const_iterator i=bins2.begin(); i!=bins2.end(); i++)
      bins.push_back(*i);

    formatLimitFile(bins,pt1,of_name.c_str());
  } else if (pt.mode==LimitPoint::lp_Muon2ECM) {
    TFile df1(df_name[0].c_str());
    TFile df2(df_name[1].c_str());

    SystematicsDB syst1,syst2;
    syst1.load(systf_name[0].c_str());
    syst1.standardSystematics("MUON"); // define the standard systematics
    syst2.load(systf_name[1].c_str());
    syst2.standardSystematics("MUON"); // define the standard systematics

    RateDB rates1,rates2;
    rates1.load(ratedb_name[0].c_str());
    rates2.load(ratedb_name[1].c_str());

    pt1.mode=LimitPoint::lp_Muon1ECM;
    pt2.mode=LimitPoint::lp_Muon1ECM;
   
    std::vector<PerBinInfo> bins1=makeLimitContent(pt1, &df1, rates1, syst1, 'a');
    std::vector<PerBinInfo> bins2=makeLimitContent(pt2, &df2, rates2, syst2, 'b');

    std::vector<PerBinInfo> bins;
    for (std::vector<PerBinInfo>::const_iterator i=bins1.begin(); i!=bins1.end(); i++)
      bins.push_back(*i);
    for (std::vector<PerBinInfo>::const_iterator i=bins2.begin(); i!=bins2.end(); i++)
      bins.push_back(*i);

    formatLimitFile(bins,pt1,of_name.c_str());
  } else {

    TFile df(df_name[0].c_str());
    SystematicsDB syst;
    syst.load(systf_name[0].c_str());
    syst.standardSystematics(whichSyst()); // define the standard systematics
    /*
      if (argc==10)
      int special_syst_mode=atoi(argv[9]);
    */
    
    RateDB rates;
    rates.load(ratedb_name[0].c_str());
    
    if (!interpolate) {
      makeLimitFile(pt,&df,rates,of_name.c_str(),syst);
    } else {
      sscanf(istring.c_str(),"%d,%d:%d,%d",
	     &pt1.mwr,&pt1.mnr,
	     &pt2.mwr,&pt2.mnr);
      printf("Interpolating %d,%d and %d,%d \n",
	   pt1.mwr,pt1.mnr,   pt2.mwr,pt2.mnr
	     );
      
      // take systematics from the first point, by convention
      pt.mnr_syst=pt1.mnr;    
      pt.mwr_syst=pt1.mwr;    
      
      makeLimitFileInterpolate(pt,&df,rates,pt1,pt2,of_name.c_str(),syst);
      
    }
  }


  return 0;
}
