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
  std::string df_name, systf_name, of_name, sigf_name, ratedb_name;
  std::string istring;
  bool interpolate=false;

  while ((opt = getopt(argc, argv, "l:w:n:x:d:o:s:i:I:r:")) != -1) {
               switch (opt) {
	       case 'l':
		 pt.lumi=atof(optarg);
		 pt1.lumi=atof(optarg);
		 pt2.lumi=atof(optarg);
		 break;
	       case 'x':
		 pt.xsec=atof(optarg);
		 pt1.xsec=atof(optarg);
		 pt2.xsec=atof(optarg);
		 break;
	       case 'w':
		 pt.mwr=atoi(optarg);
		 pt.mwr_syst=pt.mwr;
		 break;
	       case 'r':
		 ratedb_name=optarg;
		 break;
	       case 'n':
		 pt.mnr=atoi(optarg);
		 pt.mnr_syst=pt.mnr;
		 break;
	       case 'd':
		 df_name=optarg;
		 break;
	       case 'i':
		 sigf_name=optarg;
		 break;
	       case 'I':
		 istring=optarg;
		 interpolate=true;
		 break;
	       case 's':
		 systf_name=optarg;
		 break;
	       case 'o':
		 of_name=optarg;
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

  TFile df(df_name.c_str());
  SystematicsDB syst;
  syst.load(systf_name.c_str());
  syst.standardSystematics(whichSyst()); // define the standard systematics
  /*
  if (argc==10)
    int special_syst_mode=atoi(argv[9]);
  */

  RateDB rates;
  rates.load(ratedb_name.c_str());
  
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


  return 0;
}
