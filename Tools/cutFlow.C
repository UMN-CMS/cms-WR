//Hello c++!
#include "TH1.h"
#include "TFile.h"

#include <iostream>
#include <iomanip>
#include "stdio.h"
#include <vector>
#include <utility>
#include <string>

void cutFlow(bool eff = false)
{
    using namespace std;
    
    double integral, preInt = 1.0;
    
    vector<pair<string,TFile*> > fnames;
    vector<pair<string,string> > hnames;

    hnames.push_back(pair<string,string>("none","hNu/cut0_none/mWR"));
    hnames.push_back(pair<string,string>("LLJJ Pt","hNu/cut1_LLJJpt/mWR"));
    hnames.push_back(pair<string,string>("trig","hNu/cut2_TrigMatches/mWR"));
    hnames.push_back(pair<string,string>("vertex","hNu/cut3_Vertex/mWR"));
    hnames.push_back(pair<string,string>("mu1 pt","hNu/cut4_Mu1HighPt/mWR"));
    hnames.push_back(pair<string,string>("Mll","hNu/cut5_diLmass/mWR"));
    hnames.push_back(pair<string,string>("MWR","hNu/cut6_mWRmass/mWR"));

    fnames.push_back(pair<string,TFile*>("ttbar powheg",   new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_TTTo2L2Nu2B_7TeV-powheg-pythia6.root")));
    fnames.push_back(pair<string,TFile*>("ttbar madgraph", new TFile("/local/cms/user/dahmes/wr2011/bgMC/Summer11/aug30/ttbar-PFJets.root")));
    fnames.push_back(pair<string,TFile*>("Z+Jets sherpa",  new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa.root")));
    fnames.push_back(pair<string,TFile*>("W+Jets sherpa",  new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WToLNu_7TeV-sherpa.root")));
    fnames.push_back(pair<string,TFile*>("ZZ",             new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_ZZ_TuneZ2_7TeV_pythia6_tauola.root")));
    fnames.push_back(pair<string,TFile*>("WZ",             new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WZ_TuneZ2_7TeV_pythia6_tauola.root")));
    fnames.push_back(pair<string,TFile*>("WW",             new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_WW_TuneZ2_7TeV_pythia6_tauola.root")));
    fnames.push_back(pair<string,TFile*>("tW",             new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola.root")));
    fnames.push_back(pair<string,TFile*>("tbarW",          new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola.root")));
    
	printf("%8s", "cut");
	for(vector<pair<string,TFile*> >::const_iterator i = fnames.begin(); i != fnames.end(); i++)
	{
		printf(eff?" & %15s        ":" & %15s", i->first.c_str());
	}
	printf(" \\\\ \\hline\n");
	for(vector<pair<string,string> >::const_iterator j = hnames.begin(); j != hnames.end(); j++)
	{
		printf("%8s", j->first.c_str());
		for(vector<pair<string,TFile*> >::const_iterator i = fnames.begin(); i != fnames.end(); i++)
		{
			TH1* h = (TH1*)i->second->Get(j->second.c_str());
			integral = h->Integral(0, h->GetNbinsX()+1);
			
			printf((j == hnames.begin() && eff)?" & %15.0f        ":" & %15.0f", integral);
			if(j != hnames.begin() && eff) 
			{
				TH1* hpre = (TH1*)i->second->Get((j-1)->second.c_str());
				preInt = hpre->Integral(0, h->GetNbinsX()+1);
				printf("  (%4.2f)", integral/preInt);
			}
		} 
		printf(" \\\\ \\hline\n");
	}
}

void hackedbob(bool eff = false)
{
    using namespace std;
    
    double integral, preInt = 1.0;
    
    vector<pair<string,string> > fnames;
    vector<pair<string,string> > hnames;

    //TFile * tfile = new TFile("/home/ugrad/pastika/cms/HeavyNu/CMSSW_4_2_7/src/HeavyNu/AnalysisModules/analysis.root");
    TFile * tfile = new TFile("/local/cms/user/pastika/heavynu/heavynu_2011Bg_summer11_DYToLL_M-50_7TeV-sherpa_nov3.root");
    hnames.push_back(pair<string,string>("none","hNu%s/cut0_none/mWR"));
    hnames.push_back(pair<string,string>("LLJJ Pt","hNu%s/cut1_LLJJpt/mWR"));
    hnames.push_back(pair<string,string>("trig","hNu%s/cut2_TrigMatches/mWR"));
    hnames.push_back(pair<string,string>("vertex","hNu%s/cut3_Vertex/mWR"));
    hnames.push_back(pair<string,string>("mu1 pt","hNu%s/cut4_Mu1HighPt/mWR"));
    hnames.push_back(pair<string,string>("Mll","hNu%s/cut5_diLmass/mWR"));
    hnames.push_back(pair<string,string>("MWR","hNu%s/cut6_mWRmass/mWR"));

    fnames.push_back(pair<string,string>("nom Mu24", "Mu24"));
    fnames.push_back(pair<string,string>("nom Mu40", "Mu40"));
    //fnames.push_back(pair<string,string>("JES Hi Mu24", "Mu24jesHi"));
    //fnames.push_back(pair<string,string>("JES Lo Mu24", "Mu24jesLo"));
    //fnames.push_back(pair<string,string>("JES Hi Mu40", "Mu40jesHi"));
    //fnames.push_back(pair<string,string>("JES Lo Mu40", "Mu40jesLo"));
    //fnames.push_back(pair<string,string>("MES Hi Mu24", "Mu24mesHi"));
    //fnames.push_back(pair<string,string>("MES Lo Mu24", "Mu24mesLo"));
    //fnames.push_back(pair<string,string>("MES Hi Mu40", "Mu40mesHi"));
    //fnames.push_back(pair<string,string>("MES Lo Mu40", "Mu40mesLo"));
    //fnames.push_back(pair<string,string>("MuID Hi Mu24", "Mu24midHi"));
    //fnames.push_back(pair<string,string>("MuID Lo Mu24", "Mu24midLo"));
    //fnames.push_back(pair<string,string>("MuID Hi Mu40", "Mu40midHi"));
    //fnames.push_back(pair<string,string>("MuID Lo Mu40", "Mu40midLo"));
    //fnames.push_back(pair<string,string>("Trig Hi Mu24", "Mu24trigHi"));
    //fnames.push_back(pair<string,string>("Trig Lo Mu24", "Mu24trigLo"));
    //fnames.push_back(pair<string,string>("Trig Hi Mu40", "Mu40trigHi"));
    //fnames.push_back(pair<string,string>("Trig Lo Mu40", "Mu40trigLo"));
    fnames.push_back(pair<string,string>("Pileup Hi Mu24", "Mu24puHi"));
    fnames.push_back(pair<string,string>("Pileup Lo Mu24", "Mu24puLo"));
    fnames.push_back(pair<string,string>("Pileup Hi Mu40", "Mu40puHi"));
    fnames.push_back(pair<string,string>("Pileup Lo Mu40", "Mu40puLo"));
    
	printf("%8s", "cut");
	for(vector<pair<string,string> >::const_iterator i = fnames.begin(); i != fnames.end(); i++)
	{
		printf(eff?" & %15s        ":" & %15s", i->first.c_str());
	}
	printf(" \\\\ \\hline\n");
	for(vector<pair<string,string> >::const_iterator j = hnames.begin(); j != hnames.end(); j++)
	{
		printf("%8s", j->first.c_str());
		for(vector<pair<string,string> >::const_iterator i = fnames.begin(); i != fnames.end(); i++)
		{
		    char histname[256];
		    sprintf(histname, j->second.c_str(), i->second.c_str());
			TH1* h = (TH1*)tfile->Get(histname);
			integral = h->Integral(0, h->GetNbinsX()+1);
			
			if(integral > 10) printf((j == hnames.begin() && eff)?" & %15.0f        ":" & %15.0f", integral);
			else if(integral > 1) printf((j == hnames.begin() && eff)?" & %15.2f        ":" & %15.2f", integral);
			else printf((j == hnames.begin() && eff)?" & %15.3f        ":" & %15.3f", integral);
			//if(j != hnames.begin() && eff) 
			//{
			//	TH1* hpre = (TH1*)tfile->Get((j-1)->second.c_str());
			//	preInt = hpre->Integral(0, h->GetNbinsX()+1);
			//	printf("  (%4.2f)", integral/preInt);
			//}
		} 
		printf(" \\\\ \\hline\n");
	}
}
