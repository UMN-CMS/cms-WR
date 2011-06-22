#include <map>
#include <string>
#include <vector>

std::map<std::string,double> systematics10[4];
std::map<std::string,double> systematics11[4];

std::vector<std::string> systematicsList;

void loadSystematics() {
  static bool didLoad=false;
  if (didLoad) return;
  didLoad=true;
    

  systematicsList.push_back("LUMI");
  systematicsList.push_back("RECOCAL10");
  systematicsList.push_back("RECOCAL11");
  systematicsList.push_back("ISRPDF");
  systematicsList.push_back("SIGONLY");
  systematicsList.push_back("TTONLY");
  systematicsList.push_back("ZJONLY");
  systematicsList.push_back("UNCOTHER");


  // first we do common systematics

  systematics10[i_SIGNAL]["LUMI"]=1.04;
  systematics10[i_SIGNAL]["ISRPDF"]=1+sqrt(pow(0.03,2)+pow(0.04,2)+
					      pow(0.00,2));
  systematics10[i_SIGNAL]["SIGONLY"]=1+sqrt(pow(0.00,2)+pow(0.05,2));

  systematics10[i_TT]["ISRPDF"]=1+sqrt(pow(0.08,2)+pow(0.06,2)+
					  pow(0.08,2));
  systematics10[i_TT]["TTONLY"]=1+sqrt(pow(0.22,2)+pow(0.05,2));


  systematics10[i_ZJ]["ISRPDF"]=1+sqrt(pow(0.08,2)+pow(0.09,2)+
					  pow(0.14,2));
  systematics10[i_ZJ]["ZJONLY"]=1+sqrt(pow(0.13,2)+pow(0.03,2));


  systematics10[i_OTHER]["LUMI"]=1.04;
  systematics10[i_OTHER]["ISRPDF"]=1+sqrt(pow(0.08,2)+pow(0.05,2)+
					  pow(0.08,2));
  systematics10[i_OTHER]["UNCOTHER"]=1.22;

  // copy 2010 -> 2011
  systematics11[i_SIGNAL]=systematics10[i_SIGNAL];
  systematics11[i_TT]=systematics10[i_TT];
  systematics11[i_ZJ]=systematics10[i_ZJ];
  systematics11[i_OTHER]=systematics10[i_OTHER];

  systematics11[i_SIGNAL]["LUMI"]=1.06; 
  systematics11[i_OTHER]["LUMI"]=1.06; 

  systematics11[i_TT]["TTONLY"]=1+sqrt(pow(0.16,2)+pow(0.05,2));
  systematics11[i_ZJ]["ZJONLY"]=1+sqrt(pow(0.06,2)+pow(0.03,2));

  // individual systematics

  systematics10[i_SIGNAL]["RECOCAL10"]=1+sqrt(pow(0.08,2)+pow(0.02,2)+
					       pow(0.005,2)+pow(0.02,2)+
					       pow(0.02,2));
  systematics10[i_TT]["RECOCAL10"]=1+sqrt(pow(0.11,2)+pow(0.05,2)+
					   pow(0.005,2)+pow(0.02,2)+
					   pow(0.01,2));
  systematics10[i_ZJ]["RECOCAL10"]=1+sqrt(pow(0.04,2)+pow(0.02,2)+
					   pow(0.005,2)+pow(0.02,2)+
					   pow(0.005,2));
  systematics10[i_OTHER]["RECOCAL10"]=1+sqrt(pow(0.11,2)+pow(0.04,2)+
					   pow(0.005,2)+pow(0.02,2)+
					   pow(0.005,2));


  systematics11[i_SIGNAL]["RECOCAL11"]=1+sqrt(pow(0.10,2)+pow(0.02,2)+
					       pow(0.002,2)+pow(0.06,2)+
					       pow(0.002,2));
  systematics11[i_TT]["RECOCAL11"]=1+sqrt(pow(0.11,2)+pow(0.05,2)+
					   pow(0.002,2)+pow(0.008,2)+
					   pow(0.005,2));
  systematics11[i_ZJ]["RECOCAL11"]=1+sqrt(pow(0.03,2)+pow(0.03,2)+
					   pow(0.002,2)+pow(0.003,2)+
					   pow(0.005,2));
  systematics11[i_OTHER]["RECOCAL11"]=1+sqrt(pow(0.10,2)+pow(0.04,2)+
					   pow(0.002,2)+pow(0.004,2)+
					   pow(0.002,2));


}

const char* getSyst(int iChan, int iyear, const std::string& sname) {
  static char ubuff[10];

  std::map<std::string,double>::const_iterator i;

  if (iyear==2010) {
    i=systematics10[iChan].find(sname);
    if (i==systematics10[iChan].end()) return "   - ";
  } else {
    i=systematics11[iChan].find(sname);
    if (i==systematics11[iChan].end()) return "   - ";
  }
  sprintf(ubuff,"%5.2f",i->second);
  return ubuff;
}
