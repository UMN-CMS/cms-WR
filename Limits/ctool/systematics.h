#include <map>
#include <string>
#include <vector>

std::map<std::string,double> systematicsInput[4];

std::vector<std::string> systematicsList;



void loadSystematics() {
  static bool didLoad=false;
  if (didLoad) return;
  didLoad=true;
    

  systematicsList.push_back("LUMI");
  systematicsList.push_back("RECOCAL");
  systematicsList.push_back("ISRPDF");
  systematicsList.push_back("SIGONLY");
  systematicsList.push_back("TTONLY");
  systematicsList.push_back("ZJONLY");
  systematicsList.push_back("UNCOTHER");

  systematicsInput[i_SIGNAL]["LUMI"]=1.04;
  systematicsInput[i_SIGNAL]["RECOCAL"]=1+sqrt(pow(0.08,2)+pow(0.02,2)+
					       pow(0.005,2)+pow(0.02,2)+
					       pow(0.02,2));
  systematicsInput[i_SIGNAL]["ISRPDF"]=1+sqrt(pow(0.03,2)+pow(0.04,2)+
					      pow(0.00,2));
  systematicsInput[i_SIGNAL]["SIGONLY"]=1+sqrt(pow(0.00,2)+pow(0.05,2));

  systematicsInput[i_TT]["LUMI"]=1.04;
  systematicsInput[i_TT]["RECOCAL"]=1+sqrt(pow(0.11,2)+pow(0.05,2)+
					   pow(0.005,2)+pow(0.02,2)+
					   pow(0.01,2));
  systematicsInput[i_TT]["ISRPDF"]=1+sqrt(pow(0.08,2)+pow(0.06,2)+
					  pow(0.09,2));
  systematicsInput[i_TT]["TTONLY"]=1+sqrt(pow(0.15,2)+pow(0.02,2));


  systematicsInput[i_ZJ]["RECOCAL"]=1+sqrt(pow(0.04,2)+pow(0.02,2)+
					   pow(0.005,2)+pow(0.02,2)+
					   pow(0.005,2));
  systematicsInput[i_ZJ]["ISRPDF"]=1+sqrt(pow(0.00,2)+pow(0.09,2)+
					  pow(0.14,2));
  systematicsInput[i_ZJ]["ZJONLY"]=1+sqrt(pow(0.09,2)+pow(0.03,2));


  systematicsInput[i_OTHER]["LUMI"]=1.04;
  systematicsInput[i_OTHER]["UNCOTHER"]=1.22;
}

const char* getSyst(int iChan, const std::string& sname) {
  static char ubuff[10];

  std::map<std::string,double>::const_iterator i=systematicsInput[iChan].find(sname);
  if (i==systematicsInput[iChan].end()) return "   - ";
  sprintf(ubuff,"%5.2f",i->second);
  return ubuff;
}
