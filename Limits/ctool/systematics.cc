#include <map>
#include <string>
#include <vector>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "systematics.h"
#include <algorithm>

static std::string to_upper(const std::string& item) {
  std::string retval;
  for (std::string::const_iterator i=item.begin(); i!=item.end(); i++) 
    retval.push_back(toupper(*i));
  return retval;
}

SystematicsDB::SystematicsDB() : m_mode(0) {
}

static std::string makeKey(const std::string& process, const std::string& systName) {
  return to_upper(process+"-"+systName);
}

void SystematicsDB::load(const std::string& systdb) {
  FILE* f=fopen(systdb.c_str(),"r");
  if (f==0) return;

  char buffer[1024];
  float vals[10];
  char proc[100],syst[100];
  while (!feof(f)) {
    buffer[0]=0;
    fgets(buffer,1000,f);
    if (strchr(buffer,'#')!=0) *(strchr(buffer,'#'))=0;
    while (strchr(buffer,',')!=0) *(strchr(buffer,','))='\t';

    int matched=sscanf(buffer,"%s %s %f %f %f %f %f %f %f %f %f %f",proc,syst,
		       vals,vals+1,vals+2,vals+3,vals+4,
		       vals+5,vals+6,vals+7,vals+8,vals+9);
    if (matched==12) {
      std::string key=makeKey(proc,syst);
      DBitem item;
      for (int i=0; i<10; i++) item.values[i]=vals[i];
      m_db.insert(std::pair<std::string,DBitem>(key,item));
      m_processNames.insert(to_upper(proc));
    }
  }
  fclose(f);
}

void SystematicsDB::dump() const {
  for (std::map<std::string,DBitem>::const_iterator i=m_db.begin(); i!=m_db.end(); i++) {
    printf("'%s',%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",i->first.c_str(),
	   i->second.values[0],i->second.values[1],
	   i->second.values[2],i->second.values[3],
	   i->second.values[4],i->second.values[5],
	   i->second.values[6],i->second.values[7],
	   i->second.values[8],i->second.values[9]);
  }
}

double SystematicsDB::getSystematic(const std::string& systName, const std::string& process, int imassbin) const { 
  std::string key=makeKey(process,systName);

  // No PDF error for signal?
  if ((m_mode&NO_PDF_FOR_SIGNAL)!=0 && key.find("PDF")!=std::string::npos && key.find("SIGNAL")!=std::string::npos) return 0;

  std::map<std::string,DBitem>::const_iterator i=m_db.find(key);
  if (i!=m_db.end() && imassbin>=0 && imassbin<10) 
    return (i->second.values[imassbin]<1)?(i->second.values[imassbin]+1):(i->second.values[imassbin]);
  else return 0.0;
}

std::vector<std::string> SystematicsDB::getSystematicsList() const {
  return m_finalsystematics;
}
void SystematicsDB::setSimpleSystematic(const std::string& systName)  {
  m_finalsystematics.push_back(systName);
}
void SystematicsDB::defineSingleChannelSyst(const std::string& systName, const std::string& process, const std::vector<std::string>& contents) {
  std::string key=makeKey(process,systName);
  DBitem dbi;
  // clear
  for (int i=0; i<10; i++) dbi.values[i]=0;
  // iterate over inputs and sum in quadrature
  for (int i=0; i<10; i++) {
    for (std::vector<std::string>::const_iterator sourcesyst=contents.begin(); sourcesyst!=contents.end(); sourcesyst++) {
      double asyst=getSystematic(*sourcesyst,process,i);
      if (asyst<1) continue; // no contribution
      asyst-=1;
      dbi.values[i]+=pow(asyst,2);
    }
    dbi.values[i]=sqrt(dbi.values[i]);
  }
  
  m_db.insert(std::pair<std::string,DBitem>(key,dbi));
  if (std::find(m_finalsystematics.begin(),m_finalsystematics.end(),systName)==m_finalsystematics.end())
    m_finalsystematics.push_back(systName);
}
void SystematicsDB::defineCommonSyst(const std::string& systName, const std::vector<std::string>& contents) {
  for (std::set<std::string>::const_iterator i=m_processNames.begin(); i!=m_processNames.end(); i++) 
    defineSingleChannelSyst(systName,*i,contents);
}

void SystematicsDB::standardSystematics() {
  std::vector<std::string> systContents;
  systContents.push_back("SHAPE");
  systContents.push_back("NORM");

  defineSingleChannelSyst("TTONLY","TTJETS",systContents);
  defineSingleChannelSyst("ZJONLY","ZJETS",systContents);
  defineSingleChannelSyst("OTHERONLY","OTHER",systContents);

  systContents.clear();
  systContents.push_back("JES");
  systContents.push_back("MES");
  systContents.push_back("MUONID");
  systContents.push_back("PU");
  systContents.push_back("TRIG");
  defineCommonSyst("RECOID",systContents);

  systContents.clear();
  systContents.push_back("PDF");
  systContents.push_back("REN");
  systContents.push_back("FACT");
  defineCommonSyst("PDFSCALE",systContents);
}
