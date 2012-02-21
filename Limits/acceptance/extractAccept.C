#include "TFile.h"
#include "TH1.h"

double extractAccept(TFile* f,int year) {
  TH1* steps=0;
  if (year==2010) {
    steps=(TH1*)f->Get("hNuGen2010/cutProgress");
  } 
  if (year==2011) {
    steps=(TH1*)f->Get("hNuGen2011/cutProgress");
  } 
  if (steps==0) return 0.0;
  return steps->GetBinContent(7)/steps->GetBinContent(1);
}

#include <boost/regex.hpp> 
#include <string>
#include <iostream>

int main(int argc, char* argv[]) {
  TFile af(argv[1]);
  int mw=0, mn=0;
  
  std::string a1(argv[1]);
  boost::match_results<std::string::const_iterator> what; 
  boost::regex exp("_([0-9]+)_([0-9]+)");
  std::string::const_iterator start, end; 
  start = a1.begin(); 
  end = a1.end(); 

  if (regex_search(start,end,what,
		   exp,boost::match_default)) {
    std::string smw(what[1]);
    mw=atoi(smw.c_str());
    mn=atoi(std::string(what[2]).c_str());
  }
		      



  double a2010=extractAccept(&af,2010);
  double a2011=extractAccept(&af,2011);


  printf("  add(%d,%d,%f,%f);\n",mw,mn,a2010,a2011);
  return 0;
}

/*
g++ -o extractAccept.exe extractAccept.C `root-config --cflags --libs --ldflags` -lboost_regex
*/
