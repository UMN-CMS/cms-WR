#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "TKey.h"
#include "TRegexp.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TString.h"
#include "MyHistoWrapper.cc"
#include "TROOT.h"
#include "TStyle.h"
#include "mystyle.C"

using namespace std;

static bool gl_verbose=true;

const double zjet_scale = 1.57;
//const double zjet_scale = 1.46;
//const double qcdp0 =  0.06514;  // constant
//const double qcdp1 = -0.003923; // slope
const double qcdp0 =  0.08354;    // constant
const double qcdp1 = -0.009009;   // slope

const double optmwrcuts[] = {
  760, 800, 880,1000,1080,1160,1240
};
const double binwidthGeV=40.0;

//======================================================================

inline unsigned int str2int(const string& str) {
  return (unsigned int)strtoul(str.c_str(),NULL,0);
}

inline float str2flt(const string& str) {
  return (float)strtod(str.c_str(),NULL);
}

inline string int2str(int i) {
  ostringstream ss;
  ss << i;
  return ss.str();
}

//======================================================================
// Got this from
// http://www.velocityreviews.com/forums/t286357-case-insensitive-stringfind.html
//
bool ci_equal(char ch1, char ch2)
{
  return (toupper((unsigned char)ch1) ==
          toupper((unsigned char)ch2));
}

size_t ci_find(const string& str1, const string& str2)
{
  string::const_iterator pos = search(str1.begin(), str1.end(),
				      str2.begin(), str2.end(), ci_equal);
  if (pos == str1.end())
    return string::npos;
  else
    return (pos-str1.begin());
}

//======================================================================
// Got this from
// http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html

void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters = " ",
	      bool include_delimiters=false)
{
  string src=str;
  tokens.clear();

  // Skip delimiters at beginning.
  string::size_type lastPos = src.find_first_not_of(delimiters, 0);
  if (include_delimiters && lastPos>0)
    tokens.push_back(src.substr(0,lastPos));

  // Find first "non-delimiter".
  string::size_type pos = src.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(src.substr(lastPos, pos - lastPos));

    lastPos = src.find_first_not_of(delimiters, pos);

    if (include_delimiters && pos!=string::npos) {
      tokens.push_back(src.substr(pos, lastPos-pos));
    } //else skip delimiters.

    // Find next delimiter
    pos = src.find_first_of(delimiters, lastPos);
  }
}                                                            // Tokenize

//======================================================================
// Regex match a histo name in a directory
//
void regexMatchHisto( TObject    *obj,
		      TDirectory *dir,
		      TObjArray  *Args,   // list of regexes to match
		      TObjArray  *Matches)
{
  for (int i=0; i<Args->GetEntriesFast(); i++) {
    TObjString *sre = (TObjString *)(*Args)[i];
    TRegexp re(sre->GetString(),kFALSE);
    if (re.Status() != TRegexp::kOK) {
      cerr << "The regexp " << sre->GetString() << " is invalid, Status() = ";
      cerr << re.Status() << endl;
      exit(-1);
    }

    TString path( (char*)strstr( dir->GetPath(), ":" ) );
    path.Remove( 0, 2 ); // gets rid of ":/"

    TString fullspec = TString(dir->GetPath()) + "/" + obj->GetName();

    if ((fullspec.Index(re) != kNPOS) &&
	(obj->InheritsFrom("TH1"))) {
      // success, record that you read it in.
      TObjString *path = new TObjString(fullspec);
      Matches->AddLast(path);
      Matches->AddLast(obj);

      break; // don't let the object match more than one regex
    } // if we have a match
  } // Arg loop
}                                                     // regexMatchHisto

//======================================================================

void recurseDirs( TDirectory *thisdir,
		  void (*doFunc)(TObject *, TDirectory *,TObjArray *, TObjArray *),
		  TObjArray *Args,
		  TObjArray *Output)
{
  assert(doFunc);

  //thisdir->cd();

  // loop over all keys in this directory

  TIter nextkey( thisdir->GetListOfKeys() );
  TKey *key;
  while ( (key = (TKey*)nextkey())) {

    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory, recurse
      //cout << "Checking path: " << ((TDirectory *)obj)->GetPath() << endl;
      recurseDirs( (TDirectory *)obj, doFunc, Args, Output );
    } else {
      doFunc(obj, thisdir, Args, Output);
    }
  } // key loop
}                                                         // recurseDirs

//======================================================================

void getHistosFromRE(const string&   filepath,
		     const string&   sre,
		     map<string,wTH1*>&  m_wth1)
{
  if (gl_verbose)
    cout<<"Searching for regexp "<<sre<<" in "<<filepath;

  // allow for multiple regexes in OR combination
  vector<string> v_regexes;
  Tokenize(sre,v_regexes,"|");
  if (!v_regexes.size())
    v_regexes.push_back(sre);

  // Check 'em now!
  TObjArray *Args = new TObjArray();
  for (size_t i=0; i<v_regexes.size(); i++) {
    TRegexp re(v_regexes[i].c_str(),kTRUE);
    if (re.Status() != TRegexp::kOK) {
      cerr << "The regexp " << v_regexes[i] << " is invalid, Status() = ";
      cerr << re.Status() << endl;
      exit(-1);
    }
    else {
      Args->AddLast(new TObjString(v_regexes[i].c_str()));
    }
  }

  TFile *rootfile = new TFile(filepath.c_str());

  if (rootfile->IsZombie()) {
    cerr << "File failed to open, " << filepath << endl;
    Args->Delete();
    delete Args;
    return;
  }

  TObjArray *Matches = new TObjArray();
  recurseDirs(rootfile, &regexMatchHisto, Args, Matches);
  Args->Delete();
  delete Args;

  // Returns two objects per match: 
  // 1. the (string) path that was matched and
  // 2. the object whose path matched
  //
  int nx2matches = Matches->GetEntriesFast();
  if (gl_verbose) cout << "... " << nx2matches/2 << " match(es) found.";

  for (int i=0; i<nx2matches; i+=2) {
    TString fullspec = ((TObjString *)(*Matches)[i])->GetString();
    wTH1 *wth1 = new wTH1((TH1 *)((*Matches)[i+1]));
    m_wth1[string(wth1->histo()->GetName())] = wth1;
    m_wth1[string(wth1->histo()->GetName())] = wth1;

    wth1->ShutUpAlready(); // also quiets the fitting function
    wth1->histo()->UseCurrentStyle();
    wth1->histo()->SetTitle(wth1->histo()->GetName());
  }

  //Matches->Delete(); // need the histos!
  delete Matches;

  if (gl_verbose) cout << endl;
}                                                     // getHistosFromRE

//======================================================================

void fitndraw(wTH1 *wth1)
{
  gROOT->SetStyle("Plain");
  string name = wth1->histo()->GetName();
  TCanvas *c1=new TCanvas(name.c_str(),
			  name.c_str(),
			  500,400);
  c1->cd();

  gStyle->SetOptStat(110010);
  gStyle->SetOptFit(1);
  gPad->SetLogy(1);
  wth1->SetStats(true,0,.6,.7,.99,.99);
  wth1->Draw("HIST");
  gPad->Update();
  wth1->DrawFits("same");
  gPad->Update();
  wth1->DrawStats();
  gPad->Update();

  name += ".png";
  //c1->SaveAs(name.c_str());
}                                                            // fitndraw

//======================================================================

void scanCutVals(map<string,wTH1*>&  m_wth1,
		 TF1 *ttfit,
		 TF1 *zjfit,
		 TF1 *wjfit,
		 TF1 *vvfit,
		 TF1 *qcdfit)
{
  printf("%-10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
"mWR(GeV)","mNU(GeV)","optcutgev","ttbar","zjets","wjets","VV","QCD","Sum BG","Sum Sig.","maxsignif");

  for (double mwr=1000; mwr<=1600; mwr+=100) {
    cout<<"=========================================================";
    cout<<"========================================================="<<endl;
    for (double mnu=100; mnu<mwr; mnu+=100) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);

      if (!m_wth1[name]) continue;

      TH1 *sighist = m_wth1[name]->histo();

      double maxsignif=0;
      double optcutgev=0;
      double optttint   = 0;
      double optzjint   = 0;
      double optwjint   = 0;
      double optvvint   = 0;
      double optqcdint  = 0;
      double optsumback = 0;
      double optsumsig  = 0;

      for (double cutgev=600; cutgev<mwr; cutgev+=binwidthGeV) {
	double ttint   =  ttfit->Integral(cutgev,2000)/binwidthGeV;
	double zjint   =  zjfit->Integral(cutgev,2000)/(binwidthGeV*2);
	double wjint   =  wjfit->Integral(cutgev,2000)/binwidthGeV;
	double vvint   =  vvfit->Integral(cutgev,2000)/binwidthGeV;
	double qcdint  = qcdfit->Integral(cutgev,2000)/binwidthGeV;
	
	double sumback = ttint+zjint+wjint+vvint+qcdint;

	int lobin = sighist->FindBin(cutgev);
	int hibin = sighist->FindBin(2000);
      
	double sumsig = sighist->Integral(lobin,hibin-1);
	double signif = sumsig/sqrt(sumsig+sumback);

	if (signif > maxsignif) {
	  maxsignif = signif;
	  optcutgev  = cutgev;
	  optttint  = ttint;
	  optzjint  = zjint;
	  optwjint  = wjint;
	  optvvint  = vvint;
	  optqcdint = qcdint;
	  optsumback= sumback;
	  optsumsig = sumsig;
	}
      } // cutgev loop

      printf("%-10.0f%10.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n",
mwr,mnu,optcutgev,optttint,optzjint,optwjint,optvvint,optqcdint,optsumback,optsumsig,maxsignif);

    } // mnu loop
  } // mwr loop
}                                                         // scanCutVals

//======================================================================

void printDanilaTable(map<string,wTH1*>&  m_wth1,
		      TF1 *ttfit,
		      TF1 *zjfit,
		      TF1 *wjfit,
		      TF1 *vvfit,
		      TF1 *qcdfit)
{
  printf("MWR\tMNu\tL1Pt\tL2Pt\tL1Eta\tL2Eta\tMll\tMwr_cut\t");
  printf("S_nw\tS_eff\tBG_nw\tQCD_nw\tTTb_nw\tZJ_nw\tWJ_nw\tVV_nw\n");

  for (int i=0; i<7; i++) {
    double mwr = 1000. + i*100.;
    double optmwrcut = optmwrcuts[i];
    cout<<"==============================================================";
    cout<<"=============================================================="<<endl;

    for (double mnu=100; mnu<mwr; mnu+=100) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);

      if (!m_wth1[name]) continue;

      TH1 *sighist   = m_wth1[name]->histo();

      double sigeff  = sighist->GetEntries()/10000.;

      double ttint   =  ttfit->Integral(optmwrcut,2000)/binwidthGeV;
      double zjint   =  zjfit->Integral(optmwrcut,2000)/(binwidthGeV*2);
      double wjint   =  wjfit->Integral(optmwrcut,2000)/binwidthGeV;
      double vvint   =  vvfit->Integral(optmwrcut,2000)/binwidthGeV;
      double qcdint  = qcdfit->Integral(optmwrcut,2000)/binwidthGeV;
	
      double sumback = ttint+zjint+wjint+vvint+qcdint;

      int lobin = sighist->FindBin(optmwrcut);
      int hibin = sighist->FindBin(1+sighist->GetNbinsX());
      
      double sumsig = sighist->Integral(lobin,hibin);

      printf("%4.0f\t%4.0f\t60.0\t20.0\t2.1\t2.5\t200\t%4.0f\t",mwr,mnu,optmwrcut);
      printf("%6.3f\t%5.3f\t%6.3f\t%6.4f\t%6.3f\t%6.3f\t%6.4f\t%6.4f\n",
	     sumsig,sigeff,sumback,qcdint,ttint,zjint,wjint,vvint);

    } // mnu loop
  } // mwr loop
}                                                // printDanilaTable

//======================================================================

void opticut()
{
  map<string,wTH1 *> m_wth1;

  setPRDStyle();

  // preload map so I don't have to do finds and inserts
  //
  for (double mwr=1000; mwr<=1600; mwr+=1000) {
    for (double mnu=100; mnu<mwr; mnu+=100) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);
      m_wth1[name] = NULL;
    }
  }

  getHistosFromRE("optim.root",".*",m_wth1);

  wTH1 *ttbar = m_wth1["ttbar_m4"];
  wTH1 *zjets = m_wth1["zjets_m4"];
  wTH1 *wjets = m_wth1["wjets_m4"];
  wTH1 *vv    = m_wth1["vv_m4"];

  TF1 *ttfit = new TF1("ttfit","expo",650,2000);
  TF1 *zjfit = new TF1("zjfit","expo",650,2000);
  TF1 *wjfit = new TF1("wjfit","expo",650,2000);
  TF1 *vvfit = new TF1("vvfit","expo",650,2000);
  ttbar->loadFitFunction(ttfit);
  zjets->loadFitFunction(zjfit);
  wjets->loadFitFunction(wjfit);
  vv   ->loadFitFunction(vvfit);

  // Scale zjets according to Z-peak fit:
  //zjets->histo()->Scale(zjet_scale);
  zjets->histo()->Rebin(2);

  TF1 *qcdfit = new TF1("qcdfit","expo",650,2000);
  qcdfit->SetParameter(0,qcdp0);
  qcdfit->SetParameter(1,qcdp1);

  fitndraw(wjets);
  fitndraw(vv);
  fitndraw(zjets);
  fitndraw(ttbar);
  //scanCutVals(m_wth1,ttfit,zjfit,wjfit,vvfit,qcdfit);

  printDanilaTable(m_wth1,ttfit,zjfit,wjfit,vvfit,qcdfit);
}
