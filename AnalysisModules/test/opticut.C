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

#include "CLfast.C"

using namespace std;

static bool gl_verbose=true;

//const double zjet_scale = 1.57;
//const double zjet_scale = 1.46;
const double zjet_scale = 1.15;
//const double qcdp0 =  0.06514;  // constant
//const double qcdp1 = -0.003923; // slope
//const double qcdp0 =  0.08354;  // constant
//const double qcdp1 = -0.009009; // slope
const double qcdp0 =  0.1149;     // constant
const double qcdp1 = -0.008046;   // slope

const double optmwrcuts[] = {
  //760, 800, 880,1000,1080,1160,1240
  520,600,640,720,800,880,960,1040,1120,1200
};
const double binwidthGeV=40.0;

//const int rebinMWRval = 5;
const int rebinMWRval = 2;
const int rebinMLLval = 1;

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
  wth1->DrawFits("same","LL");
  gPad->Update();
  wth1->DrawStats();
  gPad->Update();

  name += ".png";
  c1->SaveAs(name.c_str());
}                                                            // fitndraw

//======================================================================

struct optVars_t {
  optVars_t():
    cutgev(0),ttint(0),zjint(0),wjint(0),vvint(0),twint(0),qcdint(0),
    sumback(0),sumsig(0),sigeff(0),signif(0),clfactor(9e99) {}
  double cutgev,ttint,zjint,wjint,vvint,twint,qcdint,sumback,sumsig,sigeff,signif,clfactor;
};

struct optFits_t {
  optFits_t():tt(0),zj(0),wj(0),vv(0),tw(0),qcd(0) {}
  TF1 *tt,*zj,*wj,*vv,*tw,*qcd;
};

//======================================================================

void loadMWRFits(map<string,wTH1*>&  m_wth1,
		 optFits_t& fits)
{
  wTH1 *ttbar = m_wth1["ttjets_m4"];
  wTH1 *zjets = m_wth1["zjets_m4"];
  wTH1 *wjets = m_wth1["wjets_m4"];
  wTH1 *vv    = m_wth1["vv_m4"];
  wTH1 *tw    = m_wth1["tw_m4"];

  fits.tt = new TF1("ttfit","expo",650,2000);
  fits.zj = new TF1("zjfit","expo",650,2000);
  fits.wj = new TF1("wjfit","expo",560,2000); // note!
  fits.vv = new TF1("vvfit","expo",650,2000);
  fits.tw = new TF1("twfit","expo",650,2000);

  ttbar->loadFitFunction(fits.tt);
  zjets->loadFitFunction(fits.zj);
  wjets->loadFitFunction(fits.wj);
  vv   ->loadFitFunction(fits.vv);
  tw   ->loadFitFunction(fits.tw);

  ttbar->histo()->Rebin(rebinMWRval);
  vv   ->histo()->Rebin(rebinMWRval);
  tw   ->histo()->Rebin(rebinMWRval);

  // Scale zjets according to Z-peak fit:
  zjets->histo()->Scale(zjet_scale);
  zjets->histo()->Rebin(rebinMWRval);

  fits.qcd = new TF1("qcdfit","expo",650,2000);
  fits.qcd->SetParameter(0,qcdp0);
  fits.qcd->SetParameter(1,qcdp1);

  fitndraw(wjets);
  fitndraw(vv);
  fitndraw(tw);
  fitndraw(zjets);
  fitndraw(ttbar);

  printf("%-9s%10s%10s%10s\n", "Sample","p0","p1", "chi2/ndof");
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","TTJets",
	 fits.tt->GetParameter(0),fits.tt->GetParameter(1),
	 fits.tt->GetChisquare(),fits.tt->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","ZJets",
	 fits.zj->GetParameter(0),fits.zj->GetParameter(1),
	 fits.zj->GetChisquare(),fits.zj->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","VV",
	 fits.vv->GetParameter(0),fits.vv->GetParameter(1),
	 fits.vv->GetChisquare(),fits.vv->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","TW",
	 fits.tw->GetParameter(0),fits.tw->GetParameter(1),
	 fits.tw->GetChisquare(),fits.tw->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","WJets",
	 fits.wj->GetParameter(0),fits.wj->GetParameter(1),
	 fits.wj->GetChisquare(),fits.wj->GetNDF());

}                                                         // loadMWRFits

//======================================================================

void loadMLLFits(map<string,wTH1*>&  m_wth1,
		 optFits_t& fits)
{
  wTH1 *ttbar = m_wth1["ttjets_m2"];
  wTH1 *zjets = m_wth1["zjets_m2"];
  wTH1 *wjets = m_wth1["wjets_m2"];
  wTH1 *vv    = m_wth1["vv_m2"];
  wTH1 *tw    = m_wth1["tw_m2"];

  fits.tt = new TF1("ttfit","expo",160,1000);
  fits.zj = new TF1("zjfit","expo",160,2000);
  fits.wj = new TF1("wjfit","expo",120,300); // note!
  fits.vv = new TF1("vvfit","expo",160,1200);
  fits.tw = new TF1("twfit","expo",160,1000);

  ttbar->loadFitFunction(fits.tt);
  zjets->loadFitFunction(fits.zj);
  wjets->loadFitFunction(fits.wj);
  vv   ->loadFitFunction(fits.vv);
  tw   ->loadFitFunction(fits.tw);

  ttbar->histo()->Rebin(rebinMLLval);
  vv   ->histo()->Rebin(rebinMLLval);
  tw   ->histo()->Rebin(rebinMLLval);

  // Scale zjets according to Z-peak fit:
  zjets->histo()->Scale(zjet_scale);
  zjets->histo()->Rebin(rebinMLLval);

  fitndraw(wjets);
  fitndraw(vv);
  fitndraw(tw);
  fitndraw(zjets);
  fitndraw(ttbar);

  printf("%-9s%10s%10s%10s\n", "Sample","p0","p1", "chi2/ndof");
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","TTJets",
	 fits.tt->GetParameter(0),fits.tt->GetParameter(1),
	 fits.tt->GetChisquare(),fits.tt->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","ZJets",
	 fits.zj->GetParameter(0),fits.zj->GetParameter(1),
	 fits.zj->GetChisquare(),fits.zj->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","VV",
	 fits.vv->GetParameter(0),fits.vv->GetParameter(1),
	 fits.vv->GetChisquare(),fits.vv->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","TW",
	 fits.tw->GetParameter(0),fits.tw->GetParameter(1),
	 fits.tw->GetChisquare(),fits.tw->GetNDF());
  printf("%-9s%10.4f%10.4f%10.4g/%d\n","WJets",
	 fits.wj->GetParameter(0),fits.wj->GetParameter(1),
	 fits.wj->GetChisquare(),fits.wj->GetNDF());

}                                                         // loadMLLFits

//======================================================================

void calcVars(TH1 *sighist,const optFits_t& fits,optVars_t& vars,double rebinval)
{
  int lobin    = sighist->FindBin(vars.cutgev);
  int hibin    = sighist->FindBin(1+sighist->GetNbinsX());
  vars.sumsig  = sighist->Integral(lobin,hibin);
  vars.sigeff  = sighist->GetEntries()/10000.;

  vars.ttint   = fits.tt ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  vars.zjint   = fits.zj ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  vars.wjint   = fits.wj ->Integral(vars.cutgev,2000)/binwidthGeV;
  vars.vvint   = fits.vv ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  vars.twint   = fits.tw ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  if (fits.qcd)
    vars.qcdint= fits.qcd->Integral(vars.cutgev,2000)/binwidthGeV;
	
  vars.sumback = vars.ttint+vars.zjint+vars.wjint+vars.vvint+vars.twint+vars.qcdint;

  vars.signif = vars.sumsig/sqrt(vars.sumsig+vars.sumback);

  vars.clfactor = CLfast_goalSeek(vars.sumsig,vars.sumback);
}                                                            // calcVars

//======================================================================

inline
void putvars(const char *fmt, double mwr, double mnu, const optVars_t& v)
{
  printf(fmt,mwr,mnu,v.cutgev,v.ttint,v.zjint,v.wjint,v.vvint,v.twint,
	 v.qcdint,v.sumback,v.sumsig,v.signif,v.clfactor);
}                                                            // putvars

//======================================================================

void scanMWRCutVals(map<string,wTH1*>&  m_wth1,
		    const optFits_t& fits)
{
  for (double mwr=700; mwr<=1600; mwr+=100) {

    cout<<"=================================================================";
    cout<<"================================================================="<<endl;

    printf("%-9s%9s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
"mWR(GeV)","mNU(GeV)","optcutgev","ttbar","zjets","wjets","VV","tW","QCD","Sum BG","Sum Sig.","signif","clfactor");

    cout<<"=================================================================";
    cout<<"================================================================="<<endl;

    for (double mnu=100; mnu<mwr; mnu+=100) {
      optVars_t minclfact;
      optVars_t maxsignif;
      optVars_t vars;

      // get signal histo for this masspoint from the input map
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);

      if (!m_wth1[name]) {
	cerr << "Didn't find " << name << endl;
	continue;
      }

      TH1 *sighist = m_wth1[name]->histo();

      for (vars.cutgev=520; vars.cutgev<mwr; vars.cutgev+=binwidthGeV) {
	calcVars(sighist,fits,vars,rebinMWRval);
	if (vars.signif   > maxsignif.signif)   maxsignif = vars;
	if (vars.clfactor < minclfact.clfactor) minclfact = vars;
      }

      putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)%7.3f\n",
	      mwr,mnu,maxsignif);
      putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)\n\n",
	      mwr,mnu,minclfact);
    } // mnu loop
  } // mwr loop
}                                                      // scanMWRCutVals

//======================================================================
// almost the same as above, but different cut optimization range
//
void scanMLLCutVals(map<string,wTH1*>&  m_wth1,
		    const optFits_t& fits)
{
  for (double mwr=1000; mwr<=1600; mwr+=100) {

    cout<<"=================================================================";
    cout<<"================================================================="<<endl;

    printf("%-9s%9s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
"mWR(GeV)","mNU(GeV)","optcutgev","ttbar","zjets","wjets","VV","tW","QCD","Sum BG","Sum Sig.","signif","clfactor");

    cout<<"=================================================================";
    cout<<"================================================================="<<endl;

    for (double mnu=100; mnu<mwr; mnu+=100) {
      optVars_t minclfact;
      optVars_t maxsignif;
      optVars_t vars;

      // get signal histo for this masspoint from the input map
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);

      if (!m_wth1[name]) {
	cerr << "Didn't find " << name << endl;
	continue;
      }

      TH1 *sighist = m_wth1[name]->histo();

      for (vars.cutgev=120; vars.cutgev<760; vars.cutgev+=binwidthGeV) {
	calcVars(sighist,fits,vars,rebinMLLval);
	if (vars.signif   > maxsignif.signif)   maxsignif = vars;
	if (vars.clfactor < minclfact.clfactor) minclfact = vars;
      }

      putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)%7.3f\n",
	      mwr,mnu,maxsignif);
      putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)\n\n",
	      mwr,mnu,minclfact);
    } // mnu loop
  } // mwr loop
}                                                      // scanMLLCutVals

//======================================================================

void printDanilaTable(map<string,wTH1*>&  m_wth1,
		      const optFits_t& fits)
{
  printf("MWR\tMNu\tL1Pt\tL2Pt\tL1Eta\tL2Eta\tMll\tMwr_cut\t");
  printf("S_nw\tS_eff\tBG_nw\tQCD_nw\tTTb_nw\tZJ_nw\tWJ_nw\tVV_nw\ttW_nw\n");

  for (int i=0; i<10; i++) {
    double mwr = 700. + i*100.;
    optVars_t v;
    v.cutgev = optmwrcuts[i];

    cout<<"==============================================================";
    cout<<"=============================================================="<<endl;

    for (double mnu=100; mnu<mwr; mnu+=100) {
      // get signal histo for this masspoint from the input map
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);

      if (!m_wth1[name]) continue;

      TH1 *sighist   = m_wth1[name]->histo();

      calcVars(sighist,fits,v,rebinMWRval);

      printf("%4.0f\t%4.0f\t60.0\t20.0\t2.1\t2.4\t200\t%4.0f\t",mwr,mnu,v.cutgev);
      printf("%6.3f\t%5.3f\t%6.3f\t%6.4f\t%6.3f\t%6.3f\t%6.4f\t%6.4f\t%6.4f\n",
	     v.sumsig,v.sigeff,v.sumback,v.qcdint,v.ttint,v.zjint,v.wjint,v.vvint,v.twint);

    } // mnu loop
  } // mwr loop
}                                                    // printDanilaTable

//======================================================================

void opticut()
{
  map<string,wTH1 *> m_wth1;
  optFits_t fits;

  setPRDStyle();

  // preload map so I don't have to do finds and inserts
  //
  for (double mwr=700; mwr<=1600; mwr+=1000) {
    for (double mnu=100; mnu<mwr; mnu+=100) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);
      m_wth1[name] = NULL;
    }
  }

#if 0
  getHistosFromRE("optimMLL.root",".*",m_wth1);
  loadMLLFits(m_wth1,fits);
  scanMLLCutVals(m_wth1,fits);
#else
  getHistosFromRE("optimMWR.root",".*",m_wth1);
  loadMWRFits(m_wth1,fits);
  scanMWRCutVals(m_wth1,fits);
  //printDanilaTable(m_wth1,fits);
#endif
}
