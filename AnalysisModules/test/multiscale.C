
#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TKey.h"
#include "TArrayD.h"

#ifdef MAIN
typedef unsigned long uint32_t;
#endif

struct FileInfo_t {
  FileInfo_t(TFile *infp,string inpath,float inxs,int innev, float inkf, int innrb,float inwt):
    fp(infp),path(inpath),xsec(inxs),kfact(inkf),nev(innev),nrebin(innrb),weight(inwt) {}
  TFile *fp;
  string path;
  float  xsec;
  float  kfact;
  int    nev;
  int    nrebin;
  float  weight;
};

typedef struct {
  TH1 *p;
  string path;
  string type;
  string descr;
}
HistInfo_t;

//======================================================================

int getFileInfo(const char *filewithpaths,
		vector<FileInfo_t>& v_rootfiles,
		float integluminvpb)
{
  char line[256];

  FILE *pathfp = fopen(filewithpaths, "r");

  if (!pathfp) {
    cerr << "File not found, " << filewithpaths << endl;
    return 0;
  }

  while (!feof(pathfp) && fgets(line, 256, pathfp)) {
    char path[256];
    int nev;
    float xsec,weight,kfactor=1.0;
    int nrebin = 0;

    if (line[0] == '#') continue;

    int nscanned = sscanf(line, "%s %f %d %f %d", path,&xsec,&nev,&kfactor,&nrebin);

    TFile *tfile =  new TFile(path);
    
    if (tfile->IsZombie()) {
      cerr << "File failed to open, " << path << endl;
      return 0;
    }

    if ((nscanned < 3) ||
	(nscanned > 5)   )  {
      cerr << "pathfile requires <pathstring> <xsec> <nevents> [kfactor] [nrebin]\n";
      return 0;
    }
    //else
    //cout << xsec << " " << integluminvpb << " " << nev << " " << nrebin << endl;

    weight = (xsec*integluminvpb*kfactor)/((float)nev);
    cout << "calculated weight for file " << path << " = ";
    cout << "("<<xsec<<"*"<<integluminvpb<<"*"<<kfactor<<")/("<<nev<<")="<<weight<<endl;
    FileInfo_t fileinfo(tfile,path,xsec,nev,kfactor,nrebin,weight);
    v_rootfiles.push_back(fileinfo);
  }
  return 1;
}                                                         // getFileInfo

//======================================================================
// most code borrowed from hadd.C tutorial

void ScaleAll1file( TDirectory *target, FileInfo_t& source, bool writeErrors ) {

  cout << "Target path: " << target->GetPath() << endl;

  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  source.fp->cd( path );

  TDirectory *current_sourcedir = gDirectory;

  // loop over all keys in this directory

  bool newdir = true;

  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;

  while ( (key = (TKey*)nextkey())) {

    // read object from first source file
    source.fp->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> scale it
#if 0
      if (newdir) {
	newdir=false;
	cout << "Scaling histograms: " << endl;
      }
      cout << obj->GetName() << " ";
#endif
      TH1 *h1 = (TH1*)obj;
      TH1 *h2 = NULL;

      if (source.nrebin > 1)
	h1->Rebin(source.nrebin);

      if (writeErrors) {
	h1->Sumw2();
	h2 = (TH1 *)h1->Clone();
      }

      //------------------------------
      h1->Scale(source.weight);
      //------------------------------

      if (writeErrors) {
	// to force weighted summing when we get to it
	// NOTE: I learned the hard way, this bit has to be set AFTER
	//      the scaling.
	// h1->SetBit(TH1::kIsAverage);
	int nbins = h1->GetNbinsX()*h1->GetNbinsY()*h1->GetNbinsZ();
	for (int ibin=1; ibin<=nbins; ibin++)
	  h1->SetBinError(ibin,source.weight*h2->GetBinError(ibin));
	delete h2;
      }

      if (strstr(obj->GetName(),"h1d_caloMet_MetQCD") &&
	  strstr(target->GetPath(),"cutNone")) {
	TArrayD *sumw2 = h1->GetSumw2();
	cout << sumw2->GetSize() << endl;
	for (int i=0; i<sumw2->GetSize(); i++)
	  printf ("%lf ", (*sumw2)[i]);
	cout << endl;
      }

    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory

      newdir = true;
#if 0
      cout << "\n=====> Found subdirectory " << obj->GetName();
      cout << "<=====\n" << endl;
#endif
      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      ScaleAll1file( newdir, source, writeErrors );

    } else {

      // object is of no type that we know or can handle
      cout << "\n======> Unknown object type, name: " 
           << obj->GetName() << " title: " << obj->GetTitle();
      cout << "<======\n" << endl;
    }

    // now write the scaledd histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();

      obj->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  //cout << endl;

  // save modifications to target file
  target->Write();

}

//======================================================================
// Example: multiscale("files2scale.txt",36.145,"NLO")
//
void multiscale(const char* filewithpaths,
		float integluminvpb,
		const char* csordersuffix="",
		bool writeErrors=false)
{
  vector<FileInfo_t> v_rootfiles;

  if (!getFileInfo(filewithpaths, v_rootfiles,integluminvpb))
    return;

  int time2sleep = (int)max(5.0,0.5*v_rootfiles.size());

  for (int i=time2sleep; i>0; i--) {
    if (!(i%30) || i==15 || i<5) cout<<i<< "..."<<flush;
    sleep(1);
  }
  cout << endl;
    
  // EXTRACT HISTOGRAMS and SCALE EACH
  for (unsigned int i=0; i<v_rootfiles.size(); i++) {

    char scalename[256];
    FileInfo_t file = v_rootfiles[i];
    char *fn = strrchr(file.path.c_str(),'/');
    if (!fn) fn = (char *)file.path.c_str();
    else fn++;
    char *ptr = strstr(fn, ".root");
    *ptr = 0;
    sprintf (scalename, "%s_%dipb%s.root", fn, (int)(integluminvpb+0.5),csordersuffix);
    cout << "Writing to " << scalename << endl;
    TFile *scaledfp = new TFile(scalename,"RECREATE");
    ScaleAll1file(scaledfp, file, writeErrors);
    scaledfp->Write();
    delete scaledfp;
    delete file.fp;
  }
}

