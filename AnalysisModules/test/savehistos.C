#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "hnuUtils.C"

using namespace std;

void savehistos(const string& filewithlist,
		const string& outputfilename,
		const string& openmode="RECREATE")
{
  FILE *fp = fopen(filewithlist.c_str(),"r");
  if (!fp) {
    cerr << "File not found, " << filewithlist << endl;
    return;
  }

  TFile *rootfile = new TFile(outputfilename.c_str(),openmode.c_str());
  if (rootfile->IsZombie()) {
    cerr << "File failed to open, " << outputfilename << endl;
    return;
  }

  string theline,newname;
  vector<string> tokens;
  while (getLine(fp,theline)) {
    if (theline[0] == '#') continue; // "comments are welcome"
    Tokenize(theline,tokens,"\t");

    TH1 *h1     = getHisto(tokens[0]);

    if (!h1) {
      cerr << "Couldn't get histo " << tokens[0] << endl;
      return;
    }
    TH1 *target = h1;

    if( tokens.size() > 1 && tokens[1].length() ) {
      newname = tokens[1];
      cout<<"Writing histo "<<newname<<" to file "<<outputfilename<<endl;
      target = (TH1 *)h1->Clone(newname.c_str());
    } else
      cout<<"Writing histo "<<h1->GetName()<<" to file "<<outputfilename<<endl;

    target->SetDirectory(rootfile);
    rootfile->cd();
    target->Write();
  }
  rootfile->Close();
}
