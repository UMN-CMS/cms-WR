#include <algorithm>
#include <assert.h>
#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <cstdlib>
#include "TClass.h"
#include "TKey.h"
#include "TRegexp.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TObjString.h"
#include "TString.h"

#include "MyHistoWrapper.cc"

//======================================================================

static bool gl_verbose=true;

inline unsigned int str2int(const std::string& str) {
  return (unsigned int)strtoul(str.c_str(),NULL,0);
}

inline float str2flt(const std::string& str) {
  return (float)strtod(str.c_str(),NULL);
}

inline std::string int2str(int i) {
  std::ostringstream ss;
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

size_t ci_find(const std::string& str1, const std::string& str2)
{
  std::string::const_iterator pos = search(str1.begin(), str1.end(),
				      str2.begin(), str2.end(), ci_equal);
  if (pos == str1.end())
    return std::string::npos;
  else
    return (pos-str1.begin());
}

//======================================================================
// Got this from
// http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html

void Tokenize(const std::string& str,
	      std::vector<std::string>& tokens,
	      const std::string& delimiters = " ",
	      bool include_delimiters=false)
{
  std::string src=str;
  tokens.clear();

  // Skip delimiters at beginning.
  std::string::size_type lastPos = src.find_first_not_of(delimiters, 0);
  if (include_delimiters && lastPos>0)
    tokens.push_back(src.substr(0,lastPos));

  // Find first "non-delimiter".
  std::string::size_type pos = src.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(src.substr(lastPos, pos - lastPos));

    lastPos = src.find_first_not_of(delimiters, pos);

    if (include_delimiters && pos!=std::string::npos) {
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
      std::cerr << "The regexp " << sre->GetString() << " is invalid, Status() = ";
      std::cerr << re.Status() << std::endl;
      exit(-1);
    }

    //TString path( (char*)strstr( dir->GetPath(), ":" ) );
    //path.Remove( 0, 2 ); // gets rid of ":/"

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
      //std::cout << "Checking path: " << ((TDirectory *)obj)->GetPath() << std::endl;
      recurseDirs( (TDirectory *)obj, doFunc, Args, Output );
    } else {
      doFunc(obj, thisdir, Args, Output);
    }
  } // key loop
}                                                         // recurseDirs

//======================================================================

void getHistosFromRE(const std::string&   filepath,
		     const std::string&   sre,
		     std::map<std::string,wTH1*>&  m_wth1)
{
  if (gl_verbose)
    std::cout<<"Searching for regexp "<<sre<<" in "<<filepath;

  // allow for multiple regexes in OR combination
  std::vector<std::string> v_regexes;
  Tokenize(sre,v_regexes,"|");
  if (!v_regexes.size())
    v_regexes.push_back(sre);

  // Check 'em now!
  TObjArray *Args = new TObjArray();
  for (size_t i=0; i<v_regexes.size(); i++) {
    TRegexp re(v_regexes[i].c_str(),kTRUE);
    if (re.Status() != TRegexp::kOK) {
      std::cerr << "The regexp " << v_regexes[i] << " is invalid, Status() = ";
      std::cerr << re.Status() << std::endl;
      exit(-1);
    }
    else {
      Args->AddLast(new TObjString(v_regexes[i].c_str()));
    }
  }

  TFile *rootfile = new TFile(filepath.c_str());

  if (rootfile->IsZombie()) {
    std::cerr << "File failed to open, " << filepath << std::endl;
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
  if (gl_verbose) std::cout << "... " << nx2matches/2 << " match(es) found.";

  for (int i=0; i<nx2matches; i+=2) {
    TString fullspec = ((TObjString *)(*Matches)[i])->GetString();
    wTH1 *wth1 = new wTH1((TH1 *)((*Matches)[i+1]));
    m_wth1[std::string(wth1->histo()->GetName())] = wth1;
    m_wth1[std::string(wth1->histo()->GetName())] = wth1;

    wth1->ShutUpAlready(); // also quiets the fitting function
    wth1->histo()->UseCurrentStyle();
    wth1->histo()->SetTitle(wth1->histo()->GetName());
  }

  //Matches->Delete(); // need the histos!
  delete Matches;

  if (gl_verbose) std::cout << std::endl;
}                                                     // getHistosFromRE
