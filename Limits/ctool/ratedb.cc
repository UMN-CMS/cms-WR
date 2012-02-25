#include "ratedb.hh"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

static std::string makeKey(const std::string& c1, const std::string& c2) {
  std::string retval;
  for (std::string::const_iterator i=c1.begin();i!=c1.end(); i++)
    retval+=toupper(*i);
  retval+="-";
  for (std::string::const_iterator i=c2.begin();i!=c2.end(); i++)
    retval+=toupper(*i);
  return retval;
}

void RateDB::load(const char* file) {
  char buffer[1025];
  FILE* f;
  f=fopen(file,"r");
  if (f==0) return;
  while (!feof(f)) {
    buffer[0]=0;
    fgets(buffer,1024,f);

    std::string proc, tag, work;
    std::vector<double> vals;

    for (int i=0; buffer[i]!=0; i++) {
      if (buffer[i]==',' || buffer[i]==' ' || buffer[i]=='\n') {
	if (work.empty()) continue;

	if (proc.empty()) proc=work;
	else if (tag.empty()) tag=work;
	else vals.push_back(atof(work.c_str()));
	
	work.clear();
      } else work+=buffer[i];
    }
    if (!work.empty()) vals.push_back(atof(work.c_str()));
    if (!vals.empty()) {
//      printf("Loading %s (%d items)\n",makeKey(proc,tag).c_str(),vals.size());
      m_db.insert(std::pair<std::string, std::vector<double> >(makeKey(proc,tag),vals));
    }
  }
  fclose(f);
}

double RateDB::get(const std::string& process, const std::string& tag, int ibin) const {
  std::map<std::string, std::vector<double> >::const_iterator i;
  std::string key=makeKey(process,tag);
  i=m_db.find(key);
  if (i==m_db.end()) {
    printf("search failure '%s'\n",key.c_str());
    return 0;
  }
  if (ibin>=int(i->second.size())) return 0;
  return i->second[ibin];
}
