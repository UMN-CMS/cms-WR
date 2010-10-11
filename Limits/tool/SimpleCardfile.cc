#include "SimpleCardfile.hh"
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

SimpleCardfile::Entry::Entry(const std::string& t) : title(t) { }

SimpleCardfile::Entry::~Entry() {
  for (unsigned int i=0; i<lines.size(); i++)
    delete (lines[i]);
}

void SimpleCardfile::Entry::addLine(const std::string& raw, const std::vector<std::string>& vals) {
  Line* l=new Line;
  l->index=int(lines.size());
  l->raw=raw;
  l->value=vals;
  lines.push_back(l);
}

SimpleCardfile::SimpleCardfile(const std::string& filename) {
  FILE* f=fopen(filename.c_str(),"r");
  if (f==NULL) return; // broken...
  char buffer[4096];
  Entry* entry=NULL;
  std::map<std::string,int> multiEntry;
  
  std::vector<std::string> tokens;
  while (!feof(f)) {
    buffer[0]=0;
    fgets(buffer,4096,f);
    
    char* c=index(buffer,'#');
    if (c!=NULL) (*c)=0;
    c=index(buffer,'\n');
    if (c!=NULL) (*c)=0;
    
    tokenize(buffer,tokens);
    if (tokens.empty()) continue;

    std::vector<std::string>::iterator q=tokens.end();
    std::string title; // just in case
    
    for (std::vector<std::string>::iterator i=tokens.begin(); i!=tokens.end() && q==tokens.end(); i++) {
      if (*i==":") q=i;
      else {
	if (!title.empty()) title+=" ";
	title+=*i;
      }
    }
    
    if (q!=tokens.end()) {
      entry=new Entry(title);
      if (m_entries.find(title)!=m_entries.end()) { // erase old one...
	Entry* temp=m_entries[title];
	m_entries.erase(m_entries.find(title));
	std::string newTitle(title);
	newTitle+="(0)";
	m_entries[newTitle]=temp;
	multiEntry[title]=1;
      } 
      if (multiEntry.find(title)!=multiEntry.end()) {
	char ext[10];
	sprintf(ext,"(%d)",multiEntry[title]);
	multiEntry[title]++;
	title+=ext;
      }
      m_entries[title]=entry;
      q++;
      tokens.erase(tokens.begin(),q); // drop all these
      //      printf("Title = '%s'\n",title.c_str());
    } 
    if (entry!=NULL) { // && !tokens.empty()) {
      entry->addLine(buffer,tokens);
      /*
      for (unsigned int i=0; i<tokens.size(); i++) 
	printf("'%s' - ",tokens[i].c_str());
      */
    }
    //    printf("\n");
  }

  fclose(f);
}

void SimpleCardfile::tokenize(const std::string& line, std::vector<std::string>& tokens) {
  tokens.clear();
  int quote=0;
  bool escape=false;
  std::string item;

  for (const char* p=line.c_str(); *p!=0; p++) {
    if (escape) {
      switch (*p) {
      case ('n') : item.push_back('\n'); break;
      case ('t') : item.push_back('\t'); break;
      case ('r') : item.push_back('\r'); break;
      case ('v') : item.push_back('\v'); break;
      case ('\\') : item.push_back('\\'); break;
      case ('\'') : item.push_back('\''); break;
      case ('\"') : item.push_back('\"'); break;
      default: item.push_back(*p); break;
      };
      escape=false;
    } else if (*p=='\\') escape=true;
    else if ((*p=='\'' || *p=='"' ) && !quote) quote=*p;
    else if ((*p=='\'' || *p=='"' ) && quote==*p) quote=0;
    else if (!quote && isspace(*p)) {
      if (item.length()>0) tokens.push_back(item);
      item.clear();
    } else item.push_back(*p);
  }
  if (item.length()>0) tokens.push_back(item);
}

SimpleCardfile::~SimpleCardfile() {
  std::map<std::string, Entry*>::iterator i;
  for (i=m_entries.begin(); i!=m_entries.end(); i++)
    delete i->second;
}

bool SimpleCardfile::hasEntry(const std::string& ename) const {
  return m_entries.find(ename)!=m_entries.end();
}

bool SimpleCardfile::requireEntry(const std::string& ename, const std::string& message) const {
  if (hasEntry(ename)) return true;
  else fprintf(stderr,"Missing required parameter: '%s' (%s)\n",ename.c_str(),message.c_str());
  return false;
}

int SimpleCardfile::getLineCount(const std::string& ename) const {
  std::map<std::string, Entry*>::const_iterator i=m_entries.find(ename);
  if (i==m_entries.end()) return -1;
  return i->second->lines.size();
}
int SimpleCardfile::getItemCount(const std::string& ename, int iline) const {
  std::map<std::string, Entry*>::const_iterator i=m_entries.find(ename);
  if (i==m_entries.end() || iline>=int(i->second->lines.size())) return -1;
  return i->second->lines[iline]->value.size();

}

static const std::string NoneSuch("NULL");

const std::string& SimpleCardfile::getItem(const std::string& ename, int iline, int iitem) const {
  std::map<std::string, Entry*>::const_iterator i=m_entries.find(ename);
  if (i==m_entries.end() || iline>=int(i->second->lines.size())) return NoneSuch;
  const Line* l=i->second->lines[iline];
  if (iitem>=int(l->value.size())) return NoneSuch;
  return l->value[iitem];
}

const std::string& SimpleCardfile::getItem(const std::string& ename, int iitem) const {
  return getItem(ename,0,iitem);
}

int SimpleCardfile::getItemInt(const std::string& ename, int iline, int iitem, int defaultval) const {
  const std::string& sval=getItem(ename,iline,iitem);
  if (sval==NoneSuch) return defaultval;
  else return strtol(sval.c_str(),0,0);
}

int SimpleCardfile::getItemInt(const std::string& ename, int item, int defaultval) const {
  return getItemInt(ename,0,item,defaultval);
}

float SimpleCardfile::getItemFloat(const std::string& ename, int iline, int iitem, float defaultval) const {
  const std::string& sval=getItem(ename,iline,iitem);
  if (sval==NoneSuch) return defaultval;
  else return atof(sval.c_str());
}
float SimpleCardfile::getItemFloat(const std::string& ename, int item, float defaultval) const {
  return getItemFloat(ename,0,item,defaultval);
}
