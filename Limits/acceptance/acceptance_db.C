#include "acceptance_db.h"

void AcceptanceDB::add(int mw, int mn, double a2010, double a2011) {
  AcceptPt pt;
  pt.mw=mw;
  pt.mn=mn;
  pt.a2010=a2010;
  pt.a2011=a2011;
  m_DB.push_back(pt);
}

static double interpol1d(double x, double x1, double y1, double x2, double y2) {
  return (x-x1)*(y2-y1)/(x2-x1)+y1;
}


double AcceptanceDB::interpol2d(int mw, int mn, 
			 int year,
			 const AcceptPt& a,const AcceptPt& b,
			 const AcceptPt& c,const AcceptPt& d) const {
  double v1=interpol1d(mn,a.mn,a.forYear(year),b.mn,b.forYear(year));
  double v2=interpol1d(mn,c.mn,c.forYear(year),d.mn,d.forYear(year));
  /*
  printf("%d %d %f %d %d %f %d %d %f  %d %d %f -- %f %f\n",
	 a.mw,a.mn,a.forYear(year),
	 b.mw,b.mn,b.forYear(year),
	 c.mw,c.mn,c.forYear(year),
	 d.mw,d.mn,d.forYear(year),
	 v1,v2);
  */
  return interpol1d(mw,a.mw,v1,c.mw,v2);
}

double AcceptanceDB::dist(const AcceptPt& a, const AcceptPt& b) const {
  return sqrt(pow(a.mw-b.mw,2)+pow(a.mn-b.mn,2));
}

double AcceptanceDB::getBestEstimate(int mw, int mn, int year) const {
  AcceptPt tool;
  tool.mw=mw;
  tool.mn=mn;

  std::vector<AcceptPt>::const_iterator match(m_DB.end());

  // same mw, closest
  std::vector<AcceptPt>::const_iterator swl(m_DB.end());
  std::vector<AcceptPt>::const_iterator swg(m_DB.end());

  // same mn, closest
  std::vector<AcceptPt>::const_iterator snl(m_DB.end());
  std::vector<AcceptPt>::const_iterator sng(m_DB.end());

  std::vector<AcceptPt>::const_iterator fpll(m_DB.end());
  std::vector<AcceptPt>::const_iterator fplg(m_DB.end());
  std::vector<AcceptPt>::const_iterator fpgl(m_DB.end());
  std::vector<AcceptPt>::const_iterator fpgg(m_DB.end());

  for (std::vector<AcceptPt>::const_iterator i=m_DB.begin(); i!=m_DB.end(); i++) {
    if (i->mw==mw && i->mn==mn) {
      match=i;
      break;
    }
    if (i->mw==mw) {
      if (i->mn>mn) {
	if (swg==m_DB.end() || abs(swg->mn-mn)>abs(i->mn-mn)) swg=i;
      } else {
	if (swl==m_DB.end() || abs(swl->mn-mn)>abs(i->mn-mn)) swl=i;
      }
    } else if (i->mn==mn) {
      if (i->mw>mw) {
	if (sng==m_DB.end() || abs(sng->mw-mw)>abs(i->mw-mw)) sng=i;
      } else {
	if (snl==m_DB.end() || abs(snl->mw-mw)>abs(i->mw-mw)) snl=i;
      }
    } else if (i->mn>mn) {
      double disti=dist(*i,tool);
      if (i->mw>mw) {
	if (fpgg==m_DB.end() || dist(*fpgg,tool)>disti) fpgg=i;
      } else { // i->mw<mw
	if (fplg==m_DB.end() || dist(*fplg,tool)>disti) fplg=i;
      }
    } else { // if (i->mn<mn) {
      double disti=dist(*i,tool);
      if (i->mw>mw) {
	if (fpgl==m_DB.end() || dist(*fpgl,tool)>disti) fpgl=i;
      } else { // i->mw<mw
	if (fpll==m_DB.end() || dist(*fpll,tool)>disti) fpll=i;
      }
    }
  }

  if (match!=m_DB.end()) return match->forYear(year);
  if (swg!=m_DB.end() && swl!=m_DB.end()) {
    return interpol1d(mn*1.0,swl->mn*1.0,swl->forYear(year),swg->mn*1.0,swg->forYear(year));
  }
  if (sng!=m_DB.end() && snl!=m_DB.end()) {
    return interpol1d(mw*1.0,snl->mw*1.0,snl->forYear(year),sng->mw*1.0,sng->forYear(year));
  }
  if (fpgl!=m_DB.end() && fpll!=m_DB.end() &&
      fpgg!=m_DB.end() && fplg!=m_DB.end()) {
    return interpol2d(mw,mn,year,*fpll,*fplg,*fpgl,*fpgg);
  }
  return -1;
}

void AcceptanceDB::loadDB() {
#include "db_jun6.C"
}
