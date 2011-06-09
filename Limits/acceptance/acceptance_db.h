#ifndef acceptance_db_included
#define acceptance_db_included

#include <vector>

class AcceptanceDB {
 public:
  AcceptanceDB() { loadDB(); }
  void add(int mw, int mn, double a2010, double a2011);
  double getBestEstimate(int mw, int mn, int year) const;

 private:
  struct AcceptPt {
    int mw, mn;
    double forYear(int year) const {
      if (year==2010) return a2010;
      if (year==2011) return a2011;
      return 0;
    }
    double a2010, a2011;
  };


  std::vector<AcceptPt> m_DB;
  double interpol2d(int mw, int mn, 
		    int year,
		    const AcceptPt& a,const AcceptPt& b,
		    const AcceptPt& c,const AcceptPt& d) const;
  double dist(const AcceptPt& a, const AcceptPt& b) const;
  void loadDB();
};

#endif // acceptance_db_included
