#ifndef ratedb_hh_included
#define ratedb_hh_included 1

#include <string>
#include <map>
#include <vector>

class RateDB {
public:
  void load(const char* file);
  double get(const std::string& process, const std::string& tag, int ibin) const;
private:
  
  std::map<std::string, std::vector<double> > m_db;
};

#endif // ratedb_hh_included
