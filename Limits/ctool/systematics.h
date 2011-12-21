#include <map>
#include <string>
#include <set>

#ifndef systematics_h_included
#define systematics_h_included 1

class SystematicsDB {
 public:
  static const int NO_PDF_FOR_SIGNAL = 0x001;

  SystematicsDB();

  double getSystematic(const std::string& systName, const std::string& process, int imassbin) const;
  std::vector<std::string> getSystematicsList() const;
  void setSimpleSystematic(const std::string& systName);
  void defineSingleChannelSyst(const std::string& systName, const std::string& process, const std::vector<std::string>& contents);
  void defineCommonSyst(const std::string& systName, const std::vector<std::string>& contents);

  void load(const std::string& systdb);
  void standardSystematics();
  void dump() const;

 private:
  std::vector<std::string> m_finalsystematics;
  std::set<std::string> m_processNames;
  struct DBitem {
    double values[10];
  };

  std::map<std::string,DBitem> m_db;

  int m_mode;

};

#endif // systematics_h_included
