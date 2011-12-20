#include <map>
#include <string>

#ifndef systematics_h_included
#define systematics_h_included 1

class SystematicsDB {
 public:
  static const int NO_PDF_FOR_SIGNAL = 0x001;

  SystematicsDB();

  double getSystematic(const std::string& systName, const std::string& process, int imassbin) const;
  std::vector<std::string> getSystematicsList() const;
  void load(const std::string& systdb);
  void dump() const;

 private:

  struct DBitem {
    double values[10];
  };

  std::map<std::string,DBitem> m_db;

  int m_mode;

};

#endif // systematics_h_included
