#ifndef SimpleCardfile_hh_included
#define SimpleCardfile_hh_included 1

#include <string>
#include <vector>
#include <map>

class SimpleCardfile {
public:
  SimpleCardfile(const std::string& textFile);  
  ~SimpleCardfile();
  
  bool hasEntry(const std::string& ename) const;
  bool requireEntry(const std::string& ename, const std::string& message) const;
  int getLineCount(const std::string& ename) const;
  int getItemCount(const std::string& ename, int iline=0) const;
  const std::string& getItem(const std::string& ename, int iline, int iitem) const;
  const std::string& getItem(const std::string& ename, int iitem=0) const;

  int getItemInt(const std::string& ename, int iline, int iitem, int defaultval) const;
  int getItemInt(const std::string& ename, int item=0, int defaultval=0) const;

  float getItemFloat(const std::string& ename, int iline, int iitem, float defaultval) const;
  float getItemFloat(const std::string& ename, int item=0, float defaultval=0) const;

private:
  void tokenize(const std::string& line, std::vector<std::string>& tokens);
  struct Line {
    int index;
    std::string raw;
    std::vector<std::string> value;
  };

  struct Entry {
    Entry(const std::string& title);
    ~Entry();
    void addLine(const std::string& raw, const std::vector<std::string>& vals);
    std::string title;
    std::vector<Line*> lines;
  };

  std::map<std::string, Entry*> m_entries;

};


#endif // SimpleCardfile_hh_included
