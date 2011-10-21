#include <map>
#include <string>
#include <vector>
#include <math.h>

#ifndef systematics_h_included
#define systematics_h_included 1


extern std::map<std::string,double> systematics10[4];
extern std::map<std::string,double> systematics11[4];

extern std::vector<std::string> systematicsList;

extern int special_syst_mode;

void loadSystematics();
const char* getSyst(int iChan, int iyear, const std::string& sname);

#endif // systematics_h_included
