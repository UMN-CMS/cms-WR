#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
/** this simple class represents a cut variable used in some selection
 *  each object of this class has a name, value, and a bool to indicate if
 *  the cut variable is an upper bound (ecal iso) or lower bound (pT)
 */

class SlimCutVar {
	public:
		SlimCutVar(std::string nameOfCut, float val, bool setAsUpperBound):
			_cutName(nameOfCut),
			_threshVal(val),
			_isUpperBound(setAsUpperBound){
			};
		
		SlimCutVar(){};
		
		///the input string must be of the form "number,sign,char array" where sign = > or <
		void setAttributesFromString(std::string initializerString){
			char sign;
			char name[50];
			sscanf(initializerString.c_str(), "%f,%c,%s", &_threshVal, &sign, name);
			_cutName = std::string(name);
			if(sign == '>') _isUpperBound=false;
			else if(sign == '<') _isUpperBound=true;

		}

		/**use this fxn to check that the minimum, maximum, and step size for a SlimCutVar are assigned reasonable
		 * values, whether the SlimCutVar is identified as an upper bound or lower bound, and the detector region
		 * (EB, tracked EE, etc) in which this SlimCutVar is relevant.
		 * the SlimCutVar name is printed first
		 */
		friend std::ostream& operator << (std::ostream& os, const SlimCutVar a){
			char c = a._isUpperBound ? '<' : '>';
			os <<  a._cutName << ",\t" << a._threshVal << ",\t" << c << std::endl;
			return os;
		}

		inline std::string printNameVal() const{
			char line[250];
			sprintf(line, "%s\t%.2f", _cutName.c_str(), _threshVal);
			return std::string(line);
		}

		std::string getCutName(){ return _cutName;}
		float getThresholdValue(){ return _threshVal;}
		bool cutIsUpperBound(){ return _isUpperBound;}

	private:
		std::string _cutName;	///< name of cut (dileptonMass, dRleptonJet, etc)
		float _threshVal;	///< threshold value of cut variable
		bool _isUpperBound;	///< indicates if this cut will be used as an upper bound (someVal < _threshVal)
};//end class SlimCutVar

