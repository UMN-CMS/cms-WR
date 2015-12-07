""" \addtogroup BackgroundFits_Group Background Fits 
@{
"""

"""
Background fits for WR
Binned at 100 GeV

Fits are in the ranges from 600 to 3000 GeV in the 4-object mass. The results of the corresponding fits are:
"""
import ROOT

### TTBar Fit \ingroup BackgroundFits_Group
TTbarFit = ROOT.TF1("TTbarFit","[0]*exp(-[1]*x**(0.5))",600,10000)
""" ttbar:
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
1  c0           6.76149e+04   3.95259e+04   1.87561e+00  -4.75975e-09
2  c1           2.69119e-01   1.95857e-02   1.52661e-06   9.11054e-03
"""
TTbarFit.SetParameters(6.76149e+04, 2.69119e-01)

DYFit = ROOT.TF1("DYFit","[0]*exp(-[1]*x**(0.5))", 600, 10000)
"""
DY Fit:
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
1  c0           5.37454e+03   4.01184e+03   5.00759e-01   9.64374e-08
2  c1           2.12876e-01   2.39904e-02   3.03108e-06  -1.85124e-02
"""
DYFit.SetParameters(5.37454e+03,2.12876e-01)

### @}
