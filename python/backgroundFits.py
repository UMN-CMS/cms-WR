import ROOT
""" \addtogroup BackgroundFits_Group Background Fits 
@{
"""

"""
Background fits for WR
Binned at 100 GeV

Fits are in the ranges from 600 to 3000 GeV in the 4-object mass. The results of the corresponding fits are:
"""

### TTBar Fit \ingroup BackgroundFits_Group
TTbar_mumu = ROOT.TF1("TTbar_mumu","[0]*exp(-[1]*x**(0.5))",600,10000)
""" ttbar:
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
1  c0           6.76149e+04   3.95259e+04   1.87561e+00  -4.75975e-09
2  c1           2.69119e-01   1.95857e-02   1.52661e-06   9.11054e-03
"""
TTbar_mumu.SetParameters(6.76149e+04/100., 2.69119e-01)

DYFit_mumu = ROOT.TF1("DYFit_mumu","[0]*exp(-[1]*x**(0.5))", 600, 10000)
"""
DY Fit:
NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
1  c0           5.37454e+03   4.01184e+03   5.00759e-01   9.64374e-08
2  c1           2.12876e-01   2.39904e-02   3.03108e-06  -1.85124e-02
"""
DYFit_mumu.SetParameters(5.37454e+03/100.,2.12876e-01)

All_ee = ROOT.TF1("All_ee","[0]*exp([1]*x)", 600, 10000)
"""
The M_EEJJ fit parameters for the sum of TTBar, DY, WZ, and ZZ MC datasets in the EEJJ channel are as follows:

1  coef         2.32541e+02   1.52408e+01   1.60990e-04  -2.26803e-01
2  expPower    -3.41355e-03   2.44022e-04   2.93939e-05   1.83283e+00
"""
All_ee.SetParameters( 2.32541e+02/50.,-3.41355e-03)

### @}
