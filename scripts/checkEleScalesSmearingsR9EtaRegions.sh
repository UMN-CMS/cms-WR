#!/bin/bash

#do not change the order of these arrays
etaMinVals=('0.' '1.0' '1.566' '2.0')
etaMaxVals=('1.0' '1.422' '2.0' '2.5')

#start in cmsWR dir
eval "cd test"
eval "sed 's@BINARYRNINE@true@g' checkEleScalesSmearingstemp.C > tempHighRNine.C"
eval "sed 's@BINARYRNINE@false@g' checkEleScalesSmearingstemp.C > tempLowRNine.C"

#loop over the pairs of eta bounds for low RNine, run the plotting macro, then high RNine, then run the macro again
for i in ${!etaMinVals[*]}
do
	eval "sed 's@ETAMAX@${etaMaxVals[$i]}@g' tempHighRNine.C > tempHighRNineHighEta.C"
	eval "sed 's@ETAMIN@${etaMinVals[$i]}@g' tempHighRNineHighEta.C > checkEleScalesSmearings.C"
	eval "root -l -b -q runCheckEleScalesSmearings.C"
	rm checkEleScalesSmearings.C tempHighRNineHighEta.C

	eval "sed 's@ETAMAX@${etaMaxVals[$i]}@g' tempLowRNine.C > tempLowRNineHighEta.C"
	eval "sed 's@ETAMIN@${etaMinVals[$i]}@g' tempLowRNineHighEta.C > checkEleScalesSmearings.C"
	eval "root -l -b -q runCheckEleScalesSmearings.C"
	rm checkEleScalesSmearings.C tempLowRNineHighEta.C

done

rm tempHighRNine.C tempLowRNine.C

eval "cd ../."
