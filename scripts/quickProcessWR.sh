#!/bin/bash
mass_n=$(seq 1 25)

for m in $mass_n
do
	eval "./bin/analysis -d analysisCppOutputRootFiles/ -m signal -c MuMu --signalN $m >& wrMuMu_${m}.txt &"
	eval "./bin/analysis -d analysisCppOutputRootFiles/ -m signal -c EE --signalN $m >& wrEE_${m}.txt &"
done

