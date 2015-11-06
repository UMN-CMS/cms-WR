#!/bin/bash

wrMass=800
nuMass=50
minNuMass=50
maxWrMass=3500
increment=50

q=1

while [ $wrMass -le $maxWrMass ]; do
	while [ $nuMass -lt $wrMass ]; do
		eval 'condor_submit request_disk=30000000 request_memory=19000 genAndRunOfflineAnalysis_WR_${wrMass}_NU_${nuMass}_${q}'

		#increment nuMass before restarting the loop
		let nuMass=nuMass+$increment
	done
	#reset the nuMass value to the lower bound value before restarting the while loop
	#with a higher wrMass value
	let nuMass=minNuMass
	let wrMass=wrMass+$increment

done

