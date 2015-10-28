#!/bin/bash

#all arrays must have the same number of elements
nuMass=(520 2080)
wrMass=(2600 2600)
#maxCount=80

for x in ${!wrMass[*]}
do

	for q in {1..80} 
	do
		#request_disk in KB, request_memory in MB
		eval 'condor_submit request_disk=30000000 request_memory=19000 inputGenToAnalyzedTreeMaker_WR_${wrMass[$x]}_NU_${nuMass[$x]}_$q'

	done

done

