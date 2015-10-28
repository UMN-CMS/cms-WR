#!/bin/bash

nuMass=(500 900 1700 2100 2500)
wrMass=(2600)
#maxCount=300


for r in ${!wrMass[*]}
do

	for w in ${!nuMass[*]}
	do

		for q in {1..300} 
		do
			#request_disk in KB, request_memory in MB
			eval 'condor_submit request_disk=30000000 request_memory=19000 genSimToAnalyzedTreeMaker_WR_${wrMass[$r]}_NU_${nuMass[$w]}_$q'

		done

	done

done


