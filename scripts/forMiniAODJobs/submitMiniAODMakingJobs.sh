#!/bin/bash

for q in {2..26} 
do
	#request_disk in KB, request_memory in MB
	eval 'condor_submit request_disk=30000000 request_memory=19000 miniAODMaker_$q'

done


