#!/bin/bash

for q in {1..200} 
do
	#request_disk in KB, request_memory in MB
	eval 'condor_submit request_disk=10000000 request_memory=5000 genSimMaker_$q'

done


