#!/bin/bash

nuMasses=(50 100 200 300)
wrMasses=(2600 2600 2600 2600)

for h in ${!nuMasses[*]}
do
	#echo 'crab submit -c WR_MWr${wrMasses[$h]}_to_ENu_MNu${nuMasses[$h]}_13TeV_GENSIM_crab.py'
	eval 'crab submit -c WR_MWr${wrMasses[$h]}_to_ENu_MNu${nuMasses[$h]}_13TeV_GENSIM_crab.py'

done

