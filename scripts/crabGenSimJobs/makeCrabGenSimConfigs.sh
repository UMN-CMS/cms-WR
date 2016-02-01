#!/bin/bash

nuMasses=(50 100 200 300)
wrMasses=(2600 2600 2600 2600)

for h in ${!nuMasses[*]}
do

	#make crab config file
	eval "sed 's/MWrUndef/MWr${wrMasses[$h]}/g' WR_MWrUndef_to_ENu_MNuUndef_13TeV_GENSIM_crab.py > WR_tempOne_crab.py"
	eval "sed 's/MNuUndef/MNu${nuMasses[$h]}/g'  WR_tempOne_crab.py > WR_MWr${wrMasses[$h]}_to_ENu_MNu${nuMasses[$h]}_13TeV_GENSIM_crab.py"
	rm WR_tempOne_crab.py

	#make python config file which governs the GEN simulation 
	eval "sed 's/MWrUndef/MWr${wrMasses[$h]}/g' WR_MWrUndef_ToENu_MNuUndef_13TeV_GEN_SIM.py > WR_tempOne.py"
	eval "sed 's/MNuUndef/MNu${nuMasses[$h]}/g'  WR_tempOne.py > WR_tempTwo.py"
	eval "sed 's/WRMASS/${wrMasses[$h]}/g'  WR_tempTwo.py > WR_tempThree.py"
	eval "sed 's/NUMASS/${nuMasses[$h]}/g'  WR_tempThree.py > WR_MWr${wrMasses[$h]}_ToENu_MNu${nuMasses[$h]}_13TeV_GEN_SIM.py"
	rm WR_tempOne.py WR_tempTwo.py WR_tempThree.py

done

