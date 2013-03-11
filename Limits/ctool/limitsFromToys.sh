#!/bin/sh

grid=$1
loc=${grid%/*}
ifile=${grid##*/}
ifile=${ifile%%.root}
ifile="${loc}/input/${ifile}"
mw=${grid##*limit_}
mw=${mw%%_*}
lf="${mw}-limit.log"

rm -f ${lf}
touch ${lf}

echo "OBS" >> $lf
combine -M HybridNew --grid=${grid} -m ${mw} ${ifile} >> $lf

echo "EXP" >> $lf
combine -M HybridNew --grid=${grid} --expectedFromGrid 0.5 -m ${mw} ${ifile} >> $lf
echo "EXP_P1S" >> $lf
combine -M HybridNew --grid=${grid} --expectedFromGrid 0.85 -m ${mw} ${ifile} >> $lf
echo "EXP_M1S" >> $lf
combine -M HybridNew --grid=${grid} --expectedFromGrid 0.15 -m ${mw} ${ifile} >> $lf
echo "EXP_P2S" >> $lf
combine -M HybridNew --grid=${grid} --expectedFromGrid 0.025 -m ${mw} ${ifile} >> $lf
echo "EXP_M2S" >> $lf
combine -M HybridNew --grid=${grid} --expectedFromGrid 0.975 -m ${mw} ${ifile} >> $lf
