ch=EE
ls -v plots/selhists_WRto${ch}*.txt | sed "s/.*JJ_\(.*\)_.*/\1/" > m
grep  global_EE $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > ee
grep  global_EB $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > eb
grep  global_BB $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > bb
grep  "global	" $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > glb
for m in $(cat m); 
	do grep WRto${ch}JJ_$m configs/datasets.dat | cut -f5; 
done > e
paste e bb eb ee glb | awk '{ print $2/$1, $3/$1, $4/$1, $5/$1}' > ee_eff

ch=MuMu
grep  global_EE $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > ee
grep  global_EB $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > eb
grep  global_BB $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > bb
grep  "global	" $(ls -v plots/selhists_WRto${ch}*.txt) | cut -f2 > glb
for m in $(cat m); 
	do grep WRto${ch}JJ_$m configs/datasets.dat | cut -f5; 
done > e
paste e bb eb ee glb | awk '{ print $2/$1, $3/$1, $4/$1, $5/$1}' > mumu_eff

echo "#mass ele_bb ele_eb ele_ee ele_global mu_bb mu_eb mu_ee mu_global" > $1.tmp
paste m ee_eff mumu_eff  >> $1.tmp

cat $1.tmp |column -t > $1
rm $1.tmp
rm m ee eb bb glb ee_eff mumu_eff

python scripts/lepton_eff_notrigger.py
