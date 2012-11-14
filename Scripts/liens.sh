#! /bin/bash

mcIn="../../MC_Cano"
mcOut="MC_Cano_L_16_Ncol_960_bunching"

ls ${mcOut}
if test $? -ne 0
then
	mkdir ${mcOut}
	cd ${mcOut}
else
	echo ${mcOut} " existe déjà."
	exit
fi

cp -r ${mcIn}/data.f90 .

listLiens="makefile mc_canonique.f90 mod_physique.f90 mod_tools.f90 bunching.f90"
listLiens=${listLiens}" notes.txt Plots Data avgRms.sh"
listLiens=${listLiens}" exec.sh plot_eTot.sh plot_activInv.sh plot_snaps.sh"

for item in ${listLiens}
do

	ln -s ${mcIn}/${item} ${item}
	echo ${item}

done

cd .. # nécessaire
