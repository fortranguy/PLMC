#! /bin/bash

outFold_i="Out"

nSimus=$(cat ../../Temp/nSimus.in)
nSimus=6 # sale

iSimu=0
until test ${iSimu} -eq ${nSimus}
do

	iSimu=$(expr ${iSimu} \+ 1)

	cd ${outFold_i}${iSimu}
	
		Ncol=$(grep "Ncol1" rapport.out | cut -d = -f 2)
		eRms=$(grep "écart-type" rapport.out | head -n1 | cut -d = -f 2)
		
		echo ${Ncol}"	"${eRms} >> ../rms_eTot_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
