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
		eMoy=$(grep "eTot.moy" rapport.out | cut -d = -f 2)
		eRms=$(grep "eTot.rms" rapport.out | cut -d = -f 2)
		
		eMoyErreur=$(awk 'NR==1 {print $3}' bunching_eTot.out)
		# premier : suffisant ?
		
		echo ${Ncol}"	"${eMoy}"	"${eMoyErreur} >> ../eTot_moy_var.out
		echo ${Ncol}"	"${eRms} >> ../eTot_rms_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
