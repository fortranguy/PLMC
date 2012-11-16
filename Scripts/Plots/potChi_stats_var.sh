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
		potChiMoy=$(grep "potChi.moy" rapport.out | cut -d = -f 2)
		potChiRms=$(grep "potChi.rms" rapport.out | cut -d = -f 2)
		
		eMoyErreur=$(cat bunching_potChi.out | head -n1 | awk '{print $NF}')
		# premier : suffisant ?
		
		echo ${Ncol}"	"${potChiMoy}"	"${eMoyErreur} >> ../potChi_moy_var.out
		echo ${Ncol}"	"${potChiRms} >> ../potChi_rms_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
