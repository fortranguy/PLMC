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
		potChiRms=$(grep "potChiEx.rms" rapport.out | cut -d = -f 2)
		
		potChiErreur=$(awk 'NR==1 {print $3/$2}' bunching_activInv.out)
		#Delta potChi = Delta activInv / activInv
			# premier : suffisant ?
		
		echo ${Ncol}"	"${potChiMoy}"	"${potChiErreur} >> ../potChi_moy_var.out
		echo ${Ncol}"	"${potChiRms} >> ../potChi_rms_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
