#! /bin/bash

outFold_i="Out"

nSimus=$(find . -name "Out*" -type d | wc -l)

iSimu=0
until test ${iSimu} -eq ${nSimus}
do

	iSimu=$(expr ${iSimu} \+ 1)

	cd ${outFold_i}${iSimu}
	
		Ncol=$(grep "Ncol1" rapport.out | cut -d = -f 2)
		potChiIdMoy=$(grep "Potentiel chimique idéal " rapport.out | \
		    cut -d = -f 2)
		potChiExMoy=$(grep "Potentiel chimique (excès) moyen" rapport.out | \
		    cut -d = -f 2)
		potChiRms=$(grep "potChiEx.rms" rapport.out | cut -d = -f 2)
		
		potChiErreur=$(awk 'NR==1 {print $3/$2}' bunching_activInv.out)
		#Delta potChi = Delta activInv / activInv
			# premier : suffisant ?
		
		echo "${Ncol}    ${potChiIdMoy} ${potChiExMoy}   ${potChiErreur}" >> \
		    ../potChi_moy_var.out
		echo "${Ncol}    ${potChiRms}" >> ../potChi_rms_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
