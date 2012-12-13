#! /bin/bash

outFold_i="Out"

nSimus=$(find . -name "Out*" -type d | wc -l)

iSimu=0
until test ${iSimu} -eq ${nSimus}
do

	iSimu=$(expr ${iSimu} \+ 1)

	cd ${outFold_i}${iSimu}
	
		vol=$(grep "Vol" rapport.out | cut -d = -f 2)
		potChiIdMoy=$(grep "Potentiel chimique idéal " rapport.out | \
		    cut -d = -f 2)
		potChiExMoy=$(grep "Potentiel chimique (excès) moyen" rapport.out | \
		    cut -d = -f 2)
		potChiRms=$(grep "potChiEx.rms" rapport.out | cut -d = -f 2)
		
		# six : assez général ?
		iBunch=6
		potChiErreur=$(awk -v i=${iBunch} 'NR==i {print $3/$2}' \
			bunching_activInv.out)
		#Delta potChi = Delta activInv / activInv
			# premier : suffisant ?
		
		echo "${vol}    ${potChiIdMoy} ${potChiExMoy}   ${potChiErreur}" >> \
		    ../potChi_moy_var.out
		echo "${vol}    ${potChiRms}" >> ../potChi_rms_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
