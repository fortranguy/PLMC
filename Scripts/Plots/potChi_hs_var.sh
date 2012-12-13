#! /bin/bash

outFold_i="Out"

nSimus=$(find . -name "Out*" -type d | wc -l)

iSimu=0
until test ${iSimu} -eq ${nSimus}
do

	iSimu=$(expr ${iSimu} \+ 1)

	cd ${outFold_i}${iSimu}
	
	    vol=$(grep "Vol" rapport.out | cut -d = -f 2)
		eta=$(grep "Compacité" rapport.out | cut -d = -f 2)
		potChiId=$(grep "Potentiel chimique idéal " rapport.out | \
		    cut -d = -f 2)
		
		echo "${vol}    ${eta}  ${potChiId}" >> ../potChi_hs_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
