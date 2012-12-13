#! /bin/bash

outFold_i="Out"

nSimus=$(find . -name "Out*" -type d | wc -l)

iSimu=0
until test ${iSimu} -eq ${nSimus}
do

	iSimu=$(expr ${iSimu} \+ 1)

	cd ${outFold_i}${iSimu}
	
		nCol=$(grep "Ncol1" rapport.out | cut -d = -f 2)
		vol=$(grep "Vol" rapport.out | cut -d = -f 2)
		eMoy=$(grep "eTot.moy" rapport.out | cut -d = -f 2)
		eMoyPart=$(grep "Energie moyenne par particule" rapport.out |\
			cut -d = -f 2)
		eRms=$(grep "eTot.rms" rapport.out | cut -d = -f 2)
		
		N=$(echo ${nCol})
		# Etrange entourloupe
		
		# six : assez général ?
		iBunch=6
		eMoyErreur=$(awk -v i=${iBunch} 'NR==i {print $3}' bunching_eTot.out)		
		eMoyErreurPart=$(awk -v i=${iBunch} -v x=${N} 'NR==i {print $3/x}' \
			bunching_eTot.out)
		
		echo ${vol}"	"${eMoy}"	"${eMoyErreur} >> ../eTot_moy_var.out
		echo ${vol}"	"${eMoyPart}"	"${eMoyErreurPart} >> \
			../epp_moy_var.out
		echo ${vol}"	"${eRms} >> ../eTot_rms_var.out
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
