#! /bin/bash

outFold_i="Out"

nSimus=$(find . -name "Out*" -type d | wc -l)

iSimu=0
until test ${iSimu} -eq ${nSimus}
do

	iSimu=$(expr ${iSimu} \+ 1)

	cd ${outFold_i}${iSimu}
	
		Ncol=$(grep "Ncol1" rapport.out | cut -d = -f 2)
              
        jour_ini=$(awk '$NF=="snapShotIni.out" {print $(NF-2)}' rapport.out)
        jour_fin=$(awk '$NF=="snapShotFin.out" {print $(NF-2)}' rapport.out)
        delta_jour=$(expr ${jour_fin} \- ${jour_ini})
        
        heure_ini=$(awk '$NF=="snapShotIni.out" {print $(NF-1)}' rapport.out |\
            cut -d : -f 1)
        heure_fin=$(awk '$NF=="snapShotFin.out" {print $(NF-1)}' rapport.out |\
            cut -d : -f 1)
        delta_heure=$(expr ${heure_fin} \- ${heure_ini})
            
        minute_ini=$(awk '$NF=="snapShotIni.out" {print $(NF-1)}' rapport.out |\
            cut -d : -f 2)
        minute_fin=$(awk '$NF=="snapShotFin.out" {print $(NF-1)}' rapport.out |\
            cut -d : -f 2)
        delta_minute=$(expr ${minute_fin} \- ${minute_ini})
        
        total_heure=$(expr ${delta_jour} \* 24 \+ ${delta_heure})
        total_minute=$(expr ${total_heure} \* 60 \+ ${delta_minute})
        
		echo "Ini :" ${jour_ini} ${heure_ini} ${minute_ini}
		echo "Fin :" ${jour_fin} ${heure_fin} ${minute_fin}
		echo "Delta" : ${delta_jour} ${delta_heure} ${delta_minute}
		echo "Total ":  "(${total_heure})" ${total_minute}
		
		echo "Rapport n°"${iSimu}" lu."
		
	cd ..

done
