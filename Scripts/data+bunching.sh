#! /bin/bash

binFold="Binaries"
outFold="Out"
plotFold="Plots"
exec2="xBunching"
statScript="avgRms.sh"

binFold_i="Bin"
outFold_i="Out"

nSimus=$(cat Temp/nSimus.in)

echo "Entrez le nombre d'itérations :"
read nBunching

# Données
cd ${outFold}

	iSimu=0
	until test ${iSimu} -eq ${nSimus}
	do
		
		iSimu=$(expr ${iSimu} \+ 1)

		ls ${outFold_i}${iSimu}
		if test $? -ne 0
		then
			mkdir ${outFold_i}${iSimu}
		else
			echo "Retraitement des données."
		fi 
		
		cd ${outFold_i}${iSimu}
		
			ls -l ../../${binFold}/${binFold_i}${iSimu} | \
				grep "snapShotIni.out" > snapIni_date.out
			ls -l ../../${binFold}/${binFold_i}${iSimu} | \
				grep "snapShotFin.out" > snapFin_date.out
		
			cp ../../${binFold}/${binFold_i}${iSimu}/*.out . # sécurité
			
			cp rapport.out rapportSave.out # sécurité bis
		
			echo "Test : présence de fichiers .out ?"
			ls *.out
			if test $? -ne 0
			then
				echo "Pas de données."
				exit
			fi
			
			# Date
			echo " Date de lancement : $(cat dateIni.out)" >> rapport.out
			rm dateIni.out
			echo " Dates des snapshots : " >> rapport.out
				cat snapIni_date.out >> rapport.out
				rm	snapIni_date.out
				cat snapFin_date.out >> rapport.out
				rm snapFin_date.out


			# Bunching
			pwd
			cp ../../${binFold}/${binFold_i}${iSimu}/${exec2} .
			./${exec2} ${nBunching}
			rm ${exec2}

			# Rapport	
			cp ../../Scripts/${statScript} .
			./${statScript}
			rm ${statScript}
			head rapport.out

		cd ..

	done

cd ..

# Plot : mauvais endroit ?

cd ${plotFold}
	
	cd Multi
	
		ls location.p
		if test $? -ne 0
		then
			cp location_save.p location.p
		fi
	
	cd ..
cd ..
