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
		
			cp ../../${binFold}/${binFold_i}${iSimu}/*.out . # sécurité ?
		
			echo "Test : présence de fichiers .out ?"
			ls *.out
			if test $? -ne 0
			then
				echo "Pas de données."
				exit
			fi
			
			# Date
			echo " Date de lancement : " $(cat dateIni.out) >> rapport.out

			# Bunching
			pwd
			cp ../../${binFold}/${binFold_i}${iSimu}/${exec2} .
			./${exec2} ${nBunching}
			rm ${exec2}

			# Rapport	
			cp ../../Scripts/${statScript} .
			./${statScript}
			rm ${statScript}
			cat rapport.out

		cd ..

	done

cd ..

# Plot : mauvais endroit ?

cd ${plotFold}
	cd Solo
		ls location.p
		if test $? -ne 0
		then
			cp SaveLoc/location.p .
		fi
	cd ..
cd ..
