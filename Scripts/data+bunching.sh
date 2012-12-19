#! /bin/bash

binFold="Binaries"
outFold="Out"
plotFold="Plots"
exec2="xBunching"
exec3="xDistrib"
statScript="avgRms.sh"

binFold_i="Bin"
outFold_i="Out"

nSimus=$(cat Temp/nSimus.in)

echo "Entrez le nombre d'itérations :"
read nBunching

# Préparation

cd ${binFold}

	iSimu=0
	until test ${iSimu} -eq ${nSimus}
	do
	
	    iSimu=$(expr ${iSimu} \+ 1)
	
	    cd ${binFold_i}${iSimu}
	    
	    	ls snapShotFin.out
	    	if test $? -ne 0
			then
				echo "Simulation pas encore finie ?"
				exit
			fi

	    cd ..
	
	done

cd ..

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
		
			cp ../../${binFold}/${binFold_i}${iSimu}/data_copy.f90 .
			cp ../../${binFold}/${binFold_i}${iSimu}/*.out .
		
			echo "Test : présence de fichiers .out ?"
			ls *.out
			if test $? -ne 0
			then
				echo "Pas de données."
				exit
			fi
			
						
			cp rapport.out rapportSave.out # sécurité
			
			ls epp_dist.out
	    	if test $? -eq 0
				cat epp_dist.out >> rapport.out
			fi

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
