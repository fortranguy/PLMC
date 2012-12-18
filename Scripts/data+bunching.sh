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
			
			ls snap1.out
	    	if test $? -eq 0
			then
				snapList=$(ls | grep "snap[0-9][0-9]*.out")
				tar cvf snap.tar ${snapList} > /dev/null
				./${exec3}
				rm ${snapList}	
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
			
			ls ../../${binFold}/${binFold_i}${iSimu}/snap.tar
	    	if test $? -eq 0
			then
				cp ../../${binFold}/${binFold_i}${iSimu}/snap.tar .
			end if
			
			cp rapport.out rapportSave.out # sécurité bis
		
			echo "Test : présence de fichiers .out ?"
			ls *.out
			if test $? -ne 0
			then
				echo "Pas de données."
				exit
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
