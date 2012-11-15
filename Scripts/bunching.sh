#! /bin/bash

binFold="Binaries"
outFold="Out"
exec2="xBunching"

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
		fi 
		
		cd ${outFold_i}${iSimu}
		
			mv ../../${binFold}/${binFold_i}${iSimu}/*.out .
		
			ls *.out > /dev/null
			if test $? -ne 0
			then
				echo "Pas de données."
				exit
			fi

			# Bunching
			cp ../../${binFold}/${binFold_i}${iSimu}/${exec2} .
			./${exec2} ${nBunching}
			rm ${exec2}

			# Rapport	
			cp ../../Scripts/avgRms.sh .
			./avgRms.sh
			rm avgRms.sh
			cat rapport.out

		cd ..

	done

cd ..
