#! /bin/bash

binFold="Binaries"
outFold="Out"

binFold_i="Bin"
outFold_i="Out"

nSimus=$(cat Temp/nSimus.in)

echo "Entrez le nombre d'itérations :"
read nBunching

cd ${outFold}

	# Bunching
	#cp ../${binFold}/xBunching .
	#./xBunching ${nBunching}
	#rm xBunching
	
cd ..

# Données
cd ${outFold}

	iSimu=0
	until test ${iSimu} -eq ${nSimus}
	do
		
		iSimu=$(expr ${iSimu} \+ 1)

		cd ${outFold_i}${iSimu}
		
			ls *.out > /dev/null
			if test $? -ne 0
			then
				echo "Pas de données."
				exit
			fi

			# Bunching
			cp ../../${binFold}/${binFold_i}${iSimu}/xBunching .
			./xBunching ${nBunching}
			rm xBunching

			# Rapport	
			cp ../../Scripts/avgRms.sh .
			./avgRms.sh
			rm avgRms.sh
			cat rapport.out
		
		cd ..
	
	done
	
cd ..
