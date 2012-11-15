#! /bin/bash

srcFold="Sources"
binFold="Binaries"
outFold="Out"
tmpFold="Temp"
exec1="xMC_Canonique"
exec2="xBunching"

condIni="cube"
binFold_i="Bin"
outFold_i="Out"
editeur="gedit"
exec1NewCore="xMC_Cano"

# ---------------------------------------------------------

# Préparation
ls ${binFold}
if test $? -ne 0
then
	mkdir ${binFold}
fi

ls ${outFold}
if test $? -ne 0
then
	mkdir ${outFold}
else
	rm -rf ${outFold}/*
fi

rm -rf ${binFold}/* # si problème

ls ${tmpFold}
if test $? -ne 0
then
	mkdir ${tmpFold}
else
	rm -rf ${tmpFold}/*
fi

# ---------------------------------------------------------

echo "Nombre de simulations : "
read nSimus
echo ${nSimus} > Temp/nSimus.in

echo "Paramètre qui varie : "
read param

# Création des dossiers pour exécutables :

cd ${binFold}

	iSimu=0
	until test ${iSimu} -eq ${nSimus}
	do

		iSimu=$(expr ${iSimu} \+ 1)
		echo "Création de Bin"${iSimu}"."
	
		mkdir ${binFold_i}${iSimu}

	done

cd ..

# Modification des données d'entrée :

cd ${srcFold}

	iSimu=0
	until test ${iSimu} -eq ${nSimus}
	do
	
		${editeur} data.f90&
		
		iSimu=$(expr ${iSimu} \+ 1)
		
		status="non"
		
		until [ ${status} = "oui" ]
		do
		
			echo "Avez-vous modifié le data n°"${iSimu}" ? (oui/non)"
			read status
		
		done
		
		make
		if test $? -ne 0
		then
			exit
		fi
		
		exec1New=${exec1NewCore}"_"${param}${iSimu}
		mv -i ${exec1} ../${binFold}/${binFold_i}${iSimu}/${exec1New}
		mv -i ${exec2} ../${binFold}/${binFold_i}${iSimu}

	done

cd ..

# Exécution

cd ${binFold}

	iSimu=0
	until test ${iSimu} -eq ${nSimus}
	do

		iSimu=$(expr ${iSimu} \+ 1)
		
		cd ${binFold_i}${iSimu}
			
			echo
			exec1New=${exec1NewCore}"_"${param}${iSimu}
			echo "Exécution de "${exec1New}
			
			date >> dateIni.out
			time ./${exec1New} ${condIni} &
			if test $? -ne 0
			then
				exit
			fi
			
		cd ..

	done

cd ..
