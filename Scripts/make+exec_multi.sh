#! /bin/bash

srcFold="Sources"
binFold="Binaries"
outFold="Out"

execFold="Bin"
editeur="gedit"

# ---------------------------------------------------------

echo "Nombre de simulations : "
read nSimus
echo ${nSimus}

echo "Paramètre qui varie : "
read param
echo ${param}

#Création des dossiers pour exécutables :

cd ${binFold}

	iSimu=0
	until test ${iSimu} -eq ${nSimus}
	do

		iSimu=$(expr ${iSimu} \+ 1)
		echo "Création de Bin"${iSimu}"."
	
		mkdir ${execFold}${iSimu}

	done

cd ..

#Modification des données d'entrée :

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

	done

cd ..
