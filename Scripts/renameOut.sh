#! /bin/bash

outBin_i="Out"

echo "Décalage : Out1 -> Outi, i entier"

cd Out

	listOut=$(ls |sort -n)
	nOut=$(ls | wc -l)
	
	echo ${listOut}
	
	# Entrée du décalage
	
	status="non"
		
	until [ ${status} = "oui" ]
	do
	
		echo "Entrez le nombre de décalage : "
		read nDecal
	
		iSimuDecal=$(expr 1 \+ ${nDecal})
		echo "Acceptez-vous le changement Out1 -> Out"${iSimuDecal}" ? (oui/non)"
		read status
		
		if test ${iSimuDecal} -le ${nOut}
		then
			echo "Le nombre est trop petit."
			ls
			status="non"
		fi
			
		
	done
	
	# Vérification

	nSimusDecal=$(expr ${iSimuDecal} \+ ${nOut})
	
	until test ${iSimuDecal} -eq ${nSimusDecal}
	do
	
		ls ${outBin_i}${iSimuDecal}
		if test $? -eq 0
		then
			echo "Le dossier n°"${iSimuDecal}" existe."
			exit
		fi
		
		iSimuDecal=$(expr ${iSimuDecal} \+ 1)
		
	done
	
	# Changement de nom
	
	iSimuDecal=$(expr 1 \+ ${nDecal})
	
	for out_i in ${listOut}
	do		
		mv ${out_i} ${outBin_i}${iSimuDecal}
		
		echo ${out_i}"->Out"${iSimuDecal}
		iSimuDecal=$(expr ${iSimuDecal} \+ 1)
	
	done

cd ..
