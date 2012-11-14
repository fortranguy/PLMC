#! /bin/bash

listOut=$(find . -name "Out" -type d)
NcolFil="Ncol.out"
eMoyFil="eMoy.out"
eMoyErreurFil="eMoyErreur.out"

rm -f ${NcolFil} ${eMoyFil} ${eMoyErreurFil} moy_eTot_var*.out

for out in ${listOut}

do
	echo ${out}

	# Energie moyenne
		
	cat ${out}"/rapport.out" | head -n4 | tail -n1 > ${NcolFil} # fragile
	cat ${out}"/rapport.out" | tail -n5 | head -n1 > ${eMoyFil} # fragile

	Ncol=$(cut -d = -f 2 ${NcolFil})
	eMoy=$(cut -d = -f 2 ${eMoyFil})
	
	# Erreur
	
	cat ${out}"/bunching_eTot.out" | head -n1 > ${eMoyErreurFil}
	
	eMoyErreur=$(awk '{print $NF}' ${eMoyErreurFil})
	
	# Ecriture
	
	echo ${Ncol}"	"${eMoy}"	"${eMoyErreur} >> moy_eTot_var.out
	

done

sort -n moy_eTot_var.out > moy_eTot_var_sorted.out
