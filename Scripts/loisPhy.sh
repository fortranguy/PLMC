#! /bin/bash

dataFold="Data"
inFold="Rho_cst_Ncol_80-400_"
outFold="Rho_cst"

plotFold="Plots"
Ndeb=80
Nfin=400

# Vérification

ls ${dataFold}/${inFold}
if test $? -ne 0
then
	echo "Le dossier ${inFold} n'existe pas dans ${dataFold}."
	exit
fi

# Statistiques

eTot_stat_script="eTot_stats_var.sh"
potChi_stat_script="potChi_stats_var.sh"

cd ${dataFold}	
	
	cd ${inFold}

		mkdir ${outFold}
		if test $? -ne 0
		then
			rm -f ${outFold}/*.out
		fi
		
		echo "Energie : statistiques"
		
		cp ../../Scripts/Plots/${eTot_stat_script} .
		./${eTot_stat_script}
		rm ./${eTot_stat_script}
		
		echo "Potentiel chimique : statistiques"
		
		cp ../../Scripts/Plots/${potChi_stat_script} .
		./${potChi_stat_script}
		rm ./${potChi_stat_script}

		mv *.out ${outFold}
		
	cd ..

cd ..

# Graphes

cd ${plotFold}
	
	cd Multi
	
		echo "dataFold='${dataFold}'" > location.p
		echo "dataFold_i='${inFold}'" >> location.p
		echo "outFold='${outFold}'" >> location.p
		echo "Ndeb="${Ndeb} >> location.p
		echo "Nfin="${Nfin} >> location.p
	
	cd ..

cd ..
