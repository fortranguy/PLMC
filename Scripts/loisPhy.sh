#! /bin/bash

dataFold="Data"

eTot_stat_script="eTot_stats_var.sh"
potChi_stat_script="potChi_stats_var.sh"

cd ${dataFold}

	inFold="Rho_cst_LNcol_var"
	outFold="Rho_cst"	
	
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
		
		echo "Potentiel chimque : statistiques"
		
		cp ../../Scripts/Plots/${potChi_stat_script} .
		./${potChi_stat_script}
		rm ./${potChi_stat_script}

		mv *.out ${outFold}
		
	cd ..

cd ..
