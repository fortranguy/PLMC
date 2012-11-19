#! /bin/bash

dataFold="Data"
plotFold="Plots"
inFold="Rho_cst_Ncol_80-400_nMove1000"
outFold="Rho_cst"

cd ${plotFold}
	
	cd Multi
	
		echo "dataFold='"${dataFold}"'" > location.p
		echo "dataFold_i='"${inFold}"'" >> location.p
		echo "outFold='"${outFold}"'" >> location.p
		echo "Ndeb=80" >> location.p
		echo "Nfin=400" >> location.p
	
	cd ..

cd ..

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
