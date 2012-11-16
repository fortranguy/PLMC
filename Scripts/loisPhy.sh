#! /bin/bash

dataFold="Data"

cd ${dataFold}

	inFold="[approx]Rho_cst_LNcol_var"
	outFold="Rho_cst"	
	
	cd ${inFold}

		mkdir ${outFold}
		if test $? -ne 0
		then
			rm -f ${outFold}/*.out
		fi
		
		echo "Energie moy"
		
		cp ../../Scripts/Plots/eTot_moy_var.sh .
		./eTot_moy_var.sh
		rm ./eTot_moy_var.sh
		
		echo "Energie rms"
		
		cp ../../Scripts/Plots/eTot_rms_var.sh .
		./eTot_rms_var.sh
		rm ./eTot_rms_var.sh

		mv *.out ${outFold}
		
	cd ..

cd ..
