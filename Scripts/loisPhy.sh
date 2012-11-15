#! /bin/bash

dataFold="Data"

cd ${dataFold}

	inFold="Rho_cst_LNcol_var"
	outFold="Rho_cst"	
	
	cd ${inFold}

		mkdir ${outFold}
		if test $? -ne 0
		then
			rm -f ${outFold}/*.out
		fi

		cp ../../Scripts/Plots/moy_eTot_var.sh .
		./moy_eTot_var.sh
		rm ./moy_eTot_var.sh

		mv *.out ${outFold}
		
	cd ..

cd ..
