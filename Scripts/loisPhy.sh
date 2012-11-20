#! /bin/bash

dataFold="Data"
inFold="RhoCst_N_80-400_"
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

cd ${dataFold}	
	
	cd ${inFold}

		mkdir ${outFold}
		if test $? -ne 0
		then
			rm -f ${outFold}/*.out
		fi
		
		echo "Energie : statistiques"
		
		cp ../../Scripts/Plots/eTot_stats_var.sh .
		./eTot_stats_var.sh
		rm ./eTot_stats_var.sh
		
		echo "Potentiel chimique : "
		
		    echo "statistiques"
		
		    cp ../../Scripts/Plots/potChi_stats_var.sh .
		    ./potChi_stats_var.sh
		    rm ./potChi_stats_var.sh
		    
		    echo "HS"
		    
		    cp ../../Scripts/Plots/potChi_hs_var.sh .
		    ./potChi_hs_var.sh
		    rm ./potChi_hs_var.sh
		
		echo "Temps de calcul"
		
		cp ../../Scripts/Plots/temps_calc_var.sh .
		./temps_calc_var.sh
		rm ./temps_calc_var.sh
		
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
