#! /bin/bash

plotFold="Plots"
editeur="gedit"

echo "		s : solo"
echo "		m : multi"


read plot_choix
	
	case ${plot_choix} in
	
		s) echo "s : solo"
		
		    # Indications 
		
			echo "	Pour y accéder : load Plots/Solo."
			echo "	Initialiser avec load Plots/Solo/init.p"
			echo "	Appuyez sur Entrée."
			read
			gnuplot ;;
			
		m) echo "m : multi"
		
            # Préparation
            
            cd ${plotFold}
	
	            cd Multi
		
		            ls potChi_moy_var_copy.p
		            if test $? -ne 0
		            then
			            cp potChi_moy_var.p potChi_moy_var_copy.p
		            fi
	
	            cd ..
            cd ..
            
            # Indications
            
			echo "	Pour y accéder : load Plots/Multi."
			gnuplot ;;
			
		*) echo "		Je n'ai pas compris." ;;
	
esac
