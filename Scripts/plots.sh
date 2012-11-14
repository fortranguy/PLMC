#! /bin/bash


echo "		e : énergie : observable et bunching"
	echo "		a : inverse de l'activité : observable et bunching"
	echo "		s : snapshot : 2D et 3D"
	read plot_choix
	
	case ${plot_choix} in
	
		e) Scripts/Plots/plot_eTot.sh ;;
		a) Scripts/Plots/plot_activInv.sh ;;
		s) Scripts/Plots/plot_snaps.sh ;;
		*) echo "		Je n'ai pas compris." ;;
	
esac
