#! /bin/bash

echo "Monte-Carlo Canonique : Menu"
echo "	go! : compilation et exécution."
echo "	bunch : bunching pour obtenir l'erreur."
echo "	plot : trace les graphes."
echo "	quit : quitter le menu."

read choix
#echo ${choix}

#until ls

#do

	case ${choix} in
		
		go!)
			echo "	go! : compilation et exécution."
			Scripts/make+exec.sh ;;
	
		bunch)
			echo "	bunch : bunching pour obtenir l'erreur."
			Scripts/bunching.sh ;;
	
		plot)
			echo "	plot : trace les graphes."
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
			;;
			
		quit) echo "	Au-revoir !";;		
	
		*) echo "	Je n'ai pas compris." ;;
		
	esac

#done
