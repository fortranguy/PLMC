#! /bin/bash

editeur="gedit"

echo "		s : solo"
echo "		m : multi"


read plot_choix
	
	case ${plot_choix} in
	
		s) echo "s : solo"
			echo "	Pour y accéder : load Plots/Solo."
			echo "	Initialiser avec load Plots/Solo/init.p"
			echo "	Appuyez sur Entrée."
			read
			gnuplot ;;
			
		m) echo "m : multi"
			echo "	Pour y accéder : load Plots/Multi."
			gnuplot ;;
			
		*) echo "		Je n'ai pas compris." ;;
	
esac
