#! /bin/bash

echo "Monte-Carlo Canonique : Menu"
echo "	n : nettoyage (binaires et données)."
echo "	go! : compilation et exécution de programme(s)."
echo "	d : traitement des données (inclut bunching)."
#echo "	p : trace les graphes."
echo "	l : lois physiques." # pompeux
echo "	q : quitter le menu."

read choix


case ${choix} in

	n)
		echo "	n : nettoyage (binaires, données)."
		Scripts/clean.sh ;;
	
	go!)
		echo "	go! : compilation et exécution de programme(s)."
		Scripts/make+exec.sh ;;
	d)
		echo "	d : traitement des données (inclut bunching)."
		Scripts/data+bunching.sh ;;

	#p)
		#echo "	p : trace les graphes."
		#Scripts/plots.sh ;;
		
	l)
		echo "	l : lois physiques."
		Scripts/loisPhy.sh ;;
		
	q) echo "	Au-revoir !";;		

	*) echo "	Je n'ai pas compris." ;;
	
esac
