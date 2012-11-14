#! /bin/bash

echo "Monte-Carlo Canonique : Menu"
echo "	go! : compilation et exécution."
echo "	b : bunching pour obtenir l'erreur."
echo "	p : trace les graphes."
echo "	l : lois physiques." # pompeux
echo "	q : quitter le menu."

read choix


case ${choix} in
	
	go!)
		echo "	go! : compilation et exécution."
		Scripts/make+exec.sh ;;

	b)
		echo "	b : bunching pour obtenir l'erreur."
		Scripts/bunching.sh ;;

	p)
		echo "	p : trace les graphes."
		Scripts/plots.sh ;;
		
	l)
		echo "	l : lois physiques."
		Scripts/loisPhy.sh ;;
		
	q) echo "	Au-revoir !";;		

	*) echo "	Je n'ai pas compris." ;;
	
esac
