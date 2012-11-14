#! /bin/bash

srcFold="Sources"
binFold="Binaries"
outFold="Out"

# Préparation
ls ${srcFold} ${binFold} ${outFold}
if test $? -ne 0
then
	echo "Il manque quelque chose..."
	exit
fi

rm -f ${outFold}/*.out
rm -f ${binFold}/*.out # si problème

# Compilation
cd ${srcFold}
	make
	if test $? -ne 0
	then
		exit
	fi
cd ..

# Exécution
cd ${binFold}
	dateIni=$(date)
	time ./xMC_Canonique cube
	if test $? -ne 0
	then
		exit
	fi
	dateFin=$(date)
cd ..
	
# Données
cd ${outFold}
	mv ../${binFold}/*.out .

	echo " Dates : " >> rapport.out
	echo "     Début : " ${dateIni} >> rapport.out
	echo "     Fin : " ${dateFin} >> rapport.out	

	# Bunching
	cp ../${binFold}/xBunching .
	./xBunching 14
	rm xBunching

	# Rapport	
	cp ../Scripts/avgRms.sh .
	./avgRms.sh
	rm avgRms.sh
	cat rapport.out
cd ..
