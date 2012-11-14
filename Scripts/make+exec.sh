#! /bin/bash

srcFold="Sources"
binFold="Binaries"
outFold="Out"

# Préparation
ls ${binFold}
if test $? -ne 0
then
	mkdir ${binFold}
fi

ls ${outFold}
if test $? -ne 0
then
	mkdir ${outFold}
else
	rm -f ${outFold}/*.out
fi

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
	echo
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
