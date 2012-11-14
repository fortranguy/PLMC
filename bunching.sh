#! /bin/bash

outFold="Out"
nBunching=14

ls ${outFold}/*.out
if test $? -ne 0
then
	echo "Pas de données."
	exit
fi

cd ${outFold}

	# Bunching
	cp ../${binFold}/xBunching .
	./xBunching ${nBunching}
	rm xBunching
	
cd ..
