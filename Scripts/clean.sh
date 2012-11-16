#! /bin/bash

binFold="Binaries"
outFold="Out"
tmpFold="Temp"

binFold_i="Bin"
outFold_i="Out"

status="non"
until [ ${status} = "oui" ]
do

	echo ${binFold} ${outFold} ${tmpFold}
	echo "Etes-vous sûr de vouloir effacer les fichiers ? (oui/non)"
	read status

	rm -rf ${binFold}/${binFold_i}*
	rm -rf ${outFold}/${outFold_i}*
	rm -rf ${tmpFold}/*.out

done

	
