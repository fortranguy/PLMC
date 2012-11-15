#! /bin/bash

binFold="Binaries"
outFold="Out"
tmpFold="Temp"

status="non"
until [ ${status} = "oui" ]
do

	echo ${binFold} ${outFold} ${tmpFold}
	echo "Etes-vous sûr de vouloir effacer les fichiers ? (oui/non)"
	read status

	rm -rf ${binFold}/*
	rm -rf ${outFold}/*
	rm -rf ${tmpFold}/*

done

	
