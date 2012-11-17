#! /bin/bash

binFold="Binaries"
outFold="Out"
tmpFold="Temp"

binFold_i="Bin"
outFold_i="Out"

echo ${binFold}":"
ls ${binFold}
echo ${outFold}":"
ls ${outFold}
echo ${tmpFold}":" 
ls ${tmpFold}

echo "Etes-vous sûr de vouloir effacer les fichiers ? (oui/non)"
read status
	
if [ ${status} = "oui" ]
then

	rm -rf ${binFold}/${binFold_i}*
	rm -rf ${outFold}/${outFold_i}*
	rm -rf ${tmpFold}/*.in
	#rm -rf ${tmpFold}/*.out
	
	echo "Les fichiers ont été effacés."
	
else
	echo "Rien n'a été effacé."
fi
