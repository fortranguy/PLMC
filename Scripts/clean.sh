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

echo "Copies :"
find . -name "*copy.f90"
find . -name "*copy.sh"

echo "Fit :"
ls fit.log

echo "Etes-vous sûr de vouloir effacer les fichiers ? (oui/non)"
read status
	
if [ ${status} = "oui" ]
then

	rm -rf ${binFold}/${binFold_i}*
	rm -rf ${outFold}/${outFold_i}*
	rm -rf ${tmpFold}/*.in
	find . -name "*copy.f90" -exec rm {} \;
	find . -name "*copy.sh" -exec rm {} \;
	rm fit.log
	
	echo "Les fichiers ont été effacés."
	
else
	echo "Rien n'a été effacé."
fi
