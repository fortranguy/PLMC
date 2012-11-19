#! /bin/bash

dataFold="Data"
scriptFold="Scripts"
editeur="gedit"

# Préparation

cd ${scriptFold}
	ls loisPhy_copy.sh
	if test $? -ne 0
	then
		echo "loisPhy.sh par défaut."
		cp loisPhy.sh loisPhy_copy.sh
	fi
cd ..

# Modification

echo Contenu de ${dataFold} :
ls -lrt ${dataFold}
echo "Modifier le script pour spécifier le chemin des données."
${editeur} ${scriptFold}/loisPhy.sh &
echo "Avez-vous spécifié le chemin et le domaine ? (oui/non)"
read answer
if [ ${answer} = "oui" ]
then
	${scriptFold}/loisPhy.sh
fi
