#! /bin/bash

dataFold="Data"
scriptFold="Scripts"
editeur="gedit"

# Préparation

cd ${scriptFold}
	cp loisPhy_save.sh loisPhy.sh
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
