#! /bin/bash

editeur="gedit"

echo "Modifier le script pour spécifier le chemin des données."
${editeur} Scripts/loisPhy.sh &
echo "Avez-vous spécifier le chemin et le domaine ? (oui/non)"
read answer
if [ ${answer} = "oui" ]
then
	Scripts/loisPhy.sh
fi
