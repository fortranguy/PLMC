#! /bin/bash

find . -name "bunching_activInv.out"
find . -name "bunching_activInv.out" -exec head -n1 {} \; > erreur_Lratio.out

