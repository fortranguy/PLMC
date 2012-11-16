load "Plots/Solo/location.p"
outIni="snapShotIni.out"
outFin="snapShotFin.out"
fileIni=outFold."/".outFold_i."/".outIni
fileFin=outFold."/".outFold_i."/".outFin

set title "Snapshot"
set xlabel "x"
set ylabel "y"
#set zlabel "z"

plot fileIni u 1:2 t "Initial", fileFin u 1:2 t "Final"
