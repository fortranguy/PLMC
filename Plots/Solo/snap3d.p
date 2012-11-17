load "Plots/Solo/location.p"
outIni="snapShotIni.out"
outFin="snapShotFin.out"
fileIni=outFold."/".outFold_i."/".outIni
fileFin=outFold."/".outFold_i."/".outFin

set title "Snapshot"
set xlabel "x"
set ylabel "y"
#set zlabel "z"

splot fileIni t "Initial", fileFin t "Final"
