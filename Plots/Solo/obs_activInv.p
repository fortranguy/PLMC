load "Plots/Solo/location.p"
out="obs.out"
file=outFold."/".outFold_i."/".out

#Observables

set title "Observable de run ".i
set xlabel "steps"
set ylabel "Inverse de l'activité (excès)"

plot file u 1:3
