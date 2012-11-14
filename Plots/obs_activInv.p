#Observables

set title "Observables"
set xlabel "steps"
set ylabel "Inverse de l'activité (excès)"

plot "Out/obs.out" u 1:3
