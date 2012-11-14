#Observables

set title "Observables"
set xlabel "steps"
set ylabel "Energie"

plot "Out/obs.out" u 1:2 #t "Energie"
