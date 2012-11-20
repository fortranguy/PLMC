load "Plots/Multi/location.p"
out="temps_calc.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Temps de calcul"
set xlabel "Ncol"
set xrange[Ndeb:Nfin]
set ylabel "Temps de calcul de minutes"
#set yrange[0.18:0.22]

plot file u 1:2 w l
