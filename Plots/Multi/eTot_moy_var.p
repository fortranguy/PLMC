load "Plots/Multi/location.p"
out="eTot_moy_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistique : Energie : moyenne"
set xlabel "Ncol"
set xrange[Ndeb:Nfin]
set ylabel "Moyenne/Ncol"
set yrange[0.18:0.22]

plot file u 1:($2/$1):(2.*$3/$1) with errorbars
