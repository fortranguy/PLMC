load "Plots/Multi/location.p"
out="moy_eTot_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistique : Energie moyenne"
set xlabel "Ncol"
set xrange[80-10:130+10]
set ylabel "Moyenne/Ncol"
set yrange[0:1]

plot file u 1:($2/$1):($3/$1) with errorbars
