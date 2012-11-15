load "location.p"
out="moy_eTot_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistiques"
set xlabel "Ncol"
set xrange[80-10:130+10]
set ylabel "Moyenne/Ncol"

plot file u 1:($2/$1):($3/$1) with errorbars
