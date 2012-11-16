load "Plots/Multi/location.p"
out="potChi_moy_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistique : Potentiel chimque : moyenne"
set xlabel "Ncol"
set xrange[80-10:130+10]
set ylabel "Moyenne/Ncol"
set yrange[1.285:1.295]

plot file u 1:2:(2.*$3/$1) with errorbars
