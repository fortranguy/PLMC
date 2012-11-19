load "Plots/Multi/location.p"
out="potChi_moy_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistique : Potentiel chimque : moyenne"
set xlabel "Ncol"
set xrange[Ndeb:Nfin]
set ylabel "Moyenne"
#set yrange[-1:0]

plot file u 1:2:(2.*$3/$1) w l, "" u 1:2:(2.*$3/$1 with errorbars
