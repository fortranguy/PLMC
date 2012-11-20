load "Plots/Multi/location.p"
out="potChi_moy_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out
outHs="potChi_hs_var.out"
fileHs=dataFold."/".dataFold_i."/".outFold."/".outHs

set title "Statistique : Potentiel chimque : moyenne"
set xlabel "Ncol"
set xrange[Ndeb:Nfin]
set ylabel "Moyenne"
#set yrange[-1:0]

mu(eta) = (eta * (8 - 9*eta + 4*eta**2))/(1-eta)**3

plot file u 1:2:(2.*$3/$1) w l, "" u 1:2:(2.*$3/$1) with errorbars, \
    fileHs u 1:(mu($2) + $3) w l t "HS"
