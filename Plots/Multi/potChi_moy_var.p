load "Plots/Multi/location.p"
out="potChi_moy_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out
outHs="potChi_hs_var.out"
fileHs=dataFold."/".dataFold_i."/".outFold."/".outHs

id = 2
ex = 3

y=ex

set title "Statistique : Potentiel chimque : moyenne"
set xlabel "Volume"
set xrange[VolDeb-100:VolFin+100]
set ylabel "Moyenne"
#set yrange[-1:0]

mu(eta) = (eta * (8 - 9*eta + 4*eta**2))/(1-eta)**3

#plot file u 1:y w l, "" u 1:y:(2.*$4/$1) with errorbars\
    , fileHs u 1:(mu($2)) w l t "HS"
plot file u 1:($2 + $3) w l t "Total", "" u 1:($2 + $3):(2.*$4/$1) with errorbars\
    , fileHs u 1:(mu($2) + ($3)) w l t "HS"
