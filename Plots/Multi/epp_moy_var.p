load "Plots/Multi/location.p"
out="epp_moy_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistique : Energie"
set xlabel "Volume"
set xrange[VolDeb-100:VolFin+100]
set ylabel "Energie par particules"
#set yrange[0.18:0.22]

plot file u 1:2:(2*$3) with errorlines
