load "Plots/Multi/location.p"
out="potChi_rms_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistiques : Potentiel chimique : écart-type"
set xlabel "Volume"
set xrange[VolDeb-100:VolFin+100]
set ylabel "rms(potChi)"
#set yrange[:]

plot file u 1:2 w l
