load "Plots/Multi/location.p"
out="potChi_rms_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistiques : Potentiel chimique : écart-type"
set xlabel "Ncol"
set xrange[Ndeb:Nfin]
set ylabel "rms(potChi)"
#set yrange[:]

plot file u 1:2 w l
