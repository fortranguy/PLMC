load "Plots/Multi/location.p"
out="potChi_rms_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistiques : Potentiel chimique : écart-type"
set xlabel "Ncol"
set xrange[Ndeb:Nfin]
set ylabel "rms(E)"
#set yrange[:]

f1(x) = a1*x + b1
a1 = -0.5
b1 = 0.

#plot file u (log($1)):(log($2)) w l, f1(x) w l
#fit f1(x) file u (log($1)):(log($2)) via a1, b1
	# pas d'a priori ?

plot file u 1:2 w l#, f1(x) w l
