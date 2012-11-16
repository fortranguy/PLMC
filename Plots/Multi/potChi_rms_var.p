load "Plots/Multi/location.p"
out="potChi_rms_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistiques : Potentiel chimique : écart-type"
set xlabel "Ncol"
set xrange[log(80-10):log(130+10)]
set ylabel "rms(E)/Ncol"
#set yrange[:]

f1(x) = a1*x + b1
a1 = -0.5
b1 = 0.

plot file u (log($1)):(log($2)) w l, f1(x) w l

fit f1(x) file u (log($1)):(log($2)) via a1, b1

plot file u (log($1)):(log($2)) w l, f1(x) w l
