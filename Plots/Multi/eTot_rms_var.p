load "Plots/Multi/location.p"
out="eTot_rms_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistiques : Energie : écart-type"
set xlabel "Ncol"
set xrange[log(Ndeb):log(Nfin)]
set ylabel "rms(E)/Ncol"
#set yrange[:]

f1(x) = a1*x + b1
a1 = -0.5
b1 = -1.

plot file u (log($1)):(log($2/$1)) w l, f1(x) w l
fit f1(x) file u (log($1)):(log($2/$1)) via a1, b1
plot file u (log($1)):(log($2/$1)) w l, f1(x) w l
