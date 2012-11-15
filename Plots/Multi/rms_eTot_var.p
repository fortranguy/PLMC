load "Plots/Multi/location.p"
out="rms_eTot_var.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Statistiques : écart-type de lénergie"
set xlabel "Ncol"
set xrange[log(80-10):log(130+10)]
set ylabel "rms(E)/Ncol"
#set yrange[:]

f1(x) = a1*x + b1
a1 = -1.
b1 = 6.

fit f1(x) file u (log($1)):(log($2/$1)) via a1, b1

plot file u (log($1)):(log($2/$1)) w l, f1(x) w l
