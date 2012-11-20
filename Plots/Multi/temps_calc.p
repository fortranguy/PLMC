load "Plots/Multi/location.p"
out="temps_calc.out"
file=dataFold."/".dataFold_i."/".outFold."/".out

set title "Temps de calcul"
set xlabel "Ncol"
set xrange[Ndeb:Nfin]
set ylabel "heures"
#set yrange[:]

f1(x) = (x/a1)**n1
a1 = 100.
n1 = 1.7

plot file u 1:($2/60.) w l, "" u 1:($2/60.) w p, f1(x) w l
fit f1(x) file u 1:($2/60.) via a1, n1
plot file u 1:($2/60.) w l, "" u 1:($2/60.) w p, f1(x) w l
