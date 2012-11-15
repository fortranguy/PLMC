load "Plots/Solo/location.p"
out="bunching_eTot.out"
file=outFold."/".outFold_i."/".out

set title "Bunching de run ".i
set xlabel "itérations"
set ylabel "erreur apparente"

plot file u 1:3 w l t "Energie"
