load "Plots/Solo/location.p"
outTherm="obsTherm.out"
fileTherm=outFold."/".outFold_i."/".outTherm
out="obs.out"
file=outFold."/".outFold_i."/".out

set title "Observable de run ".i
set xlabel "steps"
set ylabel "Energie"

plot file u 1:2, fileTherm u 1:2
