# gnuplot

set title "Bunching"
set xlabel "itérations"
set ylabel "erreur apparente"

plot "Out/bunching_eTot.out" u 1:3 w l t "Energie"
