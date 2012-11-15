# gnuplot

set title "Bunching"
set xlabel "itérations"
set ylabel "erreur apparente"

plot "bunching_eTot.out" u 1:3 w l t "Energie"
