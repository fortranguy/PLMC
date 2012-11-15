# gnuplot

set title "Bunching de l'inverse de l'activité en excès"
set xlabel "itérations"
set ylabel "erreur apparente"

plot \
"Data/Lratio_var/0.5/Out/bunching_activInv.out" u 1:3 w l t "Lratio=0.5", \
"Data/Lratio_var/0.6/Out/bunching_activInv.out" u 1:3 w l t "Lratio=0.6", \
"Data/Lratio_var/0.7/Out/bunching_activInv.out" u 1:3 w l t "Lratio=0.7", \
"Data/Lratio_var/0.8/Out/bunching_activInv.out" u 1:3 w l t "Lratio=0.8", \
"Data/Lratio_var/0.9/Out/bunching_activInv.out" u 1:3 w l t "Lratio=0.9", \
"Data/Lratio_var/1/Out/bunching_activInv.out" u 1:3 w l t "Lratio=1"
