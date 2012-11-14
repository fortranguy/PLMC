set title "Statistiques"
set xlabel "Ncol"
set ylabel "Moyenne/Ncol"

set logscale xy # why u no e ?
show logscale

plot [60:240] "Data/Ncol_var/rms_eTot_var_sorted.out" u 1:($2/$1) w l

unset logscale
