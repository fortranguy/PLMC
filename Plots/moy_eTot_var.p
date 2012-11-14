set title "Statistiques"
set xlabel "Ncol"
set ylabel "Moyenne/Ncol"

plot [60:240] "Data/Ncol_var/moy_eTot_var_sorted.out" \
	u 1:($2/$1):(100*$3/$1) with errorbars
