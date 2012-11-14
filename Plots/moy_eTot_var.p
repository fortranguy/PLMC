set title "Statistiques"
set xlabel "Ncol"
set ylabel "Moyenne/Ncol"

plot [60:240] "Data/[approx]Rho_cst_LNcol_var/Rho_cst/moy_eTot_var_sorted.out" \
	u 1:($2/$1):($3/$1) with errorbars
