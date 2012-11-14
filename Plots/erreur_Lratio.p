#gnuplot

c1 = 1.e-4
f1(x) = a1*(b1/x)**3 + c1
a1 = 4e-4
b1 = 0.2


fit f1(x) "Data/Lratio_var/erreur_Lratio_sorted.out" u 1:4 via a1, b1

plot "Data/Lratio_var/erreur_Lratio_sorted.out" u 1:4 w l, f1(x) w l
