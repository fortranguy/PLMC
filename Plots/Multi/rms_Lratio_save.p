#gnuplot

c1 = -2.3e-2
f1(x) = a1*(b1/x)**1.2 + c1
a1 = 1e-1
b1 = 4.e-1

fit f1(x) "Data/Lratio_var/rms_Lratio_sorted.out" u 1:2 via a1, b1

plot "Data/Lratio_var/rms_Lratio_sorted.out" u 1:2 w l, f1(x) w l
