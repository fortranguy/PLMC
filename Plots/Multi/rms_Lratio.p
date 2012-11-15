#gnuplot

f1(x) = a1/x**(0.9) + b1
a1 = 5e-2
b1 = -4e-2

fit f1(x) "Data/Lratio_var/rms_Lratio_sorted.out" u 1:2 via a1, b1

plot "Data/Lratio_var/rms_Lratio_sorted.out" u 1:2 w l, f1(x) w l
