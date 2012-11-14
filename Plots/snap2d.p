#snapshot

set title "Snapshot"
set xlabel "x"
set ylabel "y"
#set zlabel "z"

plot "Out/snapShotIni.out" u 1:2 t "Initial", \
	"Out/snapShotFin.out" u 1:2 t "Final"
