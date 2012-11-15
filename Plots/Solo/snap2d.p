#snapshot

set title "Snapshot"
set xlabel "x"
set ylabel "y"
#set zlabel "z"

plot "snapShotIni.out" u 1:2 t "Initial", \
	"snapShotFin.out" u 1:2 t "Final"
