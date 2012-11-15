#snapshot

set title "Snapshot"
set xlabel "x"
set ylabel "y"
set zlabel "z"

splot "snapShotIni.out" t "Initial", \
	"snapShotFin.out" t "Final"
