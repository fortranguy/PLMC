#snapshot

set title "Snapshot"
set xlabel "x"
set ylabel "y"
#set zlabel "z"

splot "Out/snapShotIni.out" t "Initial", \
	"Out/snapShotFin.out" t "Final"
