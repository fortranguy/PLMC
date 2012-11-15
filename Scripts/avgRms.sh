#! /bin/bash

echo " Statistiques :" >> rapport.out

# Energie
LC_ALL=C awk '{
	x[i++] = $2
} 
END{

	for (i=0; i<NR; i++)
		avg += x[i]
	avg /= NR
	print "	Energ.moy = ", avg
		
	for (i=0; i<NR; i++)
		deltaSqr += x[i]*x[i]
	deltaSqr = deltaSqr/NR - avg*avg
	print "	Energ.rms = " , sqrt(deltaSqr)
	
}' obs.out >> rapport.out

# Potentiel chimique / Température
LC_ALL=C awk '{
	x[i++] = $3
} 
END{

	for (i=0; i<NR; i++)
		avg += x[i]
	avg /= NR
	print "	potChi.moy = ", -log(avg)
		
	for (i=0; i<NR; i++)
		deltaSqr += x[i]*x[i]
	deltaSqr = deltaSqr/NR - avg*avg
	print "	potChi.rms = ", sqrt(deltaSqr) / (-log(avg))
	
}' obs.out >> rapport.out
