#! /bin/bash

echo " Statistiques :" >> rapport.out

echo "     Energie" >> rapport.out
LC_ALL=C awk '{
	x[i++] = $2
} 
END{

	for (i=0; i<NR; i++)
		avg += x[i]
	avg /= NR
	print "         moyenne = ", avg
		
	for (i=0; i<NR; i++)
		deltaSqr += x[i]*x[i]
	deltaSqr = deltaSqr/NR - avg*avg
	print "         écart-type = " , sqrt(deltaSqr)
	
}' obs.out >> rapport.out

echo "     Potentiel chimique / Température" >> rapport.out
LC_ALL=C awk '{
	x[i++] = $3
} 
END{

	for (i=0; i<NR; i++)
		avg += x[i]
	avg /= NR
	print "         moyenne = ", -log(avg)
		
	for (i=0; i<NR; i++)
		deltaSqr += x[i]*x[i]
	deltaSqr = deltaSqr/NR - avg*avg
	print "         écart-type = ", sqrt(deltaSqr) / (-log(avg))
	
}' obs.out >> rapport.out
