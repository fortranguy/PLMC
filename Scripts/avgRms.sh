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
	printf("	 eTot.moy = %.15g\n", avg)
		
	for (i=0; i<NR; i++)
		deltaSqr += x[i]*x[i]
	deltaSqr = deltaSqr/NR - avg*avg
	printf("	 eTot.rms = %.15g\n" , sqrt(deltaSqr))
	
}' obs.out >> rapport.out

# Potentiel chimique / Température
LC_ALL=C awk '{
	x[i++] = $3
} 
END{

	for (i=0; i<NR; i++)
		avg += x[i]
	avg /= NR
	printf("	 potChiEx.moy = %.15g\n ", -log(avg))
		
	for (i=0; i<NR; i++)
		deltaSqr += x[i]*x[i]
	deltaSqr = deltaSqr/NR - avg*avg
	printf("	 potChiEx.rms = %.15g\n ", sqrt(deltaSqr) / (-log(avg)))
	
}' obs.out >> rapport.out
