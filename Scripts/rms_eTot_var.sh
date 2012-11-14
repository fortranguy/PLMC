#! /bin/bash

listOut=$(find . -name "Out" -type d)
NcolFil="NcolFil.out"
eRmsFil="eRmsFil.out"

rm -f ${NcolFil} ${eRmsFil} rms_eTot_var*.out

for out in ${listOut}

do
	echo ${out}
	cat ${out}"/rapport.out" | head -n4 | tail -n1 > ${NcolFil} # fragile
	cat ${out}"/rapport.out" | tail -n4 | head -n1 > ${eRmsFil} # fragile

	Ncol=$(cut -d = -f 2 ${NcolFil})
	eRms=$(cut -d = -f 2 ${eRmsFil})
	
	echo ${Ncol}"	"${eRms} >> rms_eTot_var.out

done

sort -n rms_eTot_var.out > rms_eTot_var_sorted.out
