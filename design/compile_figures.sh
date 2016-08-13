#! /bin/bash


if (($# != 1))
then
    echo "Usage: compile_figures [NUM THREADS]"
    exit 1
fi

trap "wait" TERM EXIT

for tex in $(find figures/ -name "*.tex")
do
    echo $tex
    (time lualatex -output-directory=$(dirname $tex) $tex) > \
        $(dirname $tex)/$(basename $tex .tex).tmp 2>&1&
    while (($(jobs | wc -l) >= $@))
    do
        sleep 0.1
        jobs > /dev/null
    done
done
