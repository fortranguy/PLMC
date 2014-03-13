#! /bin/bash

for file in $@
do
    vi -c ":%s/\(\S\+\)\@<=\s\+$" -c "wq" ${file}
done

