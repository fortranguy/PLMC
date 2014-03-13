#! /bin/bash

for file in $@
do
    vi -es -c ":%s/\(\S\+\)\@<=\s\+$" -c "wq" ${file}
done

