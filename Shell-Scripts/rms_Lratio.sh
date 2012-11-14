#! /bin/bash

find . -name "rapport.out"
find . -name "rapport.out" -exec tail -n1 {} \; > rms_Lratio.out

