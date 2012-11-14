#! /bin/bash

gprof xMC_Canonique | Profile/gprof2dot.py | dot -Tpdf -o Profile/output.pdf