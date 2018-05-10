#!/bin/bash

resram inp.txt deltas.dat freqs.dat rpumps.dat
#resram inp_low.txt deltas_low.dat freqs_low.dat rpumps.dat


bash combine.sh
gnuplot plot.gnu

cd data
kde-open model.png
