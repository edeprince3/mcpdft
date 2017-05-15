#!/bin/bash

set term post enhanced color eps "" 20
set out "time_step.eps"

set xlabel "Energy [eV]"
set ylabel "Intensity [a.u.]"
plot "0.005.K" u 1:2 w l lw 6 lt 1 lc 0 t "100 a.u by 0.005",\
     "0.0005.K" u 1:($2/10) w l lw 4 lc 1 lt 1 t "10 a.u by 0.0005",\
     "0.00005.K" u 1:($2/10) w l lw 2 lc 2 lt 1 t "1 a.u by 0.00005"
