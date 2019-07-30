#!bin/bash

awk 'NR == 1 || NR %3 ==0' vel.dat > vel1.dat

cat vel1.dat | awk '{print $1}' > velc1.dat

seq 0 1 20000 > tmp

paste tmp velc1.dat > image.dat

xmgrace image.dat
