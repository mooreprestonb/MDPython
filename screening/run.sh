#!/bin/bash


for m in morse harmonic; do

	rm -rf $m
	mkdir $m
	cp *.py $m
	cp ${m}.in $m
	cp morse.init $m
	cd $m	

	for n in 1.2 1.3 1.4 1.5 1.6; do

		rm -rf $n
		mkdir $n
		cat morse.init | sed "s/XXX/$n/g" > $n/morse.init
		cp *.py $n
		cp ${m}.in $n
		cd $n
		python mdlammps.py ${m}.in
		python vac.py vel.dat vac.dat
		
		cd ..	
		echo ""
		echo "done with $n"
		echo ""	
	done
	cd ..
done

a=""; for x in $(ls harmonic/*/vac.dat); do a=$a" "$x; done; a=""; 
for x in $(ls harmonic/*/*eig.dat); do a=$a" "$x; done;  
  
xmgrace $a
