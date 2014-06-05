#!/bin/bash
for i in {1..500}
do
	valgrind --tool=massif --stacks=yes --massif-out-file=massifOut ./libDriver $i 16	
	echo -ne "$i\t" >> stepMemOut.txt
	python getPeak.py >> stepMemOut.txt
done

