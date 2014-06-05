#!/bin/bash
for i in {1..32}
do
	valgrind --tool=massif --stacks=yes --massif-out-file=massifOut ./libDriver 10 $i	
	echo -ne "$i\t" >> blockMemOut.txt
	python getPeak.py >> blockMemOut.txt
done

