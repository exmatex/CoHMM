#!/bin/bash
for i in {1..1500}
do
	echo -ne "$i\t"
	./libDriver $i 16
done
