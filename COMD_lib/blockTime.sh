#!/bin/bash
for i in {1..32}
do
	echo -ne "$i\t"
	./libDriver 10 $i
done
