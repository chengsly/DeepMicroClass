#!/bin/bash

for file in $1/*.npy
do
	echo $file
	python mem_measure.py $file
done
