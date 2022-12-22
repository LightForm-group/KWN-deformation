#!/bin/bash
# remove previous results
#rm -r results
# create a new results file
mkdir results
# the python code just reads the inputs in the yaml file "input.yaml" and geenrates a textfile ('input.dat') with the same inputs formatted to be readable for the main program a.out
python3 generate_input_fortran.py
# this is the main program
./a.out
# move the results
mv *.txt results
