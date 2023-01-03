#!/bin/bash
#takes as argument the name of the folder where 'input.yaml' is located
# remove previous results
rm -r $1/results
# create a new results file
mkdir $1/results
# the python code just reads the inputs in the yaml file "input.yaml" and geenrates a textfile ('input.dat') with the same inputs formatted to be readable for the main program a.out
python generate_input_fortran.py $1
# this is the main program
./a.out
# move the results
mv *.txt $1/results
