#!/bin/bash


# the python code just reads the inputs in the yaml file "input.yaml" and geenrates a textfile ('input.dat') with the same inputs formatted to be readable for the main program a.out
python generate_input_fortran.py
# this is the main program
./a.out