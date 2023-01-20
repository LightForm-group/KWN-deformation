#!/bin/bash --login

COMPILER=GNU

if [ $COMPILER == 'GNU' ]; then
export FPM_FC=gfortran
export FPM_FFLAGS="-fbounds-check -ffree-line-length-0 -fimplicit-none"
elif [ $COMPILER == 'Intel' ]; then
export FPM_FC=ifort
export FPM_FFLAGS="-traceback -check all"
fi

PREFIX_DIR=~

fpm install --prefix $PREFIX_DIR
