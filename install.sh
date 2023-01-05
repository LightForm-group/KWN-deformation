#!/bin/bash --login

export FPM_FC=gfortran
export FPM_FFLAGS="-fbounds-check -ffree-line-length-0 -fimplicit-none"
PREFIX_DIR=~

fpm install --prefix $PREFIX_DIR
