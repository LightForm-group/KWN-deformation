

# KWN model with external deformation considered via the amount of excess vacancies

KWN precipitation model including the effect of deformation via excess vacancy concentration calculation (phenomenological law for vacancy production). 

The model is described in detail in reference [3].  

## Requirements

Compiler: gfortran or ifort

Installer: [Fortran Package Manager](https://fpm.fortran-lang.org/en/index.html)

## Building with fpm

The FPM settings are defined in the file `fpm.toml`. The instructions for building and testing KWN-deform are listed below.

Tests are run using the [Test Drive](https://github.com/awvwgk/test-drive) package, as listed in the `[dev-dependencies]` section of the `fpm.toml` file. This requires the `-DWITH_QP=1` flag to be included in `FPM_FFLAGS`; if you are not running the tests this flag can be left out.

### Building with gfortran

Compilation with gfortran versions 6.4.0 to 11.3.0 has been tested, and the following options will create optimised binaries: 
```
export FPM_FC=gfortran
export FPM_FFLAGS="-fbounds-check -ffree-line-length-0 -fimplicit-none -O3 -DWITH_QP=1"
fpm build
fpm test
```

### Building with intel

Compilation with intel versions 17.0.7 to 19.1.2 have been tested, and the following options will create optimised binaries:
```
export FPM_FC=ifort
export FPM_FFLAGS="-traceback -check all -O3 -DWITH_QP=1"
fpm build
fpm test
```

## Installing with fpm

The model can be installed using the `fpm install` command. This will install the binary and libraries (in `bin` and `lib`/`include` respectively). For a standard user these will be installed with the root `~/.local/`, but to change this root you can use the `--prefix` flag:
```
fpm install --prefix <prefix directory>
```

## Running the Model

1. Fill or modify the "input.yaml" file with input corresponding to the model described in ref. [3].
2. Put the "input.yaml" file in a folder (e.g. "test_folder_name")
3. Run `./run_kwn.sh test_folder_name` in a terminal.  
4. The outputs are written in textfiles and can be visualised using the attached Jupyter notebooks.  

Some examples of input files and jupyer notebooks can be found in the "test_n" directories, with n the number of the test: 

- with deformation  in a ternary alloy (Al-Zn-Mg) containing an initial distribution
  - test_1
- without deformation in a binary alloy with no initial distribution
  - for a Cu-Co binary alloy (reproduces result of ref [7]). 

## Authors

Madeleine Bignon, University of Manchester  
Pratheek Shanthraj, University of Manchester  
Samuel Engel, University of Manchester

## Copyright & Licensing

This software has been developed as part of the [Lightform project](http://lightform.org.uk/) at the University of Manchester. 

- Copyright 2023 [University of Manchester](https://www.manchester.ac.uk/)

Licensed under the MIT license, see the file License.md for details. 

## References

[1] Robson, J. D. (2020). Metallurgical and Materials Transactions A: Physical Metallurgy and Materials Science, 51(10), 5401–5413. https://doi.org/10.1007/s11661-020-05960-5  
[2] Deschamps, A., Livet, F., & Bréchet, Y. (1998). . Acta Materialia, 47(1), 281–292. https://doi.org/10.1016/S1359-6454(98)00293-6  
[3] Bignon, M., Shanthraj, P., & Robson, J. D. (2022). Acta Materialia, 234, 118036. https://doi.org/10.1016/J.ACTAMAT.2022.118036  
[4] Deschamps, A., & De Geuser, F. (2011). Journal of Applied Crystallography, 44(2), 343–352. https://doi.org/10.1107/S0021889811003049  
[5] Perez, M. (2005). Scripta Materialia, 52(8), 709–712. https://doi.org/10.1016/j.scriptamat.2004.12.026  
[6] Nicolas, M., & Deschamps, A. (2003). Acta Materialia, 51(20), 6077–6094. https://doi.org/10.1016/S1359-6454(03)00429-4  
[7] Robson, J. D. (2004). Modelling the evolution of particle size distribution during nucleation, growth and coarsening. Materials Science and Technology, 20(4), 441–448. https://doi.org/10.1179/026708304225016725
