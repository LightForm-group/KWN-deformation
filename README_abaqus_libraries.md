## Building ABAQUS Libraries.

Abaqus can load shared libraries created by some fortran compilers.
This is a guide on how to do this (mostly) using the `fpm` package manager.

#### Building with intel

Compilation with intel versions 17.0.7 to 19.1.2 have been tested, and the following options will create optimised binaries:
```
export FPM_FC=ifort
export FPM_FFLAGS="-traceback -check all -O3 -DWITH_QP=1 -fPIC"
fpm build
fpm test
```
Note the addition of the `-fPIC` flag compared with the standard build instructions.

Now find the build directory, e.g.:
```
ls build/*/*
```
This should return something like:
```
build/ifort_F71A7146EFEE8600/kwn_data_types.mod       build/ifort_F71A7146EFEE8600/kwn_model.mod           build/ifort_F71A7146EFEE8600/testdrive.mod
build/ifort_F71A7146EFEE8600/kwn_initialise.mod       build/ifort_F71A7146EFEE8600/kwn_model_routines.mod  build/ifort_F71A7146EFEE8600/test_suite1.mod
build/ifort_F71A7146EFEE8600/kwn_io.mod               build/ifort_F71A7146EFEE8600/kwn_parameters.mod      build/ifort_F71A7146EFEE8600/test_suite2.mod
build/ifort_F71A7146EFEE8600/kwn_model_functions.mod  build/ifort_F71A7146EFEE8600/kwn_precision.mod

build/dependencies/test-drive:
CMakeLists.txt  config  fpm.toml  LICENSE-Apache  LICENSE-MIT  meson.build  meson_options.txt  README.md  requirements.txt  src  test

build/ifort_72AB886C0D171291/KWN-Deform:
libKWN-Deform.a  libKWN-Deform.a.log

build/ifort_A384C26C54167BEF/app:
KWN-Deform  KWN-Deform.log

build/ifort_A384C26C54167BEF/test:
unit-tests  unit-tests.log

build/ifort_F71A7146EFEE8600/KWN-Deform:
app_KWN_with_deformation.f90.o                            src_KWN_initialise.f90.o.digest       src_KWN_model_functions.f90.o.log    test_tester.f90.o
app_KWN_with_deformation.f90.o.digest                     src_KWN_initialise.f90.o.log          src_KWN_model_routines.f90.o         test_tester.f90.o.digest
app_KWN_with_deformation.f90.o.log                        src_KWN_io.f90.o                      src_KWN_model_routines.f90.o.digest  test_tester.f90.o.log
build_dependencies_test-drive_src_testdrive.F90.o         src_KWN_io.f90.o.digest               src_KWN_model_routines.f90.o.log     test_test_suite1.f90.o
build_dependencies_test-drive_src_testdrive.F90.o.digest  src_KWN_io.f90.o.log                  src_KWN_parameters.f90.o             test_test_suite1.f90.o.digest
build_dependencies_test-drive_src_testdrive.F90.o.log     src_KWN_model.f90.o                   src_KWN_parameters.f90.o.digest      test_test_suite1.f90.o.log
src_KWN_data_types.f90.o                                  src_KWN_model.f90.o.digest            src_KWN_parameters.f90.o.log         test_test_suite2.f90.o
src_KWN_data_types.f90.o.digest                           src_KWN_model.f90.o.log               src_KWN_precision.f90.o              test_test_suite2.f90.o.digest
src_KWN_data_types.f90.o.log                              src_KWN_model_functions.f90.o         src_KWN_precision.f90.o.digest       test_test_suite2.f90.o.log
src_KWN_initialise.f90.o                                  src_KWN_model_functions.f90.o.digest  src_KWN_precision.f90.o.log
```
The build directory will contain object files (`*.o`), in this case it is `build/ifort_F71A7146EFEE8600/KWN-Deform`.

Change to the build directory, and create the shared object file:
```
cd build/ifort_F71A7146EFEE8600/KWN-Deform/
ifort -shared src_KWN*.o -o KWN_libraries.so
```

### Using shared libraries in ABAQUS

Copy the `src_KWN*.o` and `KWN_libraries.so` files into a subdirectory of your ABAQUS running directory.

To access the subroutines use `user=subroutine.o` (or similar) for your abaqus command line call.

Create a user environment file in your home directory, called `abaqus_v6.env`. Within this create the variable `usub_lib_dir` which should contain the full path to the subdirectory in which you've copied your `*.o` and `*.so` library files. E.g.:
```
usub_lib_dir='/net/scratch2/mbessdl2/LightForm/Abaqus_KWN_Test/KWN_library'
```
Then check your abaqus environment using `abaqus information=environment` to ensure that the `usub_lib_dir` variable exists.



