cmake_minimum_required(VERSION 3.15)
enable_language(Fortran)
project(KWNDeformation VERSION 0.1)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(CMAKE_Fortran_COMPILER MATCHES "gfortran*")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0 -fimplicit-none")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-Wall -O0 -g3 -fbounds-check")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
endif()

set(TARGET KWNDeformation)
set(SOURCE_FILES src/KWN_with_deformation.f90)

add_executable(${TARGET} ${SOURCE_FILES})



enable_testing()

