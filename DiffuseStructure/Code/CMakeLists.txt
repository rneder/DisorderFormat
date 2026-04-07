# Top level for use of UNIFIED with h5fortran

cmake_minimum_required(VERSION 3.20)

project(UNIFIED)
#
enable_language(Fortran)
set (CMAKE_Fortran_COMPILER $ENV{FC})

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(HDF5 COMPONENTS Fortran)

include_directories (${UNIFIED_BINARY_DIR}/h5fortran)
include_directories (${UNIFIED_BINARY_DIR}/lib_f90)
include_directories ("${HDF5_INCLUDE_DIRS}")

set (LIBS  lib_f90 h5fortran_lib "${HDF5_LIBRARIES}" )

set (CMAKE_Fortran_FLAGS "-fPIC -cpp -fcheck=all -g -Og -fbacktrace -Wall -Wconversion -Warray-bounds -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter -Wunused-variable -Wunused-function -Wunused-value -Wno-unused-dummy-argument -fwhole-file -fimplicit-none -pedantic -Wuninitialized -ffree-line-length-400")

set (SOURCE_WRITE  unified_top_write.f90)
set (SOURCE_READ   unified_top_read.f90)

add_subdirectory(h5fortran)
add_subdirectory(lib_f90)

add_executable(unified_top_read ${SOURCE_READ})
target_link_libraries (unified_top_read ${LIBS} )

add_executable(unified_top_write ${SOURCE_WRITE})
target_link_libraries (unified_top_write ${LIBS} )
