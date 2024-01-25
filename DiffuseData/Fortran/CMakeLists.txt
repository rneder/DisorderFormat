# Toplevel for Diffuse Data

cmake_minimum_required(VERSION 3.12.0)
project(DIFFUSE_DATA)
enable_language(Fortran)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (DIFFUSE_DATA_VERSION_MAJOR 0)
set (DIFFUSE_DATA_VERSION_MINOR 0)
set (DIFFUSE_DATA_VERSION_PATCH 0)

find_package(HDF5 COMPONENTS Fortran REQUIRED)
set (LINK_LIBS ${LINK_LIBS} "${HDF5_LIBRARIES}" )
INCLUDE_DIRECTORIES ("${HDF5_INCLUDE_DIRS}")

set (CMAKE_Fortran_FLAGS "-fPIC -cpp -O3 -ffast-math -std=legacy -ffree-line-length-400")
#set (CMAKE_Fortran_FLAGS "-fPIC -cpp -fcheck=all -g -Og -fbacktrace -Warray-bounds -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wunused-parameter -Wunused-variable -Wunused-function -Wunused-value -fwhole-file -fimplicit-none -pedantic -Wuninitialized -ffree-line-length-400")

set (SOURCE_MAIN main_write.f90)
set (SOURCES_ALL precision_mod.f90 prompt_mod.f90 top_data_mod.f90 hdf_write.f90 
     hdf5_def.f90 hdf5_params.f90 hdf5_read.f90)

add_executable(write_main ${SOURCE_MAIN} ${SOURCES_ALL})
target_link_libraries(write_main ${LINK_LIBS})
