#

set (FSOURCES   error_mod.f90 lib_h5fortran_mod.f90 lib_nx_transfer_mod.f90
                lib_unified_chars_mod.f90 precision_mod.f90
                unified_read.f90 unified_write.f90
    )

add_library (lib_f90  ${FSOURCES})

add_dependencies (lib_f90 h5fortran_lib)
