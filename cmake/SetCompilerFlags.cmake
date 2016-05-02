
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    force_add_flags(CMAKE_Fortran_FLAGS_DEBUG "-W -Wall -Wextra -fbounds-check -g -fbacktrace\
        -ffpe-trap=zero,overflow,underflow")
    set(CMAKE_Fortran_FLAGS_FAST "-Ofast")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    force_add_flags(CMAKE_Fortran_FLAGS_DEBUG
        "-warn all -check all -check noarg_temp_created -fltconsistency -traceback"
    )
    set(CMAKE_Fortran_FLAGS_FAST "-fast")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

mark_as_advanced(
    CMAKE_Fortran_FLAGS_FAST
)
