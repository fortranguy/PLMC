
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    force_add_flags(CMAKE_Fortran_FLAGS_DEBUG "-W -Wall -Wextra -fbounds-check -g")
    set(CMAKE_Fortran_FLAGS_FAST "-Ofast")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    force_add_flags(CMAKE_Fortran_FLAGS_DEBUG
        "-warn all -check all -check noarg_temp_created -fltconsistency -traceback"
    )
    set(CMAKE_Fortran_FLAGS_FAST "-ipo -O3 -static-intel -no-prec-div -xHost")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

mark_as_advanced(
    CMAKE_Fortran_FLAGS_FAST
)

