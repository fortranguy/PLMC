
message("Compiler version: ${CMAKE_Fortran_COMPILER_VERSION}")

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS "4.9.1")
        message(FATAL_ERROR "Insufficient gfortran version")
    endif ()
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS "12.1.0")
        message(FATAL_ERROR "Insufficient ifort version")
    endif ()
endif ()
