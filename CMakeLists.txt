cmake_minimum_required(VERSION 3.0)

project(
    Monte-Carlo
    Fortran
    CXX
)

if (NOT CMAKE_BUILD_TYPE)
    set(
        CMAKE_BUILD_TYPE Debug CACHE STRING
        "Options are: Debug Release. Default: Debug."
        FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-W -Wall -fbounds-check -g")  
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-warn all -check all -check noarg_temp_created -fltconsistency -g -traceback")  
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2.")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

set(executable_prefix mc)

    set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

add_subdirectory(mod ${CMAKE_CURRENT_BINARY_DIR}/mod)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/mod)

find_package(
    jsonfortran-${CMAKE_Fortran_COMPILER_ID}
    1.0.0
    REQUIRED
)

include_directories(${jsonfortran_INCLUDE_DIRS})

file(
    GLOB_RECURSE
    programs_source_files
    src/program_*.f90
)

foreach(program_source ${programs_source_files})
    get_filename_component(program_name ${program_source} NAME_WE)
    string(REGEX REPLACE ^program_ ""  program_name ${program_name})
    add_executable(
        ${executable_prefix}_${program_name}
        ${program_source}
    )
    set_property(
        TARGET ${executable_prefix}_${program_name}
        PROPERTY LINKER_LANGUAGE Fortran
    )
    target_link_libraries(
        ${executable_prefix}_${program_name}
        jsonfortran
        mc_meta_module
    )
endforeach(program_source)

