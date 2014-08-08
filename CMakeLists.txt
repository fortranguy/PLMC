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
    set(CMAKE_Fortran_FLAGS_DEBUG   "-W -Wall -fbounds-check -g")  
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-warn all -check all -check noarg_temp_created -fltconsistency -g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "ipo -O3 -static-intel -no-prec-div -xHost")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG   "-W -Wall -g") 
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
elseif (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_DEBUG   "-W -Wall -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -static-intel -no-prec-div -xHost")
endif (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")

find_program(XIAR xiar)
if(XIAR)
  set(CMAKE_AR "${XIAR}")
endif(XIAR)
mark_as_advanced(XIAR)

find_program(XILD xild)
if(XILD)
  set(CMAKE_LINKER "${XILD}")
endif(XILD)
mark_as_advanced(XILD)

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
        jsonfortran-static
        mc_meta_module
    )
endforeach(program_source)

