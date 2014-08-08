cmake_minimum_required(VERSION 3.0)

project(
    Monte-Carlo
    Fortran
    CXX
)

if (NOT CMAKE_BUILD_TYPE)
    set(
        CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release Fast."
        FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)

#option(CMAKE_BUILD_TYPE)

set(
    CMAKE_INSTALL_PREFIX /usr/local CACHE STRING
    "Choose the installation directory; by default it installs in /usr/local."
    FORCE
)

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

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

