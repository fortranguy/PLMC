project(
    Monte-CarloMetaModule
    Fortran
    CXX
)

set(CMAKE_Fortran_FLAGS "-fast -no-intel-extensions")
set(CMAKE_CXX_FLAGS "-O3 -ipo")

find_package(
    jsonfortran-${CMAKE_Fortran_COMPILER_ID}
    1.0.0
    REQUIRED
)

find_package(
    Boost
    1.55.0
    REQUIRED
)

file(
    GLOB_RECURSE
    modules_source_files
    src/data_*.f90
    src/module_*.f90
    src/class_*.f90
    src/*.cpp
)

add_library(
    mc_meta_module
    STATIC
    ${modules_source_files}
)

include_directories(${Boost_INCLUDE_DIRS})

include_directories(${jsonfortran_INCLUDE_DIRS})

target_link_libraries(
    mc_meta_module
    jsonfortran-static
)
