project(
    Monte-CarloMetaModule
    Fortran
    CXX
)

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

set_property(
    TARGET mc_meta_module
    PROPERTY LINKER_LANGUAGE Fortran
)
target_link_libraries(
    mc_meta_module
    jsonfortran-static
)
