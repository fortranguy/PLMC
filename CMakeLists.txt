cmake_minimum_required(VERSION 2.8)

project(
    Monte-Carlo
    Fortran
)

set(executable_prefix mc)

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

add_subdirectory(mod ${CMAKE_CURRENT_BINARY_DIR}/mod)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/mod)

find_package(
    jsonfortran-${CMAKE_Fortran_COMPILER_ID}
    1.0.0
    REQUIRED
)

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
    target_include_directories(
        ${executable_prefix}_${program_name}
        BEFORE
        PUBLIC
        ${jsonfortran_INCLUDE_DIRS}
    )
    target_link_libraries(
        ${executable_prefix}_${program_name}
        jsonfortran
        mc_meta_module
    )
endforeach(program_source)
