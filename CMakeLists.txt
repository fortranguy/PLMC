cmake_minimum_required(VERSION 2.8)

project(
    Monte-Carlo
    Fortran
    CXX
)

set(executable_prefix "mc")

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})

find_package(
    jsonfortran-${CMAKE_Fortran_COMPILER_ID}
    1.0.0
    REQUIRED
)

file(
    GLOB_RECURSE
    programs_source_files
    "src/program_*.f90"
)

message(
    STATUS
    programs_source_files
)

file(
    GLOB_RECURSE
    modules_source_files
    "src/data_*.f90"
    "src/module_*.f90"
    "src/class_*.f90"
    "src/*.cpp"
)

foreach(program_source ${programs_source_files})
    get_filename_component(program_name ${program_source} NAME_WE)
    string(REGEX REPLACE "^program_" "" program_name ${program_name})
    add_executable(
        ${executable_prefix}_${program_name}
        ${program_source}
        ${modules_source_files}
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
    )
endforeach(program_source)

