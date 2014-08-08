
if (NOT CMAKE_BUILD_TYPE)
    set(build_type_message "Options are: Debug, Release or Fast. Default: Debug.")
    message(STATUS "Build type: ${build_type_message}")
    set(
        CMAKE_BUILD_TYPE "Debug" CACHE STRING
        ${build_type_message}
        FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)
set_property(
    CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS "Debug" "Release" "Fast"
)

