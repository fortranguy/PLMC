
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

