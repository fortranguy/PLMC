module types_component_coordinates_writer_selector

implicit none

private

    type, public :: Component_Coordinates_Writer_Selector
        logical :: write_positions = .false.
        logical :: write_orientations = .false.
    end type Component_Coordinates_Writer_Selector

end module types_component_coordinates_writer_selector
