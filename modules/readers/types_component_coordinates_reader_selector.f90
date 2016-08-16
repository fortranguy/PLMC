module types_component_coordinates_reader_selector

implicit none

private

    type, public :: Component_Coordinates_Reader_Selector
        logical :: read_positions = .false.
        logical :: read_orientations = .false.
    end type Component_Coordinates_Reader_Selector

end module types_component_coordinates_reader_selector
