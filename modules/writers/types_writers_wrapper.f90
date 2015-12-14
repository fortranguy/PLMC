module types_writers_wrapper

use class_component_coordinates_writer, only: Abstract_Coordinates_Writer
use class_components_energes_writer, only: Abstract_Components_Energies_Writer
use class_energy_writer, only: Abstract_Energy_Writer
use class_changes_writer, only: Abstract_Changes_Success_Writer

implicit none

private

    type, public :: Component_Writers_Wrapper
        class(Abstract_Coordinates_Writer), allocatable :: coordinates
        class(Abstract_Changes_Success_Writer), allocatable :: changes
    end type Component_Writers_Wrapper

    type, public :: Writers_Wrapper
        type(Component_Writers_Wrapper), allocatable :: components(:)
        class(Abstract_Components_Energies_Writer), allocatable :: short_energies, &
            long_energies_wo_reci
        class(Abstract_Energy_Writer), allocatable :: reci_energy
    end type Writers_Wrapper

end module types_writers_wrapper
