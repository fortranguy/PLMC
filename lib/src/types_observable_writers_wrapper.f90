module types_observable_writers_wrapper

use data_constants, only: num_components
use class_component_coordinates_writer, only: Abstract_Component_Coordinates_Writer
use class_particles_energy_writer, only: Abstract_Particles_Energy_Writer
use class_inter_energy_writer, only: Abstract_Inter_Energy_Writer
use class_changes_writer, only: Abstract_Changes_Success_Writer

implicit none

private

    type :: Observable_Writers_Wrapper
        class(Abstract_Component_Coordinates_Writer), allocatable :: coordinates
        class(Abstract_Particles_Energy_Writer), allocatable :: energy
        class(Abstract_Changes_Success_Writer), allocatable :: changes
    end type Observable_Writers_Wrapper

    type, public :: Mixture_Observable_Writers_Wrapper
        type(Observable_Writers_Wrapper) :: intras(num_components)
        class(Abstract_Inter_Energy_Writer), allocatable :: inter_energy
    end type Mixture_Observable_Writers_Wrapper

end module types_observable_writers_wrapper
