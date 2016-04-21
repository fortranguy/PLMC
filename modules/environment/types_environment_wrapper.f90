module types_environment_wrapper

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_temperature, only: Abstract_Temperature
use classes_external_field, only: Abstract_External_Field
use classes_permittivity, only: Abstract_Permittivity
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_walls_potential, only: Abstract_Walls_Potential
use classes_box_size_checker, only: Abstract_Box_Size_Checker

implicit none

private

    type, public :: Environment_Wrapper
        class(Abstract_Periodic_Box), allocatable :: periodic_box
        class(Abstract_Temperature), allocatable :: temperature
        class(Abstract_Permittivity), allocatable :: permittivity
        class(Abstract_External_Field), allocatable :: external_field
        class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattice
        class(Abstract_Walls_Potential), allocatable :: walls_potential
        class(Abstract_Box_Size_Checker), allocatable :: box_size_checker
    end type Environment_Wrapper

end module types_environment_wrapper
