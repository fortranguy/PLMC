module types_environment_wrapper

use class_periodic_box, only: Abstract_Periodic_Box
use class_temperature, only: Abstract_Temperature
use class_field_expression, only: Abstract_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_floor_penetration, only: Abstract_Floor_Penetration
use class_walls_potential, only: Abstract_Walls_Potential

implicit none

private

    type, public :: Environment_Wrapper
        class(Abstract_Periodic_Box), allocatable :: periodic_box
        class(Abstract_Temperature), allocatable :: temperature
        class(Abstract_External_Field), allocatable :: external_field
        class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattice
        class(Abstract_Walls_Potential), allocatable :: walls_potential
    end type Environment_Wrapper

end module types_environment_wrapper
