module types_box

use class_periodic_box, only: Abstract_Periodic_Box
use class_temperature, only: Abstract_Temperature
use class_field_expression, only: Abstract_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_floor_penetration, only: Abstract_Floor_Penetration
use class_particles_diameter, only: Abstract_Particles_Diameter
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use class_walls_potential, only: Abstract_Walls_Potential

implicit none

private

    type, public :: Box_Wrapper
        class(Abstract_Periodic_Box), allocatable :: periodic_box
        class(Abstract_Temperature), allocatable :: temperature
        class(Abstract_Field_Expression), allocatable :: field_expression
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        class(Abstract_External_Field), allocatable :: external_field
        class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattice
        class(Abstract_Floor_Penetration), allocatable :: floor_penetration
        class(Abstract_Particles_Diameter), allocatable :: wall_diameter
        class(Abstract_Potential_Expression), allocatable :: wall_expression
        class(Abstract_Pair_Potential), allocatable :: wall_pair
        class(Abstract_Walls_Potential), allocatable :: walls_potential
    end type Box_Wrapper

end module types_box
