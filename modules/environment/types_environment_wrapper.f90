module types_environment_wrapper

use classes_periodic_box, only: Abstract_Periodic_Box
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_beta_pressure, only: Abstract_Beta_Pressure
use classes_temperature, only: Abstract_Temperature
use classes_external_field, only: Abstract_External_Field
use classes_permittivity, only: Abstract_Permittivity
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use classes_visitable_walls, only: Abstract_Visitable_Walls
use classes_box_size_checker, only: Abstract_Box_Size_Checker
use classes_min_distance, only: Abstract_Min_Distance

implicit none

private

    !> @todo
    !> Remove gemc_ and the second part.
    type, public :: Environment_Wrapper
        class(Abstract_Periodic_Box), allocatable :: periodic_boxes(:)
        class(Abstract_Beta_Pressure), allocatable :: beta_pressure
        class(Abstract_Temperature), allocatable :: temperature
        class(Abstract_Permittivity), allocatable :: permittivity
        class(Abstract_External_Field), allocatable :: external_fields(:)
        class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattices(:)
        class(Abstract_Min_Distance), allocatable :: wall_min_distance
        class(Abstract_Visitable_Walls), allocatable :: gemc_visitable_walls(:)
        class(Abstract_Box_Size_Checker), allocatable :: boxes_size_checker(:)
        class(Abstract_Parallelepiped_Domain), allocatable :: accessible_domains(:)
            !!for a point particle

        class(Abstract_Periodic_Box), allocatable :: periodic_box
        class(Abstract_External_Field), allocatable :: external_field
        class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattice
        class(Abstract_Visitable_Walls), allocatable :: visitable_walls
        class(Abstract_Box_Size_Checker), allocatable :: box_size_checker
        class(Abstract_Parallelepiped_Domain), allocatable :: accessible_domain
            !!for a point particle
    end type Environment_Wrapper

end module types_environment_wrapper
