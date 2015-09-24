module types_box

use class_periodic_box, only: Abstract_Periodic_Box
use class_temperature, only: Abstract_Temperature
use class_field_expression, only: Abstract_Field_Expression
use class_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use class_external_field, only: Abstract_External_Field
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice

implicit none

private

    type, public :: Box_Wrapper
        class(Abstract_Periodic_Box), allocatable :: periodic_box
        class(Abstract_Temperature), allocatable :: temperature
        class(Abstract_Field_Expression), allocatable :: field_expression
        class(Abstract_Parallelepiped_Domain), allocatable :: parallelepiped_domain
        class(Abstract_External_Field), allocatable :: external_field
        class(Abstract_Reciprocal_Lattice), allocatable :: reciprocal_lattice
    end type Box_Wrapper

end module types_box
