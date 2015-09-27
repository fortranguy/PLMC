module types_short_potential

use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use class_particles_potential, only: Abstract_Particles_Potential
use class_visitable_list, only: Abstract_Visitable_List
use class_visitable_cells, only: Abstract_Visitable_Cells

implicit none

private

    type, public :: Short_Potential_Wrapper
        class(Abstract_Potential_Expression), allocatable :: expression
        class(Abstract_Pair_Potential), allocatable :: pair
        class(Abstract_Particles_Potential), allocatable :: particles
        class(Abstract_Visitable_List), allocatable :: list
        class(Abstract_Visitable_Cells), allocatable :: cells
    end type Short_Potential_Wrapper

    type, public :: Short_Potential_Macro_Wrapper
        class(Abstract_Particles_Potential), allocatable :: particles
        class(Abstract_Visitable_List), allocatable :: list
        class(Abstract_Visitable_Cells), allocatable :: cells
    end type Short_Potential_Macro_Wrapper

    type, public :: Short_Potential_Micro_Wrapper
        class(Abstract_Potential_Expression), allocatable :: expression
        class(Abstract_Pair_Potential), allocatable :: pair
    end type Short_Potential_Micro_Wrapper

    type, public :: Mixture_Short_Potentials_Wrapper
        type(Short_Potential_Wrapper) :: intras(2)
        type(Short_Potential_Macro_Wrapper) :: inter_macros(2)
        type(Short_Potential_Micro_Wrapper) :: inter_micro
    end type Mixture_Short_Potentials_Wrapper

end module types_short_potential
