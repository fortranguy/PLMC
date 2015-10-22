module types_short_potentials_wrapper

use data_constants, only: num_components
use class_potential_expression, only: Abstract_Potential_Expression
use class_pair_potential, only: Abstract_Pair_Potential
use class_visitable_list, only: Abstract_Visitable_List
use class_visitable_cells, only: Abstract_Visitable_Cells
use class_walls_potential_visitor, only: Abstract_Walls_Potential_Visitor
use class_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor

implicit none

private

    type, public :: Short_Potential_Wrapper
        class(Abstract_Pair_Potential), allocatable :: pair
        class(Abstract_Pair_Potential), allocatable :: wall_pair
        class(Abstract_Visitable_Cells), allocatable :: cells
    end type Short_Potential_Wrapper

    type, public :: Short_Potential_Macro_Wrapper
        class(Abstract_Visitable_Cells), allocatable :: cells
    end type Short_Potential_Macro_Wrapper

    type, public :: Mixture_Short_Potentials_Wrapper
        type(Short_Potential_Wrapper) :: intras(num_components) !+ interaction with walls
        type(Short_Potential_Macro_Wrapper) :: inters(num_components)
        class(Abstract_Pair_Potential), allocatable :: inter_pair
    end type Mixture_Short_Potentials_Wrapper

    type, public :: Pair_Potential_Wrapper
        class(Abstract_Pair_Potential), allocatable :: pair_potential
    end type Pair_Potential_Wrapper

    type, public :: Pair_Potentials_Wrapper
        type(Pair_Potential_Wrapper), allocatable :: with_components(:)
    end type Pair_Potentials_Wrapper

    type, public :: Short_Potentials_Wrapper
        class(Abstract_Short_Pairs_Visitor), allocatable :: inter_visitor
        class(Abstract_Walls_Potential_Visitor), allocatable :: walls_visitor
        type(Pair_Potentials_Wrapper), allocatable :: inter_pairs(:)
        class(Abstract_Visitable_Cells), allocatable :: inter_cells(:, :)
        type(Pair_Potential_Wrapper), allocatable :: wall_pairs(:)
    end type Short_Potentials_Wrapper

end module types_short_potentials_wrapper
