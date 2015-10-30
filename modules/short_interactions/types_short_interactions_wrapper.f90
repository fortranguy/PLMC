module types_short_interactions_wrapper

use class_pair_potential, only: Abstract_Pair_Potential
use class_visitable_cells, only: Abstract_Visitable_Cells
use class_walls_potential_visitor, only: Abstract_Walls_Potential_Visitor
use class_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor

implicit none

private

    type, public :: Pair_Potential_Wrapper
        class(Abstract_Pair_Potential), allocatable :: pair_potential
    end type Pair_Potential_Wrapper

    type, public :: Pair_Potentials_Wrapper
        type(Pair_Potential_Wrapper), allocatable :: with_components(:)
    end type Pair_Potentials_Wrapper

    type, public :: Short_Interactions_Wrapper
        class(Abstract_Short_Pairs_Visitor), allocatable :: components_visitor
        class(Abstract_Walls_Potential_Visitor), allocatable :: walls_visitor
        type(Pair_Potentials_Wrapper), allocatable :: components_pairs(:)
        class(Abstract_Visitable_Cells), allocatable :: components_cells(:, :)
        type(Pair_Potential_Wrapper), allocatable :: wall_pairs(:)
    end type Short_Interactions_Wrapper

end module types_short_interactions_wrapper