module types_short_interactions_wrapper

use classes_hard_contact, only: Abstract_Hard_Contact
use classes_walls_visitor, only: Abstract_Walls_Visitor
use types_pair_potential_wrapper, only: Pair_Potential_Wrapper, Pair_Potentials_Line
use types_neighbour_cells_wrapper, only: Neighbour_Cells_Line
use classes_visitable_cells, only: Abstract_Visitable_Cells
use classes_short_pairs_visitor, only: Abstract_Short_Pairs_Visitor

implicit none

private

    type, public :: Short_Interactions_Wrapper
        class(Abstract_Hard_Contact), allocatable :: hard_contact
        class(Abstract_Walls_Visitor), allocatable :: walls_visitor
        type(Pair_Potential_Wrapper), allocatable :: wall_pairs(:)
        class(Abstract_Short_Pairs_Visitor), allocatable :: components_visitor
        type(Pair_Potentials_Line), allocatable :: components_pairs(:)
        type(Neighbour_Cells_Line), allocatable :: neighbour_cells(:)
        class(Abstract_Visitable_Cells), allocatable :: visitable_cells(:, :)
    end type Short_Interactions_Wrapper

end module types_short_interactions_wrapper
