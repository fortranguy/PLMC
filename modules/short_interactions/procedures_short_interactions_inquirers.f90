module procedures_short_interactions_inquirers

use classes_min_distance, only: Abstract_Min_Distance, Concrete_Min_Distance
use classes_pair_potential, only: Abstract_Pair_Potential, Null_Pair_Potential
use classes_visitable_list, only: Abstract_Visitable_List, Concrete_Visitable_List, &
    Concrete_Visitable_Array

implicit none

private
public :: components_interact, component_interacts_with_wall, list_is_linked_list, list_is_array

interface components_interact
    module procedure :: components_interact_from_min_distance
    module procedure :: components_interact_from_pair_potential
end interface components_interact

contains

        pure logical function components_interact_from_min_distance(min_distance) &
        result(components_interact)
        class(Abstract_Min_Distance), intent(in) :: min_distance

        select type (min_distance)
            type is (Concrete_Min_Distance)
                components_interact = .true.
            class default
                components_interact = .false.
        end select
    end function components_interact_from_min_distance

    pure logical function components_interact_from_pair_potential(pair_potential) &
        result(components_interact)
        class(Abstract_Pair_Potential), intent(in) :: pair_potential

        select type (pair_potential)
            type is (Null_Pair_Potential)
                components_interact = .false.
            class default
                components_interact = .true.
        end select
    end function components_interact_from_pair_potential

    pure logical function component_interacts_with_wall(wall_potential) result(component_interacts)
        class(Abstract_Pair_Potential), intent(in) :: wall_potential

        select type (wall_potential)
            type is (Null_Pair_Potential)
                component_interacts = .false.
            class default
                component_interacts = .true.
        end select
    end function component_interacts_with_wall

    pure logical function list_is_linked_list(list)
        class(Abstract_Visitable_List), intent(in) :: list

        select type (list)
            type is (Concrete_Visitable_List)
                list_is_linked_list = .true.
            class default
                list_is_linked_list = .false.
        end select
    end function list_is_linked_list

    pure logical function list_is_array(list)
        class(Abstract_Visitable_List), intent(in) :: list

        select type (list)
            type is (Concrete_Visitable_Array)
                list_is_array = .true.
            class default
                list_is_array = .false.
        end select
    end function list_is_array


end module procedures_short_interactions_inquirers
