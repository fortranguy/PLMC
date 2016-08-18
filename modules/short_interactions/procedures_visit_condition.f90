module procedures_visit_condition

implicit none

    abstract interface

        pure logical function abstract_visit_condition(i_inside, i_outside)
            integer, intent(in) :: i_inside, i_outside
        end function abstract_visit_condition

    end interface

private
public :: abstract_visit_condition, visit_different, visit_lower, visit_all

contains

    pure logical function visit_different(i_inside, i_outside)
        integer, intent(in) :: i_inside, i_outside

        visit_different = i_inside /= i_outside
    end function visit_different

    pure logical function visit_lower(i_inside, i_outside)
        integer, intent(in) :: i_inside, i_outside

        visit_lower = i_inside < i_outside
    end function visit_lower

    pure logical function visit_all(i_inside, i_outside)
        integer, intent(in) :: i_inside, i_outside

        visit_all = .true.
    end function visit_all

end module procedures_visit_condition
