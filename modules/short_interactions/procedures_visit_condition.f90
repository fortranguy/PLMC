module procedures_visit_condition

implicit none

    abstract interface
        pure logical function in_range(i_inside, i_outside)
            integer, intent(in) :: i_inside, i_outside
        end function
    end interface

private
public :: in_range, different, lower, unconditional

contains

    pure logical function different(i_inside, i_outside)
        integer, intent(in) :: i_inside, i_outside

        different = i_inside /= i_outside
    end function different

    pure logical function lower(i_inside, i_outside)
    integer, intent(in) :: i_inside, i_outside

        lower = i_inside < i_outside
    end function lower

    pure logical function unconditional(i_inside, i_outside)
        integer, intent(in) :: i_inside, i_outside

        unconditional = .true.
    end function unconditional

end module procedures_visit_condition
