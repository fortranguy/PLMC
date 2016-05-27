module classes_random_coordinates

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions

implicit none

private

    type, abstract, public :: Abstract_Random_Coordinates
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_get), deferred :: get
    end type Abstract_Random_Coordinates

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Random_Coordinates
            class(Abstract_Random_Coordinates), intent(inout) ::  this
        end subroutine Abstract_destroy

        function Abstract_get(this, i_component) result(random_coordinates)
        import :: DP, num_dimensions, Abstract_Random_Coordinates
            class(Abstract_Random_Coordinates), intent(in) :: this
            integer, intent(in) :: i_component
            real(DP), dimension(num_dimensions) :: random_coordinates
        end function Abstract_get

    end interface

    type, extends(Abstract_Random_Coordinates), public :: Null_Random_Coordinates
    contains
        procedure :: destroy => Null_destroy
        procedure :: get => Null_get
    end type Null_Random_Coordinates

contains

!implementation Null_Random_Coordinates

    subroutine Null_destroy(this)
        class(Null_Random_Coordinates), intent(inout) :: this
    end subroutine Null_destroy

    function Null_get(this, i_component) result(random_position)
        class(Null_Random_Coordinates), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), dimension(num_dimensions) :: random_position
        random_position = 0._DP
    end function Null_get

!implementation Null_Random_Coordinates

end module classes_random_coordinates
