module classes_random_position

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use classes_parallelepiped_domain, only: Abstract_Parallelepiped_Domain
use classes_random_coordinates, only: Abstract_Random_Coordinates

implicit none

private

    type, extends(Abstract_Random_Coordinates), public :: Concrete_Random_Position
    private
        class(Abstract_Parallelepiped_Domain), pointer :: parallelepiped_domain => null()
        logical, allocatable :: have_positions(:)
    contains
        procedure :: construct => Concrete_construct
        procedure :: destroy => Concrete_destroy
        procedure :: get => Concrete_get
    end type Concrete_Random_Position

contains

    subroutine Concrete_construct(this, parallelepiped_domain, have_positions)
        class(Concrete_Random_Position), intent(out) :: this
        class(Abstract_Parallelepiped_Domain), target, intent(in) :: parallelepiped_domain
        logical, intent(in) :: have_positions(:)

        this%parallelepiped_domain => parallelepiped_domain
        allocate(this%have_positions, source=have_positions)
    end subroutine Concrete_construct

    subroutine Concrete_destroy(this)
        class(Concrete_Random_Position), intent(inout) :: this

        if (allocated(this%have_positions)) deallocate(this%have_positions)
        this%parallelepiped_domain => null()
    end subroutine Concrete_destroy

    function Concrete_get(this, i_component) result(random_position)
        class(Concrete_Random_Position), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), dimension(num_dimensions) :: random_position

        real(DP) :: rand_3d(num_dimensions)

        if (this%have_positions(i_component)) then
            call random_number(rand_3d)
            random_position = this%parallelepiped_domain%get_origin() + (rand_3d - 0.5_DP) * this%&
                parallelepiped_domain%get_size()
        else
            random_position = 0._DP
        end if
    end function Concrete_get

end module classes_random_position
