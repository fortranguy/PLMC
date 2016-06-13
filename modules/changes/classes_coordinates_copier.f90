module classes_coordinates_copier

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use classes_random_coordinates, only: Abstract_Random_Coordinates

implicit none

private

    type, abstract, public :: Abstract_Coordinates_Copier
    private
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_copy), deferred :: copy
    end type Abstract_Coordinates_Copier

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Coordinates_Copier
            class(Abstract_Coordinates_Copier), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_copy(this, target, source, ij_components)
        import :: DP, Abstract_Coordinates_Copier
            class(Abstract_Coordinates_Copier), intent(in) :: this
            real(DP), intent(out) :: target(:)
            real(DP), intent(in) :: source(:)
            integer, intent(in) :: ij_components(:)
        end subroutine Abstract_copy

    end interface

    type, extends(Abstract_Coordinates_Copier), public :: Random_Filling_Coordinates_Copier
    private
        class(Abstract_Random_Coordinates), pointer :: random_coordinates => null()
        logical, allocatable :: have_coordinates(:)
    contains
        procedure :: construct => Random_construct
        procedure :: destroy => Random_destroy
        procedure :: copy => Random_copy
    end type Random_Filling_Coordinates_Copier

    type, extends(Abstract_Coordinates_Copier), public :: Null_Coordinates_Copier
    contains
        procedure :: destroy => Null_destroy
        procedure :: copy => Null_copy
    end type Null_Coordinates_Copier

contains

!implementation Random_Filling_Coordinates_Copier

    subroutine Random_construct(this, random_coordinates, have_coordinates)
        class(Random_Filling_Coordinates_Copier), intent(out) :: this
        class(Abstract_Random_Coordinates), target, intent(in) :: random_coordinates
        logical, intent(in) :: have_coordinates(:)

        this%random_coordinates => random_coordinates
        allocate(this%have_coordinates, source=have_coordinates)
    end subroutine Random_construct

    subroutine Random_destroy(this)
        class(Random_Filling_Coordinates_Copier), intent(inout) :: this

        if (allocated(this%have_coordinates)) deallocate(this%have_coordinates)
        this%random_coordinates => null()
    end subroutine Random_destroy

    subroutine Random_copy(this, target, source, ij_components)
        class(Random_Filling_Coordinates_Copier), intent(in) :: this
        real(DP), intent(out) :: target(:)
        real(DP), intent(in) :: source(:)
        integer, intent(in) :: ij_components(:)

        if (.not.this%have_coordinates(ij_components(1)) .and. &
            this%have_coordinates(ij_components(2))) then
            target = this%random_coordinates%get(ij_components(2))
        else
            target = source
        end if
    end subroutine Random_copy

!end implementation Random_Filling_Coordinates_Copier

!implementation Null_Coordinates_Copier

    subroutine Null_destroy(this)
        class(Null_Coordinates_Copier), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_copy(this, target, source, ij_components)
        class(Null_Coordinates_Copier), intent(in) :: this
        real(DP), intent(out) :: target(:)
        real(DP), intent(in) :: source(:)
        integer, intent(in) :: ij_components(:)
        target = 0._DP
    end subroutine Null_copy

!end implementation Null_Coordinates_Copier

end module classes_coordinates_copier
