module classes_hetero_couples

use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Hetero_Couples
    private
        integer :: num_indices = 0
        integer, allocatable :: couple(:, :)
    contains
        procedure(Abstract_construct), deferred :: construct
        procedure :: destroy => Abstract_destroy
        procedure :: get_num_indices => Abstract_get_num_indices
        procedure :: get => Abstract_get
    end type Abstract_Hetero_Couples

    abstract interface

        subroutine Abstract_construct(this, num_components)
        import :: Abstract_Hetero_Couples
            class(Abstract_Hetero_Couples), intent(out) :: this
            integer, intent(in) :: num_components
        end subroutine Abstract_construct

    end interface

    type, extends(Abstract_Hetero_Couples), public :: Half_Hetero_Couples
    contains
        procedure :: construct => Half_construct
    end type Half_Hetero_Couples

    type, extends(Abstract_Hetero_Couples), public :: Full_Hetero_Couples
    contains
        procedure :: construct => Full_construct
    end type Full_Hetero_Couples

    type, extends(Abstract_Hetero_Couples), public :: Null_Hetero_Couples
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: get_num_indices => Null_get_num_indices
        procedure :: get => Null_get
    end type Null_Hetero_Couples

contains

!implementation Abstract_Hetero_Couples

    subroutine Abstract_destroy(this)
        class(Abstract_Hetero_Couples), intent(inout) :: this

        if (allocated(this%couple)) deallocate(this%couple)
    end subroutine Abstract_destroy

    pure integer function Abstract_get_num_indices(this) result(num_indices)
        class(Abstract_Hetero_Couples), intent(in) :: this

        num_indices = this%num_indices
    end function Abstract_get_num_indices

    pure function Abstract_get(this, index) result(couple)
        class(Abstract_Hetero_Couples), intent(in) :: this
        integer, intent(in) :: index
        integer :: couple(2)

        couple = this%couple(:, index)
    end function Abstract_get

!implementation Abstract_Hetero_Couples

!implementation Half_Hetero_Couples

    subroutine Half_construct(this, num_components)
        class(Half_Hetero_Couples), intent(out) :: this
        integer, intent(in) :: num_components

        integer :: index, i_component, j_component

        call check_positive("Half_Hetero_Couples: construct", "num_components", num_components)
        this%num_indices = num_components * (num_components - 1) / 2
        allocate(this%couple(2, this%num_indices))
        this%couple = 0
        index = 0
        do j_component = 2, num_components
            do i_component = 1, j_component - 1
                index = index + 1
                this%couple(:, index) = [j_component, i_component]
            end do
        end do
    end subroutine Half_construct

!end implementation Half_Hetero_Couples

!implementation Full_Hetero_Couples

    subroutine Full_construct(this, num_components)
        class(Full_Hetero_Couples), intent(out) :: this
        integer, intent(in) :: num_components

        integer :: index, i_component, j_component
        call check_positive("Half_Hetero_Couples: construct", "num_components", num_components)
        this%num_indices = num_components * (num_components - 1)
        allocate(this%couple(2, this%num_indices))
        this%couple = 0
        index = 0
        do j_component = 1, num_components
            do i_component = 1, num_components
                if (i_component /= j_component) then
                    index = index + 1
                    this%couple(:, index) = [j_component, i_component]
                end if
            end do
        end do
    end subroutine Full_construct

!end implementation Full_Hetero_Couples

!implementation Null_Hetero_Couples

    subroutine Null_construct(this, num_components)
        class(Null_Hetero_Couples), intent(out) :: this
        integer, intent(in) :: num_components
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Hetero_Couples), intent(inout) :: this
    end subroutine Null_destroy

    pure integer function Null_get_num_indices(this) result(num_indices)
        class(Null_Hetero_Couples), intent(in) :: this
        num_indices = 0
    end function Null_get_num_indices

    pure function Null_get(this, index) result(couple)
        class(Null_Hetero_Couples), intent(in) :: this
        integer, intent(in) :: index
        integer :: couple(2)
        couple = 0
    end function Null_get

!implementation Null_Hetero_Couples

end module classes_hetero_couples
