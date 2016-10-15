module classes_tower_sampler

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Tower_Sampler
    private
        logical :: unique_candidate = .true.
        integer :: i_unique_candidate = 0
        integer :: num_candidates = 0
        integer :: num_choices = 0
        real(DP), allocatable :: limits(:)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_reset
        procedure :: get_num_choices => Abstract_get_num_choices
        procedure :: get => Abstract_get
    end type Abstract_Tower_Sampler

    type, extends(Abstract_Tower_Sampler), public :: Concrete_Tower_Sampler

    end type Concrete_Tower_Sampler

    type, extends(Abstract_Tower_Sampler), public :: Null_Tower_Sampler
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_reset
        procedure :: get_num_choices => Null_get_num_choices
        procedure :: get => Null_get
    end type Null_Tower_Sampler

contains

!implementation Abstract_Tower_Sampler

    subroutine Abstract_construct(this, num_candidates)
        class(Abstract_Tower_Sampler), intent(out) :: this
        integer, intent(in) :: num_candidates

        call check_positive("Abstract_Tower_Sampler: construct", "num_candidates", num_candidates)
        this%num_candidates = num_candidates
        allocate(this%limits(this%num_candidates + 1))
        this%limits = 0
    end subroutine Abstract_construct

    !> Sets the ``stories'' of the tower.
    subroutine Abstract_reset(this, nums_candidates)
        class(Abstract_Tower_Sampler), intent(inout) :: this
        integer, intent(in) :: nums_candidates(:)

        real(DP) :: cumulative_weight(this%num_candidates)
        integer :: i_candidate

        if (size(nums_candidates) /= this%num_candidates) then
            call error_exit("Abstract_Tower_Sampler: set: size(nums_candidates) and "//&
                "this%num_candidates are different.")
        end if
        call check_positive("Abstract_Tower_Sampler: reset", "sum(nums_candidates)", &
            sum(nums_candidates))
        this%num_choices = sum(nums_candidates)
        if (count(nums_candidates > 0) == 1) then
            this%unique_candidate = .true.
            this%i_unique_candidate = maxloc(nums_candidates, 1)
        else
            this%unique_candidate = .false.
            do i_candidate = 1, this%num_candidates
                cumulative_weight(i_candidate) = real(sum(nums_candidates(1:i_candidate)), DP) / &
                    real(this%num_choices, DP)
            end do
            this%limits(1) = 0._DP
            this%limits(2:size(this%limits)) = cumulative_weight
        end if
    end subroutine Abstract_reset

    subroutine Abstract_destroy(this)
        class(Abstract_Tower_Sampler), intent(inout) :: this

        if (allocated(this%limits)) deallocate(this%limits)
    end subroutine Abstract_destroy

    pure integer function Abstract_get_num_choices(this) result(num_choices)
        class(Abstract_Tower_Sampler), intent(in) :: this

        num_choices = this%num_choices
    end function Abstract_get_num_choices

    integer function Abstract_get(this) result(i_random_candidate)
        class(Abstract_Tower_Sampler), intent(in) :: this

        integer :: i_candidate
        real(DP) :: rand

        if (this%unique_candidate) then
            i_random_candidate = this%i_unique_candidate
            return
        end if
        call random_number(rand)
        do i_candidate = 1, this%num_candidates
            if (this%limits(i_candidate) <= rand .and. rand < this%limits(i_candidate + 1)) then
                i_random_candidate = i_candidate
                return
            end if
        end do
    end function Abstract_get

!end implementation Abstract_Tower_Sampler

!implementation Null_Tower_Sampler

    subroutine Null_construct(this, num_candidates)
        class(Null_Tower_Sampler), intent(out) :: this
        integer, intent(in) :: num_candidates
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Tower_Sampler), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_reset(this, nums_candidates)
        class(Null_Tower_Sampler), intent(inout) :: this
        integer, intent(in) :: nums_candidates(:)
    end subroutine Null_reset

    pure integer function Null_get_num_choices(this) result(num_choices)
        class(Null_Tower_Sampler), intent(in) :: this
        num_choices = 0
    end function Null_get_num_choices

    integer function Null_get(this) result(i_random_candidate)
        class(Null_Tower_Sampler), intent(in) :: this
        i_random_candidate = 0
    end function Null_get

!end implementation Null_Tower_Sampler

end module classes_tower_sampler

