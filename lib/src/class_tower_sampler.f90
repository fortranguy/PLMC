module class_tower_sampler

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Tower_Sampler
    private
        logical :: unique_candidate
        integer :: i_unique_candidate
        integer :: num_candidates
        real(DP), allocatable :: limits(:)
    contains
        procedure :: construct => Abstract_Tower_Sampler_construct
        procedure :: destroy => Abstract_Tower_Sampler_destroy
        procedure :: get => Abstract_Tower_Sampler_get
    end type Abstract_Tower_Sampler

    type, extends(Abstract_Tower_Sampler), public :: Concrete_Tower_Sampler

    end type Concrete_Tower_Sampler

    type, extends(Abstract_Tower_Sampler), public :: Null_Tower_Sampler
    contains
        procedure :: construct => Null_Tower_Sampler_construct
        procedure :: destroy => Null_Tower_Sampler_destroy
        procedure :: get => Null_Tower_Sampler_get
    end type Null_Tower_Sampler

contains

!implementation Abstract_Tower_Sampler

    subroutine Abstract_Tower_Sampler_construct(this, nums_candidates)
        class(Abstract_Tower_Sampler), intent(out) :: this
        integer, intent(in) :: nums_candidates(:)

        real(DP), allocatable :: cumulative_weight(:)
        integer :: i_candidate

        call check_positive("Abstract_Tower_Sampler_construct", "sum(nums_candidates)", &
            sum(nums_candidates))
        this%num_candidates = size(nums_candidates)
        if (count(nums_candidates > 0) == 1) then
            this%unique_candidate = .true.
            this%i_unique_candidate = maxloc(nums_candidates, 1)
        else
            this%unique_candidate = .false.
            allocate(cumulative_weight(this%num_candidates))
            do i_candidate = 1, this%num_candidates
                cumulative_weight(i_candidate) = real(sum(nums_candidates(1:i_candidate)), DP) / &
                    real(sum(nums_candidates), DP)
            end do
            this%limits = [0._DP, cumulative_weight]
            deallocate(cumulative_weight)
        end if
    end subroutine Abstract_Tower_Sampler_construct

    subroutine Abstract_Tower_Sampler_destroy(this)
        class(Abstract_Tower_Sampler), intent(inout) :: this

        if (allocated(this%limits)) deallocate(this%limits)
    end subroutine Abstract_Tower_Sampler_destroy

    integer function Abstract_Tower_Sampler_get(this) result(i_random_candidate)
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
    end function Abstract_Tower_Sampler_get

!end implementation Abstract_Tower_Sampler

!implementation Null_Tower_Sampler

    subroutine Null_Tower_Sampler_construct(this, nums_candidates)
        class(Null_Tower_Sampler), intent(out) :: this
        integer, intent(in) :: nums_candidates(:)
    end subroutine Null_Tower_Sampler_construct

    subroutine Null_Tower_Sampler_destroy(this)
        class(Null_Tower_Sampler), intent(inout) :: this
    end subroutine Null_Tower_Sampler_destroy

    integer function Null_Tower_Sampler_get(this) result(i_random_candidate)
        class(Null_Tower_Sampler), intent(in) :: this
        i_random_candidate = 0
    end function Null_Tower_Sampler_get

!end implementation Null_Tower_Sampler

end module class_tower_sampler

