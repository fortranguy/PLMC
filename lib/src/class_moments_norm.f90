module class_moments_norm

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use procedures_errors, only: warning_continue
use procedures_checks, only: check_in_range, check_positive
use class_number, only: Abstract_Number

implicit none

private

    type, abstract, public :: Abstract_Moments_Norm
    contains
        procedure(Abstract_Moments_Norm_construct), deferred :: construct
        procedure(Abstract_Moments_Norm_destroy), deferred :: destroy
        procedure(Abstract_Moments_Norm_set), deferred :: set
        procedure(Abstract_Moments_Norm_get_num), deferred :: get_num
        procedure(Abstract_Moments_Norm_get), deferred :: get
        procedure(Abstract_Moments_Norm_add), deferred :: add
        procedure(Abstract_Moments_Norm_remove), deferred :: remove
    end type Abstract_Moments_Norm

    abstract interface

        subroutine Abstract_Moments_Norm_construct(this, number)
        import :: Abstract_Number, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(out) :: this
            class(Abstract_Number), target, intent(in) :: number
        end subroutine Abstract_Moments_Norm_construct

        subroutine Abstract_Moments_Norm_destroy(this)
        import :: Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
        end subroutine Abstract_Moments_Norm_destroy

        subroutine Abstract_Moments_Norm_set(this, i_particle, norm)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
            integer, intent(in) :: i_particle
            real(DP), intent(in) :: norm
        end subroutine Abstract_Moments_Norm_set

        pure function Abstract_Moments_Norm_get_num(this) result(num_norms)
        import :: Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(in) :: this
            integer :: num_norms
        end function Abstract_Moments_Norm_get_num

        pure function Abstract_Moments_Norm_get(this, i_particle) result(norm)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(in) :: this
            integer, intent(in) :: i_particle
            real(DP) :: norm
        end function Abstract_Moments_Norm_get

        subroutine Abstract_Moments_Norm_add(this, norm)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
            real(DP), intent(in) :: norm
        end subroutine Abstract_Moments_Norm_add

        subroutine Abstract_Moments_Norm_remove(this, i_particle)
        import :: DP, Abstract_Moments_Norm
            class(Abstract_Moments_Norm), intent(inout) :: this
            integer, intent(in) :: i_particle
        end subroutine Abstract_Moments_Norm_remove

    end interface

    type, extends(Abstract_Moments_Norm), public :: Null_Moments_Norm
    contains
        procedure :: construct => Null_Moments_Norm_construct
        procedure :: destroy => Null_Moments_Norm_destroy
        procedure :: set => Null_Moments_Norm_set
        procedure :: get_num => Null_Moments_Norm_get_num
        procedure :: get => Null_Moments_Norm_get
        procedure :: add => Null_Moments_Norm_add
        procedure :: remove => Null_Moments_Norm_remove
    end type Null_Moments_Norm

    type, extends(Abstract_Moments_Norm), public :: Uniform_Moments_Norm
    private
        real(DP) :: norm
        logical :: is_set = .false.
        class(Abstract_Number), pointer :: number
    contains
        procedure :: construct => Uniform_Moments_Norm_construct
        procedure :: destroy => Uniform_Moments_Norm_destroy
        procedure :: set => Uniform_Moments_Norm_set
        procedure :: get_num => Uniform_Moments_Norm_get_num
        procedure :: get => Uniform_Moments_Norm_get
        procedure :: add => Uniform_Moments_Norm_add
        procedure :: remove => Uniform_Moments_Norm_remove
    end type Uniform_Moments_Norm

contains

!implementation Null_Moments_Norm

    subroutine Null_Moments_Norm_construct(this, number)
        class(Null_Moments_Norm), intent(out) :: this
        class(Abstract_Number), target, intent(in) :: number
    end subroutine Null_Moments_Norm_construct

    subroutine Null_Moments_Norm_destroy(this)
        class(Null_Moments_Norm), intent(inout) :: this
    end subroutine Null_Moments_Norm_destroy

    subroutine Null_Moments_Norm_set(this, i_particle, norm)
        class(Null_Moments_Norm), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm
    end subroutine Null_Moments_Norm_set

    pure function Null_Moments_Norm_get_num(this) result(num_norms)
        class(Null_Moments_Norm), intent(in) :: this
        integer :: num_norms
        num_norms = 0
    end function Null_Moments_Norm_get_num

    pure function Null_Moments_Norm_get(this, i_particle) result(norm)
        class(Null_Moments_Norm), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: norm
        norm = 0._DP
    end function Null_Moments_Norm_get

    subroutine Null_Moments_Norm_add(this, norm)
        class(Null_Moments_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm
    end subroutine Null_Moments_Norm_add

    subroutine Null_Moments_Norm_remove(this, i_particle)
        class(Null_Moments_Norm), intent(inout) :: this
        integer, intent(in) :: i_particle
    end subroutine Null_Moments_Norm_remove

!end implementation Null_Moments_Norm

!implementation Uniform_Moments_Norm

    subroutine Uniform_Moments_Norm_construct(this, number)
        class(Uniform_Moments_Norm), intent(out) :: this
        class(Abstract_Number), target, intent(in) :: number

        this%number => number
    end subroutine Uniform_Moments_Norm_construct

    subroutine Uniform_Moments_Norm_destroy(this)
        class(Uniform_Moments_Norm), intent(inout) :: this

        this%number => null()
    end subroutine Uniform_Moments_Norm_destroy

    subroutine Uniform_Moments_Norm_set(this, i_particle, norm)
        class(Uniform_Moments_Norm), intent(inout) :: this
        integer, intent(in) :: i_particle
        real(DP), intent(in) :: norm

        call check_in_range("Uniform_Moments_Norm", this%number%get(), "i_particle", i_particle)
        call check_positive("Uniform_Moments_Norm", "norm", norm)
        if (.not. this%is_set) then
            this%is_set = .true.
        else if (abs(norm - this%norm) > real_zero) then
            call warning_continue("Uniform_Moments_Norm: setting norm is different.")
        end if
        this%norm = norm
    end subroutine Uniform_Moments_Norm_set

    pure function Uniform_Moments_Norm_get_num(this) result(num_norms)
        class(Uniform_Moments_Norm), intent(in) :: this
        integer  :: num_norms

        num_norms = this%number%get()
    end function Uniform_Moments_Norm_get_num

    pure function Uniform_Moments_Norm_get(this, i_particle) result(norm)
        class(Uniform_Moments_Norm), intent(in) :: this
        integer, intent(in) :: i_particle
        real(DP) :: norm

        norm = this%norm
    end function Uniform_Moments_Norm_get

    subroutine Uniform_Moments_Norm_add(this, norm)
        class(Uniform_Moments_Norm), intent(inout) :: this
        real(DP), intent(in) :: norm

        call this%set(this%number%get(), norm)
    end subroutine Uniform_Moments_Norm_add

    subroutine Uniform_Moments_Norm_remove(this, i_particle)
         class(Uniform_Moments_Norm), intent(inout) :: this
         integer, intent(in) :: i_particle

        call check_in_range("Uniform_Moments_Norm", this%number%get(), "i_particle", i_particle)
    end subroutine Uniform_Moments_Norm_remove

!end implementation Uniform_Moments_Norm

end module class_moments_norm
