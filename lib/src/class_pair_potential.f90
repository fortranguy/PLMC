module class_pair_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_precisions, only: real_zero
use procedures_errors, only: error_exit, warning_continue
use procedures_checks, only: check_positive
use class_potential_expression, only: Abstract_Potential_Expression
use types_potential_domain, only: Concrete_Potential_Domain

implicit none

private

    type, abstract, public :: Abstract_Pair_Potential
    private
        type(Concrete_Potential_Domain) :: domain
        class(Abstract_Potential_Expression), pointer :: potential_expression
        real(DP), allocatable :: tabulation(:)
    contains
        procedure :: construct => Abstract_Pair_Potential_construct
        procedure, private :: set_domain => Abstract_Pair_Potential_set_domain
        procedure, private :: set_tabulation => Abstract_Pair_Potential_set_tabulation
        procedure :: destroy => Abstract_Pair_Potential_destroy
        procedure :: get_max_distance => Abstract_Pair_Potential_get_max_distance
        procedure :: meet => Abstract_Pair_Potential_meet
    end type Abstract_Pair_Potential

    type, extends(Abstract_Pair_Potential), public :: Concrete_Pair_Potential

    end type Concrete_Pair_Potential

    type, extends(Abstract_Pair_Potential), public :: Hard_Pair_Potential
    contains
        procedure :: construct => Hard_Pair_Potential_construct
        procedure, private :: set_domain => Hard_Pair_Potential_set_domain
        procedure :: destroy => Hard_Pair_Potential_destroy
        procedure :: meet => Hard_Pair_Potential_meet
    end type Hard_Pair_Potential

contains

!implementation Abstract_Pair_Potential

    subroutine Abstract_Pair_Potential_construct(this, domain, potential_expression)
        class(Abstract_Pair_Potential), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        class(Abstract_Potential_Expression), target, intent(in) :: potential_expression

        this%potential_expression => potential_expression
        call this%set_domain(domain)
        call this%set_tabulation()
    end subroutine Abstract_Pair_Potential_construct

    subroutine Abstract_Pair_Potential_set_domain(this, domain)
        class(Abstract_Pair_Potential), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        real(DP) :: distance_range

        call check_positive("Abstract_Pair_Potential", "domain%min", domain%min)
        this%domain%min = domain%min
        call check_positive("Abstract_Pair_Potential", "domain%max", domain%max)
        this%domain%max = domain%max
        if (domain%min > domain%max) then
            call error_exit("Abstract_Pair_Potential: "//"domain%min > domain%max.")
        end if
        call check_positive("Abstract_Pair_Potential", "domain%delta", domain%delta)
        distance_range = domain%max - domain%min
        if (distance_range < real_zero) then
            call warning_continue("Abstract_Pair_Potential: "//"distance_range may be too small.")
        end if
        if (distance_range / domain%delta < 1._DP) then
            call warning_continue("Abstract_Pair_Potential: "//"domain%delta may be too big.")
        end if
        this%domain%delta = domain%delta
    end subroutine Abstract_Pair_Potential_set_domain

    subroutine Abstract_Pair_Potential_set_tabulation(this)
        class(Abstract_Pair_Potential), intent(inout) :: this

        real(DP) :: distance_i
        integer :: i_min, i_max, i_distance

        i_min = int(this%domain%min/this%domain%delta)
        i_max = int(this%domain%max/this%domain%delta) + 1
        allocate(this%tabulation(i_min:i_max))
        do i_distance = i_min, i_max
            distance_i = real(i_distance, DP) * this%domain%delta
            this%tabulation(i_distance) = this%potential_expression%get(distance_i)
        end do
        this%tabulation = this%tabulation - this%tabulation(i_max)
    end subroutine Abstract_Pair_Potential_set_tabulation

    subroutine Abstract_Pair_Potential_destroy(this)
        class(Abstract_Pair_Potential), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)
        this%potential_expression => null()
    end subroutine Abstract_Pair_Potential_destroy

    pure function Abstract_Pair_Potential_get_max_distance(this) result(max)
        class(Abstract_Pair_Potential), intent(in) :: this
        real(DP) :: max

        max = this%domain%max
    end function Abstract_Pair_Potential_get_max_distance

    pure subroutine Abstract_Pair_Potential_meet(this, distance, overlap, energy)
        class(Abstract_Pair_Potential), intent(in) :: this
        real(DP), intent(in) :: distance
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        real(DP) :: distance_i
        integer :: i_distance

        overlap = .false.
        energy = 0._DP
        if (distance < this%domain%min) then
            overlap = .true.
        else if (distance < this%domain%max) then
            i_distance = int(distance/this%domain%delta)
            distance_i = real(i_distance, DP) * this%domain%delta
            energy = this%tabulation(i_distance) + &
                (distance - distance_i) * &
                (this%tabulation(i_distance + 1) - this%tabulation(i_distance)) / &
                this%domain%delta
        end if
    end subroutine Abstract_Pair_Potential_meet

!end implementation Abstract_Pair_Potential

!implementation Hard_Pair_Potential

    subroutine Hard_Pair_Potential_construct(this, domain, potential_expression)
        class(Hard_Pair_Potential), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        class(Abstract_Potential_Expression), target, intent(in) :: potential_expression

        call this%set_domain(domain)
    end subroutine Hard_Pair_Potential_construct

    subroutine Hard_Pair_Potential_set_domain(this, domain)
        class(Hard_Pair_Potential), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        call check_positive("Hard_Pair_Potential", "domain%min", domain%min)
        this%domain%min = domain%min
        this%domain%max = this%domain%min
        this%domain%delta = 0._DP
    end subroutine Hard_Pair_Potential_set_domain

    subroutine Hard_Pair_Potential_destroy(this)
        class(Hard_Pair_Potential), intent(inout) :: this

    end subroutine Hard_Pair_Potential_destroy

    pure subroutine Hard_Pair_Potential_meet(this, distance, overlap, energy)
        class(Hard_Pair_Potential), intent(in) :: this
        real(DP), intent(in) :: distance
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy

        overlap = .false.
        energy = 0._DP
        if (distance < this%domain%min) then
            overlap = .true.
        end if
    end subroutine Hard_Pair_Potential_meet

!end implementation Hard_Pair_Potential

end module class_pair_potential