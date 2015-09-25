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
        class(Abstract_Potential_Expression), pointer :: expression
    contains
        procedure(Abstract_Pair_Potential_construct), deferred :: construct
        procedure(Abstract_Pair_Potential_destroy), deferred :: destroy
        procedure :: get_max_distance => Abstract_Pair_Potential_get_max_distance
        procedure(Abstract_Pair_Potential_meet), deferred :: meet
    end type Abstract_Pair_Potential

    abstract interface

        subroutine Abstract_Pair_Potential_construct(this, domain, expression)
        import :: Concrete_Potential_Domain, Abstract_Potential_Expression, Abstract_Pair_Potential
            class(Abstract_Pair_Potential), intent(out) :: this
            type(Concrete_Potential_Domain), intent(in) :: domain
            class(Abstract_Potential_Expression), target, intent(in) :: expression
        end subroutine Abstract_Pair_Potential_construct

        subroutine Abstract_Pair_Potential_destroy(this)
        import :: Abstract_Pair_Potential
            class(Abstract_Pair_Potential), intent(inout) :: this
        end subroutine Abstract_Pair_Potential_destroy

        pure subroutine Abstract_Pair_Potential_meet(this, overlap, energy, distance)
        import :: DP, Abstract_Pair_Potential
            class(Abstract_Pair_Potential), intent(in) :: this
            logical, intent(out) :: overlap
            real(DP), intent(out) :: energy
            real(DP), intent(in) :: distance
        end subroutine Abstract_Pair_Potential_meet

    end interface

    type, extends(Abstract_Pair_Potential), public :: Tabulated_Pair_Potential
    private
        real(DP), allocatable :: tabulation(:)
    contains
        procedure :: construct => Tabulated_Pair_Potential_construct
        procedure, private :: set_domain => Tabulated_Pair_Potential_set_domain
        procedure, private :: set_tabulation => Tabulated_Pair_Potential_set_tabulation
        procedure :: destroy => Tabulated_Pair_Potential_destroy
        procedure :: meet => Tabulated_Pair_Potential_meet
    end type Tabulated_Pair_Potential

    type, extends(Abstract_Pair_Potential), public :: Null_Pair_Potential
    contains
        procedure :: construct => Null_Pair_Potential_construct
        procedure :: destroy => Null_Pair_Potential_destroy
        procedure :: get_max_distance => Null_Pair_Potential_get_max_distance
        procedure :: meet => Null_Pair_Potential_meet
    end type Null_Pair_Potential

    type, extends(Abstract_Pair_Potential), public :: Raw_Pair_Potential
    private
        real(DP) :: energy_domain_max
    contains
        procedure :: construct => Raw_Pair_Potential_construct
        procedure, private :: set_domain => Raw_Pair_Potential_set_domain
        procedure :: destroy => Raw_Pair_Potential_destroy
        procedure :: meet => Raw_Pair_Potential_meet
    end type Raw_Pair_Potential

contains

!implementation Abstract_Pair_Potential

    pure function Abstract_Pair_Potential_get_max_distance(this) result(max_distance)
        class(Abstract_Pair_Potential), intent(in) :: this
        real(DP) :: max_distance

        max_distance = this%domain%max
    end function Abstract_Pair_Potential_get_max_distance

!end implementation Abstract_Pair_Potential

!implementation Tabulated_Pair_Potential

    subroutine Tabulated_Pair_Potential_construct(this, domain, expression)
        class(Tabulated_Pair_Potential), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        class(Abstract_Potential_Expression), target, intent(in) :: expression

        this%expression => expression
        call this%set_domain(domain)
        call this%set_tabulation()
    end subroutine Tabulated_Pair_Potential_construct

    subroutine Tabulated_Pair_Potential_set_domain(this, domain)
        class(Tabulated_Pair_Potential), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        real(DP) :: distance_range

        call check_positive("Tabulated_Pair_Potential", "domain%min", domain%min)
        call check_positive("Tabulated_Pair_Potential", "domain%max", domain%max)
        if (domain%min > domain%max) then
            call error_exit("Tabulated_Pair_Potential: domain%min > domain%max.")
        end if
        this%domain%min = domain%min
        this%domain%max = domain%max
        call check_positive("Tabulated_Pair_Potential", "domain%delta", domain%delta)
        distance_range = domain%max - domain%min
        if (distance_range < real_zero) then
            call warning_continue("Tabulated_Pair_Potential: "//"distance_range may be too small.")
        end if
        if (distance_range / domain%delta < 1._DP) then
            call warning_continue("Tabulated_Pair_Potential: "//"domain%delta may be too big.")
        end if
        this%domain%delta = domain%delta
    end subroutine Tabulated_Pair_Potential_set_domain

    subroutine Tabulated_Pair_Potential_set_tabulation(this)
        class(Tabulated_Pair_Potential), intent(inout) :: this

        real(DP) :: distance_i
        integer :: i_min, i_max, i_distance

        i_min = int(this%domain%min/this%domain%delta)
        i_max = int(this%domain%max/this%domain%delta) + 1
        allocate(this%tabulation(i_min:i_max))
        do i_distance = i_min, i_max
            distance_i = real(i_distance, DP) * this%domain%delta
            this%tabulation(i_distance) = this%expression%get(distance_i)
        end do
        this%tabulation = this%tabulation - this%tabulation(i_max)
    end subroutine Tabulated_Pair_Potential_set_tabulation

    subroutine Tabulated_Pair_Potential_destroy(this)
        class(Tabulated_Pair_Potential), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)
        this%expression => null()
    end subroutine Tabulated_Pair_Potential_destroy

    pure subroutine Tabulated_Pair_Potential_meet(this, overlap, energy, distance)
        class(Tabulated_Pair_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: distance

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
    end subroutine Tabulated_Pair_Potential_meet

!end implementation Tabulated_Pair_Potential

!implementation Null_Pair_Potential

    subroutine Null_Pair_Potential_construct(this, domain, expression)
        class(Null_Pair_Potential), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        class(Abstract_Potential_Expression), target, intent(in) :: expression
    end subroutine Null_Pair_Potential_construct

    subroutine Null_Pair_Potential_destroy(this)
        class(Null_Pair_Potential), intent(inout) :: this
    end subroutine Null_Pair_Potential_destroy

    pure function Null_Pair_Potential_get_max_distance(this) result(max_distance)
        class(Null_Pair_Potential), intent(in) :: this
        real(DP) :: max_distance
        max_distance = 0._DP
    end function Null_Pair_Potential_get_max_distance

    pure subroutine Null_Pair_Potential_meet(this, overlap, energy, distance)
        class(Null_Pair_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: distance
        overlap = .false.
        energy = 0._DP
    end subroutine Null_Pair_Potential_meet

!end implementation Null_Pair_Potential

!implementation Raw_Pair_Potential

    subroutine Raw_Pair_Potential_construct(this, domain, expression)
        class(Raw_Pair_Potential), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        class(Abstract_Potential_Expression), target, intent(in) :: expression

        this%expression => expression
        call this%set_domain(domain)
        this%energy_domain_max = this%expression%get(this%domain%max)
    end subroutine Raw_Pair_Potential_construct

    subroutine Raw_Pair_Potential_set_domain(this, domain)
        class(Raw_Pair_Potential), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        call check_positive("Raw_Pair_Potential", "domain%min", domain%min)
        call check_positive("Raw_Pair_Potential", "domain%max", domain%max)
        if (domain%min > domain%max) then
            call error_exit("Raw_Pair_Potential: domain%min > domain%max.")
        end if
        this%domain%min = domain%min
        this%domain%max = domain%max
        this%domain%delta = 0._DP
    end subroutine Raw_Pair_Potential_set_domain

    subroutine Raw_Pair_Potential_destroy(this)
        class(Raw_Pair_Potential), intent(inout) :: this

        this%expression => null()
    end subroutine Raw_Pair_Potential_destroy

    pure subroutine Raw_Pair_Potential_meet(this, overlap, energy, distance)
        class(Raw_Pair_Potential), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: energy
        real(DP), intent(in) :: distance

        overlap = .false.
        energy = 0._DP
        if (distance < this%domain%min) then
            overlap = .true.
        else if (distance < this%domain%max) then
            energy = this%expression%get(distance) - this%energy_domain_max
        end if
    end subroutine Raw_Pair_Potential_meet

!end implementation Raw_Pair_Potential

end module class_pair_potential
