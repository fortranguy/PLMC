module class_ewald_real_pair

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_errors, only: error_exit
use procedures_checks, only: check_positive, check_potential_domain
use procedures_ewald_micro, only: ewald_real_B, ewald_real_C
use types_potential_domain, only: Concrete_Potential_Domain

implicit none

private

    type, abstract, public :: Abstract_Ewald_Real_Pair
    private
        type(Concrete_Potential_Domain) :: domain
        real(DP) :: alpha
    contains
        procedure(Abstract_Ewald_Real_Pair_construct), deferred :: construct
        procedure(Abstract_Ewald_Real_Pair_destroy), deferred :: destroy
        generic :: meet => meet_energy, meet_field
        procedure(Abstract_Ewald_Real_Pair_expression), private, deferred :: expression
        procedure, private :: meet_energy => Abstract_Ewald_Real_Pair_meet_energy
        procedure, private :: meet_field => Abstract_Ewald_Real_Pair_meet_field
    end type Abstract_Ewald_Real_Pair

    abstract interface

        subroutine Abstract_Ewald_Real_Pair_construct(this, domain, alpha)
        import :: DP, Concrete_Potential_Domain, Abstract_Ewald_Real_Pair
            class(Abstract_Ewald_Real_Pair), intent(out) :: this
            type(Concrete_Potential_Domain), intent(in) :: domain
            real(DP), intent(in) :: alpha
        end subroutine Abstract_Ewald_Real_Pair_construct

        subroutine Abstract_Ewald_Real_Pair_destroy(this)
        import :: Abstract_Ewald_Real_Pair
            class(Abstract_Ewald_Real_Pair), intent(inout) :: this
        end subroutine Abstract_Ewald_Real_Pair_destroy

        pure function Abstract_Ewald_Real_Pair_expression(this, distance) result(expression)
        import :: DP, Abstract_Ewald_Real_Pair
            real(DP), dimension(2) :: expression
            class(Abstract_Ewald_Real_Pair), intent(in) :: this
            real(DP), intent(in) :: distance
        end function Abstract_Ewald_Real_Pair_expression

    end interface

    type, extends(Abstract_Ewald_Real_Pair), public :: Tabulated_Ewald_Real_Pair
    private
        real(DP), dimension(:, :), allocatable :: tabulation
    contains
        procedure :: construct => Tabulated_Ewald_Real_Pair_construct
        procedure :: destroy => Tabulated_Ewald_Real_Pair_destroy
        procedure, private :: set_domain => Tabulated_Ewald_Real_Pair_set_domain
        procedure, private :: set_tabulation => Tabulated_Ewald_Real_Pair_set_tabulation
        procedure, private :: expression => Tabulated_Ewald_Real_Pair_expression
    end type Tabulated_Ewald_Real_Pair

    type, extends(Abstract_Ewald_Real_Pair), public :: Raw_Ewald_Real_Pair
    private
        real(DP) :: expression_domain_max(2)
    contains
        procedure :: construct => Raw_Ewald_Real_Pair_construct
        procedure :: destroy => Raw_Ewald_Real_Pair_destroy
        procedure, private :: set_domain => Raw_Ewald_Real_Pair_set_domain
        procedure, private :: expression => Raw_Ewald_Real_Pair_expression
    end type Raw_Ewald_Real_Pair

    type, extends(Abstract_Ewald_Real_Pair), public :: Null_Ewald_Real_Pair
    contains
        procedure :: construct => Null_Ewald_Real_Pair_construct
        procedure :: destroy => Null_Ewald_Real_Pair_destroy
        procedure, private :: expression => Null_Ewald_Real_Pair_expression
        procedure, private :: meet_energy => Null_Ewald_Real_Pair_meet_energy
        procedure, private :: meet_field => Null_Ewald_Real_Pair_meet_field
    end type Null_Ewald_Real_Pair

contains

!implementation Abstract_Ewald_Real_Pair

    !> Between 2 particles
    !> \f[ (\vec{\mu}_i\cdot\vec{\mu}_j) B_\alpha(r_{ij}) -
    !>     (\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij}) C_\alpha(r_{ij}) \f]
    pure real(DP) function Abstract_Ewald_Real_Pair_meet_energy(this, vector_ij, moment_i, &
        moment_j) result(energy)
        class(Abstract_Ewald_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP), dimension(:), intent(in) :: moment_i, moment_j

        real(DP), dimension(2) :: coefficient

        coefficient(1) = dot_product(moment_i, moment_j)
        coefficient(2) =-dot_product(moment_i, vector_ij) * dot_product(moment_j, vector_ij)
        energy = dot_product(coefficient, this%expression(norm2(vector_ij)))
    end function Abstract_Ewald_Real_Pair_meet_energy

    !> Field: to check
    !> \f[
    !>      \vec{E}(\vec{r}_i) = -B_\alpha(r_{ij}) |\vec{\mu}_j) +
    !>                           (\vec{r}_{ij}\cdot\vec{\mu}_j) C_\alpha(r_{ij}) |\vec{r}_{ij})
    !> \f]
    pure function Abstract_Ewald_Real_Pair_meet_field(this, vector_ij, moment_j) result(field)
        real(DP) :: field(num_dimensions)
        class(Abstract_Ewald_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij, moment_j

        real(DP), dimension(2) :: expression

        expression = this%expression(norm2(vector_ij))
        field = -expression(1) * moment_j + &
            dot_product(vector_ij, moment_j) * expression(2) * vector_ij
    end function Abstract_Ewald_Real_Pair_meet_field

!end implementation Abstract_Ewald_Real_Pair

!implementation Tabulated_Ewald_Real_Pair

    subroutine Tabulated_Ewald_Real_Pair_construct(this, domain, alpha)
        class(Tabulated_Ewald_Real_Pair), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        real(DP), intent(in) :: alpha

        call this%set_domain(domain)
        call check_positive("Tabulated_Ewald_Real_Pair_construct", "alpha", alpha)
        this%alpha = alpha
        call this%set_tabulation()
    end subroutine Tabulated_Ewald_Real_Pair_construct

    subroutine Tabulated_Ewald_Real_Pair_set_domain(this, domain)
        class(Tabulated_Ewald_Real_Pair), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        call check_potential_domain("Tabulated_Ewald_Real_Pair_set_domain", "domain", domain)
        this%domain%min = domain%min
        this%domain%max = domain%max
        this%domain%delta = domain%delta
    end subroutine Tabulated_Ewald_Real_Pair_set_domain

    pure subroutine Tabulated_Ewald_Real_Pair_set_tabulation(this)
        class(Tabulated_Ewald_Real_Pair), intent(inout) :: this

        real(DP) :: distance_i
        integer :: i_min, i_max, i_distance

        i_min = int(this%domain%min/this%domain%delta)
        i_max = int(this%domain%max/this%domain%delta) + 1
        allocate(this%tabulation(i_min:i_max, 2))
        do i_distance = i_min, i_max
            distance_i = real(i_distance, DP) * this%domain%delta
            this%tabulation(i_distance, 1) = ewald_real_B(this%alpha, distance_i)
            this%tabulation(i_distance, 1) = ewald_real_C(this%alpha, distance_i)
        end do
        this%tabulation(:, 1) = this%tabulation(:, 1) - this%tabulation(i_max, 1)
        this%tabulation(:, 2) = this%tabulation(:, 2) - this%tabulation(i_max, 2)
    end subroutine Tabulated_Ewald_Real_Pair_set_tabulation

    subroutine Tabulated_Ewald_Real_Pair_destroy(this)
        class(Tabulated_Ewald_Real_Pair), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)
    end subroutine Tabulated_Ewald_Real_Pair_destroy

    !> Linear interpolation
    pure function Tabulated_Ewald_Real_Pair_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Tabulated_Ewald_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance

        integer :: i_distance
        real(DP) :: distance_i

        if (distance < this%domain%max) then
            i_distance = int(distance/this%domain%delta)
            distance_i = real(i_distance, DP) * this%domain%delta
            expression = this%tabulation(i_distance, :) + &
                (distance - distance_i) * &
                (this%tabulation(i_distance + 1, :) - this%tabulation(i_distance, :)) / &
                this%domain%delta
        else
            expression = 0._DP
        end if
    end function Tabulated_Ewald_Real_Pair_expression

!end implementation Tabulated_Ewald_Real_Pair

!implementation Raw_Ewald_Real_Pair

    subroutine Raw_Ewald_Real_Pair_construct(this, domain, alpha)
        class(Raw_Ewald_Real_Pair), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        real(DP), intent(in) :: alpha

        call this%set_domain(domain)
        call check_positive("Tabulated_Ewald_Real_Pair_construct", "alpha", alpha)
        this%alpha = alpha
        this%expression_domain_max = this%expression(this%domain%max)
    end subroutine Raw_Ewald_Real_Pair_construct

    subroutine Raw_Ewald_Real_Pair_set_domain(this, domain)
        class(Raw_Ewald_Real_Pair), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        call check_positive("Raw_Ewald_Real_Pair", "domain%min", domain%min)
        call check_positive("Raw_Ewald_Real_Pair", "domain%max", domain%max)
        if (domain%min > domain%max) then
            call error_exit("Raw_Ewald_Real_Pair: domain%min > domain%max.")
        end if
        this%domain%min = domain%min
        this%domain%max = domain%max
        this%domain%delta = 0._DP
    end subroutine Raw_Ewald_Real_Pair_set_domain

    pure function Raw_Ewald_Real_Pair_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Raw_Ewald_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance

        if (distance < this%domain%max) then
            expression = [ewald_real_B(this%alpha, distance), ewald_real_C(this%alpha, distance)]
            expression = expression - this%expression_domain_max
        else
            expression = 0._DP
        end if
    end function Raw_Ewald_Real_Pair_expression

    subroutine Raw_Ewald_Real_Pair_destroy(this)
        class(Raw_Ewald_Real_Pair), intent(inout) :: this
    end subroutine Raw_Ewald_Real_Pair_destroy

!end implementation Raw_Ewald_Real_Pair

!implementation Null_Ewald_Real_Pair

    subroutine Null_Ewald_Real_Pair_construct(this, domain, alpha)
        class(Null_Ewald_Real_Pair), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        real(DP), intent(in) :: alpha
    end subroutine Null_Ewald_Real_Pair_construct

    subroutine Null_Ewald_Real_Pair_destroy(this)
        class(Null_Ewald_Real_Pair), intent(inout) :: this
    end subroutine Null_Ewald_Real_Pair_destroy

    pure function Null_Ewald_Real_Pair_expression(this, distance) result(expression)
        real(DP) :: expression(2)
        class(Null_Ewald_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance
        expression = 0._DP
    end function Null_Ewald_Real_Pair_expression

    pure real(DP) function Null_Ewald_Real_Pair_meet_energy(this, vector_ij, moment_i, &
        moment_j) result(energy)
        class(Null_Ewald_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP), dimension(:), intent(in) :: moment_i, moment_j
        energy = 0._DP
    end function Null_Ewald_Real_Pair_meet_energy

    pure function Null_Ewald_Real_Pair_meet_field(this, vector_ij, moment_j) result(field)
        real(DP) :: field(num_dimensions)
        class(Null_Ewald_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij, moment_j
        field = 0._DP
    end function Null_Ewald_Real_Pair_meet_field

!end implementation Null_Ewald_Real_Pair

end module class_ewald_real_pair
