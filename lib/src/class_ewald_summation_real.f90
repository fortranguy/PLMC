module class_ewald_summation_real

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use procedures_checks, only: check_potential_domain
use procedures_ewald_summation, only: ewald_real_B, ewald_real_C
use types_potential_domain, only: Concrete_Potential_Domain

implicit none

private

    type, abstract, public :: Abstract_Ewald_Real_Pair
        type(Concrete_Potential_Domain) :: domain
        real(DP), dimension(:, :), allocatable :: tabulation
    contains
        procedure :: construct => Abstract_Ewald_Real_Pair_construct
        procedure, private :: set_domain => Abstract_Ewald_Real_Pair_set_domain
        procedure, private :: set_tabulation => Abstract_Ewald_Real_Pair_set_tabulation
        procedure :: destroy => Abstract_Ewald_Real_Pair_destroy
        generic :: meet => meet_energy, meet_field
        procedure, private :: meet_energy => Abstract_Ewald_Real_Pair_meet_energy
        procedure, private :: meet_field => Abstract_Ewald_Real_Pair_meet_field
        procedure, private :: interpolation => Abstract_Ewald_Real_Pair_interpolation
    end type Abstract_Ewald_Real_Pair

contains

    subroutine Abstract_Ewald_Real_Pair_construct(this, domain, alpha)
        class(Abstract_Ewald_Real_Pair), intent(out) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain
        real(DP), intent(in) :: alpha

        call this%set_domain(domain)
        call this%set_tabulation(alpha)
    end subroutine Abstract_Ewald_Real_Pair_construct

    subroutine Abstract_Ewald_Real_Pair_set_domain(this, domain)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: this
        type(Concrete_Potential_Domain), intent(in) :: domain

        call check_potential_domain("Abstract_Ewald_Real_Pair_set_domain", "domain", domain)
        this%domain%min = domain%min
        this%domain%max = domain%max
        this%domain%delta = domain%delta
    end subroutine Abstract_Ewald_Real_Pair_set_domain

    pure subroutine Abstract_Ewald_Real_Pair_set_tabulation(this, alpha)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: this
        real(DP), intent(in) :: alpha

        real(DP) :: distance_i
        integer :: i_min, i_max, i_distance

        i_min = int(this%domain%min/this%domain%delta)
        i_max = int(this%domain%max/this%domain%delta) + 1
        allocate(this%tabulation(i_min:i_max, 2))
        do i_distance = i_min, i_max
            distance_i = real(i_distance, DP) * this%domain%delta
            this%tabulation(i_distance, 1) = ewald_real_B(alpha, distance_i)
            this%tabulation(i_distance, 2) = ewald_real_C(alpha, distance_i)
        end do
        this%tabulation(:, 1) = this%tabulation(:, 1) - this%tabulation(i_max, 1)
        this%tabulation(:, 2) = this%tabulation(:, 2) - this%tabulation(i_max, 2)
    end subroutine Abstract_Ewald_Real_Pair_set_tabulation

    subroutine Abstract_Ewald_Real_Pair_destroy(this)
        class(Abstract_Ewald_Real_Pair), intent(inout) :: this

        if (allocated(this%tabulation)) deallocate(this%tabulation)
    end subroutine Abstract_Ewald_Real_Pair_destroy

    !> Between 2 particles
    !> \f[ (\vec{\mu}_i\cdot\vec{\mu}_j) B_\alpha(r_{ij}) -
    !>     (\vec{\mu}_i\cdot\vec{r}_{ij}) (\vec{\mu}_j\cdot\vec{r}_{ij}) C_\alpha(r_{ij}) \f]
    pure function Abstract_Ewald_Real_Pair_meet_energy(this, vector_ij, orientation_i, &
        orientation_j) result(energy)
        real(DP) :: energy
        class(Abstract_Ewald_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij
        real(DP), dimension(:), intent(in) :: orientation_i, orientation_j

        real(DP), dimension(2) :: coefficient

        coefficient(1) = dot_product(orientation_i, orientation_j)
        coefficient(2) =-dot_product(orientation_i, vector_ij) * &
                         dot_product(orientation_j, vector_ij)
        energy = dot_product(coefficient, this%interpolation(norm2(vector_ij)))
    end function Abstract_Ewald_Real_Pair_meet_energy

    !> Field: to check
    !> \f[
    !>      \vec{E}(\vec{r}_i) = -B_\alpha(r_{ij}) |\vec{\mu}_j) +
    !>                           (\vec{r}_{ij}\cdot\vec{\mu}_j) C_\alpha(r_{ij}) |\vec{r}_{ij})
    !> \f]
    pure function Abstract_Ewald_Real_Pair_meet_field(this, vector_ij, moment_j) &
                  result(field)
        real(DP), dimension(num_dimensions) :: field
        class(Abstract_Ewald_Real_Pair), intent(in) :: this
        real(DP), dimension(:), intent(in) :: vector_ij, moment_j

        real(DP), dimension(2) :: interpolation

        interpolation = this%interpolation(norm2(vector_ij))
        field = -interpolation(1) * moment_j + &
            dot_product(vector_ij, moment_j) * interpolation(2) * vector_ij
    end function Abstract_Ewald_Real_Pair_meet_field

    !> Linear interpolation
    pure function Abstract_Ewald_Real_Pair_interpolation(this, distance) result(interpolation)
        real(DP), dimension(2) :: interpolation
        class(Abstract_Ewald_Real_Pair), intent(in) :: this
        real(DP), intent(in) :: distance

        integer :: i_distance
        real(DP) :: distance_i

        if (distance < this%domain%max) then
            i_distance = int(distance/this%domain%delta)
            distance_i = real(i_distance, DP) * this%domain%delta
            interpolation(:) = this%tabulation(i_distance, :) + &
                (distance - distance_i) * &
                (this%tabulation(i_distance + 1, :) - this%tabulation(i_distance, :)) / &
                this%domain%delta
        else
            interpolation(:) = 0._DP
        end if
    end function Abstract_Ewald_Real_Pair_interpolation

end module class_ewald_summation_real
