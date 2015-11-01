!> display: public
!>          private
module class_weighted_structure

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use procedures_ewald_micro, only: set_fourier, reciprocal_size_1_sym, reciprocal_size_2_sym

implicit none

private

    type, abstract, public :: Abstract_Weighted_Structure
        real(DP) :: alpha
        integer :: reci_numbers(num_dimensions)
        real(DP), dimension(:, :, :), allocatable :: weight
        complex(DP), dimension(:, :, :), allocatable :: structure
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Reciprocal_Lattice), pointer :: reciprocal_lattice => null()
        class(Abstract_Component_Coordinates), pointer :: component_positions => null()
        class(Abstract_Component_Dipolar_Moments), pointer :: component_dipolar_moments => null()
    contains
        procedure :: construct => Abstract_Weighted_Structure_construct
        procedure :: destroy => Abstract_Weighted_Structure_destroy
        procedure :: get_weight => Abstract_Weighted_Structure_get_weight
        procedure :: get_structure => Abstract_Weighted_Structure_get_structure
        procedure, private :: set_weight => Abstract_Weighted_Structure_set_weight
        procedure, private :: set_structure => Abstract_Weighted_Structure_set_structure
    end type Abstract_Weighted_Structure

    type, extends(Abstract_Weighted_Structure) :: Concrete_Weighted_Structure

    end type Concrete_Weighted_Structure

    type, extends(Abstract_Weighted_Structure) :: Null_Weighted_Structure

    end type Null_Weighted_Structure

contains

!implementation Abstract_Weighted_Structure

    subroutine Abstract_Weighted_Structure_construct(this, periodic_box, reciprocal_lattice, &
        component_positions, component_dipolar_moments, alpha)
        class(Abstract_Weighted_Structure), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), target, intent(in) :: reciprocal_lattice
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moments
        real(DP), intent(in) :: alpha

        this%periodic_box => periodic_box
        this%reciprocal_lattice => reciprocal_lattice
        this%reci_numbers = this%reciprocal_lattice%get_numbers()
        this%component_positions => component_positions
        this%component_dipolar_moments => component_dipolar_moments
        call check_positive("Abstract_Weighted_Structure_construct", "alpha", alpha)
        this%alpha = alpha
        allocate(this%weight(-this%reci_numbers(1):this%reci_numbers(1), &
                             -this%reci_numbers(2):this%reci_numbers(2), &
                             -this%reci_numbers(3):this%reci_numbers(3)))
        call this%set_weight()
        allocate(this%structure(-this%reci_numbers(1):this%reci_numbers(1), &
                                -this%reci_numbers(2):this%reci_numbers(2), &
                                -this%reci_numbers(3):this%reci_numbers(3)))
        call this%set_structure()
    end subroutine Abstract_Weighted_Structure_construct

    !> \[
    !>      w(\alpha, \vec{k}) = \frac{e^{-k^2/4\alpha^2}}{k^2}
    !> \]
    pure subroutine Abstract_Weighted_Structure_set_weight(this)
        class(Abstract_Weighted_Structure), intent(inout) :: this

        real(DP) :: box_size(num_dimensions)
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_squared

        box_size = this%periodic_box%get_size()
        do n_3 = -this%reci_numbers(3), this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
        do n_2 = -this%reci_numbers(2), this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
        do n_1 = -this%reci_numbers(1), this%reci_numbers(1)
            wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)
            if (n_1 /= 0 .or. n_2 /= 0 .or. n_3 /= 0) then
                wave_squared = dot_product(wave_vector, wave_vector)
                this%weight(n_1, n_2, n_3) = exp(-wave_squared/this%alpha**2/4._DP) / wave_squared
            else
                this%weight(n_1, n_2, n_3) = 0._DP
            end if
        end do
        end do
        end do
    end subroutine Abstract_Weighted_Structure_set_weight

    pure real(DP) function Abstract_Weighted_Structure_get_weight(this, n_1, n_2, n_3) &
        result(weight)
        class(Abstract_Weighted_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        weight = this%weight(n_1, n_2, n_3)
    end function Abstract_Weighted_Structure_get_weight

    !> Structure factor:
    !> \[
    !>      S(\vec{k}) = \sum_{i=1}^N (\vec{k}\cdot\vec{\mu}_i) e^{i\vec{k}\cdot\vec{x}_i}
    !> \]
    pure subroutine Abstract_Weighted_Structure_set_structure(this)
        class(Abstract_Weighted_Structure), intent(inout) :: this

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_dot_positions
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_moment
        integer :: n_1, n_2, n_3
        integer :: i_particle

        box_size = this%periodic_box%get_size()
        this%structure(:, :, :) = cmplx(0._DP, 0._DP, DP)
        do i_particle = 1, this%component_positions%get_num()
            wave_dot_positions(:) = 2._DP*PI * this%component_positions%get(i_particle) / &
                this%periodic_box%get_size()
            call set_fourier(fourier_position_1, this%reci_numbers(1), wave_dot_positions(1))
            call set_fourier(fourier_position_2, this%reci_numbers(2), wave_dot_positions(2))
            call set_fourier(fourier_position_3, this%reci_numbers(3), wave_dot_positions(3))
            do n_3 = -this%reci_numbers(3), this%reci_numbers(3)
                wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
            do n_2 = -this%reci_numbers(2), this%reci_numbers(2)
                wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
            do n_1 = -this%reci_numbers(1), this%reci_numbers(1)
                wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)
                fourier_position = fourier_position_1(n_1) * fourier_position_2(n_2) * &
                    fourier_position_3(n_3)
                wave_dot_moment = dot_product(wave_vector, this%component_dipolar_moments%&
                    get(i_particle))
                this%structure(n_1, n_2, n_3) = this%structure(n_1, n_2, n_3) + &
                    cmplx(wave_dot_moment, 0._DP, DP) * fourier_position
            end do
            end do
            end do
        end do
    end subroutine Abstract_Weighted_Structure_set_structure

    pure complex(DP) function Abstract_Weighted_Structure_get_structure(this, n_1, n_2, n_3) &
        result(structure)
        class(Abstract_Weighted_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        structure = this%structure(n_1, n_2, n_3)
    end function Abstract_Weighted_Structure_get_structure

    subroutine Abstract_Weighted_Structure_destroy(this)
        class(Abstract_Weighted_Structure), intent(inout) :: this

        if (allocated(this%structure)) deallocate(this%structure)
        if (allocated(this%weight)) deallocate(this%weight)
        this%component_dipolar_moments => null()
        this%component_positions => null()
        this%reciprocal_lattice => null()
        this%periodic_box => null()
    end subroutine Abstract_Weighted_Structure_destroy

!end implementation Abstract_Weighted_Structure

end module class_weighted_structure
