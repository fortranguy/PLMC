!> display: public
!>          private
module class_ewald_reci_structures

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use procedures_checks, only: check_positive
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_component_wrapper, only: Component_Wrapper
use procedures_ewald_micro, only: set_fourier, reciprocal_size_1_sym, reciprocal_size_2_sym

implicit none

private

    type structure_wrapper
        complex(DP), dimension(:, :, :), allocatable :: structure
    end type structure_wrapper

    type, abstract, public :: Abstract_Ewald_Reci_Structures
        real(DP) :: alpha
        integer :: reci_numbers(num_dimensions)
        real(DP), dimension(:, :, :), allocatable :: weight
        complex(DP), dimension(:, :, :), allocatable :: structure
        type(structure_wrapper), allocatable :: structures(:)
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Reciprocal_Lattice), pointer :: reciprocal_lattice => null()
        type(Component_Wrapper), pointer :: components(:) => null()
    contains
        procedure :: construct => Abstract_Ewald_Reci_Structures_construct
        procedure :: destroy => Abstract_Ewald_Reci_Structures_destroy
        procedure, private :: set_weight => Abstract_Ewald_Reci_Structures_set_weight
        procedure, private :: allocate => Abstract_Ewald_Reci_Structures_allocate
        procedure, private :: set_structure => Abstract_Ewald_Reci_Structures_set_structure
        procedure, private :: deallocate => Abstract_Ewald_Reci_Structures_deallocate
    end type Abstract_Ewald_Reci_Structures

    type, extends(Abstract_Ewald_Reci_Structures) :: Concrete_Weighted_Structure

    end type Concrete_Weighted_Structure

    type, extends(Abstract_Ewald_Reci_Structures) :: Null_Weighted_Structure

    end type Null_Weighted_Structure

contains

!implementation Abstract_Ewald_Reci_Structures

    subroutine Abstract_Ewald_Reci_Structures_construct(this, periodic_box, reciprocal_lattice, &
        components, are_dipolar, alpha)
        class(Abstract_Ewald_Reci_Structures), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), target, intent(in) :: reciprocal_lattice
        type(Component_Wrapper), target, intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
        real(DP), intent(in) :: alpha

        this%periodic_box => periodic_box
        this%reciprocal_lattice => reciprocal_lattice
        this%reci_numbers = this%reciprocal_lattice%get_numbers()
        this%components => components
        call check_positive("Abstract_Ewald_Reci_Structures_construct", "alpha", alpha)
        this%alpha = alpha
        allocate(this%weight(-this%reci_numbers(1):this%reci_numbers(1), &
                             -this%reci_numbers(2):this%reci_numbers(2), &
                             -this%reci_numbers(3):this%reci_numbers(3)))
        call this%set_weight()
        call this%allocate(are_dipolar)
    end subroutine Abstract_Ewald_Reci_Structures_construct

    !> \[
    !>      w(\alpha, \vec{k}) = \frac{e^{-k^2/4\alpha^2}}{k^2}
    !> \]
    pure subroutine Abstract_Ewald_Reci_Structures_set_weight(this)
        class(Abstract_Ewald_Reci_Structures), intent(inout) :: this

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
    end subroutine Abstract_Ewald_Reci_Structures_set_weight

    subroutine Abstract_Ewald_Reci_Structures_allocate(this, are_dipolar)
        class(Abstract_Ewald_Reci_Structures), intent(inout) :: this
        logical, intent(in) :: are_dipolar(:)

        integer :: i_component

        allocate(this%structures(size(this%components)))
        do i_component = 1, size(this%components)
            if (are_dipolar(i_component)) then
                allocate(this%structures(i_component)%&
                    structure(-this%reci_numbers(1):this%reci_numbers(1), &
                              -this%reci_numbers(2):this%reci_numbers(2), &
                              -this%reci_numbers(3):this%reci_numbers(3)))
                call this%set_structure(i_component)
            end if
        end do
    end subroutine Abstract_Ewald_Reci_Structures_allocate

    !> Structure factor:
    !> \[
    !>      S(\vec{k}) = \sum_{i=1}^N (\vec{k}\cdot\vec{\mu}_i) e^{i\vec{k}\cdot\vec{x}_i}
    !> \]
    pure subroutine Abstract_Ewald_Reci_Structures_set_structure(this, i_component)
        class(Abstract_Ewald_Reci_Structures), intent(inout) :: this
        integer, intent(in) :: i_component

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_dot_position
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_moment
        integer :: n_1, n_2, n_3
        integer :: i_particle

        box_size = this%periodic_box%get_size()
        this%structures(i_component)%structure  = cmplx(0._DP, 0._DP, DP)
        do i_particle = 1, this%components(i_component)%positions%get_num()
            wave_dot_position(:) = 2._DP*PI * this%components(i_component)%positions%&
                get(i_particle) / this%periodic_box%get_size()
            call set_fourier(fourier_position_1, this%reci_numbers(1), wave_dot_position(1))
            call set_fourier(fourier_position_2, this%reci_numbers(2), wave_dot_position(2))
            call set_fourier(fourier_position_3, this%reci_numbers(3), wave_dot_position(3))
            do n_3 = -this%reci_numbers(3), this%reci_numbers(3)
                wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
            do n_2 = -this%reci_numbers(2), this%reci_numbers(2)
                wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
            do n_1 = -this%reci_numbers(1), this%reci_numbers(1)
                wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)
                fourier_position = fourier_position_1(n_1) * fourier_position_2(n_2) * &
                    fourier_position_3(n_3)
                wave_dot_moment = dot_product(wave_vector, this%components(i_component)%&
                    dipolar_moments%get(i_particle))
                this%structures(i_component)%structure(n_1, n_2, n_3) = &
                    this%structures(i_component)%structure(n_1, n_2, n_3)  + &
                    cmplx(wave_dot_moment, 0._DP, DP) * fourier_position
            end do
            end do
            end do
        end do
    end subroutine Abstract_Ewald_Reci_Structures_set_structure

    subroutine Abstract_Ewald_Reci_Structures_destroy(this)
        class(Abstract_Ewald_Reci_Structures), intent(inout) :: this

        call this%deallocate()
        if (allocated(this%weight)) deallocate(this%weight)
        this%components => null()
        this%reciprocal_lattice => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Reci_Structures_destroy

    subroutine Abstract_Ewald_Reci_Structures_deallocate(this)
        class(Abstract_Ewald_Reci_Structures), intent(inout) :: this

        integer :: i_component

        if (allocated(this%structures)) then
            do i_component = size(this%structures), 1, -1
                if (allocated(this%structures(i_component)%structure)) then
                    deallocate(this%structures(i_component)%structure)
                end if
            end do
        end if
    end subroutine Abstract_Ewald_Reci_Structures_deallocate

!end implementation Abstract_Ewald_Reci_Structures

end module class_ewald_reci_structures
