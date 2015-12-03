module class_ewald_reci_structure

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight
use procedures_ewald_micro, only: set_fourier, reci_number_1_sym, reci_number_2_sym

implicit none

private

    type, abstract, public :: Abstract_Ewald_Reci_Structure
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        class(Abstract_Component_Coordinates), pointer :: component_positions => null()
        class(Abstract_Component_Dipolar_Moments), pointer :: component_dipolar_moments => null()
        class(Abstract_Ewald_Reci_Weight), pointer :: weight => null()
        complex(DP), dimension(:, :, :), allocatable :: structure
    contains
        procedure :: construct => Abstract_Ewald_Reci_Structure_construct
        procedure :: destroy => Abstract_Ewald_Reci_Structure_destroy
        procedure :: reset => Abstract_Ewald_Reci_Structure_set
        procedure :: get => Abstract_Ewald_Reci_Structure_get
        procedure, private :: set => Abstract_Ewald_Reci_Structure_set
        procedure :: get_coordinates_delta => Abstract_Ewald_Reci_Structure_get_coordinates_delta
        procedure :: set_coordinates_delta => Abstract_Ewald_Reci_Structure_set_coordinates_delta
    end type Abstract_Ewald_Reci_Structure

    type, extends(Abstract_Ewald_Reci_Structure), public :: Concrete_Ewald_Reci_Structure

    end type Concrete_Ewald_Reci_Structure

    type, extends(Abstract_Ewald_Reci_Structure), public :: Null_Ewald_Reci_Structure
    contains
        procedure :: construct => Null_Ewald_Reci_Structure_construct
        procedure :: destroy => Null_Ewald_Reci_Structure_destroy
        procedure :: get => Null_Ewald_Reci_Structure_get
        procedure, private :: set => Null_Ewald_Reci_Structure_set
        procedure :: get_coordinates_delta => Null_Ewald_Reci_Structure_get_coordinates_delta
        procedure :: set_coordinates_delta => Null_Ewald_Reci_Structure_set_coordinates_delta
    end type Null_Ewald_Reci_Structure

contains

!implementation Abstract_Ewald_Reci_Structure

    subroutine Abstract_Ewald_Reci_Structure_construct(this, periodic_box, reciprocal_lattice, &
        component_positions, component_dipolar_moments, weight)
        class(Abstract_Ewald_Reci_Structure), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moments
        class(Abstract_Ewald_Reci_Weight), target, intent(in) :: weight

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%component_positions => component_positions
        this%component_dipolar_moments => component_dipolar_moments
        this%weight => weight
        allocate(this%structure(-this%reci_numbers(1):this%reci_numbers(1), &
                                -this%reci_numbers(2):this%reci_numbers(2), &
                                -this%reci_numbers(3):this%reci_numbers(3)))
        call this%set()
    end subroutine Abstract_Ewald_Reci_Structure_construct

    subroutine Abstract_Ewald_Reci_Structure_destroy(this)
        class(Abstract_Ewald_Reci_Structure), intent(inout) :: this

        if (allocated(this%structure)) deallocate(this%structure)
        this%weight => null()
        this%component_dipolar_moments => null()
        this%component_positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Reci_Structure_destroy

    pure subroutine Abstract_Ewald_Reci_Structure_set(this)
        class(Abstract_Ewald_Reci_Structure), intent(inout) :: this

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_1_x_position
        real(DP), dimension(num_dimensions) :: wave_vector
        real(DP) :: wave_dot_moment
        integer :: n_1, n_2, n_3
        integer :: i_particle

        box_size = this%periodic_box%get_size()
        this%structure  = cmplx(0._DP, 0._DP, DP)
        do i_particle = 1, this%component_positions%get_num()
            wave_1_x_position = 2._DP*PI * this%component_positions%get(i_particle) / box_size
            call set_fourier(fourier_position_1, this%reci_numbers(1), wave_1_x_position(1))
            call set_fourier(fourier_position_2, this%reci_numbers(2), wave_1_x_position(2))
            call set_fourier(fourier_position_3, this%reci_numbers(3), wave_1_x_position(3))
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
    end subroutine Abstract_Ewald_Reci_Structure_set

    !> Structure factor:
    !> \[
    !>      S(\vec{k}) = \sum_{i=1}^N (\vec{k}\cdot\vec{\mu}_i) e^{i\vec{k}\cdot\vec{x}_i}
    !> \]
    pure complex(DP) function Abstract_Ewald_Reci_Structure_get(this, n_1, n_2, n_3) &
        result(structure)
        class(Abstract_Ewald_Reci_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        structure = this%structure(n_1, n_2, n_3)
    end function Abstract_Ewald_Reci_Structure_get

    !> Energy delta when a particle \( i \) of component \( I \) changes its coordinates.
    !> \[
    !>      \Delta U_{I, J} = \sum_{\vec{k}} w_\alpha(k)
    !>          \Re[(s^{\prime{}I}_i(\vec{k}) - s^I_i(\vec{k})) S_J(\vec{k})^\ast] +
    !>          [I=J] \frac{1}{2}|s^{\prime{}I}_i(\vec{k}) - s^I_i(\vec{k})|^2
    !> \]
    !> with \( s_i(\vec{k}) = (\vec{k}\cdot\vec{\mu}_i) e^{i\vec{k}\cdot\vec{x}_i} \).
    !> Though \( \vec{x} \) and \( \vec{\mu} \) don't change simultaneously,
    !> only one procedure is provided for simplicity.
    pure real(DP) function Abstract_Ewald_Reci_Structure_get_coordinates_delta(this, new, old, &
        same_component) result(delta_energy)
        class(Abstract_Ewald_Reci_Structure), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        logical, intent(in) :: same_component

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        complex(DP) :: structure_wave_new, structure_wave_old
        real(DP), dimension(num_dimensions) :: wave_1_x_position_new, wave_1_x_position_old
        real(DP) :: wave_dot_moment_new, wave_dot_moment_old

        complex(DP) :: fourier_position_new
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_new_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_new_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_new_3

        complex(DP) :: fourier_position_old
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_old_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_old_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_old_3

        box_size = this%periodic_box%get_size()

        wave_1_x_position_old = 2._DP*PI * old%position / box_size
        call set_fourier(fourier_position_old_1, this%reci_numbers(1), wave_1_x_position_old(1))
        call set_fourier(fourier_position_old_2, this%reci_numbers(2), wave_1_x_position_old(2))
        call set_fourier(fourier_position_old_3, this%reci_numbers(3), wave_1_x_position_old(3))

        wave_1_x_position_new = 2._DP*PI * new%position / box_size
        call set_fourier(fourier_position_new_1, this%reci_numbers(1), wave_1_x_position_new(1))
        call set_fourier(fourier_position_new_2, this%reci_numbers(2), wave_1_x_position_new(2))
        call set_fourier(fourier_position_new_3, this%reci_numbers(3), wave_1_x_position_new(3))

        delta_energy = 0._DP
        do n_3 = 0, this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
            do n_2 = -reci_number_2_sym(this%reci_numbers, n_3), this%reci_numbers(2)
                wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
                do n_1 = -reci_number_1_sym(this%reci_numbers, n_3, n_2), this%reci_numbers(1)
                    wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)

                    fourier_position_old = fourier_position_old_1(n_1) * &
                        fourier_position_old_2(n_2) * fourier_position_old_3(n_3)
                    fourier_position_new = fourier_position_new_1(n_1) * &
                        fourier_position_new_2(n_2) * fourier_position_new_3(n_3)
                    wave_dot_moment_old = dot_product(wave_vector, old%dipolar_moment)
                    wave_dot_moment_new = dot_product(wave_vector, new%dipolar_moment)

                    structure_wave_new = wave_dot_moment_new * fourier_position_new
                    structure_wave_old = wave_dot_moment_old * fourier_position_old
                    delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                        real((structure_wave_new - structure_wave_old) * &
                            conjg(this%structure(n_1, n_2, n_3)), DP)

                    if (same_component) then
                        delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                            0.5_DP * abs(structure_wave_new - structure_wave_old)**2
                    end if
                end do
            end do
        end do
        delta_energy = 2._DP * delta_energy ! half wave vector (symmetry) -> double energy
    end function Abstract_Ewald_Reci_Structure_get_coordinates_delta

    !> Structure factor update when a particle \( i \) changes its coordinates.
    !>  \[ \Delta S(\vec{k}) = s_i^\prime(\vec{k}) - s_i(\vec{k}) \]
    pure subroutine Abstract_Ewald_Reci_Structure_set_coordinates_delta(this, new, old)
        class(Abstract_Ewald_Reci_Structure), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        complex(DP) :: structure_wave_new, structure_wave_old
        real(DP), dimension(num_dimensions) :: wave_1_x_position_new, wave_1_x_position_old
        real(DP) :: wave_dot_moment_new, wave_dot_moment_old

        complex(DP) :: fourier_position_new
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_new_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_new_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_new_3

        complex(DP) :: fourier_position_old
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_old_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_old_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_old_3

        box_size = this%periodic_box%get_size()

        wave_1_x_position_old = 2._DP*PI * old%position / box_size
        call set_fourier(fourier_position_old_1, this%reci_numbers(1), wave_1_x_position_old(1))
        call set_fourier(fourier_position_old_2, this%reci_numbers(2), wave_1_x_position_old(2))
        call set_fourier(fourier_position_old_3, this%reci_numbers(3), wave_1_x_position_old(3))

        wave_1_x_position_new = 2._DP*PI * new%position / box_size
        call set_fourier(fourier_position_new_1, this%reci_numbers(1), wave_1_x_position_new(1))
        call set_fourier(fourier_position_new_2, this%reci_numbers(2), wave_1_x_position_new(2))
        call set_fourier(fourier_position_new_3, this%reci_numbers(3), wave_1_x_position_new(3))

        ! Warning: only half wave vectors are updated
        do n_3 = 0, this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
            do n_2 = -reci_number_2_sym(this%reci_numbers, n_3), this%reci_numbers(2)
                wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
                do n_1 = -reci_number_1_sym(this%reci_numbers, n_3, n_2), this%reci_numbers(1)
                    wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)

                    fourier_position_old = fourier_position_old_1(n_1) * &
                        fourier_position_old_2(n_2) * fourier_position_old_3(n_3)
                    fourier_position_new = fourier_position_new_1(n_1) * &
                        fourier_position_new_2(n_2) * fourier_position_new_3(n_3)
                    wave_dot_moment_old = dot_product(wave_vector, old%dipolar_moment)
                    wave_dot_moment_new = dot_product(wave_vector, new%dipolar_moment)

                    structure_wave_new = wave_dot_moment_new * fourier_position_new
                    structure_wave_old = wave_dot_moment_old * fourier_position_old
                    this%structure(n_1, n_2, n_3) = this%structure(n_1, n_2, n_3) + &
                        structure_wave_new - structure_wave_old
                end do
            end do
        end do
    end subroutine Abstract_Ewald_Reci_Structure_set_coordinates_delta

!end implementation Abstract_Ewald_Reci_Structure

!implementation Null_Ewald_Reci_Structure

    subroutine Null_Ewald_Reci_Structure_construct(this, periodic_box, reciprocal_lattice, &
        component_positions, component_dipolar_moments, weight)
        class(Null_Ewald_Reci_Structure), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moments
        class(Abstract_Ewald_Reci_Weight), target, intent(in) :: weight
    end subroutine Null_Ewald_Reci_Structure_construct

    subroutine Null_Ewald_Reci_Structure_destroy(this)
        class(Null_Ewald_Reci_Structure), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Structure_destroy

    pure subroutine Null_Ewald_Reci_Structure_set(this)
        class(Null_Ewald_Reci_Structure), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Structure_set

    pure complex(DP) function Null_Ewald_Reci_Structure_get(this, n_1, n_2, n_3) &
        result(structure)
        class(Null_Ewald_Reci_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        structure = cmplx(0._DP, 0._DP, DP)
    end function Null_Ewald_Reci_Structure_get

    pure real(DP) function Null_Ewald_Reci_Structure_get_coordinates_delta(this, new, old, &
        same_component) result(delta_energy)
        class(Null_Ewald_Reci_Structure), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        logical, intent(in) :: same_component
        delta_energy = 0._DP
    end function Null_Ewald_Reci_Structure_get_coordinates_delta

    pure subroutine Null_Ewald_Reci_Structure_set_coordinates_delta(this, new, old)
        class(Null_Ewald_Reci_Structure), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
    end subroutine Null_Ewald_Reci_Structure_set_coordinates_delta

!end implementation Null_Ewald_Reci_Structure

end module class_ewald_reci_structure
