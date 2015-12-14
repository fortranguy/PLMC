module class_ewald_reci_structure_

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use class_component_coordinates, only: Abstract_Component_Coordinates
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_temporary_particle, only: Concrete_Temporary_Particle
use procedures_ewald_micro, only: set_fourier, reci_number_1_sym, reci_number_2_sym

implicit none

private

    type, abstract, public :: Abstract_Ewald_Reci_Structure_
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        class(Abstract_Component_Coordinates), pointer :: component_positions => null()
        class(Abstract_Component_Dipolar_Moments), pointer :: component_dipolar_moments => null()
        complex(DP), dimension(:, :, :), allocatable :: structure
    contains
        procedure :: construct => Abstract_Ewald_Reci_Structure__construct
        procedure :: destroy => Abstract_Ewald_Reci_Structure__destroy
        procedure :: reset => Abstract_Ewald_Reci_Structure__set
        procedure :: get => Abstract_Ewald_Reci_Structure__get
        procedure, private :: set => Abstract_Ewald_Reci_Structure__set
        procedure :: set_move_delta => Abstract_Ewald_Reci_Structure__set_move_delta
        procedure :: set_rotation_delta => Abstract_Ewald_Reci_Structure__set_rotation_delta
    end type Abstract_Ewald_Reci_Structure_

    type, extends(Abstract_Ewald_Reci_Structure_), public :: Concrete_Ewald_Reci_Structure_

    end type Concrete_Ewald_Reci_Structure_

    type, extends(Abstract_Ewald_Reci_Structure_), public :: Null_Ewald_Reci_Structure_
    contains
        procedure :: construct => Null_Ewald_Reci_Structure__construct
        procedure :: destroy => Null_Ewald_Reci_Structure__destroy
        procedure :: reset => Null_Ewald_Reci_Structure__set
        procedure :: get => Null_Ewald_Reci_Structure__get
        procedure :: set_move_delta => Null_Ewald_Reci_Structure__set_move_delta
        procedure :: set_rotation_delta => Null_Ewald_Reci_Structure__set_rotation_delta
    end type Null_Ewald_Reci_Structure_

contains

!implementation Abstract_Ewald_Reci_Structure_

    subroutine Abstract_Ewald_Reci_Structure__construct(this, periodic_box, reciprocal_lattice, &
        component_positions, component_dipolar_moments)
        class(Abstract_Ewald_Reci_Structure_), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moments

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%component_positions => component_positions
        this%component_dipolar_moments => component_dipolar_moments
        allocate(this%structure(-this%reci_numbers(1):this%reci_numbers(1), &
                                -this%reci_numbers(2):this%reci_numbers(2), &
                                -this%reci_numbers(3):this%reci_numbers(3)))
        call this%set()
    end subroutine Abstract_Ewald_Reci_Structure__construct

    subroutine Abstract_Ewald_Reci_Structure__destroy(this)
        class(Abstract_Ewald_Reci_Structure_), intent(inout) :: this

        if (allocated(this%structure)) deallocate(this%structure)
        this%component_dipolar_moments => null()
        this%component_positions => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Reci_Structure__destroy

    pure subroutine Abstract_Ewald_Reci_Structure__set(this)
        class(Abstract_Ewald_Reci_Structure_), intent(inout) :: this

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

                if (n_1**2 + n_2**2 + n_3**2 > this%reci_numbers(1)**2) cycle

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
    end subroutine Abstract_Ewald_Reci_Structure__set

    !> Structure factor:
    !> \[
    !>      S(\vec{k}) = \sum_{\mathsf{i}=1}^N (\vec{k}\cdot\vec{\mu}_\mathsf{i})
    !>          e^{i\vec{k}\cdot\vec{x}_\mathsf{i}}
    !> \]
    pure complex(DP) function Abstract_Ewald_Reci_Structure__get(this, n_1, n_2, n_3) &
        result(structure)
        class(Abstract_Ewald_Reci_Structure_), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        structure = this%structure(n_1, n_2, n_3)
    end function Abstract_Ewald_Reci_Structure__get

    !> Structure factor update when a particle \( \mathsf{i} \) moves.
    !>  \[ \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}_\mathsf{i})
    !>      (e^{i\vec{k}\cdot\vec{x}^\prime_\mathsf{i}} - e^{i\vec{k}\cdot\vec{x}_\mathsf{i}}) \]
    pure subroutine Abstract_Ewald_Reci_Structure__set_move_delta(this, new, old)
        class(Abstract_Ewald_Reci_Structure_), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position_new, wave_1_x_position_old

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

                    if (n_1**2 + n_2**2 + n_3**2 > this%reci_numbers(1)**2) cycle

                    fourier_position_old = fourier_position_old_1(n_1) * &
                        fourier_position_old_2(n_2) * fourier_position_old_3(n_3)
                    fourier_position_new = fourier_position_new_1(n_1) * &
                        fourier_position_new_2(n_2) * fourier_position_new_3(n_3)

                    this%structure(n_1, n_2, n_3) = this%structure(n_1, n_2, n_3) + &
                        dot_product(wave_vector, old%dipolar_moment) * &
                        (fourier_position_new - fourier_position_old)
                end do
            end do
        end do
    end subroutine Abstract_Ewald_Reci_Structure__set_move_delta

    !> Structure factor update when a particle \( \mathsf{i} \) rotates.
    !>  \[ \Delta S(\vec{k}) = \vec{k}\cdot(\vec{\mu}^\prime_\mathsf{i} - \vec{\mu}_\mathsf{i})
    !>      e^{i\vec{k}\cdot\vec{x}_\mathsf{i}} \]
    pure subroutine Abstract_Ewald_Reci_Structure__set_rotation_delta(this, new, old)
        class(Abstract_Ewald_Reci_Structure_), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        box_size = this%periodic_box%get_size()

        wave_1_x_position = 2._DP*PI * old%position / box_size
        call set_fourier(fourier_position_1, this%reci_numbers(1), wave_1_x_position(1))
        call set_fourier(fourier_position_2, this%reci_numbers(2), wave_1_x_position(2))
        call set_fourier(fourier_position_3, this%reci_numbers(3), wave_1_x_position(3))

        ! Warning: only half wave vectors are updated
        do n_3 = 0, this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
            do n_2 = -reci_number_2_sym(this%reci_numbers, n_3), this%reci_numbers(2)
                wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
                do n_1 = -reci_number_1_sym(this%reci_numbers, n_3, n_2), this%reci_numbers(1)
                    wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)

                    if (n_1**2 + n_2**2 + n_3**2 > this%reci_numbers(1)**2) cycle

                    fourier_position = fourier_position_1(n_1) * fourier_position_2(n_2) * &
                        fourier_position_3(n_3)

                    this%structure(n_1, n_2, n_3) = this%structure(n_1, n_2, n_3) + &
                        dot_product(wave_vector, new%dipolar_moment - old%dipolar_moment) * &
                        fourier_position
                end do
            end do
        end do
    end subroutine Abstract_Ewald_Reci_Structure__set_rotation_delta

!end implementation Abstract_Ewald_Reci_Structure_

!implementation Null_Ewald_Reci_Structure_

    subroutine Null_Ewald_Reci_Structure__construct(this, periodic_box, reciprocal_lattice, &
        component_positions, component_dipolar_moments)
        class(Null_Ewald_Reci_Structure_), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Component_Coordinates), target, intent(in) :: component_positions
        class(Abstract_Component_Dipolar_Moments), target, intent(in) :: component_dipolar_moments
    end subroutine Null_Ewald_Reci_Structure__construct

    subroutine Null_Ewald_Reci_Structure__destroy(this)
        class(Null_Ewald_Reci_Structure_), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Structure__destroy

    pure subroutine Null_Ewald_Reci_Structure__set(this)
        class(Null_Ewald_Reci_Structure_), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Structure__set

    pure complex(DP) function Null_Ewald_Reci_Structure__get(this, n_1, n_2, n_3) &
        result(structure)
        class(Null_Ewald_Reci_Structure_), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        structure = cmplx(0._DP, 0._DP, DP)
    end function Null_Ewald_Reci_Structure__get

    pure subroutine Null_Ewald_Reci_Structure__set_move_delta(this, new, old)
        class(Null_Ewald_Reci_Structure_), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
    end subroutine Null_Ewald_Reci_Structure__set_move_delta

    pure subroutine Null_Ewald_Reci_Structure__set_rotation_delta(this, new, old)
        class(Null_Ewald_Reci_Structure_), intent(inout) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
    end subroutine Null_Ewald_Reci_Structure__set_rotation_delta

!end implementation Null_Ewald_Reci_Structure_

end module class_ewald_reci_structure_
