module class_ewald_reci_structure

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use procedures_ewald_micro, only: set_fourier, reci_number_1_sym, reci_number_2_sym

implicit none

private

    type, abstract, public :: Abstract_Ewald_Reci_Structure
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        type(Component_Wrapper), pointer :: components(:) => null()
        logical, allocatable :: are_dipolar(:)
        complex(DP), dimension(:, :, :), allocatable :: structure
    contains
        procedure :: construct => Abstract_Ewald_Reci_Structure_construct
        procedure :: destroy => Abstract_Ewald_Reci_Structure_destroy
        procedure :: reset => Abstract_Ewald_Reci_Structure_set
        procedure :: is_dipolar => Abstract_Ewald_Reci_Structure_is_dipolar
        procedure :: get => Abstract_Ewald_Reci_Structure_get
        procedure, private :: set => Abstract_Ewald_Reci_Structure_set
        procedure :: update_move => Abstract_Ewald_Reci_Structure_update_move
        procedure :: update_rotation => Abstract_Ewald_Reci_Structure_update_rotation
    end type Abstract_Ewald_Reci_Structure

    type, extends(Abstract_Ewald_Reci_Structure), public :: Concrete_Ewald_Reci_Structure

    end type Concrete_Ewald_Reci_Structure

    type, extends(Abstract_Ewald_Reci_Structure), public :: Null_Ewald_Reci_Structure
    contains
        procedure :: construct => Null_Ewald_Reci_Structure_construct
        procedure :: destroy => Null_Ewald_Reci_Structure_destroy
        procedure :: reset => Null_Ewald_Reci_Structure_set
        procedure :: is_dipolar => Null_Ewald_Reci_Structure_is_dipolar
        procedure :: get => Null_Ewald_Reci_Structure_get
        procedure :: update_move => Null_Ewald_Reci_Structure_update_move
        procedure :: update_rotation => Null_Ewald_Reci_Structure_update_rotation
    end type Null_Ewald_Reci_Structure

contains

!implementation Abstract_Ewald_Reci_Structure

    subroutine Abstract_Ewald_Reci_Structure_construct(this, periodic_box, reciprocal_lattice, &
        components, are_dipolar)
        class(Abstract_Ewald_Reci_Structure), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), target, intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%components => components
        allocate(this%are_dipolar(size(are_dipolar)))
        this%are_dipolar = are_dipolar
        allocate(this%structure(-this%reci_numbers(1):this%reci_numbers(1), &
                                -this%reci_numbers(2):this%reci_numbers(2), &
                                -this%reci_numbers(3):this%reci_numbers(3)))
        call this%set()
    end subroutine Abstract_Ewald_Reci_Structure_construct

    subroutine Abstract_Ewald_Reci_Structure_destroy(this)
        class(Abstract_Ewald_Reci_Structure), intent(inout) :: this

        if (allocated(this%structure)) deallocate(this%structure)
        if (allocated(this%are_dipolar)) deallocate(this%are_dipolar)
        this%components => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Reci_Structure_destroy

    pure logical function Abstract_Ewald_Reci_Structure_is_dipolar(this, i_component) &
        result(is_dipolar)
        class(Abstract_Ewald_Reci_Structure), intent(in) :: this
        integer, intent(in) :: i_component

        is_dipolar = this%are_dipolar(i_component)
    end function Abstract_Ewald_Reci_Structure_is_dipolar

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
        integer :: i_component, i_particle

        box_size = this%periodic_box%get_size()
        this%structure  = cmplx(0._DP, 0._DP, DP)
        do i_component = 1, size(this%components)
            do i_particle = 1, this%components(i_component)%dipolar_moments%get_num()
                wave_1_x_position = 2._DP*PI * this%components(i_component)%positions%&
                    get(i_particle) / box_size
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
                    wave_dot_moment = dot_product(wave_vector, this%components(i_component)%&
                        dipolar_moments%get(i_particle))

                    this%structure(n_1, n_2, n_3) = this%structure(n_1, n_2, n_3) + &
                        cmplx(wave_dot_moment, 0._DP, DP) * fourier_position
                end do
                end do
                end do
            end do
        end do
    end subroutine Abstract_Ewald_Reci_Structure_set

    !> Structure factor:
    !> \[
    !>      S(\vec{k}) = \sum_{\mathsf{I}, \mathsf{i}} (\vec{k}\cdot\vec{\mu}^\mathsf{I}_\mathsf{i})
    !>          e^{i\vec{k}\cdot\vec{x}^\mathsf{I}_\mathsf{i}}
    !> \]
    !> where \(\mathsf{I}\) indexes a component and \(\mathsf{i}\), a particle of \(\mathsf{I}\).
    pure complex(DP) function Abstract_Ewald_Reci_Structure_get(this, n_1, n_2, n_3) &
        result(structure)
        class(Abstract_Ewald_Reci_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        structure = this%structure(n_1, n_2, n_3)
    end function Abstract_Ewald_Reci_Structure_get

    !> Structure factor update when a particle \(\mathsf{i}\) of component \(\mathsf{I}\) moves.
    !>  \[ \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu}^\mathsf{I}_\mathsf{i})
    !>      (e^{i\vec{k}\cdot\vec{x}^{\mathsf{I}\prime}_\mathsf{i}} -
    !>       e^{i\vec{k}\cdot\vec{x}^\mathsf{I}_\mathsf{i}}) \]
    !> Warning: only half wave vectors are updated.
    pure subroutine Abstract_Ewald_Reci_Structure_update_move(this, i_component, new, old)
        class(Abstract_Ewald_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
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

        if (.not.this%are_dipolar(i_component)) return

        box_size = this%periodic_box%get_size()

        wave_1_x_position_old = 2._DP*PI * old%position / box_size
        call set_fourier(fourier_position_old_1, this%reci_numbers(1), wave_1_x_position_old(1))
        call set_fourier(fourier_position_old_2, this%reci_numbers(2), wave_1_x_position_old(2))
        call set_fourier(fourier_position_old_3, this%reci_numbers(3), wave_1_x_position_old(3))

        wave_1_x_position_new = 2._DP*PI * new%position / box_size
        call set_fourier(fourier_position_new_1, this%reci_numbers(1), wave_1_x_position_new(1))
        call set_fourier(fourier_position_new_2, this%reci_numbers(2), wave_1_x_position_new(2))
        call set_fourier(fourier_position_new_3, this%reci_numbers(3), wave_1_x_position_new(3))

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
    end subroutine Abstract_Ewald_Reci_Structure_update_move

    !> Structure factor update when a particle \(\mathsf{i}\) of component \(\mathsf{I}\) rotates.
    !>  \[ \Delta S(\vec{k}) = \vec{k}\cdot(\vec{\mu}^{\mathsf{I}\prime}_\mathsf{i} -
    !>      \vec{\mu}^\mathsf{I}_\mathsf{i}) e^{i\vec{k}\cdot\vec{x}^\mathsf{I}_\mathsf{i}} \]
    !> Warning: only half wave vectors are updated.
    pure subroutine Abstract_Ewald_Reci_Structure_update_rotation(this, i_component, new, old)
        class(Abstract_Ewald_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        if (.not.this%are_dipolar(i_component)) return

        box_size = this%periodic_box%get_size()

        wave_1_x_position = 2._DP*PI * old%position / box_size
        call set_fourier(fourier_position_1, this%reci_numbers(1), wave_1_x_position(1))
        call set_fourier(fourier_position_2, this%reci_numbers(2), wave_1_x_position(2))
        call set_fourier(fourier_position_3, this%reci_numbers(3), wave_1_x_position(3))

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
    end subroutine Abstract_Ewald_Reci_Structure_update_rotation

!end implementation Abstract_Ewald_Reci_Structure

!implementation Null_Ewald_Reci_Structure

    subroutine Null_Ewald_Reci_Structure_construct(this, periodic_box, reciprocal_lattice, &
        components, are_dipolar)
        class(Null_Ewald_Reci_Structure), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), target, intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
    end subroutine Null_Ewald_Reci_Structure_construct

    subroutine Null_Ewald_Reci_Structure_destroy(this)
        class(Null_Ewald_Reci_Structure), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Structure_destroy

    pure logical function Null_Ewald_Reci_Structure_is_dipolar(this, i_component) result(is_dipolar)
        class(Null_Ewald_Reci_Structure), intent(in) :: this
        integer, intent(in) :: i_component
        is_dipolar = .false.
    end function Null_Ewald_Reci_Structure_is_dipolar

    pure subroutine Null_Ewald_Reci_Structure_set(this)
        class(Null_Ewald_Reci_Structure), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Structure_set

    pure complex(DP) function Null_Ewald_Reci_Structure_get(this, n_1, n_2, n_3) &
        result(structure)
        class(Null_Ewald_Reci_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        structure = cmplx(0._DP, 0._DP, DP)
    end function Null_Ewald_Reci_Structure_get

    pure subroutine Null_Ewald_Reci_Structure_update_move(this, i_component, new, old)
        class(Null_Ewald_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old
    end subroutine Null_Ewald_Reci_Structure_update_move

    pure subroutine Null_Ewald_Reci_Structure_update_rotation(this, i_component, new, old)
        class(Null_Ewald_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old
    end subroutine Null_Ewald_Reci_Structure_update_rotation

!end implementation Null_Ewald_Reci_Structure

end module class_ewald_reci_structure
