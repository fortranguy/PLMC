module classes_des_reci_structure

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_box_volume_memento, only: Abstract_Box_Volume_Memento
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use classes_structure_factor, only: Abstract_Structure_Factor
use procedures_dipolar_interactions_micro, only: set_fourier, reci_number_1_sym, reci_number_2_sym

implicit none

private

    !> @note
    !> Only half structure factor is stored in [[Abstract_DES_Reci_Structure:structure]] (symmetry).
    !> Hence in [[classes_des_reci_structure:Abstract_update_translation]],
    !> [[classes_des_reci_structure:Abstract_update_transmulation]],
    !> [[classes_des_reci_structure:Abstract_update_switch]],
    !> [[classes_des_reci_structure:Abstract_update_exchange]]
    !> the updates only run over half wave vectors.
    type, extends(Abstract_Structure_Factor), abstract, public :: Abstract_DES_Reci_Structure
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        class(Abstract_Box_Volume_Memento), pointer :: box_volume_memento => null()
        integer :: reci_numbers(num_dimensions) = 0
        type(Component_Wrapper), pointer :: components(:) => null()
        logical, allocatable :: are_dipolar(:)
        complex(DP), allocatable :: structure(:, :, :)
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: reset => Abstract_set
        procedure :: is_dipolar => Abstract_is_dipolar
        procedure :: get => Abstract_get
        procedure :: update_translation => Abstract_update_translation
        procedure :: update_transmutation => Abstract_update_transmulation
        procedure :: update_rotation => Abstract_update_rotation
        procedure :: update_add => Abstract_update_add
        procedure :: update_remove => Abstract_update_remove
        procedure :: update_switch => Abstract_update_switch
        procedure, private :: set => Abstract_set
        procedure, private :: update_exchange => Abstract_update_exchange
    end type Abstract_DES_Reci_Structure

    type, extends(Abstract_DES_Reci_Structure), public :: Concrete_DES_Reci_Structure

    end type Concrete_DES_Reci_Structure

    type, extends(Abstract_DES_Reci_Structure), public :: Null_DES_Reci_Structure
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: reset => Null_set
        procedure :: is_dipolar => Null_is_dipolar
        procedure :: get => Null_get
        procedure :: update_translation => Null_update_translation
        procedure :: update_transmutation => Null_update_transmulation
        procedure :: update_switch => Null_update_switch
        procedure, private :: update_exchange => Null_update_exchange
    end type Null_DES_Reci_Structure

contains

!implementation Abstract_DES_Reci_Structure

    subroutine Abstract_construct(this, periodic_box, box_volume_memento, reciprocal_lattice, &
        components, are_dipolar)
        class(Abstract_DES_Reci_Structure), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Box_Volume_Memento), target, intent(in) :: box_volume_memento
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), target, intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        this%periodic_box => periodic_box
        this%box_volume_memento => box_volume_memento
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%components => components
        allocate(this%are_dipolar(size(are_dipolar)))
        this%are_dipolar = are_dipolar
        allocate(this%structure(-this%reci_numbers(1):this%reci_numbers(1), &
                                -this%reci_numbers(2):this%reci_numbers(2), &
                                                    0:this%reci_numbers(3)))
        call this%set()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this

        if (allocated(this%structure)) deallocate(this%structure)
        if (allocated(this%are_dipolar)) deallocate(this%are_dipolar)
        this%components => null()
        this%box_volume_memento => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    pure logical function Abstract_is_dipolar(this, i_component) result(is_dipolar)
        class(Abstract_DES_Reci_Structure), intent(in) :: this
        integer, intent(in) :: i_component

        is_dipolar = this%are_dipolar(i_component)
    end function Abstract_is_dipolar

    pure subroutine Abstract_set(this)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this

        type(Concrete_Temporary_Particle) :: particle
        integer :: i_component, i_particle

        this%structure  = cmplx(0._DP, 0._DP, DP)
        do i_component = 1, size(this%components)
            do i_particle = 1, this%components(i_component)%dipole_moments%get_num()
                particle%position = this%components(i_component)%positions%get(i_particle)
                particle%dipole_moment = this%components(i_component)%dipole_moments%&
                    get(i_particle)
                call this%update_add(i_component, particle)
            end do
        end do
    end subroutine Abstract_set

    !> Structure factor:
    !> \[
    !>      S(\vec{k}) = \sum_{\vec{x}, \vec{\mu}} (\vec{k}\cdot\vec{\mu}) e^{i\vec{k}\cdot\vec{x}}
    !> \]
    pure complex(DP) function Abstract_get(this, n_1, n_2, n_3) result(structure)
        class(Abstract_DES_Reci_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3

        structure = this%structure(n_1, n_2, n_3)
    end function Abstract_get

    !> Structure factor update when a particle is translated:
    !> \( (\vec{x}, \vec{\mu}) \to (\vec{x}^\prime, \vec{\mu}) \).
    !>  \[
    !>      \Delta S(\vec{k}) = (\vec{k}\cdot\vec{\mu})
    !>          \left( e^{i\vec{k}\cdot\vec{x}^\prime} - e^{i\vec{k}\cdot\vec{x}} \right)
    !>  \]
    pure subroutine Abstract_update_translation(this, i_component, new_position, old)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: old

        real(DP) :: box_size(num_dimensions), box_edge_ratio
        real(DP), dimension(num_dimensions) :: wave_1_x_position_new, wave_1_x_position_old, &
            wave_vector
        integer :: n_1, n_2, n_3

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
        box_edge_ratio = this%box_volume_memento%get_ratio()**(-1._DP/3._DP)

        wave_1_x_position_old = 2._DP*PI * old%position / box_size
        call set_fourier(fourier_position_old_1, this%reci_numbers(1), wave_1_x_position_old(1))
        call set_fourier(fourier_position_old_2, this%reci_numbers(2), wave_1_x_position_old(2))
        call set_fourier(fourier_position_old_3, this%reci_numbers(3), wave_1_x_position_old(3))
        wave_1_x_position_new = 2._DP*PI * new_position / box_size
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
                        box_edge_ratio * dot_product(wave_vector, old%dipole_moment) * &
                        (fourier_position_new - fourier_position_old)
                end do
            end do
        end do
    end subroutine Abstract_update_translation

    !> Structure factor update when a particle is transmuted:
    !> \( (\vec{x}, \vec{\mu}) \to (\vec{x}, \vec{\mu}^\prime) \).
    !>  \[
    !>      \Delta S(\vec{k}) = \vec{k}\cdot(\vec{\mu}^\prime - \vec{\mu}) e^{i\vec{k}\cdot\vec{x}}
    !>  \]
    pure subroutine Abstract_update_transmulation(this, ij_components, new_dipole_moment, old)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        real(DP), intent(in) :: new_dipole_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old

        real(DP) :: box_size(num_dimensions), box_edge_ratio
        real(DP), dimension(num_dimensions) :: wave_1_x_position, wave_vector
        integer :: n_1, n_2, n_3

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        if (.not.(this%are_dipolar(ij_components(1)) .or. this%are_dipolar(ij_components(2)))) &
            return

        box_size = this%periodic_box%get_size()
        box_edge_ratio = this%box_volume_memento%get_ratio()**(-1._DP/3._DP)

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
                        box_edge_ratio * dot_product(wave_vector, &
                            new_dipole_moment - old%dipole_moment) * fourier_position
                end do
            end do
        end do
    end subroutine Abstract_update_transmulation

    pure subroutine Abstract_update_rotation(this, i_component, new_dipole_moment, old)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_dipole_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old

        call this%update_transmutation([i_component, i_component], new_dipole_moment, old)
    end subroutine Abstract_update_rotation

    !> cf. [[classes_des_reci_structure:Abstract_update_exchange]]
    pure subroutine Abstract_update_add(this, i_component, particle)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        call this%update_exchange(i_component, particle, +1._DP)
    end subroutine Abstract_update_add

    !> cf. [[classes_des_reci_structure:Abstract_update_exchange]]
    pure subroutine Abstract_update_remove(this, i_component, particle)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        call this%update_exchange(i_component, particle, -1._DP)
    end subroutine Abstract_update_remove

    !> Structure factor update when a particle of coordinates \( (\vec{x}, \vec{\mu}) \) is added
    !> (\( + )\) or removed (\( - \)).
    !> \[
    !>      \Delta S(\vec{k}) = \pm (\vec{k}\cdot\vec{\mu}) e^{i\vec{k}\cdot\vec{x}}
    !> \]
    pure subroutine Abstract_update_exchange(this, i_component, particle, signed)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed

        real(DP) :: box_size(num_dimensions), box_edge_ratio
        real(DP), dimension(num_dimensions) :: wave_1_x_position, wave_vector
        integer :: n_1, n_2, n_3

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        if (.not.this%are_dipolar(i_component)) return

        box_size = this%periodic_box%get_size()
        box_edge_ratio = this%box_volume_memento%get_ratio()**(-1._DP/3._DP)

        wave_1_x_position = 2._DP*PI * particle%position / box_size
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
                        signed*box_edge_ratio * dot_product(wave_vector, particle%dipole_moment) * &
                        fourier_position
                end do
            end do
        end do
    end subroutine Abstract_update_exchange

    !> Structure factor update when 2 particles of coordinates \( (\vec{x}_1, \vec{\mu}_1) \) and
    !> \( (\vec{x}_2, \vec{\mu}_2) \) are switched.
    !> \[
    !>      \Delta S(\vec{k}) = \vec{k}\cdot(\vec{\mu}_1 - \vec{\mu}_2)
    !>          \left( e^{i\vec{k}\cdot\vec{x}_2} - e^{i\vec{k}\cdot\vec{x}_1} \right)
    !> \]
    pure subroutine Abstract_update_switch(this, ij_components, particles)
        class(Abstract_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        real(DP) :: box_size(num_dimensions), box_edge_ratio
        real(DP), dimension(num_dimensions) :: wave_1_x_position_1, wave_1_x_position_2, &
            wave_vector
        integer :: n_1, n_2, n_3

        complex(DP) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_1_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_1_3
        complex(DP) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_2_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_2_3

        if (.not.(this%are_dipolar(ij_components(1)) .or. this%are_dipolar(ij_components(2)))) &
            return

        box_size = this%periodic_box%get_size()
        box_edge_ratio = this%box_volume_memento%get_ratio()**(-1._DP/3._DP)

        wave_1_x_position_1 = 2._DP*PI * particles(1)%position / box_size
        call set_fourier(fourier_position_1_1, this%reci_numbers(1), wave_1_x_position_1(1))
        call set_fourier(fourier_position_1_2, this%reci_numbers(2), wave_1_x_position_1(2))
        call set_fourier(fourier_position_1_3, this%reci_numbers(3), wave_1_x_position_1(3))
        wave_1_x_position_2 = 2._DP*PI * particles(2)%position / box_size
        call set_fourier(fourier_position_2_1, this%reci_numbers(1), wave_1_x_position_2(1))
        call set_fourier(fourier_position_2_2, this%reci_numbers(2), wave_1_x_position_2(2))
        call set_fourier(fourier_position_2_3, this%reci_numbers(3), wave_1_x_position_2(3))

        do n_3 = 0, this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
            do n_2 = -reci_number_2_sym(this%reci_numbers, n_3), this%reci_numbers(2)
                wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
                do n_1 = -reci_number_1_sym(this%reci_numbers, n_3, n_2), this%reci_numbers(1)
                    wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)

                    if (n_1**2 + n_2**2 + n_3**2 > this%reci_numbers(1)**2) cycle

                    fourier_position_1 = fourier_position_1_1(n_1) * fourier_position_1_2(n_2) * &
                        fourier_position_1_3(n_3)
                    fourier_position_2 = fourier_position_2_1(n_1) * fourier_position_2_2(n_2) * &
                        fourier_position_2_3(n_3)

                    this%structure(n_1, n_2, n_3) = this%structure(n_1, n_2, n_3) + &
                        box_edge_ratio * dot_product(wave_vector, particles(1)%dipole_moment - &
                            particles(2)%dipole_moment) * (fourier_position_2 - fourier_position_1)
                end do
            end do
        end do
    end subroutine Abstract_update_switch

!end implementation Abstract_DES_Reci_Structure

!implementation Null_DES_Reci_Structure

    subroutine Null_construct(this, periodic_box, box_volume_memento, reciprocal_lattice, &
        components, are_dipolar)
        class(Null_DES_Reci_Structure), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Box_Volume_Memento), target, intent(in) :: box_volume_memento
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), target, intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DES_Reci_Structure), intent(inout) :: this
    end subroutine Null_destroy

    pure logical function Null_is_dipolar(this, i_component) result(is_dipolar)
        class(Null_DES_Reci_Structure), intent(in) :: this
        integer, intent(in) :: i_component
        is_dipolar = .false.
    end function Null_is_dipolar

    pure subroutine Null_set(this)
        class(Null_DES_Reci_Structure), intent(inout) :: this
    end subroutine Null_set

    pure complex(DP) function Null_get(this, n_1, n_2, n_3) result(structure)
        class(Null_DES_Reci_Structure), intent(in) :: this
        integer, intent(in) :: n_1, n_2, n_3
        structure = cmplx(0._DP, 0._DP, DP)
    end function Null_get

    pure subroutine Null_update_translation(this, i_component, new_position, old)
        class(Null_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: old
    end subroutine Null_update_translation

    pure subroutine Null_update_transmulation(this, ij_components, new_dipole_moment, old)
        class(Null_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        real(DP), intent(in) :: new_dipole_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old
    end subroutine Null_update_transmulation

    pure subroutine Null_update_rotation(this, i_component, new_dipole_moment, old)
        class(Null_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_dipole_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old
    end subroutine Null_update_rotation

    pure subroutine Null_update_exchange(this, i_component, particle, signed)
        class(Null_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed
    end subroutine Null_update_exchange

    pure subroutine Null_update_switch(this, ij_components, particles)
        class(Null_DES_Reci_Structure), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)
    end subroutine Null_update_switch

!end implementation Null_DES_Reci_Structure

end module classes_des_reci_structure
