module classes_dlc_structures

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use classes_periodic_box, only: Abstract_Periodic_Box
use classes_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_component_wrapper, only: Component_Wrapper
use types_particle_wrapper, only: Concrete_Particle
use classes_structure_factor, only: Abstract_Structure_Factor
use procedures_dipolar_interactions_micro, only: set_fourier, set_exp_kz, reci_number_1_sym

implicit none

private

    type, extends(Abstract_Structure_Factor), abstract, public :: Abstract_DLC_Structures
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(2) = 0
        type(Component_Wrapper), pointer :: components(:) => null()
        logical, allocatable :: are_dipolar(:)
        complex(DP), dimension(:, :), allocatable :: structure_p, structure_m
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: target => Abstract_target
        procedure :: reset => Abstract_set
        procedure :: is_dipolar => Abstract_is_dipolar
        procedure :: get_plus => Abstract_get_plus
        procedure :: get_minus => Abstract_get_minus
        procedure :: update_translation => Abstract_update_translation
        procedure :: update_transmutation => Abstract_update_transmutation
        procedure :: update_rotation => Abstract_update_rotation
        procedure :: update_add => Abstract_update_add
        procedure :: update_remove => Abstract_update_remove
        procedure :: update_switch => Abstract_update_switch
        procedure, private :: set => Abstract_set
        procedure, private :: update_exchange => Abstract_update_exchange
    end type Abstract_DLC_Structures

    type, extends(Abstract_DLC_Structures), public :: Concrete_DLC_Structures

    end type Concrete_DLC_Structures

    type, extends(Abstract_DLC_Structures), public :: Null_DLC_Structures
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: target => Null_target
        procedure :: reset => Null_set
        procedure :: is_dipolar => Null_is_dipolar
        procedure :: get_plus => Null_get
        procedure :: get_minus => Null_get
        procedure :: update_translation => Null_update_translation
        procedure :: update_transmutation => Null_update_transmutation
        procedure :: update_switch => Null_update_switch
        procedure, private :: update_exchange => Null_update_exchange
    end type Null_DLC_Structures

contains

!implementation Abstract_DLC_Structures

    subroutine Abstract_construct(this, periodic_box, reciprocal_lattice, components, are_dipolar)
        class(Abstract_DLC_Structures), intent(out) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)

        call this%target(periodic_box, components)
        this%reci_numbers = reshape(reciprocal_lattice%get_numbers(), [2])
        allocate(this%are_dipolar(size(are_dipolar)))
        this%are_dipolar = are_dipolar
        allocate(this%structure_p(-this%reci_numbers(1):this%reci_numbers(1), &
                                                      0:this%reci_numbers(2)))
        allocate(this%structure_m(-this%reci_numbers(1):this%reci_numbers(1), &
                                                      0:this%reci_numbers(2)))
        call this%set()
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DLC_Structures), intent(inout) :: this

        if (allocated(this%structure_m)) deallocate(this%structure_m)
        if (allocated(this%structure_p)) deallocate(this%structure_p)
        if (allocated(this%are_dipolar)) deallocate(this%are_dipolar)
        this%components => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    subroutine Abstract_target(this, periodic_box, components)
        class(Abstract_DLC_Structures), intent(inout) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        type(Component_Wrapper), target, intent(in) :: components(:)

        this%periodic_box => periodic_box
        this%components => components
    end subroutine Abstract_target

    subroutine Abstract_set(this)
        class(Abstract_DLC_Structures), intent(inout) :: this

        type(Concrete_Particle) :: particle
        integer :: i_component, i_particle

        this%structure_p = cmplx(0._DP, 0._DP, DP)
        this%structure_m = cmplx(0._DP, 0._DP, DP)
        do i_component = 1, size(this%components)
            do i_particle = 1, this%components(i_component)%dipole_moments%get_num()
                particle%position = this%components(i_component)%positions%get(i_particle)
                particle%dipole_moment = this%components(i_component)%dipole_moments%&
                    get(i_particle)
                call this%update_add(i_component, particle)
            end do
        end do
    end subroutine Abstract_set

    pure logical function Abstract_is_dipolar(this, i_component) result(is_dipolar)
        class(Abstract_DLC_Structures), intent(in) :: this
        integer, intent(in) :: i_component

        is_dipolar = this%are_dipolar(i_component)
    end function Abstract_is_dipolar

    !> Structure factors
    !> \[
    !>      S_\pm(\vec{k}_{1:2}) = \sum_{\vec{x}, \vec{\mu}}
    !>          (\pm k_{1:2} \mu_3 + i \vec{k}_{1:2} \cdot \vec{\mu}_{1:2})
    !>          e^{\pm k_{1:2} x_3} e^{i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !> \]
    pure complex(DP) function Abstract_get_plus(this, n_1, n_2) result(structure_p)
        class(Abstract_DLC_Structures), intent(in) :: this
        integer, intent(in) :: n_1, n_2

        structure_p = this%structure_p(n_1, n_2)
    end function Abstract_get_plus

    !> Cf. [[Abstract_get_plus]]
    pure complex(DP) function Abstract_get_minus(this, n_1, n_2) result(structure_m)
        class(Abstract_DLC_Structures), intent(in) :: this
        integer, intent(in) :: n_1, n_2

        structure_m = this%structure_m(n_1, n_2)
    end function Abstract_get_minus

    !> Structure factors update when a particle is translated:
    !> \( (\vec{x}, \vec{\mu}) \to (\vec{x}^\prime, \vec{\mu}) \).
    !> \[
    !>      \Delta S_\pm(\vec{k}_{1:2}) =
    !>          (\pm k_{1:2} \mu_3 + i \vec{k}_{1:2} \cdot \vec{\mu}_{1:2})
    !>          \left(
    !>              e^{\pm k_{1:2} x^\prime_3} e^{i \vec{k}_{1:2} \cdot \vec{x}^\prime_{1:2}} -
    !>              e^{\pm k_{1:2} x_3} e^{i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !>          \right)
    !> \]
    !> Warning: only half wave vectors are updated.
    pure subroutine Abstract_update_translation(this, i_component, new_position, old)
        class(Abstract_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Particle), intent(in) :: old

        real(DP) :: surface_size(2)
        real(DP), dimension(2) :: wave_1_x_position_old, wave_1_x_position_new, wave_vector
        real(DP) :: wave_dot_moment_12, wave_x_moment_3
        integer :: n_1, n_2

        complex(DP) :: fourier_position_new
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_new_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_new_2
        real(DP) :: exp_kz_new
        real(DP), dimension(0:this%reci_numbers(1), 0:this%reci_numbers(2)) :: exp_kz_new_tab
        complex(DP) :: fourier_position_old
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_old_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_old_2
        real(DP) :: exp_kz_old
        real(DP), dimension(0:this%reci_numbers(1), 0:this%reci_numbers(2)) :: exp_kz_old_tab

        if (.not.this%are_dipolar(i_component)) return

        surface_size = reshape(this%periodic_box%get_size(), [2])
        wave_1_x_position_old = 2._DP*PI * old%position(1:2) / surface_size
        call set_fourier(fourier_position_old_1, this%reci_numbers(1), wave_1_x_position_old(1))
        call set_fourier(fourier_position_old_2, this%reci_numbers(2), wave_1_x_position_old(2))
        call set_exp_kz(exp_kz_old_tab, surface_size, old%position(3))
        wave_1_x_position_new = 2._DP*PI * new_position(1:2) / surface_size
        call set_fourier(fourier_position_new_1, this%reci_numbers(1), wave_1_x_position_new(1))
        call set_fourier(fourier_position_new_2, this%reci_numbers(2), wave_1_x_position_new(2))
        call set_exp_kz(exp_kz_new_tab, surface_size, new_position(3))

        do n_2 = 0, this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / surface_size(2)
            do n_1 = -reci_number_1_sym(this%reci_numbers, 0, n_2), this%reci_numbers(1)
                wave_vector(1) = 2._DP*PI * real(n_1, DP) / surface_size(1)

                if (n_1**2 + n_2**2 > this%reci_numbers(1)**2) cycle

                fourier_position_old = fourier_position_old_1(n_1) * fourier_position_old_2(n_2)
                exp_kz_old = exp_kz_old_tab(abs(n_1), abs(n_2))
                fourier_position_new = fourier_position_new_1(n_1) * fourier_position_new_2(n_2)
                exp_kz_new = exp_kz_new_tab(abs(n_1), abs(n_2))

                wave_dot_moment_12 = dot_product(wave_vector, old%dipole_moment(1:2))
                wave_x_moment_3 = norm2(wave_vector) * old%dipole_moment(3)

                this%structure_p(n_1, n_2) = this%structure_p(n_1, n_2) + &
                    cmplx(+wave_x_moment_3, wave_dot_moment_12, DP) * &
                    (fourier_position_new * exp_kz_new - fourier_position_old * exp_kz_old)
                this%structure_m(n_1, n_2) = this%structure_m(n_1, n_2) + &
                    cmplx(-wave_x_moment_3, wave_dot_moment_12, DP) * &
                    (fourier_position_new / exp_kz_new - fourier_position_old / exp_kz_old)
            end do
        end do
    end subroutine Abstract_update_translation

    !> Structure factors update when a particle is transmuted:
    !> \( (\vec{x}, \vec{\mu}) \to (\vec{x}, \vec{\mu}^\prime) \).
    !>  \[
    !>      \Delta S_\pm(\vec{k}_{1:2}) =
    !>          \left[
    !>              \pm k_{1:2} (\mu^\prime_3 - \mu_3) +
    !>              i \vec{k}_{1:2} \cdot (\vec{\mu}^\prime_{1:2} - \vec{\mu}_{1:2})
    !>          \right]
    !>          e^{\pm k_{1:2} x_3} e^{i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !>  \]
    !> Warning: only half wave vectors are updated.
    pure subroutine Abstract_update_transmutation(this, ij_components, new_dipole_moment, old)
        class(Abstract_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        real(DP), intent(in) :: new_dipole_moment(:)
        type(Concrete_Particle), intent(in) :: old

        real(DP) :: surface_size(2)
        real(DP), dimension(2) :: wave_1_x_position, wave_vector
        real(DP) :: wave_dot_delta_moment_12, wave_delta_moment_3
        integer :: n_1, n_2

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        real(DP) :: exp_kz
        real(DP), dimension(0:this%reci_numbers(1), 0:this%reci_numbers(2)) :: exp_kz_tab

        if (.not.(this%are_dipolar(ij_components(1)) .or. this%are_dipolar(ij_components(2)))) &
            return

        surface_size = reshape(this%periodic_box%get_size(), [2])
        wave_1_x_position = 2._DP*PI * old%position(1:2) / surface_size
        call set_fourier(fourier_position_1, this%reci_numbers(1), wave_1_x_position(1))
        call set_fourier(fourier_position_2, this%reci_numbers(2), wave_1_x_position(2))
        call set_exp_kz(exp_kz_tab, surface_size, old%position(3))

        do n_2 = 0, this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / surface_size(2)
            do n_1 = -reci_number_1_sym(this%reci_numbers, 0, n_2), this%reci_numbers(1)
                wave_vector(1) = 2._DP*PI * real(n_1, DP) / surface_size(1)

                if (n_1**2 + n_2**2 > this%reci_numbers(1)**2) cycle

                fourier_position = fourier_position_1(n_1) * fourier_position_2(n_2)
                exp_kz = exp_kz_tab(abs(n_1), abs(n_2))

                wave_dot_delta_moment_12 = dot_product(wave_vector, new_dipole_moment(1:2) - old%&
                    dipole_moment(1:2))
                wave_delta_moment_3 = norm2(wave_vector) * (new_dipole_moment(3) - old%&
                    dipole_moment(3))

                this%structure_p(n_1, n_2) = this%structure_p(n_1, n_2) + &
                    cmplx(+wave_delta_moment_3, wave_dot_delta_moment_12, DP) * &
                    (fourier_position * exp_kz)
                this%structure_m(n_1, n_2) = this%structure_m(n_1, n_2) + &
                    cmplx(-wave_delta_moment_3, wave_dot_delta_moment_12, DP) * &
                    (fourier_position / exp_kz)
            end do
        end do
    end subroutine Abstract_update_transmutation

    pure subroutine Abstract_update_rotation(this, i_component, new_dipole_moment, old)
        class(Abstract_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_dipole_moment(:)
        type(Concrete_Particle), intent(in) :: old

        call this%update_transmutation([i_component, i_component], new_dipole_moment, old)
    end subroutine Abstract_update_rotation

    !> cf. [[classes_dlc_structures:Abstract_update_exchange]]
    pure subroutine Abstract_update_add(this, i_component, particle)
        class(Abstract_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Particle), intent(in) :: particle

        call this%update_exchange(i_component, particle, +1._DP)
    end subroutine Abstract_update_add

    !> cf. [[classes_dlc_structures:Abstract_update_exchange]]
    pure subroutine Abstract_update_remove(this, i_component, particle)
        class(Abstract_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Particle), intent(in) :: particle

        call this%update_exchange(i_component, particle, -1._DP)
    end subroutine Abstract_update_remove

    !> Structure factors update when a particle of coordinates \( (\vec{x}, \vec{\mu}) \) is added
    !> (\( + )\) or removed (\( - \)):
    !> \[
    !>      \Delta S_\pm(\vec{k}_{1:2}) =
    !>          {\bf\pm} (\pm k_{1:2} \mu_3 + i \vec{k}_{1:2} \cdot \vec{\mu}_{1:2})
    !>          e^{\pm k_{1:2} x_3} e^{i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}.
    !> \]
    pure subroutine Abstract_update_exchange(this, i_component, particle, signed)
        class(Abstract_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed

        real(DP) :: surface_size(2)
        real(DP), dimension(2) :: wave_1_x_position, wave_vector
        real(DP) :: wave_dot_moment_12, wave_x_moment_3
        integer :: n_1, n_2

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        real(DP) :: exp_kz
        real(DP), dimension(0:this%reci_numbers(1), 0:this%reci_numbers(2)) :: exp_kz_tab

        if (.not.this%are_dipolar(i_component)) return

        surface_size = reshape(this%periodic_box%get_size(), [2])
        wave_1_x_position = 2._DP*PI * particle%position(1:2) / surface_size
        call set_fourier(fourier_position_1, this%reci_numbers(1), wave_1_x_position(1))
        call set_fourier(fourier_position_2, this%reci_numbers(2), wave_1_x_position(2))
        call set_exp_kz(exp_kz_tab, surface_size, particle%position(3))

        do n_2 = 0, this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / surface_size(2)
            do n_1 = -reci_number_1_sym(this%reci_numbers, 0, n_2), this%reci_numbers(1)
                wave_vector(1) = 2._DP*PI * real(n_1, DP) / surface_size(1)

                if (n_1**2 + n_2**2 > this%reci_numbers(1)**2) cycle

                fourier_position = fourier_position_1(n_1) * fourier_position_2(n_2)
                exp_kz = exp_kz_tab(abs(n_1), abs(n_2))
                wave_dot_moment_12 = dot_product(wave_vector, particle%dipole_moment(1:2))
                wave_x_moment_3 = norm2(wave_vector) * particle%dipole_moment(3)

                this%structure_p(n_1, n_2) = this%structure_p(n_1, n_2) + &
                    signed * cmplx(+wave_x_moment_3, wave_dot_moment_12, DP) * &
                    fourier_position * exp_kz
                this%structure_m(n_1, n_2) = this%structure_m(n_1, n_2) + &
                    signed * cmplx(-wave_x_moment_3, wave_dot_moment_12, DP) * &
                    fourier_position / exp_kz
            end do
        end do
    end subroutine Abstract_update_exchange

    !> Structure factors update when 2 particles of coordinates \( (\vec{x}_1, \vec{\mu}_1) \) and
    !> \( (\vec{x}_2, \vec{\mu}_2) \) are switched.
    !> \[
    !>      \Delta S_\pm(\vec{k}_{1:2}) =
    !>          \left[
    !>              \pm k_{1:2} (\mu_{1, 3} - \mu_{2, 3}) +
    !>              i \vec{k}_{1:2} \cdot (\vec{\mu}_{1, 1:2} - \vec{\mu}_{2, 1:2})
    !>          \right] \\
    !>          \left(
    !>              e^{\pm k_{1:2} x_{2, 3}} e^{i \vec{k}_{1:2} \cdot \vec{x}_{2, 1:2}} -
    !>              e^{\pm k_{1:2} x_{1, 3}} e^{i \vec{k}_{1:2} \cdot \vec{x}_{1, 1:2}}
    !>          \right)
    !> \]
    pure subroutine Abstract_update_switch(this, ij_components, particles)
        class(Abstract_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Particle), intent(in) :: particles(:)

        real(DP) :: surface_size(2)
        real(DP), dimension(2) :: wave_1_x_position_1, wave_1_x_position_2, wave_vector
        real(DP) :: wave_dot_delta_moment_12, wave_delta_moment_3
        integer :: n_1, n_2

        complex(DP) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_1_2
        real(DP) :: exp_kz_1
        real(DP), dimension(0:this%reci_numbers(1), 0:this%reci_numbers(2)) :: exp_kz_1_tab
        complex(DP) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_2_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2_2
        real(DP) :: exp_kz_2
        real(DP), dimension(0:this%reci_numbers(1), 0:this%reci_numbers(2)) :: exp_kz_2_tab

        if (.not.(this%are_dipolar(ij_components(1)) .or. this%are_dipolar(ij_components(2)))) &
            return

        surface_size = reshape(this%periodic_box%get_size(), [2])
        wave_1_x_position_1 = 2._DP*PI * particles(1)%position(1:2) / surface_size
        call set_fourier(fourier_position_1_1, this%reci_numbers(1), wave_1_x_position_1(1))
        call set_fourier(fourier_position_1_2, this%reci_numbers(2), wave_1_x_position_1(2))
        call set_exp_kz(exp_kz_1_tab, surface_size, particles(1)%position(3))
        wave_1_x_position_2 = 2._DP*PI * particles(2)%position(1:2) / surface_size
        call set_fourier(fourier_position_2_1, this%reci_numbers(1), wave_1_x_position_2(1))
        call set_fourier(fourier_position_2_2, this%reci_numbers(2), wave_1_x_position_2(2))
        call set_exp_kz(exp_kz_2_tab, surface_size, particles(2)%position(3))

        do n_2 = 0, this%reci_numbers(2)
            wave_vector(2) = 2._DP*PI * real(n_2, DP) / surface_size(2)
            do n_1 = -reci_number_1_sym(this%reci_numbers, 0, n_2), this%reci_numbers(1)
                wave_vector(1) = 2._DP*PI * real(n_1, DP) / surface_size(1)

                if (n_1**2 + n_2**2 > this%reci_numbers(1)**2) cycle

                fourier_position_1 = fourier_position_1_1(n_1) * fourier_position_1_2(n_2)
                exp_kz_1 = exp_kz_1_tab(abs(n_1), abs(n_2))
                fourier_position_2 = fourier_position_2_1(n_1) * fourier_position_2_2(n_2)
                exp_kz_2 = exp_kz_2_tab(abs(n_1), abs(n_2))

                wave_dot_delta_moment_12 = dot_product(wave_vector, &
                    particles(1)%dipole_moment(1:2) - particles(2)%dipole_moment(1:2))
                wave_delta_moment_3 = norm2(wave_vector) * &
                    (particles(1)%dipole_moment(3) - particles(2)%dipole_moment(3))

                this%structure_p(n_1, n_2) = this%structure_p(n_1, n_2) + &
                    cmplx(+wave_delta_moment_3, wave_dot_delta_moment_12, DP) * &
                    (fourier_position_2 * exp_kz_2 - fourier_position_1 * exp_kz_1)
                this%structure_m(n_1, n_2) = this%structure_m(n_1, n_2) + &
                    cmplx(-wave_delta_moment_3, wave_dot_delta_moment_12, DP) * &
                    (fourier_position_2 / exp_kz_2 - fourier_position_1 / exp_kz_1)
            end do
        end do
    end subroutine Abstract_update_switch

!end implementation Abstract_DLC_Structures

!implementation Null_DLC_Structures

    subroutine Null_construct(this, periodic_box, reciprocal_lattice, components, are_dipolar)
        class(Null_DLC_Structures), intent(out) :: this
        class(Abstract_Periodic_Box), intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DLC_Structures), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_target(this, periodic_box, components)
        class(Null_DLC_Structures), intent(inout) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        type(Component_Wrapper), target, intent(in) :: components(:)
    end subroutine Null_target

    subroutine Null_set(this)
        class(Null_DLC_Structures), intent(inout) :: this
    end subroutine Null_set

    pure logical function Null_is_dipolar(this, i_component) result(is_dipolar)
        class(Null_DLC_Structures), intent(in) :: this
        integer, intent(in) :: i_component
        is_dipolar = .false.
    end function Null_is_dipolar

    pure complex(DP) function Null_get(this, n_1, n_2) result(structure_pm)
        class(Null_DLC_Structures), intent(in) :: this
        integer, intent(in) :: n_1, n_2
        structure_pm = cmplx(0._DP, 0._DP, DP)
    end function Null_get

    pure subroutine Null_update_translation(this, i_component, new_position, old)
        class(Null_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Particle), intent(in) :: old
    end subroutine Null_update_translation

    pure subroutine Null_update_transmutation(this, ij_components, new_dipole_moment, old)
        class(Null_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        real(DP), intent(in) :: new_dipole_moment(:)
        type(Concrete_Particle), intent(in) :: old
    end subroutine Null_update_transmutation

    pure subroutine Null_update_exchange(this, i_component, particle, signed)
        class(Null_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: i_component
        type(Concrete_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed
    end subroutine Null_update_exchange

    pure subroutine Null_update_switch(this, ij_components, particles)
        class(Null_DLC_Structures), intent(inout) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Particle), intent(in) :: particles(:)
    end subroutine Null_update_switch

!end implementation Null_DLC_Structures

end module classes_dlc_structures
