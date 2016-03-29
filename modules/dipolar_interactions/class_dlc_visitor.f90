module class_dlc_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_dlc_weight, only: Abstract_DLC_Weight
use class_dlc_structures, only: Abstract_DLC_Structures
use procedures_dipolar_interactions_micro, only: set_fourier, set_exp_kz, reci_number_1_sym

implicit none

private

    type, abstract, public :: Abstract_DLC_Visitor
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        class(Abstract_DLC_Weight), pointer :: weight => null()
        class(Abstract_DLC_Structures), pointer :: structures => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: visit => Abstract_visit
        procedure :: visit_move => Abstract_visit_move
        procedure :: visit_rotation => Abstract_visit_rotation
        procedure :: visit_switch => Abstract_visit_switch
    end type Abstract_DLC_Visitor

    type, extends(Abstract_DLC_Visitor), public :: Concrete_DLC_Visitor

    end type Concrete_DLC_Visitor

    type, extends(Abstract_DLC_Visitor), public :: Null_DLC_Visitor
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
        procedure :: visit_move => Null_visit_move
        procedure :: visit_rotation => Null_visit_rotation
        procedure :: visit_switch => Null_visit_switch
    end type Null_DLC_Visitor

contains

!implementation Abstract_DLC_Visitor

    subroutine Abstract_construct(this, periodic_box, reciprocal_lattice, weight, structures)
        class(Abstract_DLC_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_DLC_Weight), target, intent(in) :: weight
        class(Abstract_DLC_Structures), target, intent(in) :: structures

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%weight => weight
        this%structures => structures
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_DLC_Visitor), intent(inout) :: this

        this%structures => null()
        this%weight => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    !> \[
    !>      U = \sum_{\vec{k}_{1:2}} w(\vec{k}_{1:2})
    !>          \Re[S_+(\vec{k}_{1:2}) S_-(\vec{k}_{1:2})^\ast]
    !> \]
    !> where \( w(\vec{k}_{1:2}) \) is [[Abstract_DLC_Weight]] and
    !> \( S_\pm(\vec{k}_{1:2}) \) is [[Abstract_DLC_Structures]].
    pure real(DP) function Abstract_visit(this) result(energy)
        class(Abstract_DLC_Visitor), intent(in) :: this

        integer :: n_1, n_2

        energy = 0._DP
        do n_2 = 0, this%reci_numbers(2)
            do n_1 = -reci_number_1_sym(this%reci_numbers, 0, n_2), this%reci_numbers(1)
                if (n_1**2 + n_2**2 > this%reci_numbers(1)**2) cycle
                energy = energy + this%weight%get(n_1, n_2) * &
                    real(this%structures%get_plus(n_1, n_2) * &
                        conjg(this%structures%get_minus(n_1, n_2)), DP)
            end do
        end do
        energy = 2._DP * energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit

    !> Energy delta when a particle of coordinates \( (\vec{x}, \vec{\mu}) \) moves.
    !> \[
    !>      \Delta U = \sum_{\vec{k}_{1:2}} w(\vec{k}_{1:2}) \Re\left[
    !>          S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2}) +
    !>          S_-^\ast(\vec{k}_{1:2}) \Delta S_+(\vec{k}_{1:2}) +
    !>          \Delta S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2})
    !>      \right]
    !> \]
    !> where:
    !> \[
    !>      \Delta S_-^\ast(\vec{k}_{1:2}) =
    !>          (-k_{1:2} \mu_3 - i \vec{k}_{1:2} \cdot \vec{\mu}_{1:2})
    !>          \left(
    !>              e^{-k_{1:2} x^\prime_3} e^{-i \vec{k}_{1:2} \cdot \vec{x}^\prime_{1:2}} -
    !>              e^{-k_{1:2} x_3} e^{-i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !>          \right)
    !> \]
    !> \[
    !>      \Delta S_+(\vec{k}_{1:2}) = (+k_{1:2} \mu_3 + i \vec{k}_{1:2} \cdot \vec{\mu}_{1:2})
    !>          \left(
    !>              e^{+k_{1:2} x^\prime_3} e^{+i \vec{k}_{1:2} \cdot \vec{x}^\prime_{1:2}} -
    !>              e^{+k_{1:2} x_3} e^{+i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !>          \right)
    !> \]
    !> \[
    !>      \Delta S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2}) =
    !>          \left[
    !>              -(k_{1:2} \mu_3)^2 + (\vec{k}_{1:2} \cdot \vec{\mu}_{1:2})^2 -
    !>              2 i k_{1:2} \mu_3 \vec{k}_{1:2} \cdot \vec{\mu}_{1:2}
    !>          \right]
    !>          \left(
    !>              2 - \\
    !>              e^{+k_{1:2} x^\prime_3} e^{+i \vec{k}_{1:2} \cdot \vec{x}^\prime_{1:2}}
    !>              e^{-k_{1:2} x_3} e^{-i \vec{k}_{1:2} \cdot \vec{x}_{1:2}} - \\
    !>              e^{-k_{1:2} x^\prime_3} e^{-i \vec{k}_{1:2} \cdot \vec{x}^\prime_{1:2}}
    !>              e^{+k_{1:2} x_3} e^{+i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !>          \right)
    !> \]
    pure real(DP) function Abstract_visit_move(this, i_component, new_position, old) &
        result(delta_energy)
        class(Abstract_DLC_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: old

        real(DP) :: real_part_1, real_part_2, real_part_3

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

        delta_energy = 0._DP
        if (.not.this%structures%is_dipolar(i_component)) return

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

                wave_dot_moment_12 = dot_product(wave_vector, old%dipolar_moment(1:2))
                wave_x_moment_3 = norm2(wave_vector) * old%dipolar_moment(3)

                real_part_1 = real(this%structures%get_plus(n_1, n_2) * &
                    cmplx(-wave_x_moment_3, -wave_dot_moment_12, DP) * &
                    conjg(fourier_position_new/exp_kz_new - fourier_position_old/exp_kz_old), DP)
                real_part_2 = real(conjg(this%structures%get_minus(n_1, n_2)) * &
                    cmplx(+wave_x_moment_3, +wave_dot_moment_12, DP) * &
                    (fourier_position_new * exp_kz_new - fourier_position_old * exp_kz_old), DP)
                real_part_3 = real(cmplx(-wave_x_moment_3**2 + wave_dot_moment_12**2, &
                    -2._DP*wave_x_moment_3 * wave_dot_moment_12, DP) * (2._DP - &
                    fourier_position_new*exp_kz_new * conjg(fourier_position_old)/exp_kz_old - &
                    fourier_position_old*exp_kz_old * conjg(fourier_position_new)/exp_kz_new), DP)

                delta_energy = delta_energy + this%weight%get(n_1, n_2) * &
                    (real_part_1 + real_part_2 + real_part_3)
            end do
        end do
        delta_energy = 2._DP * delta_energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit_move

    !> Energy delta when a particle of coordinates \( (\vec{x}, \vec{\mu}) \) rotates.
    !> \[
    !>      \Delta U = \sum_{\vec{k}_{1:2}} w(\vec{k}_{1:2}) \Re\left[
    !>          S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2}) +
    !>          S_-^\ast(\vec{k}_{1:2}) \Delta S_+(\vec{k}_{1:2}) +
    !>          \Delta S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2})
    !>      \right]
    !> \]
    !> where:
    !> \[
    !>      \Delta S_-^\ast(\vec{k}_{1:2}) =
    !>          \left[
    !>              -k_{1:2} (\mu^\prime_3 - \mu_3) -
    !>              i \vec{k}_{1:2} \cdot (\vec{\mu}^\prime_{1:2} - \vec{\mu}_{1:2})
    !>          \right]
    !>          e^{-k_{1:2} x_3} e^{-i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !> \]
    !> \[
    !>      \Delta S_+(\vec{k}_{1:2}) =
    !>          \left[
    !>              +k_{1:2} (\mu^\prime_3 - \mu_3) +
    !>              i \vec{k}_{1:2} \cdot (\vec{\mu}^\prime_{1:2} - \vec{\mu}_{1:2})
    !>          \right]
    !>          e^{+k_{1:2} x_3} e^{+i \vec{k}_{1:2} \cdot \vec{x}_{1:2}}
    !> \]
    !> \[
    !>      \Re[\Delta S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2})] =
    !>          -[k_{1:2} (\mu^\prime_3 - \mu_3)]^2 +
    !>          [\vec{k}_{1:2} \cdot (\vec{\mu}^\prime_{1:2} - \vec{\mu}_{1:2})]^2
    !> \]
    pure real(DP) function Abstract_visit_rotation(this, i_component, new_dipolar_moment, old) &
        result(delta_energy)
        class(Abstract_DLC_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_dipolar_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old

        real(DP) :: real_part_1, real_part_2, real_part_3

        real(DP) :: surface_size(2)
        real(DP), dimension(2) :: wave_1_x_position, wave_vector
        real(DP) :: wave_dot_delta_moment_12, wave_delta_moment_3
        integer :: n_1, n_2

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        real(DP) :: exp_kz
        real(DP), dimension(0:this%reci_numbers(1), 0:this%reci_numbers(2)) :: exp_kz_tab

        delta_energy = 0._DP
        if (.not.this%structures%is_dipolar(i_component)) return

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

                wave_dot_delta_moment_12 = dot_product(wave_vector, new_dipolar_moment(1:2) - old%&
                    dipolar_moment(1:2))
                wave_delta_moment_3 = norm2(wave_vector) * (new_dipolar_moment(3) - old%&
                    dipolar_moment(3))

                real_part_1 = real(this%structures%get_plus(n_1, n_2) * &
                    cmplx(-wave_delta_moment_3, -wave_dot_delta_moment_12, DP) * &
                    conjg(fourier_position) / exp_kz, DP)
                real_part_2 = real(conjg(this%structures%get_minus(n_1, n_2)) * &
                    cmplx(+wave_delta_moment_3, +wave_dot_delta_moment_12, DP) * &
                    fourier_position * exp_kz, DP)
                real_part_3 = -wave_delta_moment_3**2 + wave_dot_delta_moment_12**2

                delta_energy = delta_energy + this%weight%get(n_1, n_2) * &
                    (real_part_1 + real_part_2 + real_part_3)
            end do
        end do
        delta_energy = 2._DP * delta_energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit_rotation

    !> Energy delta when 2 particles of coordinates \( (\vec{x}_1, \vec{\mu}_1) \) and
    !> \( (\vec{x}_2, \vec{\mu}_2) \) are switched.
    !> \[
    !>      \Delta U = \sum_{\vec{k}_{1:2}} w(\vec{k}_{1:2}) \Re\left[
    !>          S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2}) +
    !>          S_-^\ast(\vec{k}_{1:2}) \Delta S_+(\vec{k}_{1:2}) +
    !>          \Delta S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2})
    !>      \right]
    !> \]
    !> where:
    !> \[
    !>      \Delta S_-^\ast(\vec{k}_{1:2}) =
    !>          \left[
    !>              -k_{1:2} (\mu_{1, 3} - \mu_{2, 3}) -
    !>              i \vec{k}_{1:2} \cdot (\vec{\mu}_{1, 1:2} - \vec{\mu}_{2, 1:2})
    !>          \right] \\
    !>          \left(
    !>              e^{-k_{1:2} x_{2, 3}} e^{-i \vec{k}_{1:2} \cdot \vec{x}_{2, 1:2}} -
    !>              e^{-k_{1:2} x_{1, 3}} e^{-i \vec{k}_{1:2} \cdot \vec{x}_{1, 1:2}}
    !>          \right)
    !> \]
    !> \[
    !>      \Delta S_+(\vec{k}_{1:2}) =
    !>          \left[
    !>              +k_{1:2} (\mu_{1, 3} - \mu_{2, 3}) +
    !>              i \vec{k}_{1:2} \cdot (\vec{\mu}_{1, 1:2} - \vec{\mu}_{2, 1:2})
    !>          \right] \\
    !>          \left(
    !>              e^{+k_{1:2} x_{2, 3}} e^{+i \vec{k}_{1:2} \cdot \vec{x}_{2, 1:2}} -
    !>              e^{+k_{1:2} x_{1, 3}} e^{+i \vec{k}_{1:2} \cdot \vec{x}_{1, 1:2}}
    !>          \right)
    !> \]
    !> \[
    !>      \Delta S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2}) =
    !>          \left\{
    !>              -[k_{1:2} (\mu_{1, 3} - \mu_{2, 3})]^2 +
    !>              [\vec{k}_{1:2} \cdot (\vec{\mu}_{1, 1:2} - \vec{\mu}_{2, 1:2})]^2 - \\
    !>              2 i k_{1:2} (\mu_{1, 3} - \mu_{2, 3})
    !>              \vec{k}_{1:2} \cdot (\vec{\mu}_{1, 1:2} - \vec{\mu}_{2, 1:2})
    !>          \right\}
    !>          \left(
    !>              2 - \\
    !>              e^{+k_{1:2} x_{2, 3}} e^{+i \vec{k}_{1:2} \cdot \vec{x}_{2, 1:2}}
    !>              e^{-k_{1:2} x_{1, 3}} e^{-i \vec{k}_{1:2} \cdot \vec{x}_{1, 1:2}} - \\
    !>              e^{+k_{1:2} x_{1, 3}} e^{+i \vec{k}_{1:2} \cdot \vec{x}_{1, 1:2}}
    !>              e^{-k_{1:2} x_{2, 3}} e^{-i \vec{k}_{1:2} \cdot \vec{x}_{2, 1:2}}
    !>          \right)
    !> \]
    pure real(DP) function Abstract_visit_switch(this, ij_components, particles) &
        result(delta_energy)
        class(Abstract_DLC_Visitor), intent(in) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        real(DP) :: real_part_1, real_part_2, real_part_3

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

        delta_energy = 0._DP
        if (.not.(this%structures%is_dipolar(ij_components(1)) .or. &
            this%structures%is_dipolar(ij_components(2)))) return

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
                    particles(1)%dipolar_moment(1:2) - particles(2)%dipolar_moment(1:2))
                wave_delta_moment_3 = norm2(wave_vector) * &
                    (particles(1)%dipolar_moment(3) - particles(2)%dipolar_moment(3))

                real_part_1 = real(this%structures%get_plus(n_1, n_2) * &
                    cmplx(-wave_delta_moment_3, -wave_dot_delta_moment_12, DP) * &
                    conjg(fourier_position_2 / exp_kz_2 - fourier_position_1 / exp_kz_1), DP)
                real_part_2 = real(conjg(this%structures%get_minus(n_1, n_2)) * &
                    cmplx(+wave_delta_moment_3, +wave_dot_delta_moment_12, DP) * &
                    (fourier_position_2 * exp_kz_2 - fourier_position_1 * exp_kz_1), DP)
                real_part_3 = real(cmplx(-wave_delta_moment_3**2 + wave_dot_delta_moment_12**2, &
                    -2._DP*wave_delta_moment_3 * wave_dot_delta_moment_12, DP) * (2._DP - &
                    fourier_position_2 * exp_kz_2 * conjg(fourier_position_1) / exp_kz_1 - &
                    fourier_position_1 * exp_kz_1 * conjg(fourier_position_2) / exp_kz_2), DP)

                delta_energy = delta_energy + this%weight%get(n_1, n_2) * &
                    (real_part_1 + real_part_2 + real_part_3)
            end do
        end do
        delta_energy = 2._DP * delta_energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit_switch

!end implementation Abstract_DLC_Visitor

!implementation Null_DLC_Visitor

    subroutine Null_construct(this, periodic_box, reciprocal_lattice, weight, structures)
        class(Null_DLC_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_DLC_Weight), target, intent(in) :: weight
        class(Abstract_DLC_Structures), target, intent(in) :: structures
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_DLC_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    pure real(DP) function Null_visit(this) result(energy)
        class(Null_DLC_Visitor), intent(in) :: this
        energy = 0._DP
    end function Null_visit

    pure real(DP) function Null_visit_move(this, i_component, new_position, old) &
        result(delta_energy)
        class(Null_DLC_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: old
        delta_energy = 0._DP
    end function Null_visit_move

    pure real(DP) function Null_visit_rotation(this, i_component, new_dipolar_moment, old) &
        result(delta_energy)
        class(Null_DLC_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_dipolar_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old
        delta_energy = 0._DP
    end function Null_visit_rotation

    pure real(DP) function Null_visit_switch(this, ij_components, particles) result(delta_energy)
        class(Null_DLC_Visitor), intent(in) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        delta_energy = 0._DP
    end function Null_visit_switch

!end implementation Null_DLC_Visitor

end module class_dlc_visitor
