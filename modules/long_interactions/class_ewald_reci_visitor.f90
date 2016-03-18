module class_ewald_reci_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight
use class_ewald_reci_structure, only: Abstract_Ewald_Reci_Structure
use procedures_long_interactions_micro, only: set_fourier, reci_number_1_sym, reci_number_2_sym

implicit none

private

    type, abstract, public :: Abstract_Ewald_Reci_Visitor
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        class(Abstract_Ewald_Reci_Weight), pointer :: weight => null()
        class(Abstract_Ewald_Reci_Structure), pointer :: structure => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: visit => Abstract_visit
        procedure :: visit_move => Abstract_visit_move
        procedure :: visit_rotation => Abstract_visit_rotation
        procedure :: visit_add => Abstract_visit_add
        procedure :: visit_remove => Abstract_visit_remove
        procedure :: visit_switch => Abstract_visit_switch
        procedure, private :: visit_exchange => Abstract_visit_exchange
    end type Abstract_Ewald_Reci_Visitor

    type, extends(Abstract_Ewald_Reci_Visitor), public :: Concrete_Ewald_Reci_Visitor

    end type Concrete_Ewald_Reci_Visitor

    type, extends(Abstract_Ewald_Reci_Visitor), public :: Null_Ewald_Reci_Visitor
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
        procedure :: visit_move => Null_visit_move
        procedure :: visit_rotation => Null_visit_rotation
        procedure :: visit_switch => Null_visit_switch
        procedure, private :: visit_exchange => Null_visit_exchange
    end type Null_Ewald_Reci_Visitor

contains

!implementation Abstract_Ewald_Reci_Visitor

    subroutine Abstract_construct(this, periodic_box, reciprocal_lattice, weight, structure)
        class(Abstract_Ewald_Reci_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Ewald_Reci_Weight), target, intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), target, intent(in) :: structure

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%weight => weight
        this%structure => structure
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Ewald_Reci_Visitor), intent(inout) :: this

        this%structure => null()
        this%weight => null()
        this%periodic_box => null()
    end subroutine Abstract_destroy

    !> \[
    !>      U = \sum_{\vec{k}} w_\alpha(\vec{k}) |S(\vec{k})|^2
    !> \]
    !> where \( w_\alpha(\vec{k}) \) is [[Abstract_Ewald_Reci_Weight]]
    !> and \( S(\vec{k}) \) is [[Abstract_Ewald_Reci_Structure]].
    pure real(DP) function Abstract_visit(this) result(energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this

        integer :: n_1, n_2, n_3

        energy = 0._DP
        do n_3 = 0, this%reci_numbers(3)
            do n_2 = -reci_number_2_sym(this%reci_numbers, n_3), this%reci_numbers(2)
                do n_1 = -reci_number_1_sym(this%reci_numbers, n_3, n_2), this%reci_numbers(1)
                    if (n_1**2 + n_2**2 + n_3**2 > this%reci_numbers(1)**2) cycle
                    energy = energy + this%weight%get(n_1, n_2, n_3) * real(this%structure%&
                        get(n_1, n_2, n_3) * conjg(this%structure%get(n_1, n_2, n_3)), DP)
                end do
            end do
        end do
        energy = 2._DP * energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit

    !> Energy delta when a particle of coordinates \( (\vec{x}, \vec{\mu}) \) moves.
    !> \[
    !>      \Delta U = 2 \sum_{\vec{k}} w_\alpha(\vec{k}) (\vec{k}\cdot\vec{\mu}) \left\{
    !>              \Re\left[
    !>                  \left( e^{i\vec{k}\cdot\vec{x}^\prime} - e^{i\vec{k}\cdot\vec{x}} \right)
    !>                  S^\ast(\vec{k})
    !>              \right] + (\vec{k}\cdot\vec{\mu})\left[
    !>                  1 - \Re\left(
    !>                      e^{i\vec{k}\cdot\vec{x}^\prime} e^{-i\vec{k}\cdot\vec{x}}
    !>                  \right)
    !>              \right]
    !>          \right\}
    !> \]
    pure real(DP) function Abstract_visit_move(this, i_component, new_position, old) &
        result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: old

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position_new, wave_1_x_position_old
        real(DP) :: real_delta_fourier_x_conjg_structure, wave_dot_moment

        complex(DP) :: fourier_position_new
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_new_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_new_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_new_3

        complex(DP) :: fourier_position_old
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_old_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_old_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_old_3

        delta_energy = 0._DP
        if (.not.this%structure%is_dipolar(i_component)) return

        box_size = this%periodic_box%get_size()
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
                    wave_dot_moment = dot_product(wave_vector, old%dipolar_moment)

                    real_delta_fourier_x_conjg_structure = &
                        real((fourier_position_new - fourier_position_old) * &
                            conjg(this%structure%get(n_1, n_2, n_3)), DP)

                    delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                        wave_dot_moment * (real_delta_fourier_x_conjg_structure + &
                        wave_dot_moment * (1._DP - real(fourier_position_new * &
                            conjg(fourier_position_old), DP)))
                end do
            end do
        end do
        delta_energy = 4._DP * delta_energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit_move

    !> Energy delta when a particle of coordinates \( (\vec{x}, \vec{\mu}) \) rotates.
    !> \[
    !>      \Delta U = \sum_{\vec{k}} w_\alpha(\vec{k}) \vec{k}\cdot(\vec{\mu}^\prime - \vec{\mu})
    !>          \left[
    !>              2\Re\left( e^{i\vec{k}\cdot\vec{x}} S^\ast(\vec{k}) \right) +
    !>              \vec{k}\cdot(\vec{\mu}^\prime - \vec{\mu})
    !>          \right]
    !> \]
    pure real(DP) function Abstract_visit_rotation(this, i_component, new_dipolar_moment, old) &
        result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_dipolar_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position
        real(DP) :: real_fourier_x_conjg_structure, wave_dot_delta_moment

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        delta_energy = 0._DP
        if (.not.this%structure%is_dipolar(i_component)) return

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
                    wave_dot_delta_moment = dot_product(wave_vector, &
                        new_dipolar_moment - old%dipolar_moment)

                    real_fourier_x_conjg_structure = real(fourier_position * &
                        conjg(this%structure%get(n_1, n_2, n_3)), DP)

                    delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                        wave_dot_delta_moment * (2._DP*real_fourier_x_conjg_structure + &
                            wave_dot_delta_moment)
                end do
            end do
        end do
        delta_energy = 2._DP * delta_energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit_rotation

    pure real(DP) function Abstract_visit_add(this, i_component, particle) result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        delta_energy = this%visit_exchange(i_component, particle, 1._DP)
    end function Abstract_visit_add

    pure real(DP) function Abstract_visit_remove(this, i_component, particle) result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle

        delta_energy = this%visit_exchange(i_component, particle, -1._DP)
    end function Abstract_visit_remove

    !> Energy delta when a particle of coordinates \( (\vec{x}, \vec{\mu}) \) is added (\( + \)) or
    !> removed (\( - \)).
    !> \[
    !>      \Delta U = \sum_{\vec{k}} w_\alpha(\vec{k}) (\vec{k} \cdot \vec{\mu}) \left[
    !>              2 \Re\left( e^{i\vec{k}\cdot\vec{x}} S^\ast(\vec{k}) \right) \pm
    !>              (\vec{k} \cdot \vec{\mu})
    !>          \right]
    !> \]
    pure real(DP) function Abstract_visit_exchange(this, i_component, particle, signed) &
        result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position
        real(DP) :: real_fourier_x_conjg_structure, wave_dot_moment

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        delta_energy = 0._DP
        if (.not.this%structure%is_dipolar(i_component)) return

        box_size = this%periodic_box%get_size()
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
                    wave_dot_moment = dot_product(wave_vector, particle%dipolar_moment)

                    real_fourier_x_conjg_structure = real(fourier_position * &
                        conjg(this%structure%get(n_1, n_2, n_3)), DP)

                    delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                        wave_dot_moment * (2._DP*real_fourier_x_conjg_structure + &
                        signed*wave_dot_moment)
                end do
            end do
        end do
        delta_energy = 2._DP * delta_energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit_exchange

    !> Energy delta when 2 particles of coordinates \( (\vec{x}_1, \vec{\mu}_1) \) and
    !> \( (\vec{x}_2, \vec{\mu}_2) \) are switched.
    !> \[
    !>      \Delta U = 2 \sum_{\vec{k}} w_\alpha(\vec{k}) \vec{k}\cdot(\vec{\mu}_1 - \vec{\mu}_2)
    !>      \left\{
    !>              \Re\left[
    !>                  \left(e^{i\vec{k}\cdot\vec{x}_2} - e^{i\vec{k}\cdot\vec{x}_1} \right)
    !>                      S^\ast(\vec{k})
    !>              \right] +
    !>              \vec{k}\cdot(\vec{\mu}_1 - \vec{\mu}_2)
    !>              \left[
    !>                  1 - \Re\left( e^{i\vec{k}\cdot\vec{x}_1} e^{-i\vec{k}\cdot\vec{x}_2} \right)
    !>              \right]
    !>      \right\}
    !> \]
    pure real(DP) function Abstract_visit_switch(this, ij_components, particles) &
        result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position_1, wave_1_x_position_2
        real(DP) :: real_delta_fourier_x_conjg_structure, wave_dot_delta_moment

        complex(DP) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_1_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_1_3
        complex(DP) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_2_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_2_3

        delta_energy = 0._DP
        if (.not.(this%structure%is_dipolar(ij_components(1)) .or. &
            this%structure%is_dipolar(ij_components(2)))) return

        box_size = this%periodic_box%get_size()
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
                    wave_dot_delta_moment = dot_product(wave_vector, particles(1)%dipolar_moment - &
                        particles(2)%dipolar_moment)

                    real_delta_fourier_x_conjg_structure = &
                        real((fourier_position_2 - fourier_position_1) * &
                            conjg(this%structure%get(n_1, n_2, n_3)), DP)

                    delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                        wave_dot_delta_moment * (real_delta_fourier_x_conjg_structure + &
                            wave_dot_delta_moment * (1._DP - real(fourier_position_1 * &
                                conjg(fourier_position_2), DP)))
                end do
            end do
        end do
        delta_energy = 4._DP * delta_energy ! symmetry: half wave vectors -> double energy
    end function Abstract_visit_switch
!end implementation Abstract_Ewald_Reci_Visitor

!implementation Null_Ewald_Reci_Visitor

    subroutine Null_construct(this, periodic_box, reciprocal_lattice, weight, structure)
        class(Null_Ewald_Reci_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Ewald_Reci_Weight), target, intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), target, intent(in) :: structure
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Ewald_Reci_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    pure real(DP) function Null_visit(this) result(energy)
        class(Null_Ewald_Reci_Visitor), intent(in) :: this
        energy = 0._DP
    end function Null_visit

    pure real(DP) function Null_visit_move(this, i_component, new_position, old) &
        result(delta_energy)
        class(Null_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_position(:)
        type(Concrete_Temporary_Particle), intent(in) :: old
        delta_energy = 0._DP
    end function Null_visit_move

    pure real(DP) function Null_visit_rotation(this, i_component, new_dipolar_moment, old) &
        result(delta_energy)
        class(Null_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        real(DP), intent(in) :: new_dipolar_moment(:)
        type(Concrete_Temporary_Particle), intent(in) :: old
        delta_energy = 0._DP
    end function Null_visit_rotation

    pure real(DP) function Null_visit_exchange(this, i_component, particle, signed) &
        result(delta_energy)
        class(Null_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: particle
        real(DP), intent(in) :: signed
        delta_energy = 0._DP
    end function Null_visit_exchange

    pure real(DP) function Null_visit_switch(this, ij_components, particles) result(delta_energy)
        class(Null_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: ij_components(:)
        type(Concrete_Temporary_Particle), intent(in) :: particles(:)
        delta_energy = 0._DP
    end function Null_visit_switch

!end implementation Null_Ewald_Reci_Visitor

end module class_ewald_reci_visitor
