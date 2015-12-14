module class_ewald_reci_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions, PI
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice
use types_temporary_particle, only: Concrete_Temporary_Particle
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight
use class_ewald_reci_structure, only: Abstract_Ewald_Reci_Structure
use procedures_ewald_micro, only: set_fourier, reci_number_1_sym, reci_number_2_sym

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
    !>      U = \frac{1}{2} \sum_{\vec{k}} w_\alpha(\vec{k}) |S(\vec{k})|^2
    !> \]
    pure real(DP) function Abstract_visit(this) result(energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this

        integer :: n_1, n_2, n_3

        energy = 0._DP
        do n_3 = -this%reci_numbers(3), this%reci_numbers(3)
        do n_2 = -this%reci_numbers(2), this%reci_numbers(2)
        do n_1 = -this%reci_numbers(1), this%reci_numbers(1)
            if (n_1**2 + n_2**2 + n_3**2 > this%reci_numbers(1)**2) cycle
            energy = energy + this%weight%get(n_1, n_2, n_3) * real(this%structure%&
                get(n_1, n_2, n_3) * conjg(this%structure%get(n_1, n_2, n_3)), DP)
        end do
        end do
        end do
        energy = energy / 2._DP
    end function Abstract_visit

    !> Energy delta when a particle \( \mathsf{i} \) of component \( \mathsf{I} \) moves.
    !> \[
    !>      \Delta U_{\mathsf{I}, \mathsf{J}} = \sum_{\vec{k}} w_\alpha(\vec{k})
    !>          (\vec{k}\cdot\vec{\mu}_\mathsf{i}) \{
    !>          \Re[(e^{i\vec{k}\cdot\vec{x}^\prime_\mathsf{i}} -
    !>              e^{i\vec{k}\cdot\vec{x}_\mathsf{i}}) S_\mathsf{J}^\ast(\vec{k})] +
    !>          (\vec{k}\cdot\vec{\mu}_\mathsf{i})
    !>              [1 - \Re(e^{i\vec{k}\cdot\vec{x}^\prime_\mathsf{i}}
    !>                  e^{-i\vec{k}\cdot\vec{x}_\mathsf{i}})]
    !>          \}
    !> \]
    pure real(DP) function Abstract_visit_move(this, i_component, new, old) result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

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
        delta_energy = 2._DP * delta_energy ! half wave vector (symmetry) -> double energy
    end function Abstract_visit_move

    !> Energy delta when a particle \( \mathsf{i} \) of component \( \mathsf{\mathsf{I}} \) rotates.
    !> \[
    !>      \Delta U_{\mathsf{\mathsf{I}}, \mathsf{J}} = \sum_{\vec{k}} w_\alpha(\vec{k})
    !>          \vec{k}\cdot(\vec{\mu}^\prime_\mathsf{i} - \vec{\mu}_\mathsf{i}) [
    !>              \Re(e^{i\vec{k}\cdot\vec{x}_\mathsf{i}} S_\mathsf{J}^\ast(\vec{k})) +
    !>              \vec{k}\cdot(\vec{\mu}^\prime_\mathsf{i} - \vec{\mu}_\mathsf{i})
    !>      ]
    !> \]
    pure real(DP) function Abstract_visit_rotation(this, i_component, new, old) result(delta_energy)
        class(Abstract_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old

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
                        new%dipolar_moment - old%dipolar_moment)

                    real_fourier_x_conjg_structure = real(fourier_position * &
                        conjg(this%structure%get(n_1, n_2, n_3)), DP)

                    delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                        wave_dot_delta_moment * (real_fourier_x_conjg_structure + &
                            0.5_DP * wave_dot_delta_moment)
                end do
            end do
        end do
        delta_energy = 2._DP * delta_energy ! half wave vector (symmetry) -> double energy
    end function Abstract_visit_rotation

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

    pure real(DP) function Null_visit_move(this, i_component, new, old) result(delta_energy)
        class(Null_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        delta_energy = 0._DP
    end function Null_visit_move

    pure real(DP) function Null_visit_rotation(this, i_component, new, old) result(delta_energy)
        class(Null_Ewald_Reci_Visitor), intent(in) :: this
        integer, intent(in) :: i_component
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        delta_energy = 0._DP
    end function Null_visit_rotation

!end implementation Null_Ewald_Reci_Visitor

end module class_ewald_reci_visitor
