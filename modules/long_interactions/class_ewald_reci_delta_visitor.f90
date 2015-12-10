module class_ewald_reci_delta_visitor

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

    type, abstract, public :: Abstract_Ewald_Reci_Delta_Visitor
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        class(Abstract_Ewald_Reci_Weight), pointer :: weight => null()
        class(Abstract_Ewald_Reci_Structure), pointer :: structure => null()
    contains
        procedure :: construct => Abstract_Ewald_Reci_Delta_Visitor_construct
        procedure :: destroy => Abstract_Ewald_Reci_Delta_Visitor_destroy
        procedure :: visit_move => Abstract_Ewald_Reci_Delta_Visitor_visit_move
        procedure :: visit_rotation => Abstract_Ewald_Reci_Delta_Visitor_visit_rotation
    end type Abstract_Ewald_Reci_Delta_Visitor

    type, extends(Abstract_Ewald_Reci_Delta_Visitor), public :: Concrete_Ewald_Reci_Delta_Visitor

    end type Concrete_Ewald_Reci_Delta_Visitor

    type, extends(Abstract_Ewald_Reci_Delta_Visitor), public :: Null_Ewald_Reci_Delta_Visitor
    contains
        procedure :: construct => Null_Ewald_Reci_Delta_Visitor_construct
        procedure :: destroy => Null_Ewald_Reci_Delta_Visitor_destroy
        procedure :: visit_move => Null_Ewald_Reci_Delta_Visitor_visit_move
        procedure :: visit_rotation => Null_Ewald_Reci_Delta_Visitor_visit_rotation
    end type Null_Ewald_Reci_Delta_Visitor

contains

!implementation Abstract_Ewald_Reci_Delta_Visitor

    subroutine Abstract_Ewald_Reci_Delta_Visitor_construct(this, periodic_box, reciprocal_lattice, &
        weight, structure)
        class(Abstract_Ewald_Reci_Delta_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Ewald_Reci_Weight), target, intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), target, intent(in) :: structure

        this%periodic_box => periodic_box
        this%reci_numbers = reciprocal_lattice%get_numbers()
        this%weight => weight
        this%structure => structure
    end subroutine Abstract_Ewald_Reci_Delta_Visitor_construct

    subroutine Abstract_Ewald_Reci_Delta_Visitor_destroy(this)
        class(Abstract_Ewald_Reci_Delta_Visitor), intent(inout) :: this

        this%structure => null()
        this%weight => null()
        this%periodic_box => null()
    end subroutine Abstract_Ewald_Reci_Delta_Visitor_destroy

    !> Energy delta when a particle \( \mathsf{i} \) of component \( \mathsf{I} \) moves.
    !> \[
    !>      \Delta U_{\mathsf{I}, \mathsf{J}} = \sum_{\vec{k}} w_\alpha(\vec{k})
    !>          (\vec{k}\cdot\vec{\mu}_\mathsf{i}) \{
    !>          \Re[(e^{i\vec{k}\cdot\vec{x}^\prime_\mathsf{i}} -
    !>              e^{i\vec{k}\cdot\vec{x}_\mathsf{i}}) S_\mathsf{J}^\ast(\vec{k})] +
    !>          [\mathsf{I}=\mathsf{J}] (\vec{k}\cdot\vec{\mu}_\mathsf{i})
    !>              [1 - \Re(e^{i\vec{k}\cdot\vec{x}^\prime_\mathsf{i}}
    !>                  e^{-i\vec{k}\cdot\vec{x}_\mathsf{i}})]
    !>          \}
    !> \]
    !> with \( s_\mathsf{i}(\vec{k}) = (\vec{k}\cdot\vec{\mu}_\mathsf{i})
    !>      e^{i\vec{k}\cdot\vec{x}_\mathsf{i}} \).
    pure real(DP) function Abstract_Ewald_Reci_Delta_Visitor_visit_move(this, new, old, &
        same_component) result(delta_energy)
        class(Abstract_Ewald_Reci_Delta_Visitor), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        logical, intent(in) :: same_component

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
                    wave_dot_moment = dot_product(wave_vector, old%dipolar_moment)

                    real_delta_fourier_x_conjg_structure = &
                        real((fourier_position_new - fourier_position_old) * &
                            conjg(this%structure%get(n_1, n_2, n_3)), DP)

                    if (same_component) then
                        delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                            wave_dot_moment * (real_delta_fourier_x_conjg_structure + &
                            wave_dot_moment * (1._DP - real(fourier_position_new * &
                                conjg(fourier_position_old), DP)))
                    else
                        delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                            wave_dot_moment * real_delta_fourier_x_conjg_structure
                    end if
                end do
            end do
        end do
        delta_energy = 2._DP * delta_energy ! half wave vector (symmetry) -> double energy
    end function Abstract_Ewald_Reci_Delta_Visitor_visit_move

    !> Energy delta when a particle \( \mathsf{i} \) of component \( \mathsf{\mathsf{I}} \) rotates.
    !> \[
    !>      \Delta U_{\mathsf{\mathsf{I}}, \mathsf{J}} = \sum_{\vec{k}} w_\alpha(\vec{k})
    !>          \vec{k}\cdot(\vec{\mu}^\prime_\mathsf{i} - \vec{\mu}_\mathsf{i}) [
    !>              \Re(e^{i\vec{k}\cdot\vec{x}_\mathsf{i}} S_\mathsf{J}^\ast(\vec{k})) +
    !>              [\mathsf{\mathsf{I}}=\mathsf{J}] \vec{k}\cdot(\vec{\mu}^\prime_\mathsf{i} -
    !>                  \vec{\mu}_\mathsf{i})
    !>      ]
    !> \]
    !> with \( s_\mathsf{i}(\vec{k}) = (\vec{k}\cdot\vec{\mu}_\mathsf{i})
    !>      e^{i\vec{k}\cdot\vec{x}_\mathsf{i}} \).
    pure real(DP) function Abstract_Ewald_Reci_Delta_Visitor_visit_rotation(this, new, old, &
        same_component) result(delta_energy)
        class(Abstract_Ewald_Reci_Delta_Visitor), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        logical, intent(in) :: same_component

        real(DP) :: box_size(num_dimensions)
        real(DP), dimension(num_dimensions) :: wave_vector
        integer :: n_1, n_2, n_3
        real(DP), dimension(num_dimensions) :: wave_1_x_position
        real(DP) :: real_fourier_x_conjg_structure, wave_dot_delta_moment

        complex(DP) :: fourier_position
        complex(DP), dimension(-this%reci_numbers(1):this%reci_numbers(1)) :: fourier_position_1
        complex(DP), dimension(-this%reci_numbers(2):this%reci_numbers(2)) :: fourier_position_2
        complex(DP), dimension(-this%reci_numbers(3):this%reci_numbers(3)) :: fourier_position_3

        box_size = this%periodic_box%get_size()

        wave_1_x_position = 2._DP*PI * old%position / box_size
        call set_fourier(fourier_position_1, this%reci_numbers(1), wave_1_x_position(1))
        call set_fourier(fourier_position_2, this%reci_numbers(2), wave_1_x_position(2))
        call set_fourier(fourier_position_3, this%reci_numbers(3), wave_1_x_position(3))

        delta_energy = 0._DP
        do n_3 = 0, this%reci_numbers(3)
            wave_vector(3) = 2._DP*PI * real(n_3, DP) / box_size(3)
            do n_2 = -reci_number_2_sym(this%reci_numbers, n_3), this%reci_numbers(2)
                wave_vector(2) = 2._DP*PI * real(n_2, DP) / box_size(2)
                do n_1 = -reci_number_1_sym(this%reci_numbers, n_3, n_2), this%reci_numbers(1)
                    wave_vector(1) = 2._DP*PI * real(n_1, DP) / box_size(1)

                    fourier_position = fourier_position_1(n_1) * fourier_position_2(n_2) * &
                        fourier_position_3(n_3)
                    wave_dot_delta_moment = dot_product(wave_vector, &
                        new%dipolar_moment - old%dipolar_moment)

                    real_fourier_x_conjg_structure = real(fourier_position * &
                        conjg(this%structure%get(n_1, n_2, n_3)), DP)

                    if (same_component) then
                        delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                            wave_dot_delta_moment * (real_fourier_x_conjg_structure + &
                                0.5_DP * wave_dot_delta_moment)
                    else
                        delta_energy = delta_energy + this%weight%get(n_1, n_2, n_3) * &
                            wave_dot_delta_moment * real_fourier_x_conjg_structure
                    end if
                end do
            end do
        end do
        delta_energy = 2._DP * delta_energy ! half wave vector (symmetry) -> double energy
    end function Abstract_Ewald_Reci_Delta_Visitor_visit_rotation

!end implementation Abstract_Ewald_Reci_Delta_Visitor

!implementation Null_Ewald_Reci_Delta_Visitor

    subroutine Null_Ewald_Reci_Delta_Visitor_construct(this, periodic_box, reciprocal_lattice, &
        weight, structure)
        class(Null_Ewald_Reci_Delta_Visitor), intent(out) :: this
        class(Abstract_Periodic_Box), target, intent(in) :: periodic_box
        class(Abstract_Reciprocal_Lattice), intent(in) :: reciprocal_lattice
        class(Abstract_Ewald_Reci_Weight), target, intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), target, intent(in) :: structure
    end subroutine Null_Ewald_Reci_Delta_Visitor_construct

    subroutine Null_Ewald_Reci_Delta_Visitor_destroy(this)
        class(Null_Ewald_Reci_Delta_Visitor), intent(inout) :: this
    end subroutine Null_Ewald_Reci_Delta_Visitor_destroy

    pure real(DP) function Null_Ewald_Reci_Delta_Visitor_visit_move(this, new, old, &
        same_component) result(delta_energy)
        class(Null_Ewald_Reci_Delta_Visitor), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        logical, intent(in) :: same_component
        delta_energy = 0._DP
    end function Null_Ewald_Reci_Delta_Visitor_visit_move

    pure real(DP) function Null_Ewald_Reci_Delta_Visitor_visit_rotation(this, new, old, &
        same_component) result(delta_energy)
        class(Null_Ewald_Reci_Delta_Visitor), intent(in) :: this
        type(Concrete_Temporary_Particle), intent(in) :: new, old
        logical, intent(in) :: same_component
        delta_energy = 0._DP
    end function Null_Ewald_Reci_Delta_Visitor_visit_rotation

!end implementation Null_Ewald_Reci_Delta_Visitor

end module class_ewald_reci_delta_visitor
