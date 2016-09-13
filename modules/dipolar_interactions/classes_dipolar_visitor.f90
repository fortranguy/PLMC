module classes_dipolar_visitor

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_component_wrapper, only: Component_Wrapper
use types_dipolar_interactions_wrapper, only: Dipolar_Interactions_Wrapper
use types_reals_line, only: Reals_Line
use procedures_plmc_visit, only: visit_dipolar

implicit none

private

    type, abstract, public :: Abstract_Dipolar_Visitor
    contains
        procedure(Abstract_destroy), deferred :: destroy
        procedure(Abstract_visit), deferred :: visit
    end type Abstract_Dipolar_Visitor

    abstract interface

        subroutine Abstract_destroy(this)
        import :: Abstract_Dipolar_Visitor
            class(Abstract_Dipolar_Visitor), intent(inout) :: this
        end subroutine Abstract_destroy

        subroutine Abstract_visit(this, new_energies, new_mixture_energy, box_volume_ratio, &
            energies, mixture_energy)
        import :: DP, Reals_Line, Abstract_Dipolar_Visitor
            class(Abstract_Dipolar_Visitor), intent(in) :: this
            type(Reals_Line), intent(inout) :: new_energies(:)
            real(DP), intent(out) :: new_mixture_energy
            real(DP), intent(in) :: box_volume_ratio !! \( \frac{V}{V_0} \)
            type(Reals_Line), intent(in) :: energies(:)
            real(DP), intent(in) :: mixture_energy
        end subroutine Abstract_visit

    end interface

    type, extends(Abstract_Dipolar_Visitor), public :: Scalable_Dipolar_Visitor
    contains
        procedure :: destroy => Scalable_destroy
        procedure :: visit => Scalable_visit
    end type Scalable_Dipolar_Visitor

    type, extends(Abstract_Dipolar_Visitor), public :: Unscalable_Dipolar_Visitor
    private
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Dipolar_Interactions_Wrapper), pointer :: dipolar_interactions => null()
    contains
        procedure :: construct => Unscalable_construct
        procedure :: destroy => Unscalable_destroy
        procedure :: visit => Unscalable_visit
    end type Unscalable_Dipolar_Visitor

    type, extends(Abstract_Dipolar_Visitor), public :: Null_Dipolar_Visitor
    contains
        procedure :: destroy => Null_destroy
        procedure :: visit => Null_visit
    end type Null_Dipolar_Visitor

contains

!implementation Scalable_Dipolar_Visitor

    subroutine Scalable_destroy(this)
        class(Scalable_Dipolar_Visitor), intent(inout) :: this
    end subroutine Scalable_destroy

    subroutine Scalable_visit(this, new_energies, new_mixture_energy, box_volume_ratio, &
        energies, mixture_energy)
        class(Scalable_Dipolar_Visitor), intent(in) :: this
        type(Reals_Line), intent(inout) :: new_energies(:)
        real(DP), intent(out) :: new_mixture_energy
        real(DP), intent(in) :: box_volume_ratio
        type(Reals_Line), intent(in) :: energies(:)
        real(DP), intent(in) :: mixture_energy

        integer :: i_component

        do i_component = 1, size(new_energies)
            new_energies(i_component)%line = energies(i_component)%line / box_volume_ratio
        end do
        new_mixture_energy = mixture_energy / box_volume_ratio
    end subroutine Scalable_visit

!implementation Scalable_Dipolar_Visitor

!implementation Unscalable_Dipolar_Visitor

    subroutine Unscalable_construct(this, components, dipolar_interactions)
        class(Unscalable_Dipolar_Visitor), intent(out) :: this
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Dipolar_Interactions_Wrapper), target, intent(in) :: dipolar_interactions

        this%components => components
        this%dipolar_interactions => dipolar_interactions
    end subroutine Unscalable_construct

    subroutine Unscalable_visit(this, new_energies, new_mixture_energy, box_volume_ratio, &
        energies, mixture_energy)
        class(Unscalable_Dipolar_Visitor), intent(in) :: this
        type(Reals_Line), intent(inout) :: new_energies(:)
        real(DP), intent(out) :: new_mixture_energy
        real(DP), intent(in) :: box_volume_ratio
        type(Reals_Line), intent(in) :: energies(:)
        real(DP), intent(in) :: mixture_energy

        call visit_dipolar(new_energies, new_mixture_energy, this%components, this%&
            dipolar_interactions)
    end subroutine Unscalable_visit

    subroutine Unscalable_destroy(this)
        class(Unscalable_Dipolar_Visitor), intent(inout) :: this

        this%dipolar_interactions => null()
        this%components => null()
    end subroutine Unscalable_destroy

!end implementation Unscalable_Dipolar_Visitor

!implementation Null_Dipolar_Visitor

    subroutine Null_destroy(this)
        class(Null_Dipolar_Visitor), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_visit(this, new_energies, new_mixture_energy, box_volume_ratio, &
        energies, mixture_energy)
        class(Null_Dipolar_Visitor), intent(in) :: this
        type(Reals_Line), intent(inout) :: new_energies(:)
        real(DP), intent(out) :: new_mixture_energy
        real(DP), intent(in) :: box_volume_ratio
        type(Reals_Line), intent(in) :: energies(:)
        real(DP), intent(in) :: mixture_energy

        integer :: i_component

        do i_component = 1, size(new_energies)
            new_energies(i_component)%line = 0._DP
        end do
        new_mixture_energy = 0._DP
    end subroutine Null_visit

!end implementation Null_Dipolar_Visitor

end module classes_dipolar_visitor
