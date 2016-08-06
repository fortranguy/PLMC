module classes_box_volume_change

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use types_environment_wrapper, only: Environment_Wrapper
use types_mixture_wrapper, only: Mixture_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use classes_changed_box_size, only: Abstract_Changed_Box_Size
use types_temporary_observables, only: Concrete_Single_Delta_Energies
use types_generating_observables_wrapper, only: Generating_Observables_Wrapper
use classes_metropolis_algorithm, only: Abstract_Metropolis_Algorithm

implicit none

private

    type, extends(Abstract_Metropolis_Algorithm), abstract, public :: Abstract_Box_Volume_Change
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Mixture_Wrapper), pointer :: mixture => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
        class(Abstract_Changed_Box_Size), pointer :: changed_box_size => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
    end type Abstract_Box_Volume_Change

contains

!implementation Abstract_Metropolis_Algorithm

    subroutine Abstract_construct(this, environment, mixture, short_interactions, changed_box_size)
        class(Abstract_Box_Volume_Change), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Mixture_Wrapper), target, intent(in) :: mixture
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
        class(Abstract_Changed_Box_Size), target, intent(in) :: changed_box_size

        this%environment => environment
        this%mixture => mixture
        this%short_interactions => short_interactions
        this%changed_box_size => changed_box_size
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Box_Volume_Change), intent(inout) :: this

        this%changed_box_size => null()
        this%short_interactions => null()
        this%mixture => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, observables)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        type(Generating_Observables_Wrapper), intent(inout) :: observables
    end subroutine Abstract_try

    subroutine Abstract_test_metropolis(this, success, deltas)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: success
        type(Concrete_Single_Delta_Energies), intent(inout) :: deltas

        real(DP) :: box_size_ratio(num_dimensions)
        integer :: i_component

        success = .false.
        box_size_ratio = this%changed_box_size%get_ratio()
        do i_component = 1, size(this%mixture%components)
            call this%mixture%components(i_component)%positions%rescale_all(box_size_ratio)
        end do
    end subroutine Abstract_test_metropolis

    subroutine Abstrac_visit_walls(this, overlap, new_energy)
        class(Abstract_Box_Volume_Change), intent(in) :: this
        logical, intent(out) :: overlap
        real(DP), intent(out) :: new_energy

    end subroutine Abstrac_visit_walls

!end implementation Abstract_Metropolis_Algorithm

end module classes_box_volume_change
