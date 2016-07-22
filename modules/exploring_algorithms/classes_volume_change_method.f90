module classes_volume_change_method

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use classes_number_to_string, only: Concrete_Number_to_String
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_temporary_particle, only: Concrete_Temporary_Particle
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use procedures_visit_condition, only: visit_condition_lower => lower
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper

implicit none

private

    type, abstract, public :: Abstract_Volume_Change_Method
    private
        type(Environment_Wrapper), pointer :: environment => null()
        type(Component_Wrapper), pointer :: components(:) => null()
        type(Short_Interactions_Wrapper), pointer :: short_interactions => null()
    contains
        procedure :: construct => Abstract_construct
        procedure :: destroy => Abstract_destroy
        procedure :: try => Abstract_try
    end type Abstract_Volume_Change_Method

    type, extends(Abstract_Volume_Change_Method), public :: Concrete_Volume_Change_Method

    end type Concrete_Volume_Change_Method

    type, extends(Abstract_Volume_Change_Method), public :: Null_Volume_Change_Method
    contains
        procedure :: construct => Null_construct
        procedure :: destroy => Null_destroy
        procedure :: try => Null_try
    end type Null_Volume_Change_Method

contains

!implementation Abstract_Volume_Change_Method

    subroutine Abstract_construct(this, environment, components, short_interactions)
        class(Abstract_Volume_Change_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions

        this%environment => environment
        this%components => components
        this%short_interactions => short_interactions
    end subroutine Abstract_construct

    subroutine Abstract_destroy(this)
        class(Abstract_Volume_Change_Method), intent(inout) :: this

        this%short_interactions => null()
        this%components => null()
        this%environment => null()
    end subroutine Abstract_destroy

    subroutine Abstract_try(this, observables)
        class(Abstract_Volume_Change_Method), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables

        real(DP) :: contacts, conctacts_j
        integer :: i_component, j_component, num_particles_max, i_particle, i_exclude
        logical :: overlap
        type(Concrete_Temporary_Particle) :: particle
        type(Concrete_Number_to_String) :: string

        contacts = 0._DP
        do j_component = 1, size(this%components)
            do i_component = 1, j_component
                num_particles_max = max(this%components(j_component)%positions%get_num(), this%&
                    components(i_component)%positions%get_num())
                do i_particle = 1, this%components(j_component)%positions%get_num()
                    particle%i = i_particle
                    particle%position = this%components(j_component)%positions%get(particle%i)
                    i_exclude = merge(particle%i, num_particles_max + 1, i_component == j_component)
                    call this%short_interactions%visitable_cells(i_component, j_component)%&
                        visit_contacts(overlap, conctacts_j, particle, visit_condition_lower, &
                            i_exclude)
                    if (overlap) then
                        call error_exit("Abstract_Volume_Change_Method: try: components "//string%&
                            get(i_component)//" and "//string%get(j_component)//" overlap.")
                    end if
                    contacts = contacts + conctacts_j
                end do
            end do
        end do
        observables%beta_pressure_excess = this%short_interactions%hard_contact%&
            get_beta_pressure_excess(contacts)
    end subroutine Abstract_try

!end implementation Abstract_Volume_Change_Method

!implementation Null_Volume_Change_Method

    subroutine Null_construct(this, environment, components, short_interactions)
        class(Null_Volume_Change_Method), intent(out) :: this
        type(Environment_Wrapper), target, intent(in) :: environment
        type(Component_Wrapper), target, intent(in) :: components(:)
        type(Short_Interactions_Wrapper), target, intent(in) :: short_interactions
    end subroutine Null_construct

    subroutine Null_destroy(this)
        class(Null_Volume_Change_Method), intent(inout) :: this
    end subroutine Null_destroy

    subroutine Null_try(this, observables)
        class(Null_Volume_Change_Method), intent(in) :: this
        type(Exploring_Observables_Wrapper), intent(inout) :: observables
    end subroutine Null_try

!end implementation Null_Volume_Change_Method

end module classes_volume_change_method
