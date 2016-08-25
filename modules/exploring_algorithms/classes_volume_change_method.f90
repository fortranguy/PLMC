module classes_volume_change_method

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_errors, only: error_exit
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper
use types_exploring_observables_wrapper, only: Exploring_Observables_Wrapper
use classes_exploring_algorithm, only: Abstract_Exploring_Algorithm
use procedures_plmc_visit, only: visit_short

implicit none

private

    type, extends(Abstract_Exploring_Algorithm), abstract, public :: Abstract_Volume_Change_Method
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
        type(Exploring_Observables_Wrapper), intent(inout) :: observables !too much?

        logical :: overlap
        real(DP) :: contacts

        call visit_short(overlap, contacts, this%components, this%short_interactions)
        if (overlap) call error_exit("Abstract_Volume_Change_Method: try: visit_short: overlap")
        observables%beta_pressure_excess = this%short_interactions%beta_pressure_excess%&
            get(contacts)
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
