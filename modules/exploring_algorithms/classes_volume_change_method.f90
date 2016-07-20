module classes_volume_change_method

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_environment_wrapper, only: Environment_Wrapper
use types_component_wrapper, only: Component_Wrapper
use types_short_interactions_wrapper, only: Short_Interactions_Wrapper

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
    end type Abstract_Volume_Change_Method

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

!end implementation Abstract_Volume_Change_Method

end module classes_volume_change_method
