module procedures_ewald_self_factory

use class_permittivity, only: Abstract_Permittivity
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use types_component_wrapper, only: Component_Wrapper
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_self_component, only: Abstract_Ewald_Self_Component, &
    Concrete_Ewald_Self_Component, Null_Ewald_Self_Component
use types_long_interactions_wrapper, only: Ewald_Self_Component_Wrapper

implicit none

private
public :: ewald_self_create, ewald_self_destroy

interface ewald_self_create
    module procedure :: create_self_components
    module procedure :: create_self
end interface ewald_self_create

interface ewald_self_destroy
    module procedure :: destroy_self
    module procedure :: destroy_self_components
end interface ewald_self_destroy

contains

    subroutine create_self_components(self_components, permittivity, components, are_dipolar, alpha)
        type(Ewald_Self_Component_Wrapper), allocatable, intent(out) :: self_components(:)
        class(Abstract_Permittivity), intent(in) :: permittivity
        type(Component_Wrapper), intent(in) :: components(:)
        logical, intent(in) :: are_dipolar(:)
        class(Abstract_Ewald_Convergence_Parameter), intent(in) :: alpha
        integer :: i_component

        allocate(self_components(size(components)))
        do i_component = 1, size(self_components)
            call ewald_self_create(self_components(i_component)%self, permittivity, &
                components(i_component)%dipolar_moments, are_dipolar(i_component), alpha)
        end do
    end subroutine create_self_components

    subroutine destroy_self_components(self_components)
        type(Ewald_Self_Component_Wrapper), allocatable, intent(inout) :: self_components(:)

        integer :: i_component

        if (allocated(self_components)) then
            do i_component = size(self_components), 1, -1
                call ewald_self_destroy(self_components(i_component)%self)
            end do
            deallocate(self_components)
        end if
    end subroutine destroy_self_components

    subroutine create_self(self, permittivity, dipolar_moments, dipoles_exist, alpha)
        class(Abstract_Ewald_Self_Component), allocatable, intent(out) :: self
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Component_Dipolar_Moments), intent(in) :: dipolar_moments
        logical, intent(in) :: dipoles_exist
        class(Abstract_Ewald_Convergence_Parameter), intent(in) :: alpha

        if (dipoles_exist) then
            allocate(Concrete_Ewald_Self_Component :: self)
        else
            allocate(Null_Ewald_Self_Component :: self)
        end if
        call self%construct(permittivity, dipolar_moments, alpha)
    end subroutine create_self

    subroutine destroy_self(self)
        class(Abstract_Ewald_Self_Component), allocatable, intent(inout) :: self

        if (allocated(self)) then
            call self%destroy()
            deallocate(self)
        end if
    end subroutine destroy_self

end module procedures_ewald_self_factory
