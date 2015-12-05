module procedures_ewald_self_factory

use class_permittivity, only: Abstract_Permittivity
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_self, only: Abstract_Ewald_Self, Concrete_Ewald_Self, Null_Ewald_Self
use types_long_interactions_wrapper, only: Ewald_Self_Wrapper

implicit none

private
public :: ewald_self_create, ewald_self_destroy

interface ewald_self_create
    module procedure :: create_selves
    module procedure :: create_self
end interface ewald_self_create

interface ewald_self_destroy
    module procedure :: destroy_self
    module procedure :: destroy_selves
end interface ewald_self_destroy

contains

    subroutine create_selves(selves, permittivity, alpha, are_dipolar)
        type(Ewald_Self_Wrapper), allocatable, intent(out) :: selves(:)
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Ewald_Convergence_Parameter), intent(in) :: alpha
        logical, intent(in) :: are_dipolar(:)

        integer :: i_component

        allocate(selves(size(are_dipolar)))
        do i_component = 1, size(selves)
            call ewald_self_create(selves(i_component)%self, permittivity, alpha, &
                are_dipolar(i_component))
        end do
    end subroutine create_selves

    subroutine destroy_selves(selves)
        type(Ewald_Self_Wrapper), allocatable, intent(inout) :: selves(:)

        integer :: i_component

        if (allocated(selves)) then
            do i_component = size(selves), 1, -1
                call ewald_self_destroy(selves(i_component)%self)
            end do
            deallocate(selves)
        end if
    end subroutine destroy_selves

    subroutine create_self(self, permittivity, alpha, dipoles_exist)
        class(Abstract_Ewald_Self), allocatable, intent(out) :: self
        class(Abstract_Permittivity), intent(in) :: permittivity
        class(Abstract_Ewald_Convergence_Parameter), intent(in) :: alpha
        logical, intent(in) :: dipoles_exist

        if (dipoles_exist) then
            allocate(Concrete_Ewald_Self :: self)
        else
            allocate(Null_Ewald_Self :: self)
        end if
        call self%construct(permittivity, alpha)
    end subroutine create_self

    subroutine destroy_self(self)
        class(Abstract_Ewald_Self), allocatable, intent(inout) :: self

        if (allocated(self)) then
            call self%destroy()
            deallocate(self)
        end if
    end subroutine destroy_self

end module procedures_ewald_self_factory
