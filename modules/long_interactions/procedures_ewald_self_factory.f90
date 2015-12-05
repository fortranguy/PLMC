module procedures_ewald_self_factory

use class_permittivity, only: Abstract_Permittivity
use class_ewald_convergence_parameter, only: Abstract_Ewald_Convergence_Parameter
use class_ewald_self, only: Abstract_Ewald_Self, Concrete_Ewald_Self, Null_Ewald_Self

implicit none

private
public :: ewald_self_create, ewald_self_destroy

contains

    subroutine ewald_self_create(self, permittivity, alpha, dipoles_exist)
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
    end subroutine ewald_self_create

    subroutine ewald_self_destroy(self)
        class(Abstract_Ewald_Self), allocatable, intent(inout) :: self

        if (allocated(self)) then
            call self%destroy()
            deallocate(self)
        end if
    end subroutine ewald_self_destroy

end module procedures_ewald_self_factory
