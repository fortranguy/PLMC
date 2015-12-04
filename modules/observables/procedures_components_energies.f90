module procedures_components_energies

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_observables_wrapper, only: Concrete_Components_Energies

implicit none

private
public :: operator(+)

interface assignment(=)
    module procedure :: Concrete_Components_Energies_assignment
end interface

interface operator(+)
    module procedure :: Concrete_Components_Energies_addition
end interface

contains

    pure subroutine Concrete_Components_Energies_assignment(energies_target, energies_value)
        type(Concrete_Components_Energies), intent(out) :: energies_target(:)
        type(Concrete_Components_Energies), intent(in) :: energies_value(:)

        integer :: i_component
        do i_component = 1, size(energies_value)
            energies_target(i_component)%with_components = &
                energies_value(i_component)%with_components
        end do
    end subroutine Concrete_Components_Energies_assignment

    pure function Concrete_Components_Energies_addition(energies_1, energies_2) &
        result(additive_energies)
        type(Concrete_Components_Energies), intent(in) :: energies_1(:), energies_2(:)
        type(Concrete_Components_Energies) :: additive_energies(size(energies_1))

        integer :: i_component
        do i_component = 1, size(energies_1)
            allocate(additive_energies(i_component)%with_components(size(energies_1(i_component)%&
                with_components)))
            additive_energies(i_component)%with_components = energies_1(i_component)%&
                with_components + energies_2(i_component)%with_components
        end do
    end function Concrete_Components_Energies_addition

end module procedures_components_energies
