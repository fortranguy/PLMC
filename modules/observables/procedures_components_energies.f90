module procedures_components_energies

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_observables_wrapper, only: Concrete_Components_Energies

implicit none

private
public :: Concrete_Components_Energies_init, Concrete_Components_Energies_add

interface Concrete_Components_Energies_add
    module procedure :: Concrete_Components_Energies_add_line
    module procedure :: Concrete_Components_Energies_add_triangle
end interface Concrete_Components_Energies_add

contains

    pure subroutine Concrete_Components_Energies_init(energies)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)

        integer :: i_component
        do i_component = 1, size(energies)
            energies(i_component)%with_components = 0._DP
        end do
    end subroutine Concrete_Components_Energies_init

    pure subroutine Concrete_Components_Energies_add_line(energies, energies_i)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        real(DP), intent(in) :: energies_i(:)

        integer :: i_component
        do i_component = 1, size(energies)
            energies(i_component)%with_components(i_component) = energies(i_component)%&
                with_components(i_component) + energies_i(i_component)
        end do
    end subroutine Concrete_Components_Energies_add_line

    pure subroutine Concrete_Components_Energies_add_triangle(energies, energies_i)
        type(Concrete_Components_Energies), intent(inout) :: energies(:)
        type(Concrete_Components_Energies), intent(in) :: energies_i(:)

        integer :: i_component
        do i_component = 1, size(energies)
            energies(i_component)%with_components = energies(i_component)%with_components + &
                energies_i(i_component)%with_components
        end do
    end subroutine Concrete_Components_Energies_add_triangle

end module procedures_components_energies
