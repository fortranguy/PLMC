module procedures_triangle_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line

implicit none

private
public :: triangle_observables_init, triangle_observables_add

interface triangle_observables_add
    module procedure :: add_line
    module procedure :: add_triangle
end interface triangle_observables_add

contains

    elemental pure subroutine triangle_observables_init(energies)
        type(Reals_Line), intent(inout) :: energies

        energies%line = 0._DP
    end subroutine triangle_observables_init

    pure subroutine add_line(energies, energies_i)
        type(Reals_Line), intent(inout) :: energies(:)
        real(DP), intent(in) :: energies_i(:)

        integer :: i_component
        do i_component = 1, size(energies)
            energies(i_component)%line(i_component) = energies(i_component)%line(i_component) + &
                energies_i(i_component)
        end do
    end subroutine add_line

    elemental pure subroutine add_triangle(energies, energies_i)
        type(Reals_Line), intent(inout) :: energies
        type(Reals_Line), intent(in) :: energies_i

        energies%line = energies%line + energies_i%line
    end subroutine add_triangle

end module procedures_triangle_observables
