module procedures_metropolis_micro

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_line_observables, only: Concrete_Line_Observables

implicit none

private
public :: update_energies

contains

    subroutine update_energies(energies, deltas, i_actor)
        type(Concrete_Line_Observables), intent(inout) :: energies(:)
        real(DP), intent(in) :: deltas(:)
        integer, intent(in) :: i_actor

        integer :: i_component, j_observable, i_observable

        do i_component = 1, size(deltas)
            j_observable = maxval([i_actor, i_component])
            i_observable = minval([i_actor, i_component])
            energies(j_observable)%with_components(i_observable) = energies(j_observable)%&
                with_components(i_observable) + deltas(i_component)
        end do
    end subroutine update_energies

end module procedures_metropolis_micro
