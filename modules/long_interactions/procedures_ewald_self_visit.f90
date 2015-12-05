module procedures_ewald_self_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_component_dipolar_moments, only: Abstract_Component_Dipolar_Moments
use class_ewald_self, only: Abstract_Ewald_Self

implicit none

private
public :: ewald_self_visit

contains

    pure real(DP) function ewald_self_visit(moments, self) result(energy)
        class(Abstract_Component_Dipolar_Moments), intent(in) :: moments
        class(Abstract_Ewald_Self), intent(in) :: self

        integer :: i_particle

        energy = 0._DP
        do i_particle = 1, moments%get_num()
            energy = energy + self%get(moments%get(i_particle))
        end do
    end function ewald_self_visit

end module procedures_ewald_self_visit
