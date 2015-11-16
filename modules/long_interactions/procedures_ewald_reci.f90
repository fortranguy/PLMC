module procedures_ewald_reci

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_ewald_reci_weight, only: Abstract_Ewald_Reci_Weight
use class_ewald_reci_structure, only: Abstract_Ewald_Reci_Structure

implicit none

private
public :: ewald_reci_visit

interface ewald_reci_visit
    module procedure :: visit_intra
    module procedure :: visit_inter
end interface ewald_reci_visit

contains

    pure real(DP) function visit_intra(reci_numbers, weight, structure) result(energy)
        integer, intent(in) :: reci_numbers(:)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), intent(in) :: structure

        integer :: n_1, n_2, n_3

        energy = 0._DP
        do n_3 = 1, reci_numbers(3)
        do n_2 = 1, reci_numbers(2)
        do n_1 = 1, reci_numbers(1)
            energy = energy + weight%get(n_1, n_2, n_3) * real(structure%get(n_1, n_2, n_3) * &
                conjg(structure%get(n_1, n_2, n_3)), DP)
        end do
        end do
        end do
        energy = energy / 2._DP
    end function visit_intra

    pure real(DP) function visit_inter(reci_numbers, weight, structure_i, structure_j) &
        result(energy)
        integer, intent(in) :: reci_numbers(:)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), intent(in) :: structure_i, structure_j

        integer :: n_1, n_2, n_3

        energy = 0._DP
        do n_3 = 1, reci_numbers(3)
        do n_2 = 1, reci_numbers(2)
        do n_1 = 1, reci_numbers(1)
            energy = energy + weight%get(n_1, n_2, n_3) * real(structure_i%get(n_1, n_2, n_3) * &
                conjg(structure_j%get(n_1, n_2, n_3)), DP)
        end do
        end do
        end do
    end function visit_inter

end module procedures_ewald_reci
