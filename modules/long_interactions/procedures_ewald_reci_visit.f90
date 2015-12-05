module procedures_ewald_reci_visit

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
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

    !> \[
    !>      U_I = \sum_{\vec{k}} w_\alpha(\vec{k}) |S_I(\vec{k})|^2
    !> \]
    pure real(DP) function visit_intra(weight, structure) result(energy)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), intent(in) :: structure

        integer :: reci_numbers(num_dimensions)
        integer :: n_1, n_2, n_3

        reci_numbers = weight%get_reci_numbers()
        energy = 0._DP
        do n_3 = -reci_numbers(3), reci_numbers(3)
        do n_2 = -reci_numbers(2), reci_numbers(2)
        do n_1 = -reci_numbers(1), reci_numbers(1)
            energy = energy + weight%get(n_1, n_2, n_3) * real(structure%get(n_1, n_2, n_3) * &
                conjg(structure%get(n_1, n_2, n_3)), DP)
        end do
        end do
        end do
        energy = energy / 2._DP
    end function visit_intra

    !> \[
    !>      U_{IJ} = \sum_{\vec{k}} w_\alpha(\vec{k}) S_I(\vec{k}) S_J(\vec{k})^\ast
    !> \]
    pure real(DP) function visit_inter(weight, structure_i, structure_j) result(energy)
        class(Abstract_Ewald_Reci_Weight), intent(in) :: weight
        class(Abstract_Ewald_Reci_Structure), intent(in) :: structure_i, structure_j

        integer :: reci_numbers(num_dimensions)
        integer :: n_1, n_2, n_3

        reci_numbers = weight%get_reci_numbers()
        energy = 0._DP
        do n_3 = -reci_numbers(3), reci_numbers(3)
        do n_2 = -reci_numbers(2), reci_numbers(2)
        do n_1 = -reci_numbers(1), reci_numbers(1)
            energy = energy + weight%get(n_1, n_2, n_3) * real(structure_i%get(n_1, n_2, n_3) * &
                conjg(structure_j%get(n_1, n_2, n_3)), DP)
        end do
        end do
        end do
    end function visit_inter

end module procedures_ewald_reci_visit
