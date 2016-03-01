module class_dlc_weight

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use data_constants, only: num_dimensions
use class_periodic_box, only: Abstract_Periodic_Box
use class_reciprocal_lattice, only: Abstract_Reciprocal_Lattice

implicit none

private

    type, abstract, public :: Abstract_DLC_Weight
    private
        class(Abstract_Periodic_Box), pointer :: periodic_box => null()
        integer :: reci_numbers(num_dimensions)
        real(DP) :: permittivity
        real(DP), dimension(:, :), allocatable :: weight
    end type Abstract_DLC_Weight

end module class_dlc_weight
