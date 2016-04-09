module types_radial_distribution

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Radial_Distribution_Component
        integer :: num_particles_sum, num_particles
        real(DP) :: density
        real(DP), allocatable :: positions(:, :)
    end type Concrete_Radial_Distribution_Component

end module types_radial_distribution
