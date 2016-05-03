module types_radial_distribution_component

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Radial_Distribution_Component
        integer :: num_particles_sum = 0, num_particles = 0
        real(DP) :: density = 0
        real(DP), allocatable :: positions(:, :)
    end type Concrete_Radial_Distribution_Component

end module types_radial_distribution_component
