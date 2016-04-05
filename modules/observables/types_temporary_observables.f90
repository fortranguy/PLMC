module types_temporary_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64

implicit none

private

    type, public :: Concrete_Single_Delta_Energies
        real(DP) :: field, walls
        real(DP), allocatable :: short(:), dipolar(:)
        real(DP) :: dipolar_mixture
    end type Concrete_Single_Delta_Energies

    type, public :: Concrete_Double_Delta_Energies
        real(DP), dimension(2) :: field, walls
        real(DP), allocatable :: short(:, :), dipolar(:, :)
        real(DP) :: dipolar_mixture
    end type Concrete_Double_Delta_Energies

end module types_temporary_observables
