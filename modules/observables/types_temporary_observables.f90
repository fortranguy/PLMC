module types_temporary_observables

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line

implicit none

private

    type, public :: Concrete_Single_Energies
        real(DP) :: field = 0._DP, walls = 0._DP
        real(DP), allocatable :: short(:), dipolar(:)
        real(DP) :: dipolar_mixture = 0._DP
    end type Concrete_Single_Energies

    type, public :: Concrete_Double_Energies
        real(DP), dimension(2) :: field = 0._DP, walls = 0._DP
        real(DP), allocatable :: short(:, :), dipolar(:, :)
        real(DP) :: dipolar_mixture = 0._DP
    end type Concrete_Double_Energies

    type, public :: Concrete_Energies !in Observables?
        real(DP), allocatable :: field_energies(:), walls_energies(:)
        type(Reals_Line), allocatable :: short_energies(:), dipolar_energies(:)
        real(DP) :: dipolar_mixture_energy = 0._DP
    end type Concrete_Energies

end module types_temporary_observables
