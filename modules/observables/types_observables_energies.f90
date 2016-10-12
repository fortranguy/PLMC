module types_observables_energies

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use types_reals_line, only: Reals_Line

implicit none

private

    type, public :: Concrete_Single_Energies
        real(DP) :: walls_energy = 0._DP, field_energy = 0._DP
        real(DP), allocatable :: short_energies(:), dipolar_energies(:)
        real(DP) :: dipolar_shared_energy = 0._DP
    end type Concrete_Single_Energies

    type, public :: Concrete_Double_Energies
        real(DP), dimension(2) :: walls_energies = 0._DP, field_energies = 0._DP
        real(DP), allocatable :: short_energies(:, :), dipolar_energies(:, :)
        real(DP) :: dipolar_shared_energy = 0._DP
    end type Concrete_Double_Energies

    type, public :: Concrete_Observables_Energies
        real(DP), allocatable :: walls_energies(:), field_energies(:)
        type(Reals_Line), allocatable :: short_energies(:), dipolar_energies(:)
        real(DP) :: dipolar_shared_energy = 0._DP
    end type Concrete_Observables_Energies

end module types_observables_energies
