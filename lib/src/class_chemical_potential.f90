module class_chemical_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Chemical_Potential
    private
        real(DP) :: density, excess
    contains
        procedure :: set => Abstract_Chemical_Potential_set
        procedure :: get_density => Abstract_Chemical_Potential_get_density
        procedure :: get_excess => Abstract_Chemical_Potential_get_excess
    end type Abstract_Chemical_Potential

    type, extends(Abstract_Chemical_Potential), public :: Concrete_Chemical_Potential

    end type Concrete_Chemical_Potential

contains

    subroutine Abstract_Chemical_Potential_set(this, density, excess)
        class(Abstract_Chemical_Potential), intent(inout) :: this
        real(DP), intent(in) :: density, excess

        call check_positive("Abstract_Chemical_Potential", "density", density)
        this%density = density
        call check_positive("Abstract_Chemical_Potential", "excess", excess)
        this%excess = excess
    end subroutine Abstract_Chemical_Potential_set

    pure function Abstract_Chemical_Potential_get_density(this) result(density)
        class(Abstract_Chemical_Potential), intent(in) :: this
        real(DP) :: density

        density = this%density
    end function Abstract_Chemical_Potential_get_density

    pure function Abstract_Chemical_Potential_get_excess(this) result(excess)
        class(Abstract_Chemical_Potential), intent(in) :: this
        real(DP) :: excess

        excess = this%excess
    end function Abstract_Chemical_Potential_get_excess

end module class_chemical_potential
