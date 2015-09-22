module class_particles_chemical_potential

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use procedures_checks, only: check_positive

implicit none

private

    type, abstract, public :: Abstract_Particles_Chemical_Potential
    private
        real(DP) :: density, excess
    contains
        procedure :: set => Abstract_Particles_Chemical_Potential_set
        procedure :: get_density => Abstract_Particles_Chemical_Potential_get_density
        procedure :: get_excess => Abstract_Particles_Chemical_Potential_get_excess
    end type Abstract_Particles_Chemical_Potential

    type, extends(Abstract_Particles_Chemical_Potential), public :: &
        Concrete_Particles_Chemical_Potential

    end type Concrete_Particles_Chemical_Potential

    type, extends(Abstract_Particles_Chemical_Potential), public :: &
        Null_Particles_Chemical_Potential
    contains
        procedure :: set => Null_Particles_Chemical_Potential_set
        procedure :: get_density => Null_Particles_Chemical_Potential_get_density
        procedure :: get_excess => Null_Particles_Chemical_Potential_get_excess
    end type Null_Particles_Chemical_Potential

contains

!implementation Abstract_Particles_Chemical_Potential

    subroutine Abstract_Particles_Chemical_Potential_set(this, density, excess)
        class(Abstract_Particles_Chemical_Potential), intent(inout) :: this
        real(DP), intent(in) :: density, excess

        call check_positive("Abstract_Particles_Chemical_Potential", "density", density)
        this%density = density
        call check_positive("Abstract_Particles_Chemical_Potential", "excess", excess)
        this%excess = excess
    end subroutine Abstract_Particles_Chemical_Potential_set

    pure function Abstract_Particles_Chemical_Potential_get_density(this) result(density)
        class(Abstract_Particles_Chemical_Potential), intent(in) :: this
        real(DP) :: density

        density = this%density
    end function Abstract_Particles_Chemical_Potential_get_density

    pure function Abstract_Particles_Chemical_Potential_get_excess(this) result(excess)
        class(Abstract_Particles_Chemical_Potential), intent(in) :: this
        real(DP) :: excess

        excess = this%excess
    end function Abstract_Particles_Chemical_Potential_get_excess

!end implementation Abstract_Particles_Chemical_Potential

!implementation Null_Particles_Chemical_Potential

    subroutine Null_Particles_Chemical_Potential_set(this, density, excess)
        class(Null_Particles_Chemical_Potential), intent(inout) :: this
        real(DP), intent(in) :: density, excess
    end subroutine Null_Particles_Chemical_Potential_set

    pure function Null_Particles_Chemical_Potential_get_density(this) result(density)
        class(Null_Particles_Chemical_Potential), intent(in) :: this
        real(DP) :: density
        density = 0._DP
    end function Null_Particles_Chemical_Potential_get_density

    pure function Null_Particles_Chemical_Potential_get_excess(this) result(excess)
        class(Null_Particles_Chemical_Potential), intent(in) :: this
        real(DP) :: excess
        excess = 0._DP
    end function Null_Particles_Chemical_Potential_get_excess

!end implementation Null_Particles_Chemical_Potential

end module class_particles_chemical_potential
