module class_inter_particles_diameter

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use class_particles_diameter, only: Abstract_Particles_Diameter

implicit none

private

    type, abstract, public :: Abstract_Inter_Particles_Diameter
    private
        real(DP) :: offset
        class(Abstract_Particles_Diameter), pointer :: diameter_1, diameter_2
    contains
        procedure :: construct => Abstract_Inter_Particles_Diameter_construct
        procedure :: destroy => Abstract_Inter_Particles_Diameter_destroy
        procedure :: get => Abstract_Inter_Particles_Diameter_get
        procedure :: get_min => Abstract_Inter_Particles_Diameter_get_min
    end type Abstract_Inter_Particles_Diameter

    type, extends(Abstract_Inter_Particles_Diameter), public :: Concrete_Inter_Particles_Diameter

    end type Concrete_Inter_Particles_Diameter

contains

    subroutine Abstract_Inter_Particles_Diameter_construct(this, diameter_1, diameter_2, offset)
        class(Abstract_Inter_Particles_Diameter), intent(out) :: this
        class(Abstract_Particles_Diameter), target, intent(in) :: diameter_1, diameter_2
        real(DP), intent(in) :: offset

        this%diameter_1 => diameter_1
        this%diameter_2 => diameter_2
        this%offset = offset
    end subroutine Abstract_Inter_Particles_Diameter_construct

    subroutine Abstract_Inter_Particles_Diameter_destroy(this)
        class(Abstract_Inter_Particles_Diameter), intent(inout) :: this

        this%diameter_2 => null()
        this%diameter_1 => null()
    end subroutine Abstract_Inter_Particles_Diameter_destroy

    pure function Abstract_Inter_Particles_Diameter_get(this) result(diameter)
        class(Abstract_Inter_Particles_Diameter), intent(in) :: this
        real(DP) :: diameter

        diameter = (this%diameter_1%get() + this%diameter_2%get()) / 2._DP + this%offset
    end function Abstract_Inter_Particles_Diameter_get

    pure function Abstract_Inter_Particles_Diameter_get_min(this) result(min_diameter)
        class(Abstract_Inter_Particles_Diameter), intent(in) :: this
        real(DP) :: min_diameter

        min_diameter = (this%diameter_1%get_min() + this%diameter_2%get_min()) / 2._DP + this%offset
    end function Abstract_Inter_Particles_Diameter_get_min

end module class_inter_particles_diameter
