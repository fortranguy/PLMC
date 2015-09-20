module module_changes

use class_moved_positions, only: Abstract_Moved_Positions
use class_rotated_orientations, only: Abstract_Rotated_Orientations
use class_particles_exchange, only: Particles_Exchange_Facade

implicit none

private
public Concrete_Changes_construct, Concrete_Changes_destroy

    type, public :: Concrete_Changes
        class(Abstract_Moved_Positions), pointer :: moved_positions
        class(Abstract_Rotated_Orientations), pointer :: rotated_orientations
        type(Particles_Exchange_Facade), pointer :: particles_exchange
    end type Concrete_Changes

contains

    subroutine Concrete_Changes_construct(changes, moved_positions, rotated_orientations, &
            particles_exchange)
        type(Concrete_Changes), intent(out) :: changes
        class(Abstract_Moved_Positions), target, intent(in) :: moved_positions
        class(Abstract_Rotated_Orientations), target, intent(in) :: rotated_orientations
        type(Particles_Exchange_Facade), target, intent(in) :: particles_exchange

        changes%moved_positions => moved_positions
        changes%rotated_orientations => rotated_orientations
        changes%particles_exchange => particles_exchange
    end subroutine Concrete_Changes_construct

    subroutine Concrete_Changes_destroy(changes)
        type(Concrete_Changes), intent(inout) :: changes

        changes%particles_exchange => null()
        changes%rotated_orientations => null()
        changes%moved_positions => null()
    end subroutine Concrete_Changes_destroy

end module module_changes
