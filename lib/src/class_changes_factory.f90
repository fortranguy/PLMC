module class_changes_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64
use json_module, only: json_file
use module_data, only: test_data_found
use class_moved_particles_positions, only: Abstract_Moved_Particles_Positions
use class_rotated_particles_orientations, only: Abstract_Rotated_Particles_Orientations
use class_particles_exchange, only: Particles_Exchange_Facade
use types_changes, only: Changes_Wrapper

implicit none

private
public :: Changes_Wrapper_construct, Changes_Wrapper_destroy



contains

    subroutine Changes_Wrapper_construct(changes, moved_positions, rotated_orientations, &
            particles_exchange)
        type(Changes_Wrapper), intent(out) :: changes
        class(Abstract_Moved_Particles_Positions), target, intent(in) :: moved_positions
        class(Abstract_Rotated_Particles_Orientations), target, intent(in) :: rotated_orientations
        type(Particles_Exchange_Facade), target, intent(in) :: particles_exchange

        changes%moved_positions => moved_positions
        changes%rotated_orientations => rotated_orientations
        changes%particles_exchange => particles_exchange
    end subroutine Changes_Wrapper_construct

    subroutine Changes_Wrapper_destroy(changes)
        type(Changes_Wrapper), intent(inout) :: changes

        changes%particles_exchange => null()
        changes%rotated_orientations => null()
        changes%moved_positions => null()
    end subroutine Changes_Wrapper_destroy

end module class_changes_factory
