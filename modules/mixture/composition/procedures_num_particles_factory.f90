module procedures_num_particles_factory

use classes_num_particles, only: Abstract_Num_Particles, Concrete_Num_Particles, &
    Null_Num_Particles

implicit none

private
public :: create, destroy

contains

    !> num_particles will be set with coordinates, cf.
    !> [[classes_component_coordinates_reader:Abstract_Component_Coordinates_Reader]].
    !> @warning Is it too fragile?
    subroutine create(num_particles, exists)
        class(Abstract_Num_Particles), allocatable, intent(out) :: num_particles
        logical, intent(in) :: exists

        if (exists) then
            allocate(Concrete_Num_Particles :: num_particles)
        else
            allocate(Null_Num_Particles :: num_particles)
        end if
    end subroutine create

    subroutine destroy(num_particles)
        class(Abstract_Num_Particles), allocatable, intent(inout) :: num_particles

        if (allocated(num_particles)) deallocate(num_particles)
    end subroutine destroy

end module procedures_num_particles_factory
