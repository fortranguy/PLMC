module procedures_num_particles_factory

use classes_num_particles, only: Abstract_Num_Particles, Concrete_Num_Particles, &
    Null_Num_Particles

implicit none

private
public :: create_line, create_element, destroy_line, destroy_element

contains

    subroutine create_line(nums_particles, num_components, needed)
        class(Abstract_Num_Particles), allocatable, intent(out) :: nums_particles(:)
        integer, intent(in) :: num_components
        logical, intent(in) :: needed

        if (needed) then
            allocate(Concrete_Num_Particles :: nums_particles(num_components))
        else
            allocate(Null_Num_Particles :: nums_particles(num_components))
        end if
    end subroutine create_line

    !> num_particles will be set with coordinates, cf.
    !> [[classes_component_coordinates_reader:Abstract_Component_Coordinates_Reader]].
    !> @warning Is it too fragile?
    subroutine create_element(num_particles, exists)
        class(Abstract_Num_Particles), allocatable, intent(out) :: num_particles
        logical, intent(in) :: exists

        if (exists) then
            allocate(Concrete_Num_Particles :: num_particles)
        else
            allocate(Null_Num_Particles :: num_particles)
        end if
    end subroutine create_element

    subroutine destroy_line(nums_particles)
        class(Abstract_Num_Particles), allocatable, intent(inout) :: nums_particles(:)

        if (allocated(nums_particles)) deallocate(nums_particles)
    end subroutine destroy_line

    subroutine destroy_element(num_particles)
        class(Abstract_Num_Particles), allocatable, intent(inout) :: num_particles

        if (allocated(num_particles)) deallocate(num_particles)
    end subroutine destroy_element

end module procedures_num_particles_factory
