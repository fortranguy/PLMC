module procedures_particle_insertion_method_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use classes_num_particles, only: Abstract_Num_Particles
use procedures_composition_factory, only: composition_create => create, composition_destroy => &
    destroy
use types_physical_model_wrapper, only: Physical_Model_Wrapper
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_particle_insertion_method, only: Abstract_Particle_Insertion_Method, &
    Concrete_Particle_Insertion_Method, Null_Particle_Insertion_Method

implicit none

private
public :: create, destroy

contains

    subroutine create(particle_insertion_method, physical_model, random_position, &
        random_orientation, measure_inv_pow_activities, exploring_data, prefix)
        class(Abstract_Particle_Insertion_Method), allocatable, intent(out) :: &
            particle_insertion_method
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Random_Coordinates), intent(in) :: random_position, random_orientation
        logical, intent(in) :: measure_inv_pow_activities
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Num_Particles), allocatable :: nums_particles(:)
        integer :: num_particles
        type(Concrete_Number_to_String) :: string
        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: i_component

        call composition_create(nums_particles, size(physical_model%mixture%components), &
            measure_inv_pow_activities)
        if (measure_inv_pow_activities) then
            allocate(Concrete_Particle_Insertion_Method :: particle_insertion_method)
            do i_component = 1, size(nums_particles)
                data_field = prefix//"Component "//string%get(i_component)//".number of particles"
                call exploring_data%get(data_field, num_particles, data_found)
                call check_data_found(data_field, data_found)
                call nums_particles(i_component)%set(num_particles)
            end do
        else
            allocate(Null_Particle_Insertion_Method :: particle_insertion_method)
        end if

        call particle_insertion_method%construct(physical_model%environment, nums_particles, &
            physical_model%mixture%components, physical_model%short_interactions, physical_model%&
            dipolar_interactions_dynamic, random_position, random_orientation)
        call composition_destroy(nums_particles)
    end subroutine create

    subroutine destroy(particle_insertion_method)
        class(Abstract_Particle_Insertion_Method), allocatable, intent(inout) :: &
            particle_insertion_method

        if (allocated(particle_insertion_method)) then
            call particle_insertion_method%destroy()
            deallocate(particle_insertion_method)
        end if
    end subroutine destroy

end module procedures_particle_insertion_method_factory
