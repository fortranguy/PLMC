module procedures_widom_method_factory

use json_module, only: json_file
use procedures_checks, only: check_data_found
use classes_number_to_string, only: Concrete_Number_to_String
use classes_component_number, only: Abstract_Component_Number, Concrete_Component_Number, &
    Null_Component_Number
use classes_random_coordinates, only: Abstract_Random_Coordinates
use classes_widom_method, only: Abstract_Widom_Method, Concrete_Widom_Method, Null_Widom_Method
use types_physical_model_wrapper, only: Physical_Model_Wrapper

implicit none

private
public :: create, destroy

contains

    subroutine create(widom_method, physical_model, random_position, random_orientation, &
        num_components, measure_inv_pow_activities, exploring_data, prefix)
        class(Abstract_Widom_Method), allocatable, intent(out) :: widom_method
        type(Physical_Model_Wrapper), intent(in) :: physical_model
        class(Abstract_Random_Coordinates), intent(in) :: random_position, random_orientation
        integer, intent(in) :: num_components
        logical, intent(in) :: measure_inv_pow_activities
        type(json_file), intent(inout) :: exploring_data
        character(len=*), intent(in) :: prefix

        class(Abstract_Component_Number), allocatable :: numbers(:)
        integer :: num_particles
        type(Concrete_Number_to_String) :: string
        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: i_component

        if (measure_inv_pow_activities) then
            allocate(Concrete_Widom_Method :: widom_method)
            allocate(Concrete_Component_Number :: numbers(num_components))
            do i_component = 1, size(numbers)
                data_field = prefix//"Component "//string%get(i_component)//".number"
                call exploring_data%get(data_field, num_particles, data_found)
                call check_data_found(data_field, data_found)
                call numbers(i_component)%set(num_particles)
            end do
        else
            allocate(Null_Widom_Method :: widom_method)
            allocate(Null_Component_Number :: numbers(1))
        end if
        call widom_method%construct(physical_model%environment, physical_model%mixture%components, &
            physical_model%short_interactions, physical_model%dipolar_interactions, numbers, &
            random_position, random_orientation)
    end subroutine create

    subroutine destroy(widom_method)
        class(Abstract_Widom_Method), allocatable, intent(inout) :: widom_method

        if (allocated(widom_method)) then
            call widom_method%destroy()
            deallocate(widom_method)
        end if
    end subroutine destroy

end module procedures_widom_method_factory
