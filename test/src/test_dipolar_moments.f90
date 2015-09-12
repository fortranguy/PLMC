module procedures_dipolar_moments_print

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file
use module_data, only: test_data_found
use class_moments_norm, only: Abstract_Moments_Norm
use class_orientations, only: Abstract_Orientations
use class_dipolar_moments, only: Dipolar_Moments_Facade

implicit none

private
public set_objects, print_dipolar_moments

contains

    subroutine set_objects(input_data, moments_norm_name, moments_norm, &
                           orientations_name, orientations)
        type(json_file), intent(inout) :: input_data
        class(Abstract_Moments_Norm), intent(inout) :: moments_norm
        class(Abstract_Orientations), intent(inout) :: orientations
        character(len=*), intent(in) :: moments_norm_name, orientations_name

        character(len=:), allocatable :: data_field
        logical :: data_found
        integer :: i_particle
        character(len=1024) :: string_i
        real(DP) :: moments_norm_value
        real(DP), allocatable :: orientation(:)

        write(output_unit, *) moments_norm_name//" + "//orientations_name

        do i_particle = 1, moments_norm%get_num()
            write(string_i, *) i_particle
            data_field = "Dipoles."//trim(adjustl(string_i))//".moment norm"
            call input_data%get(data_field, moments_norm_value, data_found)
            call test_data_found(data_field, data_found)
            call moments_norm%set(i_particle, moments_norm_value)
        end do

        do i_particle = 1, orientations%get_num()
            write(string_i, *) i_particle
            data_field = "Dipoles."//trim(adjustl(string_i))//".orientation"
            call input_data%get(data_field, orientation, data_found)
            call test_data_found(data_field, data_found)
            call orientations%set(i_particle, orientation)
            deallocate(orientation)
        end do
    end subroutine set_objects

    subroutine print_dipolar_moments(dipolar_moments)
        type(Dipolar_Moments_Facade), intent(in) :: dipolar_moments

        integer :: i_particle

        do i_particle = 1, dipolar_moments%get_num()
            write(output_unit, *) dipolar_moments%get(i_particle)
        end do
    end subroutine print_dipolar_moments

end module procedures_dipolar_moments_print

program test_dipolar_moments

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_number, only: Abstract_Number, Concrete_Number
use class_moments_norm, only: Abstract_Moments_Norm, Null_Moments_Norm, Uniform_Moments_Norm
use class_orientations, only: Abstract_Orientations, Null_Orientations, Concrete_Orientations
use class_dipolar_moments, only: Dipolar_Moments_Facade
use procedures_dipolar_moments_print, only: set_objects, print_dipolar_moments

implicit none

    type(Dipolar_Moments_Facade) :: dipolar_moments
    class(Abstract_Orientations), allocatable :: orientations
    class(Abstract_Moments_Norm), allocatable :: moments_norm
    class(Abstract_Number), allocatable :: number
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found
    integer :: num_particles

    call json_initialize()
    data_filename = "dipolar_moments.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Number :: number)
    data_field = "Dipoles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call number%set(num_particles)

    allocate(Null_Moments_Norm :: moments_norm)
    allocate(Null_Orientations :: orientations)
    call moments_norm%construct(number)
    call orientations%construct(number)
    call set_objects(input_data, "Null_Moments_Norm", moments_norm, &
                     "Null_Orientations", orientations)
    call dipolar_moments%construct(moments_norm, orientations)
    call print_dipolar_moments(dipolar_moments)
    call dipolar_moments%destroy()
    call orientations%destroy()
    call moments_norm%destroy()
    deallocate(orientations)
    deallocate(moments_norm)

    allocate(Uniform_Moments_Norm :: moments_norm)
    allocate(Concrete_Orientations :: orientations)
    call moments_norm%construct(number)
    call orientations%construct(number)
    call set_objects(input_data, "Uniform_Moments", moments_norm, &
                     "Concrete_Orientations", orientations)
    call dipolar_moments%construct(moments_norm, orientations)
    call print_dipolar_moments(dipolar_moments)
    call dipolar_moments%destroy()
    call orientations%destroy()
    call moments_norm%destroy()
    deallocate(orientations)
    deallocate(moments_norm)

    deallocate(number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_dipolar_moments