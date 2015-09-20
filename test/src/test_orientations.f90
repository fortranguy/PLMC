module procedures_orientations_manipulate

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file
use module_data, only: test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_particles_orientations, only: Abstract_Particles_Orientations

implicit none

private
public manipulate_orientations

contains

    subroutine manipulate_orientations(orientations, input_data)

        class(Abstract_Particles_Orientations), intent(inout) :: orientations
        type(json_file), intent(inout) :: input_data

        class(Abstract_Particles_Number), allocatable :: number
        character(len=:), allocatable :: data_field
        logical :: found
        integer :: orientations_num, i_particle
        real(DP), allocatable :: orientation(:)
        character(len=1024) :: string_i


        data_field = "Orientations.number"
        call input_data%get(data_field, orientations_num, found)
        call test_data_found(data_field, found)
        allocate(Concrete_Particles_Number :: number)
        call number%set(orientations_num)

        call orientations%construct(number)

        do i_particle = 1, orientations%get_num()
            write(string_i, *) i_particle
            data_field = "Orientations."//trim(adjustl(string_i))
            call input_data%get(data_field, orientation, found)
            call test_data_found(data_field, found)
            call orientations%set(i_particle, orientation)
            write(output_unit, *) "orientation", i_particle, "=", orientations%get(i_particle)
        end do

        data_field = "Orientations.add"
        call input_data%get(data_field, orientation, found)
        call test_data_found(data_field, found)
        call number%set(number%get() + 1)
        call orientations%add(orientation)
        write(output_unit, *) "orientation", orientations%get_num(), "=", &
            orientations%get(orientations%get_num())

        data_field = "Orientations.remove"
        call input_data%get(data_field, i_particle, found)
        call test_data_found(data_field, found)
        call orientations%remove(i_particle)
        call number%set(number%get() - 1)
        do i_particle = 1, orientations%get_num()
            write(output_unit, *) "orientation", i_particle, "=", orientations%get(i_particle)
        end do

        deallocate(number)
        deallocate(data_field)

    end subroutine manipulate_orientations

end module procedures_orientations_manipulate

program test_orientations

use, intrinsic :: iso_fortran_env, only: output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_orientations, only: Abstract_Particles_Orientations, &
    Concrete_Particles_Orientations, Null_Particles_Orientations
use procedures_orientations_manipulate, only: manipulate_orientations

implicit none

    class(Abstract_Particles_Orientations), allocatable :: orientations
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()

    data_filename = "orientations.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    write(output_unit, *) "Null"
    allocate(Null_Particles_Orientations :: orientations)
    call manipulate_orientations(orientations, input_data)
    deallocate(orientations)

    write(output_unit, *) "Concrete"
    allocate(Concrete_Particles_Orientations :: orientations)
    call manipulate_orientations(orientations, input_data)
    deallocate(orientations)

    deallocate(data_filename)
    call input_data%destroy()

end program test_orientations
