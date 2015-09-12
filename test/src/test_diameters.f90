program test_diameters

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_number, only: Abstract_Number, Concrete_Number
use class_diameters, only: Abstract_Diameters, Uniform_Diameters

implicit none

    class(Abstract_Number), allocatable :: number
    class(Abstract_Diameters), allocatable :: diameters
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: diameters_num, i_particle
    real(DP) :: diameters_diameter

    call json_initialize()

    data_filename = "diameters.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(Concrete_Number :: number)
    data_field = "Diameters.number"
    call input_data%get(data_field, diameters_num, found)
    call test_data_found(data_field, found)
    call number%set(diameters_num)

    allocate(Uniform_Diameters :: diameters)
    call diameters%construct(number)
    data_field = "Diameters.diameter"
    call input_data%get(data_field, diameters_diameter, found)
    call test_data_found(data_field, found)
    if (diameters%get_num() > 0) then
        call diameters%set(1, diameters_diameter)
    end if

    do i_particle = 1, diameters%get_num()
        write(output_unit, *) "diameter", i_particle, "=", diameters%get(i_particle)
    end do

    data_field = "Diameters.add"
    call input_data%get(data_field, diameters_diameter, found)
    call test_data_found(data_field, found)
    call number%set(number%get() + 1)
    call diameters%add(diameters_diameter)
    write(output_unit, *) "new diameter", diameters%get_num(), "=", &
        diameters%get(diameters%get_num())

    data_field = "Diameters.remove"
    call input_data%get(data_field, i_particle, found)
    call test_data_found(data_field, found)
    call diameters%remove(i_particle)
    call number%set(number%get() - 1)
    do i_particle = 1, diameters%get_num()
        write(output_unit, *) "diameter", i_particle, "=", diameters%get(i_particle)
    end do

    call diameters%destroy()
    deallocate(diameters)
    deallocate(number)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_diameters
