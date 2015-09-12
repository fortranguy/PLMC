program test_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_number, only: Abstract_Number, Concrete_Number
use class_positions, only: Abstract_Positions, Concrete_Positions

implicit none

    class(Abstract_Periodic_Box), allocatable :: periodic_box
    class(Abstract_Number), allocatable :: number
    class(Abstract_Positions), allocatable :: positions
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: positions_num, i_particle
    real(DP), allocatable :: periodic_box_size(:), position(:)
    character(len=1024) :: string_i

    call json_initialize()

    data_filename = "positions.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(XYZ_Periodic_Box :: periodic_box)
    write(output_unit, *) "XYZ_Periodic_Box"
    data_field = "Periodic Box.size"
    call input_data%get(data_field, periodic_box_size, found)
    call test_data_found(data_field, found)
    call periodic_box%set_size(periodic_box_size)

    allocate(Concrete_Number :: number)
    data_field = "Positions.number"
    call input_data%get(data_field, positions_num, found)
    call test_data_found(data_field, found)
    call number%set(positions_num)

    allocate(Concrete_Positions :: positions)
    call positions%construct(periodic_box, number)

    do i_particle = 1, positions%get_num()
        write(string_i, *) i_particle
        data_field = "Positions."//trim(adjustl(string_i))
        call input_data%get(data_field, position, found)
        call test_data_found(data_field, found)
        call positions%set(i_particle, position)
        write(output_unit, *) "position", i_particle, "=", positions%get(i_particle)
    end do

    data_field = "Positions.add"
    call input_data%get(data_field, position, found)
    call test_data_found(data_field, found)
    call number%set(number%get() + 1)
    call positions%add(position)
    write(output_unit, *) "position", positions%get_num(), "=", &
        positions%get(positions%get_num())

    data_field = "Positions.remove"
    call input_data%get(data_field, i_particle, found)
    call test_data_found(data_field, found)
    call positions%remove(i_particle)
    call number%set(number%get() - 1)
    do i_particle = 1, positions%get_num()
        write(output_unit, *) "position", i_particle, "=", positions%get(i_particle)
    end do

    call positions%destroy()
    deallocate(positions)
    deallocate(number)
    deallocate(periodic_box)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_positions
