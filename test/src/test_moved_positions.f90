program test_moved_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use data_geometry, only: num_dimensions
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_number, only: Abstract_Number, Concrete_Number
use class_positions, only: Abstract_Positions, Concrete_Positions
use class_moved_positions, only: Abstract_Moved_Positions, Concrete_Moved_Positions

implicit none

    class(Abstract_Moved_Positions), allocatable :: moved_positions
    class(Abstract_Positions), allocatable :: positions
    class(Abstract_Number), allocatable :: number
    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    integer :: num_steps, i_step, num_particles, i_particle
    character(len=1024) :: string_i
    real(DP), allocatable :: periodic_box_size(:), position(:), moved_positions_delta(:)
    real(DP), dimension(num_dimensions) :: old_position, new_position
    integer, allocatable :: positions_units(:)

    call json_initialize()
    data_filename = "moved_positions.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(XYZ_Periodic_Box :: periodic_box)
    write(output_unit, *) "XYZ_Periodic_Box"
    data_field = "Periodic Box.size"
    call input_data%get(data_field, periodic_box_size, data_found)
    call test_data_found(data_field, data_found)
    call periodic_box%set(periodic_box_size)

    allocate(Concrete_Number :: number)
    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call number%set(num_particles)

    allocate(Concrete_Positions :: positions)
    call positions%construct(periodic_box, number)

    allocate(positions_units(positions%get_num()))
    do i_particle = 1, positions%get_num()
        write(string_i, *) i_particle
        data_field = "Particles."//trim(adjustl(string_i))//".initial position"
        call input_data%get(data_field, position, data_found)
        call test_data_found(data_field, data_found)
        call positions%set(i_particle, position)
        open(newunit=positions_units(i_particle), recl=4096, &
             file="positions_"//trim(adjustl(string_i))//".out", action="write")
    end do

    allocate(Concrete_Moved_Positions :: moved_positions)
    call moved_positions%construct(positions)
    data_field = "Small Move.delta"
    call input_data%get(data_field, moved_positions_delta, data_found)
    call test_data_found(data_field, data_found)
    call moved_positions%set(moved_positions_delta)

    data_field = "Number of Steps"
    call input_data%get(data_field, num_steps, data_found)
    call test_data_found(data_field, data_found)

    do i_step = 1, num_steps
        do i_particle = 1, positions%get_num()
            old_position = positions%get(i_particle)
            call positions%set(i_particle, moved_positions%get(i_particle))
            new_position = positions%get(i_particle)
            write(positions_units(i_particle), *) i_step, old_position, new_position - old_position
        end do
    end do

    do i_particle = positions%get_num(), 1, -1
        close(positions_units(i_particle))
    end do
    deallocate(positions_units)

    call moved_positions%destroy()
    deallocate(moved_positions)
    call positions%destroy()
    deallocate(positions)
    deallocate(number)
    deallocate(periodic_box)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_moved_positions
