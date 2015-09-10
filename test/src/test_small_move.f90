program test_small_move

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_periodic_box, only: Abstract_Periodic_Box, XYZ_Periodic_Box
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_positions, only: Abstract_Positions, Concrete_Positions
use class_small_move, only: Abstract_Small_Move, Concrete_Small_Move

implicit none

    class(Abstract_Small_Move), allocatable :: small_move
    class(Abstract_Positions), allocatable :: positions
    class(Abstract_Particles_Number), allocatable :: particles_number
    class(Abstract_Periodic_Box), allocatable :: periodic_box
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: data_found

    integer :: num_particles, i_particle
    character(len=1024) :: string_i
    real(DP), allocatable :: periodic_box_size(:), position(:), small_move_delta(:)

    integer :: positions_unit

    call json_initialize()
    data_filename = "small_move.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)

    allocate(XYZ_Periodic_Box :: periodic_box)
    write(output_unit, *) "XYZ_Periodic_Box"
    data_field = "Periodic Box.size"
    call input_data%get(data_field, periodic_box_size, data_found)
    call test_data_found(data_field, data_found)
    call periodic_box%set_size(periodic_box_size)

    allocate(Concrete_Particles_Number :: particles_number)
    data_field = "Particles.number"
    call input_data%get(data_field, num_particles, data_found)
    call test_data_found(data_field, data_found)
    call particles_number%set(num_particles)

    allocate(Concrete_Positions :: positions)
    call positions%construct(periodic_box, particles_number)

    do i_particle = 1, particles_number%get()
        write(string_i, *) i_particle
        data_field = "Particles."//trim(adjustl(string_i))//".position"
        call input_data%get(data_field, position, data_found)
        call test_data_found(data_field, data_found)
        call positions%set(i_particle, position)
    end do

    allocate(Concrete_Small_Move :: small_move)
    call small_move%construct(positions)
    data_field = "Small Move.delta"
    call input_data%get(data_field, small_move_delta, data_found)
    call test_data_found(data_field, data_found)
    call small_move%set(small_move_delta)

    open(newunit=positions_unit, recl=4096, file="moved_positions.out", action="write")
    do i_particle = 1, particles_number%get()
        write(positions_unit, *) positions%get(i_particle), &
            small_move%get(i_particle) - positions%get(i_particle)
        write(positions_unit, *)
        write(positions_unit, *)
    end do
    close(positions_unit)

    call small_move%destroy()
    deallocate(small_move)
    call positions%destroy()
    deallocate(positions)
    deallocate(particles_number)
    deallocate(periodic_box)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()
    
end program test_small_move
