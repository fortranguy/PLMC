program test_positions

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists, test_data_found
use class_particles_number, only: Abstract_Particles_Number, Concrete_Particles_Number
use class_positions, only: Abstract_Positions, Concrete_Positions

implicit none
    
    class(Abstract_Particles_Number), allocatable :: particles_num
    class(Abstract_Positions), allocatable :: positions
    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename, data_field
    logical :: found
    integer :: positions_num, i_particle
    real(DP), allocatable :: position(:)
    character(len=1024) :: string_i

    call json_initialize()

    data_filename = "positions.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    
    allocate(Concrete_Particles_Number :: particles_num)
    data_field = "Positions.number"
    call input_data%get(data_field, positions_num, found)
    call test_data_found(data_field, found)
    call particles_num%set(positions_num)
    
    allocate(Concrete_Positions :: positions)
    call positions%construct(particles_num)
    
    do i_particle = 1, particles_num%get()
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
    call particles_num%set(particles_num%get() + 1)
    call positions%add(position)
    write(output_unit, *) "position", particles_num%get(), "=", positions%get(particles_num%get())
    
    data_field = "Positions.remove"
    call input_data%get(data_field, i_particle, found)
    call test_data_found(data_field, found)
    call positions%remove(i_particle)
    call particles_num%set(particles_num%get() - 1)
    do i_particle = 1, particles_num%get()
        write(output_unit, *) "position", i_particle, "=", positions%get(i_particle)
    end do

    call positions%destroy()
    deallocate(positions)
    deallocate(particles_num)
    deallocate(data_field)
    deallocate(data_filename)
    call input_data%destroy()

end program test_positions
