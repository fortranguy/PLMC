program test_box_factory

use, intrinsic :: iso_fortran_env, only: DP => REAL64, output_unit
use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use types_box, only: Box_Wrapper
use procedures_box_factory, only: box_factory_construct, box_factory_destroy

implicit none

    type(Box_Wrapper) :: box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "box_factory.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call box_factory_construct(box, input_data, "System")
    write(output_unit, *) "box size", box%periodic_box%get_size()
    write(output_unit, *) "temperature", box%temperature%get()
    write(output_unit, *) "external field at origin", box%external_field%get([0._DP, 0._DP, 0._DP])
    write(output_unit, *) "reciprocal lattice numbers", box%reciprocal_lattice%get_numbers()
    write(output_unit, *) "walls gap", box%walls_potential%get_gap()

    call box_factory_destroy(box)
    call input_data%destroy()

end program test_box_factory
