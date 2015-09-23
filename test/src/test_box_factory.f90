program test_box_factory

use json_module, only: json_file, json_initialize
use module_data, only: test_file_exists
use types_box, only: Box_Wrapper
use class_box_factory, only: Concrete_Box_Factory

implicit none

    type(Concrete_Box_Factory) :: box_factory
    type(Box_Wrapper) :: box

    type(json_file) :: input_data
    character(len=:), allocatable :: data_filename

    call json_initialize()
    data_filename = "box_factory.json"
    call test_file_exists(data_filename)
    call input_data%load_file(filename = data_filename)
    deallocate(data_filename)

    call box_factory%allocate(box, input_data, "PLMC")
    call box_factory%construct(box)
    call box_factory%destroy(box)

    call input_data%destroy()

end program test_box_factory
