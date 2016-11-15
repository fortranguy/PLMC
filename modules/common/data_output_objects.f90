module data_output_objects

implicit none

private
public :: generating_report_filename, exploring_report_filename, random_number_generator_object, &
    generating_algorithms_weights_object

    character(len=*), parameter :: generating_report_filename = "generating_report.json"
    character(len=*), parameter :: exploring_report_filename = "exploring_report.json"
    character(len=*), parameter :: random_number_generator_object = "Random Number"
    character(len=*), parameter :: generating_algorithms_weights_object = &
        "Algorithms Weight"

end module data_output_objects
