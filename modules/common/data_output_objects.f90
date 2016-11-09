module data_output_objects

implicit none

private
public :: report_filename, random_number_generator_object, generating_algorithms_object

    character(len=*), parameter :: report_filename = "report.json"
    character(len=*), parameter :: random_number_generator_object = "Random Number Generator"
    character(len=*), parameter :: generating_algorithms_object = "Generating Algorithms"

end module data_output_objects
