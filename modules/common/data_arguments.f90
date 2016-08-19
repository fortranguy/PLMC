module data_arguments

implicit none

private
public :: i_generating, i_exploring, num_json_arguments

    integer, parameter :: i_generating = 1
    integer, parameter :: i_exploring = 2
    integer, parameter :: num_json_arguments = 2

end module data_arguments
