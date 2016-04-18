module data_strings

implicit none

private
public :: max_line_length, max_word_length

    integer, parameter :: max_line_length = 4096
    integer, parameter :: max_word_length = 1024

end module data_strings
