module data_strings

implicit none

private
public :: max_line_length, max_word_length

    integer, parameter :: max_line_length = 1024
    integer, parameter :: max_word_length = 128

end module data_strings
