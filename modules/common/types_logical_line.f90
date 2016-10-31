module types_logical_line

implicit none

private

    type, public :: Logical_Line
        logical, allocatable :: line(:)
    end type Logical_Line

    type, public :: Logical_Triangle
        type(Logical_Line), allocatable :: triangle(:)
    end type Logical_Triangle

end module types_logical_line
