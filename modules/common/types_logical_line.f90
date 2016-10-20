module types_logical_line

implicit none

private

    type, public :: Logical_Line
        logical, allocatable :: line(:)
    end type Logical_Line

end module types_logical_line
