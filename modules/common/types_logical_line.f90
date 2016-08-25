module types_logical_line

implicit none

private

    type, public :: Concrete_Logical_Line
        logical, allocatable :: line(:)
    end type Concrete_Logical_Line

end module types_logical_line
